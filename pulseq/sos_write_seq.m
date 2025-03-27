function sos_write_seq(varargin)
% creates the pulseq files for stack-of-(variable density) spirals sequence
%
% by David Frey (djfrey@umich.edu)
%
% inputs:
% te - echo time (s) ('min' to shrink to minimum)
% tr - repetition time (s) ('min' to shrink to minimum)
% fov - field of view (cm)
% F - vds FOV coefficients (cm)
% dt - raster time (s)
% N - 2D matrix size
% nint - number of interleaves (2D in-plane rotations)
% nprj - number of projections (kz blips)
% gmax - max gradient amplitude (G/cm)
% smax - max slew rate (G/cm/s)
% plotseq - option to plot the sequence
% pislquant - number of TRs to use for prescan
% writepge - option to convert seq to pge file
%
% output files:
% sos.seq file - seq file for pulseq
% sos.pge file - pge file for pge2 interpreter
% seq_args.mat - .mat file containing copy of input arguments
%

    % set default arguments
    arg.te = 'min'; % extra delay for TE (s)
    arg.tr = 200e-3; % extra delay for TR (s)
    arg.fov = 16; % fov (cm)
    arg.F = []; % vds fov coefficients (frac of nominal fov)
    arg.dt = 4e-6; % raster time (s)
    arg.N = 128; % 3D matrix size
    arg.oversamp = 1.5;
    arg.nprj = 131; % number of projections (kz blips)
    arg.nint = 4; % number of interleaves
    arg.slthick = 16; % slab thickness (cm)
    arg.gmax = 4; % max gradient amplitude (G/cm)
    arg.smax = 12000; % max slew rate (G/cm/s)
    arg.writepge = true;
    arg.pislquant = 1;
    arg.plotseq = false; % option to plot the sequence
    
    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % set system limits
    sys = mr.opts('MaxGrad', arg.gmax*10, 'GradUnit', 'mT/m', ...
        'MaxSlew', arg.smax*1e-2, 'SlewUnit', 'mT/m/ms', ...
        'rfDeadTime', 100e-6, ...
        'rfRingdownTime', 60e-6, ...
        'adcRasterTime', arg.dt, ...
        'gradRasterTime', arg.dt, ...
        'rfRasterTime', arg.dt, ...
        'blockDurationRaster', 4e-6, ...
        'B0', 3, ...
        'adcDeadTime', 0e-6);
    
    % initialize sequence object
    seq = mr.Sequence(sys);
    warning('OFF', 'mr:restoreShape');
    
    % Create 90 degree non-selective excitation pulse
    [rf,g_slicesel] = mr.makeSincPulse(pi/2, ...
        'system', sys, ...
        'Duration', 3e-3,...
        'use', 'excitation', ...
        'SliceThickness', arg.slthick*1e-2, ...
        'apodization', 0.5, ...
        'timeBwProduct', 4, ...
        'system',sys);
    g_refoc = mr.makeTrapezoid('z',sys, ...
        'Area', -g_slicesel.area/2, ...
        'system', sys);

    % create gz blip (for kz max)
    gz_blip = mr.makeTrapezoid('z', sys, ...
        'Area', arg.N/(arg.fov*1e-2)/2, ...
        'system', sys);

    % determine echo time delay
    te_min = mr.calcDuration(rf)/2 + mr.calcDuration(g_refoc) + ...
        mr.calcDuration(gz_blip);
    if strcmpi(arg.te,'min')
        arg.te = te_min;
        te_delay = 0;
    elseif arg.te >= te_min
        te_delay = arg.te - te_min;
    else
        error('echo time >= %.3fms', te_min*1e3);
    end
    
    % form initial 2D spiral gradients
    if isempty(arg.F)
        arg.F = [1,0]; % default - archimedian spiral
    end
    arg.F = [arg.F(:).',0];
    [~,g] = sos.vds(arg.smax, ...
        arg.gmax, ...
        arg.dt, ...
        arg.nint, ...
        arg.fov*arg.oversamp*arg.F, ...
        arg.N/arg.fov/2);
    g_wav = sys.gamma * 1e-2 * g; % kHz/m
    acq_len = sys.adcSamplesDivisor*ceil(length(g_wav)/sys.adcSamplesDivisor);
    
    % append ramp
    nramp = ceil(arg.gmax / arg.smax / arg.dt);
    ramp_down = (1 - linspace(0,1,nramp));
    g_wav = [g_wav, g_wav(end)*ramp_down];
    G0 = [real(g_wav); imag(g_wav)];
    
    % create ADC
    adc = mr.makeAdc(acq_len, ...
        'Duration', sys.adcRasterTime*acq_len, ...
        'Delay', te_delay, ...
        'system', sys);
    
    % create spoiler
    gz_spoil = mr.makeTrapezoid('z', sys, ...
        'Area', arg.N/(arg.fov*1e-2)*4, ...
        'system', sys);

    % determine repetition time delay
    tr_min = mr.calcDuration(rf) + mr.calcDuration(g_refoc) + ...
        mr.calcDuration(gz_blip) + te_delay + acq_len*arg.dt + ...
        mr.calcDuration(gz_spoil);
    if strcmpi(arg.tr,'min')
        arg.tr = tr_min;
        tr_delay = 0;
    elseif arg.tr >= tr_min
        tr_delay = arg.tr - tr_min;
    else
        error('repetition time >= %.3fms', tr_min*1e3);
    end
    
    % define sequence blocks
    for iint = 1:arg.nint
        for iprj = 1:arg.nprj

            % write the excitation to sequence
            seq.addBlock(rf, g_slicesel, mr.makeLabel('SET', 'TRID', 1));
            seq.addBlock(g_refoc);

            % get kz fraction for partition
            kz_frac = (-1)^(iprj)*2*floor(iprj/2)/(arg.nprj-1);
            seq.addBlock(mr.scaleGrad(gz_blip,kz_frac));

            % rotate the gradients by golden angle per interleaf
            R = eul2rotm((iint-1)*pi*(3-sqrt(5))*[1,0,0],"ZYX");
            iG = R(1:2,1:2) * G0;

            % write gradients to sequence
            gx_sp = mr.makeArbitraryGrad('x', 0.99*iG(1,:), ...
                'Delay', te_delay, 'system', sys, 'first', 0, 'last', 0);
            gy_sp = mr.makeArbitraryGrad('y', 0.99*iG(2,:), ...
                'Delay', te_delay, 'system', sys, 'first', 0, 'last', 0);
            seq.addBlock(gx_sp, gy_sp, adc);

            % add spoiler and delay
            seq.addBlock(gz_spoil, mr.makeDelay(tr_delay));

        end
    end
    
    % check whether the timing of the sequence is correct
    [ok, error_report] = seq.checkTiming;
    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end
    
    % write out sequence and save args
    seq.setDefinition('FOV', arg.fov*1e-2*ones(1,3));
    seq.setDefinition('Name', 'sos');
    seq.write('sos.seq');
    if arg.writepge
        ceq = seq2ceq('sos.seq');
        writeceq(ceq, 'sos.pge', 'pislquant', arg.pislquant);
    end
    save seq_args.mat -struct arg
    
    % the sequence is ready, so let's see what we got
    if arg.plotseq
        seq.plot();
    end

end