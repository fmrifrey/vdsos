function [kdata,ktraj,seq_args] = sos_convert_data(safile, h5file)
% gets data from scanarchive file, and formats the kspace trajectories and
% sequence arguments based on seq_args.mat, then writes formatted data to
% h5 file for external recon if specified
% by David Frey (djfrey@umich.edu)
%
% inputs:
% safile - name of scanarchive file to read in
% h5file - name of output h5 file to write formatted data to (leave empty
% to not write to file)
%
% outputs:
% kdata - kspace data (ndat x nprj x nint x ncoil)
% ktraj - 3D kspace locations (ndat x nprj x nint x 3)
% seq_args - struct containing pulse sequence arguments
%

    % get directory and file names
    d = dir(safile);
    sadir = d(1).folder;
    safile = d(1).name;

    % load in sequence arguments
    seq_args = load([sadir,'/seq_args.mat']);
    
    % load first shot and get data size
    archive = GERecon('Archive.Load', [sadir,'/',safile]);
    shot = GERecon('Archive.Next', archive);
    [ndat,nc] = size(shot.Data);
    
    % load data
    kdata = zeros(ndat, nc, seq_args.nint*seq_args.nprj);
    kdata(:, :, 1) = shot.Data;
    for l = 2:seq_args.nint*seq_args.nprj
        shot = GERecon('Archive.Next', archive);
        kdata(:, :, l) = shot.Data;
    end
    kdata = reshape(kdata,[ndat,nc,seq_args.nprj,seq_args.nint]);
    kdata = permute(kdata,[1,3,4,2]); % ndat x nprj x nint x nc

    % generate kspace trajectory
    k0 = sos.vds(seq_args.smax, ...
        seq_args.gmax, ...
        seq_args.dt, ...
        seq_args.nint, ...
        seq_args.fov*seq_args.oversamp*seq_args.F, ...
        seq_args.N/seq_args.fov/2);
    k0 = padarray([real(k0(:)), imag(k0(:))], [0,1], ...
        seq_args.N/seq_args.fov/2, 'post');
    
    % transform kspace trajectory
    ktraj = zeros(ndat,3,seq_args.nprj,seq_args.nint);
    for iint = 1:seq_args.nint
        for iprj = 1:seq_args.nprj
            % get kz fraction for partition
            kz_frac = (-1)^(iprj)*2*floor(iprj/2)/(seq_args.nprj-1);

            % rotate the gradients by golden angle per interleaf
            R = eul2rotm((iint-1)*pi*(3-sqrt(5))*[1,0,0],"ZYX");

            ktraj(:,:,iprj,iint) = (k0 * R') .* [1,1,kz_frac];
        end
    end

    % save h5 file
    if nargin > 1 && ~isempty(h5file)

        if isfile(h5file)
            system(sprintf('rm %s',h5file));
        end

        % save kspace data
        h5create(h5file, '/kdata/real', size(kdata), ...
            'Datatype', class(real(kdata)));
        h5write(h5file, '/kdata/real', real(kdata));
        h5create(h5file, '/kdata/imag', size(kdata), ...
            'Datatype', class(imag(kdata)));
        h5write(h5file, '/kdata/imag', imag(kdata));

        % save sampling locations
        h5create(h5file, '/ktraj', size(ktraj), 'Datatype', class(ktraj));
        h5write(h5file, '/ktraj', ktraj);

        % save sequence arguments
        seq_args_fields = fieldnames(seq_args);
        for i = 1:numel(seq_args_fields)
            field = seq_args_fields{i};
            val = seq_args.(field);
            if islogical(val)
                val = 1*val;
            end
            h5create(h5file, sprintf('/seq_args/%s',field), size(val), ...
                'Datatype', class(val));
            h5write(h5file, sprintf('/seq_args/%s',field), val)
        end

    end

end