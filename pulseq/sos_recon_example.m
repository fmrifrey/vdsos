% basic reconstruction example for stack of spirals data
% by David Frey

%% load the data
safile = './scanarc.h5'; % scan archive file name
[kdata,ktraj,seq_args] = sos_convert_data(safile);
nc = size(kdata,4);

%% set up the 2D NUFFT operator
omega2d = 2*pi*seq_args.fov/seq_args.N*squeeze(ktraj(:,1:2,1,:));
omega2d = permute(omega2d,[1,3,2]);
omega2d = reshape(omega2d,[],2);
nufft_args = {seq_args.N*ones(1,2), 6*ones(1,2), 2*seq_args.N*ones(1,2), ...
        seq_args.N/2*ones(1,2), 'table', 2^10, 'minmax:kb'};
F2d = Gnufft( ...
    true(seq_args.N*ones(1,2)), ...
    [omega2d, nufft_args]);

%% calculate density compensation
W2d = sos.dcf_pipe(F2d);
W = Gdiag(repmat(W2d.arg.diag,[size(ktraj,3),1]));

%% set up the 3D NUFFT operator
omega = 2*pi*seq_args.fov/seq_args.N*ktraj;
omega = permute(omega,[1,3,4,2]);
omega = reshape(omega,[],3);
nufft_args = {seq_args.N*ones(1,3), 6*ones(1,3), 2*seq_args.N*ones(1,3), ...
        seq_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};
F = Gnufft( ...
    true(seq_args.N*ones(1,3)), ...
    [omega, nufft_args]);
b = reshape(kdata,[],nc);

%% calculate density compensation
W = sos.dcf_pipe(F2d);

%% reconstruct w/ coil by coil dc-NUFFT
niter = 10;
par_coils = true;
WA = kronI(nc, W*F, 'parfor', par_coils);
x0 = WA' * b;
A = kronI(nc, F, 'parfor', par_coils);
qp = Reg1(true([seq_args.N*ones(1,3),nc]),'beta',2^12);
x_star = qpwls_pcg1(x0,A,1,b(:),qp.C,'niter',niter);
img_coilwise = reshape(x_star,[seq_args.N*ones(1,3),nc]);