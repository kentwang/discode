
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to do list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - set upper limit of K_inf and L_inf
% - clarify the goals: what matrices to compare
% - improve the initialization of feature matrices
% - rebuild synthetic.m to mimic the buffet process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load package
% pkg load statistics;
pkg load all;

% fix random seed
% test random seed trand.m and bibetarnd.m are OK

seed = 20160410;
rand("seed", seed); randn("seed", seed); randg("seed", seed);


% trandn(0, Inf)
% [X, Y] = bibetarnd(3, 1, 1, 1)
% X
% Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_true = [[1, 1, 0, 0, 0, 0]; 
    [0, 0, 1, 0, 0, 1]; 
    [0, 0, 0, 1, 1, 1]; 
    [0, 1, 1, 0, 0, 0]];

I = 15; J = 20;
K = 4; L = 6;

[U_orig, V_orig, Z_orig, X] = synthetic(W_true, I, J, K, L, -2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize the chain
E = 1000; % length of MCMC sample
BURN_IN = 0;
SAMPLE_SIZE = 1000; % number of Monte Carlo sample

sigma_w = 1;
nuep = 1;
a = 1; b = 1;
a_sigw = 1; b_sigw = 1;
K_inf = 10; L_inf = 10;

% use sampledIBP for U and V
[U, V, K_plus, L_plus] = sampleDIBP(a, b, I, J);
% use normal for W
W = W_true + randn(size(W_true)) / 10;
Z = Z_orig + randn(size(Z_orig)) / 10;
mu_u = sort(mean(U, 1)', 'descend');
mu_v = sort(mean(V, 1)', 'descend');

chain.U = zeros(SAMPLE_SIZE, I, K_inf);
chain.V = zeros(SAMPLE_SIZE, J, L_inf);
chain.Z = zeros(SAMPLE_SIZE, I, J);
chain.W = zeros(SAMPLE_SIZE, K, L); 
chain.K = zeros(SAMPLE_SIZE, 1);
chain.L = zeros(SAMPLE_SIZE, 1);
chain.sigma_w = zeros(SAMPLE_SIZE, 1);
chain.nuep = zeros(SAMPLE_SIZE, 1);
chain.a = zeros(SAMPLE_SIZE, 1);
chain.b = zeros(SAMPLE_SIZE, 1);

%% MCMC
s_counter = 0;
for e=1:20

	[mu_u, mu_v] = samplemu(U, V, mu_u, mu_v, a, b);
	[U, V, W] = sampleUV(Z, U, V, W, nuep, a, b, sigma_w);
	W = sampleW(Z, U, V, nuep, sigma_w);
	[K, L] = size(W);
	Z = sampleZ(X, U, V, W, nuep);
	[a, b] = updateAB(a, b, mu_u, mu_v);
	% sigma_w = updateSigw(sigma_w, U, V, W, Z, nuep);
	[sigma_w, a_sigw, b_sigw] = sampleSigw(a_sigw, b_sigw, K, L, W);
	nuep = sampleNuep(U, V, W, Z);

	printf('iter %d: a = %f, b = %f, nuep = %f, sigma_w = %f \n', e, a, b, nuep, sigma_w);
	printf('K and L: %d, %d\n', K, L);
end













