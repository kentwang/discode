
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to do list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - set upper limit of K_inf and L_inf
% - clarify the goals: what matrices to compare
% - improve the initialization of feature matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load package
pkg load statistics;

% fix random seed
% test random seed trand.m and bibetarnd.m are OK

seed = 20160209;
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

[U_true, V_true, Z_true, X] = synthetic(W_true, I, J, K, L, -.2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize the chain
E = 100; % length of MCMC sample
BURN_IN = 0;
SAMPLE_SIZE = 100; % number of Monte Carlo sample

sigma_w = .1;
nuep = -.1;
a = 2; b = 5;
K_inf = 15; L_inf = 15;

U = U_true + randn(size(U_true)) / 10;
V = V_true + randn(size(V_true)) / 10;
Z = Z_true + randn(size(Z_true)) / 10;
W = W_true + randn(size(W_true)) / 10;

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
for e=1:E
	
end




