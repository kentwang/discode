%% This script tests all functions related to dIBP sampling
tic();

%% load package
pkg load statistics;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fix random seed
seed = 20160209;
rand ("seed", seed)

% initialization (tell the difference between initialization and ground truth)
sigma_w = .1;
nuep = -.2;
alpha = 1; beta = 1;
K_inf = 15; L_inf = 15;
a = 2; b = 5; c = 1;

% synthetic data
I = 40; J = 60;
K = 4; L = 6;
M = max(K, L);
W_true = [[1, 1, 0, 0, 0, 0]; 
    [0, 0, 1, 0, 0, 1]; 
    [0, 0, 0, 1, 1, 1]; 
    [0, 1, 1, 0, 0, 0]];
[U_true, V_true, Z_true, X] = synthetic(W_true, I, J, K, L, nuep);

% initialization
U = U_true; V = V_true; W = W_true; Z = Z_true;

[lambda_u, lambda_v] = bibetarnd(M, a, b, 1);
mu_u= cumprod(lambda_u);
mu_v = cumprod(lambda_v);
mu = [mu_u, mu_v];

% open gnuplot before plotting
% initiate graph

% x = 1:10;
% plot(sin(x));
% clf;
% subplot(2,2,1); imagesc(U_true); colormap(gray); title('True feature U'); axis off
% subplot(2,2,2); imagesc(V_true); colormap(gray); title('True feature V'); axis off
% subplot(2,2,3); imagesc(W_true); colormap(gray); title('True weight W'); axis off
% subplot(2,2,4); imagesc(X); colormap(gray); title('True bipartite network'); axis off
% print("./output/figure/data_synthetic.jpg", "-djpg");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test sampling of MU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test samplemu.m 
printf('Sampling MU ...\n');
[mu_u, mu_v] = samplemu(U, V, mu_u, mu_v, a, b);
mu = match_with_zero(mu_u, mu_v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test sampling of U and V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test sampleUV.m
printf('Updating features for U and V ...\n');
[U, V, W] = sampleUV(Z, U, V, W, nuep, a, b, sigma_w); 
% TEST OK. BUT A LITTLE. BUT A LITTLE TOO LONG

% clf;
% subplot(3, 1, 1); imagesc(U); colormap(gray); title('Updated U'); axis off
% subplot(3, 1, 2); imagesc(V); colormap(gray); title('Updated V'); axis off
% subplot(3, 1, 3); imagesc(W); colormap(gray); title('Expanded W'); axis off
% print("./output/figure/UVW.jpg", "-djpg");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test sampling of U and V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf('Updating weight matrix W ...\n');
W = sampleW(Z, U, V, nuep, sigma_w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test sampling of U and V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf('Updating auxiliary data Z ...\n');
Z = sampleZ(X, U, V, W, nuep);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test updating of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a, b] = updateAB(a, b, mu);
% warning: division by zero
% warning: called from
%     updateAB at line 17 column 9
%     dIBP_test.m at line 96 column 6
sigma_w = updateSigw(sigma_w, U, V, W, Z, nuep);
samplenuep(U, V, W, Z);

toc();







