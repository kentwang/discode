% MCMC inference

% fix random seed
seed = 20160209;

% synthetic data
I = 40; J = 60;
K = 4; L = 6;
[U_true, V_true, Z_true, X] = synthetic(W, I, J, K, L);
W_true = [[1, 1, 0, 0, 0, 0]; 
    [0, 0, 1, 0, 0, 1]; 
    [0, 0, 0, 1, 1, 1]; 
    [0, 1, 1, 0, 0, 0]];

subplot(2,2,1); imagesc(U_true); colormap(gray); axis off
subplot(2,2,2); imagesc(V_true); colormap(gray); axis off
subplot(2,2,3); imagesc(W_true); colormap(gray); axis off
subplot(2,2,4); imagesc(X); colormap(gray); axis off

% MCMC sampling
E = 1000;
Burn_in =100;
SAMPLE_SIZE = 1000;

% initialization
sigma_w = .1;
nuep = -.2;
alpha = 1; beta = 1;
K_inf = 15; L_inf = 15;





