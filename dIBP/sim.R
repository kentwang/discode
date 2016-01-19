source("function.R")
set.seed(01192016)

I = 10
J = 15
K = 5
L = 10
M = max(K, L)
sigma.w = 1 # standard Normal

# hyper pars for biBeta dist. Need some tuning.
a = rgamma(1, 2)
b = rgamma(1, 1)
c = 1

# ------------------------------
# dIBP Feature Generation
# ------------------------------

# bivariate beta dist. [Olkin and Liu, 2003]
U.bb = rgamma(M, a)
V.bb = rgamma(M, b)
W.bb = rgamma(M, c)

lambda.u = U.bb / (U.bb + W.bb)
lambda.v = V.bb / (V.bb + W.bb) 

# generate row and column features U and V
U = matrix(rep(0, I * K), ncol = K)
V = matrix(rep(0, J * L), ncol = L)

for(k in 1:K) {
  mu.u = prod(lambda.u[1:k])
  cat("mu.u ", mu.u, "\n")
  U[, k] = rbinom(I, 1, mu.u)
}

for(l in 1:L) {
  mu.v = prod(lambda.v[1:l])
  cat("mu.v", mu.v, "\n")
  V[, l] = rbinom(J, 1, mu.v)
}

# ------------------------------
# Likelihood Function
# ------------------------------

# generate the interaction weight W using standar normal
W = matrix(rnorm(K * L), ncol = L)

# generate bipartite observation X_{IJ}
X = sigmoid(U %*% W %*% t(V)) > 0.5

dev.new(width=20, height=25)
plot(rep(1:I, each = J), rep(-(1:J), I), axes = FALSE, ann = FALSE, 
pch = ifelse(X, 15, 0), cex = 6, asp = 1, xpd = NA)