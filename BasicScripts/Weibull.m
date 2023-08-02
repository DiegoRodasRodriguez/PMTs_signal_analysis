function W = Weibull(x, N, b, C, d)

W = (N*x.^b).*exp(-C*x.^d);