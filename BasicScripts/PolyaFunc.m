% Polya function

function Y = PolyaFunc(par, N)

theta = par(1);
Nmean = par(2);

Y = 1/Nmean*((theta+1)^(theta+1))/gamma(theta+1)...
    .* ((N/Nmean).^theta) .* exp(-(theta+1)*N/Nmean);
Y = Y';
end