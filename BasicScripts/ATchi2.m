function chi2=ATchi2(x,A,T)

A0=x(1);
T0=x(2);
tau=x(3);

chi2=sum((A-ATfun(x,T)).^2);
