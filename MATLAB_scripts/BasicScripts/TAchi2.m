function chi2=TAchi2(x,A,T)

A0=x(1);
T0=x(2);
tau=x(3);

chi2=sum((T-TAfun(x,A)).^2);
