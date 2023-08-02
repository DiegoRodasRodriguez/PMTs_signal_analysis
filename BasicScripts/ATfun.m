function A=ATfun(x,T)

A0=x(1);
T0=x(2);
tau=x(3);

A=A0+exp((T-T0)/tau);
