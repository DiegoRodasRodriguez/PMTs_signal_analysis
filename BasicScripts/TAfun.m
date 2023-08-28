function T=TAfun(x,A)

A0=x(1);
T0=x(2);
tau=x(3);

I1=find(A>A0);
I2=find(A==A0);
I3=find(A<A0);

T=A*0;

T(I1)=T0+tau*log(A(I1)-A0);
T(I2)=T0+I2*0;
T(I3)=T0-tau*log(A0-A(I3));

