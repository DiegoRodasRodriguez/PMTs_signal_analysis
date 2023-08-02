function [A,B]=complexRaiz(a,b,n)

k=0;
I=0+1*i;
ro=sqrt(a*a+b*b);
if(a<0)
   fi=atan(b/a)+pi;
else
   fi=atan(b/a);
end

if((a<0 & b<0))
   if(n==2)k=1;end
   if(n==3)k=2;end
end


A=real((exp(log(ro)/n)*(cos((fi+2*k*pi)/n)+I*sin((fi+2*k*pi)/n))));
B=imag((exp(log(ro)/n)*(cos((fi+2*k*pi)/n)+I*sin((fi+2*k*pi)/n))));


return


