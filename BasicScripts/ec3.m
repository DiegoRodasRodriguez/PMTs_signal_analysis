function [x1,x2,x3]=ec3(A,B,C,D)
%
%Resuelve ecuaciones de trercer grado donde la forma canonica de la 
%ecuacion es: Ax^3+Bx^2+CX+D

q=((B^3)/(27*(A^3)))-((B*C)/(6*(A^2)))+(D/(2*A));
p=((3*A*C)-(B^2))/(9*(A^2));

Dverification=q^2+p^3;

u=(-q+sqrt(q^2+p^3))^(1/3);
v=(-q-sqrt(q^2+p^3))^(1/3);
e1=-0.5+0.86602540378444i;
e2=-0.5-0.86602540378444i;

y1=u+v;
y2=e1*u+e2*v;
y3=e2*u+e1*v;

x1=y1-(B/(3*A));
x2=y2-(B/(3*A));
x3=y3-(B/(3*A));

return