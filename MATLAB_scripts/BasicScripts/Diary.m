Trpc2=(U(:,19)+U(:,20))/2; Arpc2=(U(:,7)+U(:,8))/2;  DT=U(:,19)-U(:,20);  Iok=find(U(:,19)<(2048) & U(:,20)<2048);% strip de baixo
??? Undefined function or variable U.

Trpc1=(U(:,17)+U(:,18))/2; Arpc1=(U(:,5)+U(:,6))/2;  DT=U(:,17)-U(:,18); Iok=find(U(:,17)<(2048) & U(:,18)<2048);% strip de cima
??? Undefined function or variable U.

figure,histf(Arpc2,0:1000)
??? Undefined function or variable Arpc2.

Ibad=find(U(:,17)>=(2048) | U(:,18)>=2048);
??? Undefined function or variable U.

Igood=find(U(:,19)<(2048) & U(:,20)<2048);
??? Undefined function or variable U.

Icross=find(Arpc2>60);
??? Undefined function or variable Arpc2.

IAtail=Igood(find(Arpc(Igood)<55));
??? Undefined function or variable Igood.

diary off
