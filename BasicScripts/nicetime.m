function nicetime(Tdist);
% nice time distribution plot
[Mean, Sigma, Max]=stdfit(Tdist)
SigmaEffective=sqrt(Sigma^2-35^2)
x=-1000:20:1000;
y=[Tdist; x'];
[N,X]=histf(y,x);
figure;stairs(X,N); hold on;lastline('color','k')
x2=(-1.5*Sigma:5:1.5*Sigma)+Mean;
plot(x2,g(x2,Mean,Sigma)*Max*4,'k')
lastline('linewidth',2);logy;
%plot(-3*Sigma:5:3*Sigma,g(-3*Sigma:5:3*Sigma,Mean,Sigma)*Max*4,'--b')
plot(x,g(x,Mean,Sigma)*Max*4,'--k')
yaxis(1,2*max(N));
set(gca,'Xtick',[-1000 -500 -300 0 300 500 1000])
set(gca,'FontSize',[14])

text(min(xaxis)*.9,max(yaxis)/2,['Sigma=' num2str(round(SigmaEffective)) 'ps'])
lastline('FontSize',[14])
