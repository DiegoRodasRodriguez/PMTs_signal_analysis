function Tuncorr=uncorr(Tmain, Tsecond,plotit)
%function Tuncorr=uncorr(Tmain, Tsecond,plotit)

if nargin==2, plotit=0; end


[x,n]=histf(Tmain,-1000:1000);

[Ts,Tm]=cut(find(Tmain-(n(x==max(x)))<50),Tsecond,Tmain);
CoefCorr=corrcoef(Ts,Tm); CoefCorr=CoefCorr(1,2)
if abs(CoefCorr)<0.01
	Tuncorr=Tmain;
	return;
end

[Mean, Sigma, Max]=stdfit(Tm);Sigma0=Sigma
p=polyfit(Ts,Tm,2);
%p=polyfit(Ts,Tm,1);
if plotit
figure;plot(Ts,Tm,'.');hold on;plot(Ts,polyval(p,Ts),'ro');drawnow
end
Tuncorr=Tmain-polyval(p,Tsecond);
Tm=Tm-polyval(p,Ts);
CoefCorr=corrcoef(Ts,Tm);CoefCorr=CoefCorr(1,2)
[Mean, Sigma, Max]=stdfit(Tm);Sigma=Sigma

if Sigma0<Sigma
	Tuncorr=Tmain;
end


