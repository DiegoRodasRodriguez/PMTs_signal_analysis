function Tcorr=betacorr(Tmain,Tsecond,box)
% Tmain, Tsecond in bins, centered on 0
% box=[x1 x2 y1 y2];
% correct for particle velocity

T_=Tmain-mean(Tmain);
Tsecond=Tsecond-mean(Tsecond);
%I=find(abs(T_-5)<15);
I=find(T_>box(1) & T_<box(2) & Tsecond>box(3) & & Tsecond<box(4));
T_=T_(I); Tsecond_=Tsecond(I);
[t,p]=uncorr2(T_,Tsecond_,'betacorr');
Tcorr=Tmain-polyval(p,Tsecond);

return


T_=T-mean(T);
I=find(abs(T_-5)<15);
T_=T_(I); TsM_=TsM(I);
[t,p]=uncorr2(T_,TsM_,'betacorr');
T_=T-polyval(p,TsM);

