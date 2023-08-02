OKflag(Icut)=Icut*0;
Iok=find(OKflag);
Events=length(Iok);
IinEff=find(OKflag & NodetMain);% & NodetI & NodetII);
IinEff2=find(OKflag & NodetMain2 & NodetI2 & NodetII2);
Ineighbors=find(OKflag & NodetMain & (DetI | DetII) );
Doubles=length(Ineighbors)/Events*100;
Eff=(1-length(IinEff)/Events)*100;
Eff2=(1-length(IinEff2)/Events)*100;

LocalOKflag=ones(length(Trpc),1);
LocalOKflag(Icut)=Icut*0;
LocalIok=find(LocalOKflag);
LocalEvents=length(LocalIok);
%LocalIinEff=find(U(LocalIok,TchLeft)>=(2048) & U(LocalIok,TchRight)>=2048);
LocalIinEff=find(LocalOKflag & U(:,TchLeft)>=(2048) & U(:,TchRight)>=2048);
LocalEff=(1-length(LocalIinEff)/Events)*100;

disp('Events	Eff	LocalEvents	LocalEff');
s=sprintf('%d	%4.1f 	%d 		%4.1f',Events,Eff,LocalEvents,LocalEff);
disp(s);
disp(' ')