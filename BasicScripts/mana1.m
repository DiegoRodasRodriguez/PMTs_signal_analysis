% Analise multiple runs

global ParList RunList runindex fptr TChList AChList OutFile

% dont forget to comment the necessary lines in anarun
%OutFile='556-564.txt',RunList=[556 557 558 559  560  561  562 563 564 556 557 558 559 560 561 562 563 564];ParList=[100 200 500 1000 1000 1000 100 100 100 100 100 50  20  20  20  10  10  10 ];TChList=[ones(1,9)*33 ones(1,9)*34] ;AChList=TChList*0+22;
%OutFile='569-583.txt';RunList=[569 570 571 572 573 574 575 576 577 578 579 580 581 583]; RunList=[RunList RunList];ParList=[56  56  56  54  54  52  50  48  46  58  60  62  64  66]*100;  ParList=[ParList ParList];TChList=[ones(1,14)*33 ones(1,14)*34];AChList=TChList*0+22;
%OutFile='634-656.txt';RunList=634:656; ParList=1:length(RunList); TChList=RunList*0+34;AChList=TChList*0+22;
OutFile='600-606.txt';RunList=[600 601 602    604    605  605];ParList=[200 400 3000/4 7500/4 3600 7200];  TChList=ParList*0+34;AChList=TChList*0+22;

fptr = fopen(OutFile,'a'); if fptr<0, error('File exists'),end

for runindex=1:length(RunList)

%readrun(['Paulo' num2str(RunList(runindex))])

close all
clear variables
clear functions

global ParList RunList runindex fptr TChList AChList OutFile

eval(['load data\Paulo' num2str(RunList(runindex)) ';']); Tstart=(Ch25+Ch26)/2; 
eval(['Trpc1=Ch' num2str(TChList(runindex)) ';']); 
eval(['Arpc1=Ch' num2str(AChList(runindex)) ';']); 
Nevents=length(Trpc1);
anarun
fprintf(fptr,'\n run%d (T=Ch%d, A=Ch%d) Events=%d/%d, Par=%d : <A>=%d, <T>=%2.3e ps, Eff=%d%%, Res=%d ps',...
	RunList(runindex),TChList(runindex),AChList(runindex),...
	length(Tcorr),Nevents, ParList(runindex),...
	round(MeanA),round(RawMeanT),round(Eff_TOF_Hardware),round(SigmaRpc));

end

fclose(fptr);

OutFile
