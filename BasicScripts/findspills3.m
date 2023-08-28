function [InewSpill, LIspillNr]=findspills2(Clock,plotit)
% function [InewSpill, LIspillNr]=findspills(Clock)
% numbers the spills.

LI=diff(Clock)>1;
InewSpill=find(LI)+1;
LIspillNr=cumsum(LI)+1;

if exist('plotit')
   figure
   plot(Clock,'.'); hold on; 
   plot([InewSpill'; InewSpill'],[(InewSpill'*0); (InewSpill'*0+max(yaxis))]);
   plot(LIspillNr);
   %hold off; plot(InewSpill,'.');
end