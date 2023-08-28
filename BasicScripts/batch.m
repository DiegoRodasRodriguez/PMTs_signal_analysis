Data=[];save collect Data % reset collect
%for i=52:62
for i=53:57	% só para
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(4)
   end
end
csvwrite('V-scan.csv',Data);
return

Data=[];save collect Data % reset collect
for i=45:50
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(4)
   end
end
csvwrite('rate-scan.csv',Data);
return


% must turn on 2nd channel cut
Data=[];save collect Data % reset collect
for i=149:170
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(4)
   end
end
csvwrite('X2-scan.csv',Data);
return


Data=[];save collect Data % reset collect
for i=0:8
   filename=['AllX' sprintf('%d',i) 'dwn'];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(10)
   end
end
csvwrite('AllX.csv',Data);
return


Data=[];save collect Data % reset collect
for i=62:102
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(4)
   end
end
csvwrite('X1-scan.csv',Data);
return



Data=[];save collect Data % reset collect
for i=[134 139 147 141 140 142 143]				% full analysis
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(4)
   end
end
csvwrite('Y2-scan.csv',Data); % select lower strip
return


Data=[];save collect Data % reset collect
for i=[137 138 135 134 139 147 141] % full analysis. Sorted by Y
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(4)
   end
end
csvwrite('Y1-scan.csv',Data); % select upper strip
return

Data=[];save collect Data % reset collect
for i=[136 137 138 135 134 139 147 141 140 142 143] 	% efficiency only. Remove crosstalk cuts
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      Yeff, collect
      %mosaic, pause(4)
   end
end
csvwrite('Y-eff.csv',Data); % select lower strip
return


% must change channels in anabig
Data=[];save collect Data % reset collect
for i=350:353
   filename=['dat_' sprintf('%04d',i)];
   if exist(['data\' filename '.mat'])
      anabig, collect
      mosaic, pause(4)
   end
end
csvwrite('2strips.csv',Data);

return

