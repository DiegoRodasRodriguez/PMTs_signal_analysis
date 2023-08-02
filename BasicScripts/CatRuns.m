clear all
close all

clear all; BigU=[];
for i=164:168						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
   s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:4000);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX0dwn.mat U;
return

clear all; BigU=[];
for i=170:-1:166						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
   s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:4000);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX1dwn.mat U;

clear all; BigU=[];
for i=165:-1:161						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:4000);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX2dwn.mat U;

clear all; BigU=[];
for i=[149 151:153]						% Time here different from rest
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:4000);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX3dwn.mat U;

clear all; BigU=[];
for i=154:156						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:4000);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX4dwn.mat U;

clear all; BigU=[];
for i=157:160						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:4000);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX5dwn.mat U;

clear all; BigU=[];
for i=[170:-1:161]						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:1800);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX6dwn.mat U;

clear all; BigU=[];
for i=[154:160]						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:1800);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX7dwn.mat U;


clear all; BigU=[];
for i=[170:-1:161 154:160]						% strip B -  left to right
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   I=randperm(length(U));I=I(1:1000);
   BigU=[BigU;U(I,:)];   
end
U=BigU;
save data\AllX8dwn.mat U;
return




for i=[134 139 147 141 140 142 143]		% Sorted by Y	
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   BigU=[BigU;U(1:5000,:)];   
end
U=BigU;
save data\AllYdown.mat U;
return


for i=[137 138 135 134 139 147 141] % Sorted by Y
   filename=['dat_' sprintf('%04d',i)];
	s=['load data\' filename], eval(s);
   BigU=[BigU;U(1:5000,:)];   
end
U=BigU;
save data\AllYup.mat U;
return


