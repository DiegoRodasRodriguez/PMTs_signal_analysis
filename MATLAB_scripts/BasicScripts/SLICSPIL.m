% select the first events from spill

Events=length(Arpc)

timeslice=timeslice

if timeslice==1, return,end

TEvs=sort(Tevent);

TimeSlice=find(Tevent<TEvs(timeslice*Events));
[TrpcT1,Arpc,Tstart]=cut(TimeSlice,TrpcT1,Arpc,Tstart);
clear TEvs TimeSlice

Events=length(Arpc)

