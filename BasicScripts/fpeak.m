function  [centro]=fpeak(dat);


[N,X]=histf(dat,10:10:1980);peaks=find(N==max(N));centro=mean(peaks);
