%Function for creating average in N steps

function [Yav, Xav]=average(Y, X, N)

Icount = 1;
i      = 1;
while ((Icount+N)<length(X))
    Xav_ = 0;
    Yav_ = 0;
    for IcountTmp = Icount:(Icount + N - 1)
        Xav_ = Xav_ + X(IcountTmp)/N; %#ok<*AGROW>
        Yav_ = Yav_ + Y(IcountTmp)/N;
    end
    Xav(i) = Xav_;
    Yav(i) = Yav_;
    Icount     = Icount + N;
    i          = i      + 1;
end
Xav = [Xav mean(X(Icount:length(X)))];
Yav = [Yav mean(Y(Icount:length(X)))];
    
return
