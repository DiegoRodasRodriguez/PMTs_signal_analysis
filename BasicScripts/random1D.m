%disp('                                                                       ')
%disp('           *---------------------------------------------------------* ')
%disp('           *          Function for generating random numbers from    * ')
%disp('           *          an arbitrary 1-D distribution (this is done    * ')
%disp('           *          through 2-D MC sampling, therefore a range     * ')
%disp('           *          in the random variable must be specified)      * ')
%disp('           *---------------------------------------------------------* ')
%disp('                                                                       ')

% -> X is a vector that samples the -X- variable. It is recomended to use at least 100 points/bins.
% -> Y is a vector that samples the -Y- variable.
% -> Nentries is the size of X_out.
% -> X_out is an array containing as many random numbers from the given
% distribution as indicated by Nentries.

function [X_out] = random1D(X,Y,Nentries,a_title)

if nargin==2, Nentries=1; end

Xmin  = min(X);
Xmax  = max(X);

X_adim = (X-Xmin)/(Xmax-Xmin);              % Make the X coordinate adimensional between 0 and 1

Y_adim = Y/max(Y);                          % Set  the Y coordinate to be between 0 and 1

i=0;

while(i<Nentries)

    X_rnd = rand;
    Y_rnd = rand;
      
    if Y_rnd<=interp1(X_adim, Y_adim, X_rnd)  
        i=i+1;
        X_out(i) = X_rnd*(Xmax-Xmin) + Xmin;       
    end
       
end

if nargin==2, return; end

%Comparison between the parent distribution and the performance of the
%random generator. 

if nargin==4
figure;

[N_,X_]=histf(X_out,X);
stairs(X_,N_/max(N_)); hold on;
plot(X, Y_adim);

title(a_title);
xlabel('random variable');
end

return