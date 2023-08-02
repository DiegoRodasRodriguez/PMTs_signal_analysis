%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        Code for reading data from Tektronix oscilloscope        %%%%%
%%%%%                    (DGD 10/Dec/2020)                            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chain the Nfiles associated to the calibration and squeeze the data matrix:
% data[iPM, jTime, kEvts]   : iPM == PM index  //  jTime == time bin index;  //  kEvts == event ID;

for i=1:Nfiles
    load([DIR,FILE,'-', num2str(i-1)]);
    eval(['data', num2str(i-1), '=(squeeze(data(1,:,:)));']);
end
data = data0;
for i=2:Nfiles, eval(['data = [data, data', num2str(i-1),'];']); end 

% Retrieve timebin per point in x-axis and amplitude factor for V-scaling
eval(['load(''',DIR,FILE, '.hdr.mat'')']);
TBin      = Scales(1).xincr / 1e-9;  %ns
Amp2mV    = 1000;                    %scope indeed seems to always save in 'V'
% Obtain amplitude bin
dataSize  = size(data);
dataArray = reshape(data, 1, numel(data));
DiffArray = abs(diff(dataArray));
ABin      = min(DiffArray(DiffArray>0));
% Define amplitude as positive and convert to mV
data      = -data*Amp2mV;     % [V] to [mV] (and sign inverted)
ABin      = ABin*Amp2mV;      % will be useful later for smearing
% Define time range and convert to ns
time = 0:(dataSize(1))-1;
time = time*TBin;
% Number of events
Nevts = dataSize(2);
