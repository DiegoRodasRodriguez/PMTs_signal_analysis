%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%            Code for reading data from CAEN card                 %%%%%
%%%%%                     (DGD 10/Dec/2020)                           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = [DIR, FILE, '.mat'];
fid      = fopen(filename,'r');

if (fid ~= -1)
    eval(['load ', filename]);    
    disp(' loading workspace...');
    fclose(fid);
else    
    disp(' ');
    disp(['NOTE!: CAEN ''', filename, '']);
    disp(' workspace not found, reading from ASCII file and creating workspace...');
    disp(' ');
    
    for i=1:nCh                                                 
        filename = [DIR, FILE, num2str(i-1), '.txt'];
        fid      = fopen(filename,'r');
        
        Nevt = 1;
        
        if (fid ~= -1)
            while(~feof(fid))
                for k=1:7, fgets(fid); end
                data(i, :, Nevt) = fscanf(fid,'%d');
                if(rem(Nevt, 100)==0 && i<=4)
                    disp(['PMT ', num2str(i), ' // evt ', num2str(Nevt)]);
                elseif(rem(Nevt, 100)==0 && i==5)
                    disp(['MWPC ', num2str(i), ' // evt ', num2str(Nevt)]);
                end
                if (Nevt==Nsave), break; end               
                Nevt = Nevt + 1;
            end
        end
        if(fid ~= -1), fclose(fid); end %Do not really understand why is needed...
    end
    if (Nevt==Nsave && i==1)
        disp(['WARNING: Workspace limit exceeded, number of events trimmed down to ', num2str(Nsave)]);
    end
    disp('saving workspace...');
    eval(['save ', DIR, FILE, '.mat data']);
end

dataSize  = size(data);
ABin = 2/(2^14-1) * 1000;
data = - data * ABin;                    %[mV] (convert to mV and flip sign)
TBin = 4;
time = (0:(length(data(1,:,1))-1))*TBin; %[ns] (sampling time of the CAEN card)
Nevts = dataSize(3);


