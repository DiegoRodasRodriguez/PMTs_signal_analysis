%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%               Code for aligning waveforms                       %%%%%
%%%%%                    (DGD 14/Mar/2021)                            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wvf_f = TimeAlign(time, wvf_i, ti, tf)

if(~isempty(find(time==ti)) && ~isempty(find(time==tf))) %#ok<EFIND>        %If times coincide with time bins: SHIFT WITHOUT MANIPULATION
    Ii = find(time == ti);
    If = find(time == tf);
    
    if(tf<ti)                                                               %Left-shift
        Ishift = Ii - If;
        Imax  = length(time);
        wvf_f(1 : Imax-Ishift) = wvf_i(Ishift+1 : Imax);
        wvf_f(Imax)      = 0;
    elseif(tf>ti)                                                           %Right-shift
        Ishift = If - Ii;
        Imax = length(time);
        wvf_f(1:Ishift)    = 0;
        wvf_f(1+Ishift : Imax) = wvf_i(1 : Imax-Ishift);
        wvf_f(Imax)          = 0;
    else
        wvf_f = wvf_i;
    end
    
else                                                                        %If times coincide with time bins: SHIFT WITH LINEAR INTERPOLATION
    timeShift = time - (ti - tf);
    wvf_f = interp1(timeShift, wvf_i, time, 'linear', 0);
end

end
