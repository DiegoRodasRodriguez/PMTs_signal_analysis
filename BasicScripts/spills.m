
disp('==================')

counts=diff(InewSpill);

if (size(counts)>=1)

    if      (RPC==1) 
        rate=0.000004*counts.*counts + 0.2057*counts;
    
    else if (RPC==3)  
       rate=0.00003*counts.*counts + 0.2207*counts; 

    end
end 
 

plot(rate,'.','MarkerSize',15);hold on;

 
plot(rate'*0 + mean(rate));
plot(rate'*0 + mean(rate) + 0.1*mean(rate),'r');
plot(rate'*0 + mean(rate) - 0.1*mean(rate),'r');
ylabel('rate (Hz/cm2)');
xlabel('#spill');


devrate=std(counts);

%plot(mean(round(diff(InewSpill)))

end

end