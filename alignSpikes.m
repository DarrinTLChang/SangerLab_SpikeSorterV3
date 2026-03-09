
function [spikes1,index] = alignSpikes(data,index,par,biasIndex)
    if ~isempty(index)
        %% Get spikes
        samplesBefore = par.w_pre; % how much time to take into the window before the peak (0.5)
        samplesAfter = par.w_post; % how much time to take into the window after the peak (1.5)
        index(index<=samplesBefore+par.int_factor) = [];
        index(index>=length(data)-samplesAfter-par.int_factor) = [];
        spikes = zeros(length(index),samplesBefore+samplesAfter+par.int_factor*2+1); %? why add int
        for i=1:length(index)                          
             spikes(i,:) = detrend(data(index(i)-samplesBefore-par.int_factor:index(i)+samplesAfter+par.int_factor)); % detrend to remove the average
        end       
        %% Interpolate  spikes
        ls = samplesBefore + samplesAfter;
        nspk = size(spikes,1);
        extra = (size(spikes,2)-ls)/2;
        s = 1:size(spikes,2);
        ints = 1/par.int_factor:1/par.int_factor:size(spikes,2);
        intspikes = zeros(1,length(ints));
        spikes1 = zeros(nspk,ls);
        for i=1:nspk
            intspikes(:) = spline(s,spikes(i,:),ints);
            [~,iaux] = max(abs(intspikes((samplesBefore+extra-1)*par.int_factor:(samplesBefore+extra+1)*par.int_factor))); 
%             [~,iaux] = min((intspikes((samplesBefore+extra-1)*par.int_factor:(samplesBefore+extra+1)*par.int_factor))); 
            iaux = iaux + (samplesBefore+extra-1)*par.int_factor -1;
            spikes1(i,1:ls)= intspikes(iaux-samplesBefore*par.int_factor+par.int_factor:par.int_factor:iaux+samplesAfter*par.int_factor);
        end
    else
        spikes1 = [];
    end
end


