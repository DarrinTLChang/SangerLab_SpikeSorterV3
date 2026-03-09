

function [index,thresh_spk,thersh_artifact] = getSpikeIndex(data,fs,par,mode,ichannel) % Detect spike index by adaptative threshold

    % Check data size. Data must be a row vector
    [r,c] = size(data);
    if r>1 && r*c == length(data)
        data = data';
    end
    % Tresholding on the peaks
    [peakValues, peakTimes, peakWidths, peakProminence] = findpeaks(data);  %normalize the heights for plotting

    % Fit a gaussian distribution on the peaks
    if mode == "FitAndApply" 
        pd = fitdist(peakValues', 'Normal'); 
        thresh_spk = pd.mu + par.stdmin*pd.sigma;
        thersh_artifact = pd.mu + par.stdmax*pd.sigma;
    elseif mode =="OnlyApply" 
        thresh_spk = par.thresholdAll{ichannel};
        thersh_artifact = par.thersoldArtifacthAll{ichannel};
    end

    iPeak = peakTimes(find(abs(peakValues) > thresh_spk)); 
    
    
    % Eliminate indices larger than max threshold and assign the output 
    if ~isempty(iPeak)
            %% Eliminate indices larger than maximum threshold
            for i=1:length(iPeak)                          
                if abs(data(iPeak(i))) >= thersh_artifact
                    iPeak(i) = 0;
                end
            end
            iPeak(~iPeak) = [];
            index = iPeak;        
    else
        index = [];
    end
end
