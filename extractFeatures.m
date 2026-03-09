function par = extractFeatures(spkFolder,microLabels,par,mode)
    
    % start parallel pooling
%     delete(gcp('nocreate'));
%     parpool(par.parpoolworkers); % the number of cores we use
    if isempty(gcp('nocreate'))
        parpool(par.parpoolworkers);
    end
    
    
    %% configuration for features extraction
    numChannels = length(microLabels);
    criteria = par.feature_select_criteria;
    good_features = cell(numChannels,1);

    fprintf("\nSelected settings: \n(Feature type = Wavelet, FamilyType: " + string(par.wave_type) + ...
        ", DecompositionLevel= " + num2str(par.scales) + ", NumberOfFeatures= " + num2str(par.inputs) + ")\n");
   
    %% looping through the channels
    parfor ichannel=1:numChannels 
        % getting file name for this channel
        filename = [ spkFolder filesep microLabels{ichannel} '_spikes.mat'];
        fprintf("Extracting features of channel " +string(microLabels{ichannel})+ "\n");

        % check for existance of file name
        if isfile(filename)
            % loading apikes for this channel
            spikeFile = load(filename);
            spikes = spikeFile.spikes;
            index = spikeFile.index;

            % check if any spikes are available
            if length(spikeFile.index)==0
                continue
            end

            % select method for feature selection
            %% method 1: based on custom feature prioritizing
            if (criteria=="custom") 
%                 fprintf("Using custom method \n");

                % assigning parameters for feature extraction
                coefScores= [];
                waveType= par.wave_type; %best sym3 - bior3.3 interesting
                waveLevel = par.scales;
                numFeatures = par.inputs;
                cAll=[];
                
                % loop over all spikes to find the most relevant features
                for i=1:size(spikes,1)

                    % load individual spike
                    spike= spikes(i,:);

                    % wavelet decomposition
                    [c,l] = wavedec(spike,waveLevel,waveType);
                    if isempty(cAll)
                        cAll=c;
                    else
                        cAll=[cAll;c];
                    end
                    
                    if mode == "FitAndApply" 
                        % separating certain coefficients for prioritization 
                        firstcoefs=[];
                        secondcoefs=[];
                        currentpoint=0;
                        secondcoefHigh=[];
                        for i=1:length(l)-1
                            firsthalf = round(l(i)/2);
                            secondhalf = l(i)-firsthalf;
                            firstcoefs=[firstcoefs, currentpoint+1:currentpoint+firsthalf];
                            currentpoint = currentpoint+firsthalf;
                            secondcoefs= [secondcoefs, currentpoint+1:currentpoint+secondhalf];
                            if i>length(l)/2
                                secondcoefHigh=[secondcoefHigh, currentpoint+1:currentpoint+secondhalf];
                            end
                            currentpoint = currentpoint+secondhalf;
                        
                        end
                        
                        % selecting coefficients that are most contributing
                        % based priority 
                        % first half
                        firstHalf=c(firstcoefs);
                        m=mean(abs(firstHalf));
                        s=std(abs(firstHalf));
                        firstHalf(abs(firstHalf)<m)=0;
                        c(firstcoefs) =firstHalf;
                        
                        % second half
                        secondHalf=c(secondcoefs);
                        m=mean(abs(secondHalf));
                        s=std(abs(secondHalf));
                        secondHalf(abs(secondHalf)<m+s)=0;
                        c(secondcoefs) =secondHalf;
                        c(secondcoefHigh) =0;
                        
                        if isempty(coefScores)
                            coefScores=[c~=0];
                        else
                            coefScores= coefScores+[c~=0];
                        end
                    end
                end
                if mode == "FitAndApply"
                    % sort most contributing coeffs
                    [b goodCoefsIndeces]=sort(coefScores,'descend');
                    % choose final coeff to be extracted
                    goodCoefsIndeces = goodCoefsIndeces(1:numFeatures);
                elseif mode == "OnlyApply"
                    goodCoefsIndeces = par.good_features_all{ichannel};
                end

                good_features{ichannel} = goodCoefsIndeces;
                
                % assigning featrues for all spikes
                features = zeros(size(spikes,1),numFeatures);
                for i=1:size(spikes,1)
                    spike= spikes(i,:);
                    c = cAll(i,:)
                    features(i,:) = c(goodCoefsIndeces);
                end

                % saving features
                parsave([spkFolder filesep microLabels{ichannel} '_features.mat'],spikes ,features,index);

            %% method 2: based on distribution of feature values (take most un-normal ones)
            elseif(criteria ==  "normalityTest") 
%                 fprintf("Using normality method \n");

                % assigning parameters for feature extraction
                scales = par.scales;
                numFeatures = par.inputs;
                waveType=par.wave_type;
                nspk = size(spikes,1);
                ls = size(spikes,2);
                cc=zeros(nspk,ls);
                % check for wavelet package
                if ~~exist('wavedec','file')   
                    for i=1:nspk                
                        % wavelet decomposition
                        [c,l]=wavedec(spikes(i,:),scales,waveType);
                        cc(i,1:ls)=c(1:ls);
                    end
                else
                    fprintf("Install wavelets toolbox!! \n");
                end
                sd=[];
                
                if mode == "FitAndApply"
                    % KS test for coefficient selection  
                    for i=1:ls                                   
                        thr_dist = std(cc(:,i)) * 3;
                        thr_dist_min = mean(cc(:,i)) - thr_dist;
                        thr_dist_max = mean(cc(:,i)) + thr_dist;
                        
                        aux = cc(find(cc(:,i)>thr_dist_min & cc(:,i)<thr_dist_max),i);
                        
                        if length(aux) > 10;
                            [ksstat]=test_ks(aux);
                            sd(i)=ksstat;
                        else
                            sd(i)=0;
                        end
                    end
                    
                    % choose final coeff to be extracted
                    coeff=[];
                    [max_sd ind]=sort(sd);
                    coeff(1:numFeatures)=ind(ls:-1:ls-numFeatures+1);
                    good_features{ichannel} = coeff(1:numFeatures);
                elseif mode == "OnlyApply"
                    good_features{ichannel} = par.good_features_all{ichannel};
                end

                % assigning featrues for all spikes
                features = zeros(size(spikes,1),numFeatures);
                for i=1:nspk 
                    features(i,:) = cc(i,good_features{ichannel})
                end

                % saving features
                parsave([spkFolder '\' microLabels{ichannel} '_features.mat'], spikes, features, index);
            end
        else
            fprintf("file does not exist\n");
        end
    end
    par.good_features_all = good_features;
end



%% %%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving function v7.3 works with parallel loops
function parsave(fname, spikes, features, index)
    save(fname, 'spikes', 'features','index','-v7.3')
end

% Normality test function
function [KSmax] = test_ks(x)
    % 
    % Calculates the CDF (expcdf)
    %[y_expcdf,x_expcdf]=cdfcalc(x);
    
    yCDF = [];
    xCDF = [];
    x = x(~isnan(x));
    n = length(x);
    x = sort(x(:));
    % Get cumulative sums
    yCDF = (1:n)' / n;
    % Remove duplicates; only need final one with total count
    notdup = ([diff(x(:)); 1] > 0);
    x_expcdf = x(notdup);
    y_expcdf = [0; yCDF(notdup)];
    
    %
    % The theoretical CDF (theocdf) is assumed to be normal  
    % with unknown mean and sigma
    
    zScores  =  (x_expcdf - mean(x))./std(x);
    
    %theocdf  =  normcdf(zScores , 0 , 1);
    mu = 0; 
    sigma = 1; 
    theocdf = 0.5 * erfc(-(zScores-mu)./(sqrt(2)*sigma));
    
    
    %
    % Compute the Maximum distance: max|S(x) - theocdf(x)|.
    %
    
    delta1    =  y_expcdf(1:end-1) - theocdf;   % Vertical difference at jumps approaching from the LEFT.
    delta2    =  y_expcdf(2:end)   - theocdf;   % Vertical difference at jumps approaching from the RIGHT.
    deltacdf  =  abs([delta1 ; delta2]);
    
    KSmax =  max(deltacdf);
end



