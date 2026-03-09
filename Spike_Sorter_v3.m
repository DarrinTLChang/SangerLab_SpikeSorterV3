function par = Spike_Sorter_v3(directoryRead,directorySave,microLabels,par,mode,dataID,period,considered_window)


%% Creating destination folders if not already existing and navigating into them
spikeFolderName = ['Spikes_' num2str(par.stdmin) 'std'];
clusterFolderName = ['SpikeClusters_' num2str(par.stdmin) 'std_' par.features];
spikeTime = struct; spikeWaveform = struct; spikesNotClassified = struct;

spkFolder = [directorySave filesep spikeFolderName];
if ~exist(spkFolder)
    mkdir(spkFolder)
end
clusterFolder = [directorySave filesep clusterFolderName];
if ~exist(clusterFolder)
    mkdir(clusterFolder)
end
diary([directorySave '\log.txt'])
diary on


%% Detect spikes
tic
fprintf('\nDetecting spikes from all micro contacts in %s period %i\n',dataID,period) 
[~,par] = detectSpikes(directoryRead,microLabels,par,spkFolder,mode,considered_window); %% options = {"FitAndApply", "OnlyApply"}
fprintf('\nDone with spike detection\n')
toc
%% Extract features 
tic
fprintf('\nExtracting features for %s period %i\n',dataID,period)
par = extractFeatures(spkFolder,microLabels,par,mode);
fprintf("Done with extracting features \n ")
toc
%% Generate clusters 
tic
fprintf('\nBuilding clusters for %s period %i\n',dataID,period)
par = clusterSpikes(spkFolder,microLabels, par,mode);
fprintf("Done with clustering \n ")
toc

%% Create data structures
tic
% Calculates fs and segmentlength
fprintf('\nGenerating spike structures: %s\n',dataID)
sampleFile = dir([directoryRead filesep microLabels{1} '*.mat']);
load([sampleFile.folder filesep sampleFile.name]) %loads fs and data
lengthData = size(data,2)/fs; clear data;

for i = 1:numel(microLabels)
    % Initialize structures
    spikeTime(i).electrode = microLabels{i};      spikeWaveform(i).electrode = microLabels{i};   
    spikeTime(i).fs = fs;                         spikeWaveform(i).fs = fs;                        
    spikeTime(i).time = [];                       spikeWaveform(i).time = [];                       
                                                  spikeWaveform(i).waveform = [];                 
    spikeTime(i).dataSegmentLength=lengthData;    spikeWaveform(i).dataSegmentLength=lengthData;    
    fname = [spkFolder filesep microLabels{i} '_clusters.mat'];
    if exist(fname,'file') == 2 % If there were enough spikes to generate 'times_' file:
        load(fname);
%         numNeurons = unique(clusterID);
%         spikeTime(i).numberNeurons = numel(numNeurons); spikeWaveform(i).numberNeurons = numel(numNeurons);
         numNeurons = par.gm_All{i}.NumComponents;
         spikeTime(i).numberNeurons = numNeurons; spikeWaveform(i).numberNeurons = numNeurons;
        if numNeurons~=0
            for cn = 1:numNeurons
                counts = length(clusterID(clusterID==cn));
                times = index(clusterID==cn);
                spikeTime(i).time{cn} = times;
                spikeTime(i).counts(cn) = counts;
                spikeWaveform(i).time{cn} = times;
                spikeWaveform(i).counts(cn) = counts;
                spikeWaveform(i).waveform{cn} = spikes(clusterID==cn,:);
                % SNR
                thisWaveform = mean(spikeWaveform(i).waveform{cn});
                noise = spikeWaveform(i).waveform{cn}-repmat(thisWaveform,size(spikeWaveform(i).waveform{cn},1),1);    
                powerSignal = sum(thisWaveform.^2)/length(thisWaveform);
                powerNoise = sum(noise(:).^2)/length(noise(:));           
                snr = (powerSignal./powerNoise);  
                spikeTime(i).snr(cn) = snr;
                spikeWaveform(i).snr(cn) = snr; 
            end
        else
           spikeTime(i).counts = 0; 
           spikeWaveform(i).counts = 0; 
        end        
    else
        spikeTime(i).numberNeurons = 0; spikeWaveform(i).numberNeurons = 0;
        spikeTime(i).counts = 0; spikeWaveform(i).counts = 0;
    end
end

fprintf("Done with generating spike structures \n ")
toc
%% Save data
fprintf('Saving: %s\n',dataID)  
save([clusterFolder filesep 'sortingParameters.mat'],'par')
save([clusterFolder filesep 'spikeTime.mat'],'spikeTime')
save([clusterFolder filesep 'spikeWaveform.mat'],'spikeWaveform','-v7.3')
fprintf('\nDone.\n\n\n')
diary off


end


