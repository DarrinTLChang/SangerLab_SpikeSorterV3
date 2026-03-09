clear all
close all
clc
gcp('nocreate');
%% Instructions to run the code
% 1- insert 3 paths mentioned in this code corresponding the "header path",
% "Reading directory" which contains .mat channel recordings and "Saving
% directory" to save all results
% 2- Choose how many parallel cores you want to use in set_parameters
% script
% 3- choose the periods you want to run the code for
% 4- if needed change parameters for different stages of the pipeline in
% set_parameters
% Run!

%% User Inputs
dataID = 'test_subject'; % String: s + 3 digits patient ID + _ + YYYYMM (recording date)
nPeriod = [1];%5 8 9 10 11 12 13]; % Vector with the period number(s) to be clustered
nmode = ["OnlyApply"]; %% options for each period = {"FitAndApply", "OnlyApply"}
% par = set_parameter(); % All modifications of parameters should be done in this function and not outside!
par_path = 'E:\test_subject\Mat Data\test\sorting_v2_autcluster_3std_4_fixed\period1\SpikeClusters_3std_wav\sortingParameters.mat';
 
fit_period.start = 1; % in index
fit_period.end = 16190000; % in index 16190000

%% Running header files to get the names of electrodes ** Insert header file path here**
% for testing : 'E:\test_subject\test_subject_header_bipolar'
header = ['E:\test_subject\test_subject_header_bipolar']; run(header) % subject specific header file (contact names)
fprintf('\nStarting spike sorting for %s\n',dataID)

for iperiod = 1:length(nPeriod)
    period = nPeriod(iperiod);
    directoryRead = ['E:\test_subject\Mat Data\test\micro_bipolar\period' num2str(period)];
    directorySave = ['E:\test_subject\Mat Data\test\sorting_v3_autcluster_3std_12_fixed\period' num2str(period)];

    % check if spike sorter should fit or only predict. 
    % if needs to be fitted, fit on the specified period
    % 
    mode = nmode(iperiod);
    if mode== "FitAndApply" % fit with 
        fprintf("\n fitting parameters ... \n");
        par = set_parameter();
        par = Spike_Sorter_v3(directoryRead,directorySave,microLabels,par,"FitAndApply",dataID,period,fit_period);

    elseif mode== "OnlyApply"
        load(par_path); % loads par
    end
    
        
    % Apply spike sorter using the pre-fit parameters
    % check if file needs to get sliced 
    queryfile=dir([directoryRead filesep microLabels{1} '*.mat']);
    % find the maximum array size to see if split is needed
    info = whos('-file',[queryfile.folder filesep queryfile.name]);
    vars = numel(info);
    for ivar=1:vars
        m(ivar)= max(info(ivar).size);
    end
    max_array_size = max(m);
    

    if max_array_size > 10000000 % 50000000 is around 30 mins
        max_acc_arr_size=2000000; % best 10000000 with 1 core (test 1000000)
        loops = floor(max_array_size/max_acc_arr_size)-1;
        fprintf("\n file too big ... deviding the data into " + num2str(loops) + "... \n");
        for index = 0:loops-1
            fprintf("   processing part " + num2str(index) + "... \n");
            % index*max_acc_arr_size+1:(index+1)*max_acc_arr_size
            considered_period.start = index*max_acc_arr_size+1;
            considered_period.end = (index+1)*max_acc_arr_size;
            directorySave_part = [directorySave 'Part' num2str(index)];
            Spike_Sorter_v3(directoryRead,directorySave_part,microLabels,par,"OnlyApply",dataID,period,considered_period);
        end
        % loops*max_acc_arr_size+1:max_array_size
        fprintf("   processing part " + num2str(loops) + "... \n");
        directorySave_part = [directorySave 'Part' num2str(loops)];
        considered_period.start = loops*max_acc_arr_size+1;
        considered_period.end = max_array_size;
        Spike_Sorter_v3(directoryRead,directorySave_part,microLabels,par,"OnlyApply",dataID,period,considered_period);

        %% stitching the spike times
        fprintf("\n Stitching the spike times for " + num2str(loops) + " files ... \n");
        clusterFolderName = ['SpikeClusters_' num2str(par.stdmin) 'std_' par.features];
        directorySave_part = [directorySave 'Part' num2str(0)  filesep clusterFolderName filesep 'spikeTime.mat'];
        load(directorySave_part)
        for index = 1:loops
            directorySave_part = [directorySave 'Part' num2str(index)  filesep clusterFolderName filesep 'spikeTime.mat'];
            part=load(directorySave_part);
            for ichannel = 1:numel(microLabels)
                for ineuron=1:size(spikeTime(ichannel).time,2)
                    spikeTime(ichannel).counts(ineuron)=spikeTime(ichannel).counts(ineuron) + part.spikeTime(ichannel).counts(ineuron);
                    spikeTime(ichannel).time{1,ineuron}(end+1: end+part.spikeTime(ichannel).counts(ineuron)) = part.spikeTime(ichannel).time{1,ineuron};
                    
                end
                spikeTime(ichannel).snr=[spikeTime(ichannel).snr;part.spikeTime(ichannel).snr];
            end
        end
        FinalFolder = [directorySave filesep];
        if ~exist(FinalFolder)
            mkdir(FinalFolder)
        end
        save([FinalFolder 'spikeTime.mat'],'spikeTime','-v7.3')
        save([FinalFolder 'parameters.mat'],'par','-v7.3')
        fprintf("\n Done stitching! \n");

    else
        considered_period.start = 1;
        considered_period.end = max_array_size;
        Spike_Sorter_v3(directoryRead,directorySave,microLabels,par,"OnlyApply",dataID,period,considered_period);
    end
end

