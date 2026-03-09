function par = set_parameter() 
    fprintf("\nSetting parameters!! \n");

    % Load params
%     par.segments_length = 60;             % length (in minutes) of segments in which the data is cutted (default 5min).   % diff 5min
    par.sr = 24414;                      % sampling rate (in Hz). This parameter will be only used if the data file don't have a sr. 

    % Parpool parameters
    par.parpoolworkers = 4;
    
    % Event detection parameters
    par.neo_order = 1;                   % NEO filter order 

    par.stdmin = 3;                      % minimum threshold for detection
    par.stdmax = 70;                     % maximum threshold for detection (bigger than this is disregarded)

    par.detect_fmin = 350;               % high pass filter for detection 
    par.detect_fmax = 3000;              % low pass filter for detection 
    par.detect_order = 8;                % filter order for detection 

    par.sort_fmin = 3e2;                 % high pass filter for sorting was 300
    par.sort_fmax = 3500;                % low pass filter for sorting (default 3000) was 3.5k
    par.sort_order = 8;                  % filter order for sorting 

    par.ms_repeated_index = 0.1;         % ms window to ignore repeated indices when multi electrode array is used 
    par.max_num_coincedences = 5;        % maximum number of coincedences very close to each other (coincedences more than this will be interpreted as artifact)
                                         % if want to skip this part put max_num_coincedences to veryhugh number 
    % Feature extraction parameters
    par.inputs = 15;                     % number of inputs to the clustering 
    par.scales = 6;                      % level of decomposition for wavelet
    par.wave_type = ['sym3'];            % wavelet family type could be ('haar', 'sym3' , etc.)
    par.feature_select_criteria = "custom"; % this is the criteria with which few features are selected from many (could be "custom" or "normalityTest")
    par.features = 'wav';                % type of feature not an option now!!!! only 'wav'


    % Interpolation parameters for peak alignment 
    par.int_factor = 6;                  % interpolation factor 
    par.w_pre = round(0.5*par.sr/1000);  % number of pre-event data points stored  
    par.w_post = round(1.5*par.sr/1000); % number of post-event data points stored 

    % Clustering parameters
    par.initNumCluster=15;                   % number of clusters for initial round of clustering for cleaning up spikes (remove spike too far from clusters) 
    par.minNumCluster=5;                     % min number of cluster when searching for appropriate number of clusters (mode auto)
    par.maxNumCluster=15;                    % max number of cluster when searching for appropriate number of clusters (mode auto)
    par.fixedNumCluster=5;                   % the number of cluster when fixed number of clusters is used (mode fixed)
    par.numClusterMode = "fixed";                  % could be ("auto" or "fixed")
    par.clusterMethod = "GMM";               % Currerntly only "GMM"
    par.minClsuterSize =50;                  % min number of spikes allowed in a cluster

    par.initProbThreshold = 0.05;            % (percentage of mean pdf) spikes with probability value lower than this based on GM will be removed as bad spikes in the cleaning stage
    par.finalProbThreshold = 0.01;           % (percentage of mean pdf) spikes with probability value lower than this based on GM will be removed as bad spikes in the final stage



end
