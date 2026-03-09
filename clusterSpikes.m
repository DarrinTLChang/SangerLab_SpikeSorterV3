function par = clusterSpikes(spkFolder,microLabels,par,mode)
% thngs to do :
%some channels find one spike group because of awful noises:
% maybe redo cluster again after removing bad spikes and bad clusters
% maybe with eval. ororrrrr use false in shared weights which makes it able
% to make the noise one bigger and put them into one cluster 
% important: minClsuterSize ahould be  FUNCTION OF TIME OR TOTAL SPIKES

    % parameter initiation
    numChannels = length(microLabels);
    gm_All = cell(numChannels,1);
    C_norm_All = cell(numChannels,1);
    S_norm_All = cell(numChannels,1);
    pca_coeff_ALL = cell(numChannels,1);

    % start parallel pooling
%     delete(gcp('nocreate'));
%     parpool(par.parpoolworkers);
    if isempty(gcp('nocreate'))
        parpool(par.parpoolworkers);
    end

    % configuration for features extraction
    initNumCluster = par.initNumCluster; 
    minNumCluster = par.minNumCluster; 
    maxNumCluster = par.maxNumCluster; 
    clusterMethod = par.clusterMethod; 
    minClsuterSize = par.minClsuterSize; 
    initProbThreshold = par.initProbThreshold; 
    finalProbThreshold = par.finalProbThreshold; 
    numClusterMode = par.numClusterMode; 
    fixedNumCluster = par.fixedNumCluster;  

    % stats of spike removal
    allSpikes=0;
    goodSpikes=0;

    % looping through the channels
%     for ichannel=2:2  %%% UNCOMMENT FOR TESTING ... %%%
    parfor ichannel=1:numChannels
        % Warning off
        warning('off');

        % creating reading path
        filename = [ spkFolder filesep microLabels{ichannel} '_features.mat'];
        fprintf("Clustering channel " +string(microLabels{ichannel})+ "\n");

        
        % check if there are enough events for clustering then read them
        if isfile(filename)
            featureFile = load(filename);
            features = featureFile.features;
            spikes = featureFile.spikes;
            index = featureFile.index;
            if length(index)<minClsuterSize
                continue
            end
            
            % add to the total spikes numbers
            allSpikes = allSpikes+length(index);
            
            if mode=="FitAndApply" 
                %% Initial clustering for cleaning up spikes
                % 1. cluster with the selected method
                switch clusterMethod
                    % method 1 gaussian mixture method
                    case "GMM"
                        % need to use orthogonal features so use pca
                        [coeff,orthfeatures,latent] = pca(normalize(features)); 
                        options = statset('MaxIter', 5000);
                        gm=[];
                        try
                            gm = fitgmdist(orthfeatures(:,1:5),initNumCluster, 'CovarianceType', 'full', 'SharedCovariance', true,  ...
                        'Options', options,'replicates', 10); %full better / shared=true better false not sure
                            features2 = orthfeatures(:,1:5);
                            % shared does not make much sense unless we
                            % take the bad spikes out and rerun spikes again
                            % false gives more clusters while true less
                            % clusters
                        catch 
                            fprintf(" could not fit gm!\n")
                            gm = fitgmdist(orthfeatures(:,1:5),initNumCluster, 'CovarianceType', 'full', 'SharedCovariance', true,  ...
                        'Options', options,'replicates', 10); %full better / shared=true better false not sure
                            features2 = orthfeatures(:,1:5); 
                        end
    
                        if isempty(gm)
                            fprintf("gm did not work! skipping channel!\n")
                            continue
                        end
                    
                        
                        clusterID = cluster(gm,features2);
    
    
    
                    % method 2 k means
                    case "Kmeans"
                        [clusterID, clusterCenter] = kmeans(features,initNumCluster);
    
                    % method 3 lin kage
                    case "Linkage"
                        Z = linkage(features,'ward');
                        clusterID = cluster(Z,initNumCluster);
    
                end
    %             plotClusteredSpikes(clusterID,spikes);  %%% UNCOMMENT FOR TESTING ... %%%
    
                % Calculate clustering scores
    %                 silhouetteScore = [silhouetteScore, evalclusters(features,clusterID,'silhouette')];
    %                 CalinskiHarabaszScoer = [CalinskiHarabaszScore, evalclusters(features,clusterID,'CalinskiHarabasz')];
    %                 DaviesBouldinScore = [DaviesBouldinScore, evalclusters(features,clusterID,'DaviesBouldin')];
    
    
    
                % 2. Remove spikes based on their probablities
                meanpdf = mean(pdf(gm, features2));
                badindices= pdf(gm, features2)<(initProbThreshold*meanpdf);
                spikes(badindices,:)=[];
                index(badindices)=[];
                features(badindices,:)=[];
                features2(badindices,:)=[];
                clusterID(badindices)=[];
    %             plotClusteredSpikes(clusterID,spikes);  %%% UNCOMMENT FOR TESTING ... %%%
    
                % 3. Remove clusters with low number of spikes
                clusterSizes =zeros(1,max(clusterID));
                newIDCounter=1;
                for i =1:max(clusterID)
                    if length(clusterID(clusterID==i))>minClsuterSize
                        clusterID(clusterID==i) = newIDCounter;
                        newIDCounter = newIDCounter+1;
                    else
                        spikes(clusterID==i,:) = [];
                        index(clusterID==i)=[];
                        features(clusterID==i,:)=[];
                        features2(clusterID==i,:)=[];
                        clusterID(clusterID==i) = [];         
                    end
                end
    %             plotClusteredSpikes(clusterID,spikes);  %%% UNCOMMENT FOR TESTING ... %%%
                if length(index)<minClsuterSize
                    continue
                end
            end
            %% Final clustering 
            if mode=="FitAndApply" 
                % figuring out the correct number of clusters
                [normalized_features,C,S] = normalize(features);
                [coeff,orthfeatures,latent] = pca(normalize(normalized_features));
                C_norm_All{ichannel} = C;
                S_norm_All{ichannel} = S;
                pca_coeff_ALL{ichannel} = coeff;
                features2 = orthfeatures(:,1:5);
                if numClusterMode=="auto"
                    if(size(features2,1)>5000)
                        features3=features2(1:5000,:);
                    else
                        features3=features2;
                    end
                    eval = evalclusters(features3,'linkage','DaviesBouldin','KList',minNumCluster:maxNumCluster);
                    scoreList = eval.CriterionValues;
                    scoreList = [max(scoreList),scoreList, max(scoreList)];
                    InspectedK = eval.InspectedK;
                    InspectedK = [0,InspectedK, 0];
                    [pn ln] = findpeaks(-scoreList);
                    optimalKs = InspectedK(ln);
                    rng(0, 'twister');
                    eval = evalclusters(features3,'gmdistribution','silhouette','KList',optimalKs);
                    OptimalK = eval.OptimalK;
                else 
                    OptimalK = fixedNumCluster;
                end
    
                % 1. cluster with the selected method on the cleaned dataset
                gm=[];
                rng(0, 'twister');
                try
                    gm = fitgmdist(orthfeatures(:,1:5),OptimalK, 'CovarianceType', 'full', 'SharedCovariance', false,  ...
                            'Options', options,'replicates', 10); %full better / shared=true better false not sure
                catch
                    fprintf(" could not fit gm2!\n")
                    gm = fitgmdist(orthfeatures(:,1:5),OptimalK, 'CovarianceType', 'full', 'SharedCovariance', true,  ...
                            'Options', options,'replicates', 10); %full better / shared=true better false not sure
                end
                if isempty(gm)
                    fprintf(" skipping channel!\n")
                    continue
                end

            elseif mode == "OnlyApply"
                C_norm_All{ichannel} = par.C_norm_All{ichannel};
                S_norm_All{ichannel} = par.S_norm_All{ichannel};
                pca_coeff_ALL{ichannel} = par.pca_coeff_ALL{ichannel};
                normalized_features = normalize(features,"center",C_norm_All{ichannel},"scale",S_norm_All{ichannel});
                orthfeatures =  normalized_features*pca_coeff_ALL{ichannel};
                gm = par.gm_All{ichannel};
            end
            gm_All{ichannel} = gm;
            clusterID = cluster(gm,orthfeatures(:,1:5));
            features2 = orthfeatures(:,1:5);

%             eval = evalclusters(features2,'gmdistribution','silhouette','KList',1:10); 
%             eval = evalclusters(features2,'gmdistribution','DaviesBouldin','KList',1:10);
%             eval = evalclusters(features2,'gmdistribution','CalinskiHarabasz','KList',1:10);
% % 
%             eval =[];
%             eval = evalclusters(features,'kmeans','silhouette','KList',1:10);
%             eval = evalclusters(features,'kmeans','DaviesBouldin','KList',1:10);
%             eval = evalclusters(features,'kmeans','CalinskiHarabasz','KList',1:10);
% 
%             eval = evalclusters(features,'linkage','silhouette','KList',1:10);
%             eval = evalclusters(features,'linkage','DaviesBouldin','KList',1:10);
%             eval = evalclusters(features,'linkage','CalinskiHarabasz','KList',1:10);
%             % show clusters
%             clusterID = eval.OptimalY;


%             plotClusteredSpikes(clusterID,spikes);  %%% UNCOMMENT FOR TESTING ... %%%

            % 2. Remove spikes based on their distance to the centroids
            meanpdf = mean(pdf(gm, features2));
            badindices= pdf(gm, features2)<(finalProbThreshold*meanpdf);
            spikes(badindices,:)=[];
            index(badindices)=[];
            features(badindices,:)=[];
            features2(badindices,:)=[];
            clusterID(badindices)=[];
%             plotClusteredSpikes(clusterID,spikes);  %%% UNCOMMENT FOR TESTING ... %%%

            % 3. Remove clusters with low number of spikes
            clusterSizes =zeros(1,max(clusterID));
%             newIDCounter=1;
            for i =1:max(clusterID)
                if length(clusterID(clusterID==i))<minClsuterSize
%                     clusterID(clusterID==i) = newIDCounter;
%                     newIDCounter = newIDCounter+1;
%                 else
                    spikes(clusterID==i,:) = [];
                    index(clusterID==i)=[];
                    features(clusterID==i,:)=[];
                    features2(clusterID==i,:)=[];
                    clusterID(clusterID==i) = [];         
                end
            end

%             plotClusteredSpikes(clusterID,spikes);  %%% UNCOMMENT FOR TESTING ... %%%

            goodSpikes = goodSpikes+length(index);

            % Save cluster information
            parsave([spkFolder filesep microLabels{ichannel} '_clusters.mat'], clusterID, spikes, features2, index);
    
        else
            fprintf("file does not exist\n");
        end

    end
    par.gm_All = gm_All;
    par.C_norm_All = C_norm_All;
    par.S_norm_All = S_norm_All;
    par.pca_coeff_ALL = pca_coeff_ALL;

    fprintf('Removed: %i (%2.2f%%) |  Accepted: %i (%2.2f%%) \n',(allSpikes-goodSpikes),100*(1-goodSpikes/allSpikes),goodSpikes,100*goodSpikes/allSpikes)
    fprintf('Done.\n')

end

%% %%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving function v7.3 works with parallel loops
function parsave(fname, clusterID, spikes, features, index)
    save(fname, 'clusterID', 'spikes', 'features','index','-v7.3')
end

%  Visulaization functions for testing purposes
function plotClusteredSpikes(clusterID,spikes)
    clusterSizes =zeros(1,max(clusterID));
    figure,
    for i =1:max(clusterID)
        subplot(4,5,i)
        clusterSizes(i) = length(clusterID(clusterID==i));
        plot(spikes(clusterID==i,:)','color','blue') 
        title("cluster "+string(i) + "(" + string(clusterSizes(i)) + ")");
    end
end