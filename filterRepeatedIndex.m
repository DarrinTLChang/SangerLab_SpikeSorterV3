
function [nonRepeatedIndex] = filterRepeatedIndex(indexAll,fs,par)
    
    % Getting event times in ms
    allSpikesIndex = cell2mat(indexAll');
    allSpikesIndex_ms = allSpikesIndex/fs*1e3;
    
    % Create a histogram to find very close events
    edges = min(allSpikesIndex_ms):par.ms_repeated_index:max(allSpikesIndex_ms);
    [N,~,bin] = histcounts(allSpikesIndex_ms,edges);

    % Find reapeated indices
    repeatedEvents = find(N>par.max_num_coincedences);
    repIndex = [];
    parfor i = 1:numel(repeatedEvents)
        repIndex = [repIndex find(bin == repeatedEvents(i))];
    end

    [~, nIndexPerChannel] = cellfun(@size, indexAll); 
    % Save indices
    nonRepeatedIndex_tmp = allSpikesIndex;
    nonRepeatedIndex_tmp(repIndex) = 0;
    
    % Generate output
    for i = 1:length(nIndexPerChannel)
        if numel(nIndexPerChannel(i)) ~= 0 
            index_tmp = nonRepeatedIndex_tmp(1:nIndexPerChannel(i)); 
            nonRepeatedIndex{i} = index_tmp(index_tmp~=0);
            nonRepeatedIndex_tmp(1:nIndexPerChannel(i)) = [];
        else
            nonRepeatedIndex{i} = [];
        end
    end    
end

