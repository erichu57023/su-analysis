function percentages = EvaluateSpikeFile(path, max_ISI)
    % refractory = EVALUATESPIKEFILE(path, window)
    %
    % Returns the percentage of spikes in the ntt file with an interspike 
    % interval less than the given time window.
    % 
    % path: (char[]) the absolute path of the NTT file
    % max_ISI: (float) max of the ISI to keep, in ms
    % percentages: (float[]) vector of percentages of spikes within max_ISI 
    %             of the previous spike. Each element corresponds to a
    %             single cluster. If a cluster has 0 or 1 spikes, the
    %             element is NaN.

    %% Import NTT file.
    if ~isfile(path) || ~strcmp(path(end-3:end), '.ntt')
        throw(MException('ExtractISISpikes:FileNotFound', 'The path does not exist, or is not a valid ntt file.'))
    end

    [NTTTimeStamps, ScNumbers, CellNumbers, Features, Samples] = Nlx2MatSpike(path, ones(1,5), 0, 1, []);
    
    %% Calculate percentage of spikes in refractory period of previous spike.
    totalClusters = max(CellNumbers) + 1;
    percentages = zeros(1, totalClusters);
    for cluster = 0 : totalClusters - 1
        timestamps = NTTTimeStamps(CellNumbers == cluster);
        if length(timestamps) <= 1
            percentages(cluster) = NaN;
            continue; 
        end
        
        intervals = diff(timestamps);
        in_refrac = sum(intervals < max_ISI * 1e3);
        percentages(cluster + 1) =  in_refrac / length(intervals);
    end
end
