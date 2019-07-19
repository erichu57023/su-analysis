function windowGraph = ExtractEventSpikes(root, tetrode, cluster, window, exportModifiedNTTFlag)
    % windowGraph = 
    %    EXTRACTEVENTSPIKES(root, tetrode, cluster, exportModifiedNTTFlag)
    % Extracts spikes in the assigned cluster with a timestamp within a
    % certain time window after a TTL ON event.
    %
    % root: (char[]) the path of the folder containing the NTT file
    % tetrode: (int) number of the tetrode to be analyzed
    % cluster: (int) number of the cluster to be analyzed
    % window: (float) length of the time window to keep, in ms
    % exportModifiedNTTFlag: (logical) if true, will save a copy of the
    %       extracted spike file as TT#_events.ntt.
    
    % windowGraph: a handle to the figure created

    %% Import NEV file
    eventFile = [root, '\Events.nev'];
    [EventTimeStamps, TTLs] = Nlx2MatEV(eventFile, [1 0 1 0 0], 0, 1, []);
    EventTimeStamps = EventTimeStamps(logical(TTLs));

    %% Import NTT file, extract event windows, export modified NTT
    try 
        nttFile = [root, '\TT', num2str(tetrode), '_s.ntt'];
        [NTTTimeStamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nttFile, ones(1,5), 1, 1, []);
    catch
        if cluster ~= 0
            throw(MException('EventExtractSpikes:SortedFileNotFound', 'Non-0 clusters require a sorted file (TT#_s.NTT)'))
        else
            try
                nttFile = [root, '\TT', num2str(tetrode), '.ntt'];
                [NTTTimeStamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nttFile, ones(1,5), 1, 1, []);
            catch
                throw(MException('EventExtractSpikes:FileNotFound', ['TT', num2str(tetrode), '.ntt was not found']))
            end
        end
    end

    Samples_volts = single(Samples) * 250e-6 / 32767;
    samplesToKeep = false(1, length(ScNumbers));
    
    clusterIndex = find(CellNumbers == cluster);
    if isempty(clusterIndex)
        throw(MException('EventExtractSpikes:ClusterNotFound', ['Cluster ', num2str(cluster), ' does not exist in tetrode ', num2str(tetrode)])) 
    end

    for jj = 1 : length(EventTimeStamps)
        inEvent = intersect(find(NTTTimeStamps >= EventTimeStamps(jj)), find(NTTTimeStamps <= (EventTimeStamps(jj) + window * 1e3)));
        samplesToKeep(inEvent) = 1;
    end
    
    windowGraph = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    all = [subplot(2, 4, 1), subplot(2, 4, 2), subplot(2, 4, 5), subplot(2, 4, 6)];
    win = [subplot(2, 4, 3), subplot(2, 4, 4), subplot(2, 4, 7), subplot(2, 4, 8)];
    
    linkaxes([all, win]);
    subplot(all(1));
    xlim([1, 32]);
    ylim([-250e-6, 250e-6]);
    
    for ii = [all win]
        subplot(ii);
        hold on;
    end
    
    for ii = 1 : length(all)
        for jj = clusterIndex
            line(all(ii), 1:32, Samples_volts(:, ii, jj), 'Color', 'black');
        end
        average = mean(Samples_volts(:, ii, clusterIndex), 3);
        line(all(ii), 1:32, average, 'Color', 'red');
    end        

    NTTTimeStamps = NTTTimeStamps(samplesToKeep);
    ScNumbers = ScNumbers(samplesToKeep);
    CellNumbers = CellNumbers(samplesToKeep);
    Features = Features(:, samplesToKeep);
    Samples = Samples(:, :, samplesToKeep);
    Samples_volts = Samples_volts(:, :, samplesToKeep);
    
    clusterIndex = find(CellNumbers == cluster);
    if isempty(clusterIndex)
        disp(['No spikes exist in cluster ', num2str(cluster), ' after event-extraction.']);
    else
        for ii = 1 : length(win)
            for jj = clusterIndex
                line(win(ii), 1:32, Samples_volts(:, ii, jj), 'Color', 'black');
            end
            average = mean(Samples_volts(:, ii, clusterIndex), 3);
            line(win(ii), 1:32, average, 'Color', 'red');
        end        
    end

    if exportModifiedNTTFlag
        newNTTFile = [nttFile(1:end-4), '_events.ntt'];
        Mat2NlxSpike(newNTTFile, 0, 1, [], ones(1,6), NTTTimeStamps, ScNumbers, CellNumbers, Features, Samples, Header);
        disp('Event-extracted spikes saved in \TT', num2str(tetrode), '_events.ntt')
    end
end
