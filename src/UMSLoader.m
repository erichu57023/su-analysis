clearvars;
tic;

%% Setup
root = 'C:\Users\erich\Dropbox\Research Files (Eric)\su-analysis\test\SUBLAT3-7\';
fileNumbers = [2 : 5];    % select multiple files to collate
channels = [21 : 24];   % select channels to sort

trialLength = 15;           % trial length in seconds
sampleRate = 32000;         % sampling rate in Hz

localAverageFlag = true;       % subtracts local average from all active channels
deadChannels = [];          % specify dead channels to skip when averaging

%% Averaging all channels
tic;
liveChannels = setdiff(1 : 32, deadChannels);
fileNumbers = fileNumbers + 2;
nextChannel = cell(4,1);
for ii = 1 : length(liveChannels)
    rootDir = dir(root);
    for jj = 1 : length(fileNumbers)
        try
            fileName = rootDir(fileNumbers(jj)).name;
        catch
            throw(MException('UMSTest:fileNotFound', ['File number ', num2str(fileNumbers(jj) - 2), ' does not exist.']));
        end
        file = [root, fileName, '\CSC', num2str(liveChannels(ii)), '.ncs'];
        nextChannel{jj} = reshape(single(Nlx2MatCSC(file, [0 0 0 0 1], 0, 1, [])),[],1);
    end
    try
        average = average + cell2mat(nextChannel);
    catch
        average = cell2mat(nextChannel);
    end
    disp(['Completed ', num2str(ii)]);
end
average = average / length(liveChannels);
disp(['Averaging: ', num2str(toc), ' seconds']);

%% Channel data import (variable length experiment; >1 hr not recommended)
samples = [];
fileNumbers = fileNumbers + 2;
for ii = 1 : length(channels)
    rootDir = dir(root);
    nextChannel = [];
    for jj = 1 : length(fileNumbers)
        try
            fileName = rootDir(fileNumbers(jj)).name;
        catch
            throw(MException('UMSTest:fileNotFound', ['File number ', num2str(fileNumbers(jj) - 2), ' does not exist.']));
        end
        file = [root, fileName, '\CSC', num2str(channels(ii)), '.ncs'];
        nextChannel = vertcat(nextChannel, reshape(int16(Nlx2MatCSC(file, [0 0 0 0 1], 0, 1, [])),[],1));
    end
    samples = horzcat(samples, nextChannel);
end
clear nextChannel;
disp(['Loading: ', num2str(toc), ' seconds']);

%% Sectioning into trials
tic;
N = length(samples);
trialLength = trialLength * sampleRate;
numTrials = ceil(N / trialLength);
trialSizes = ones(numTrials, 1) * trialLength;
trialSizes(end) = mod(N, trialLength);
samples = single(samples) * 2e-3 / 32767;
data = mat2cell(samples, trialSizes)';
clear file fileName i j localAverageFlag N numTrials rootDir trialLength trialSizes;
disp(['Sectioning: ', num2str(toc), ' seconds']);

%% Running UMS
tic;
spikes = ss_default_params_TL(sampleRate);
spikes = ss_detect_BI(data,spikes);
spikes = ss_align_BI(spikes);
spikes = ss_kmeans_BI(spikes);
spikes = ss_energy_parallel_BI(spikes);
spikes = ss_aggregate_BI(spikes);
splitmerge_tool(spikes) 

disp(['Total: ', num2str(toc), ' seconds']);
