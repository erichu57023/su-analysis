clearvars;

%% Setup
root = 'C:\Users\erich\Dropbox\Research Files (Eric)\su-analysis\test\SUBLAT3-7\';
fileNumber = 1;   % select file in root directory
tetrodes = 6;     % select tetrodes to import
window = 8;       % time window for post-event spikes (in ms)

%% Import NEV file
rootDir = dir(root);
try
    fileName = rootDir(fileNumber + 2).name;
catch
    throw(MException('ExtractEventSpikes:fileNotFound', ['File number ', num2str(fileNumber), ' does not exist.']));
end
eventFile = [root, fileName, '\Events.nev'];
[EventTimeStamps, TTLs] = Nlx2MatEV(eventFile, [1 0 1 0 0], 0, 1, []);
EventTimeStamps = EventTimeStamps(logical(TTLs));

%% Import NTT files, extract event windows, export modified NTT
for ii = 1 : length(tetrodes)
    nttFile = [root, fileName, '\TT', num2str(tetrodes(ii)), '.ntt'];
    [NTTTimeStamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nttFile, ones(1,5), 1, 1, []);
    samplesToKeep = false(1, length(ScNumbers));
    
    for jj = 1 : length(EventTimeStamps)
        inEvent = intersect(find(NTTTimeStamps >= EventTimeStamps(jj)), find(NTTTimeStamps <= (EventTimeStamps(jj) + window * 1e3)));
        samplesToKeep(inEvent) = 1;
    end
    
    NTTTimeStamps = NTTTimeStamps(samplesToKeep);
    ScNumbers = ScNumbers(samplesToKeep);
    CellNumbers = CellNumbers(samplesToKeep);
    Features = Features(:, samplesToKeep);
    Samples = Samples(:, :, samplesToKeep);
    
    newNTTFile = [nttFile(1:end-4), '_events.ntt'];
    Mat2NlxSpike(newNTTFile, 0, 1, [], ones(1,6), NTTTimeStamps, ScNumbers, CellNumbers, Features, Samples, Header);
end