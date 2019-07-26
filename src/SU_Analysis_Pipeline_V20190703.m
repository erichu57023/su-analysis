%% this code uses the following functions: biobserve_import, 
% clear all;clc;close all;
%% variables
for Variables=1:2
Date='2019-04-30';
MouseName='SUBLAT1-1'; %%% specify mouse name
TetrodeNumber=[3];% what Tetrodes to run [1 2 3 4 5 6 7 8]
ComputerDir='H:\'; %location of NLX files
BioObserveDir='G:\Single unit videos\Compressed';%location of Bioobserve video files my computer
ExperimentalConditionList={'1hz','2hz','20hz','Empty','Chow','Jelly-Exposure','Jelly-OFF','Jelly_Exposure','Jelly_OFF','Jelly - Exposure','Jelly - OFF'};
attribute = 'height'; % height, pca_1, etc....
cluster_attribute = 'height';
min_pulse_interval = 50000; % (µs) minimum distance between pulses (there's some noise on the TTL)
threshold_method = 'Increasing threshold';%threshold_method = '5sigma';'Increasing threshold';
interval = [1 20]; % Present the original figure and set interval
automated_clustering = false;
show_electrodes = [1 2 4];
nr_clusters = [2, 10];
end % Variables
%% code options
for Code_Options=1:1
OfflineThreshold=false ; LastThreshold=55;
SplitFiles=true ; % indicate 1 for YES 0 for NO
SortedFile=false; % if you wish to run a sorted file in Han's code
ClusteredFile=true;
CollectFiles=true;
ConcatenateFiles=true; 
SortConcatenatedFile=true;
CreateTetrodeObject=true;
OfflineFilter=false;
PlotFigures=true; % indicate if should plot average figures 0 for no 1 for Yes
opto_tagging=false;
SortEmpty=false;
SortChow=false;
SortJellyExposure=false;
SortJellyOFF=false;
Replace_Values=false;
RemoveAllButSelected=false;
SaveAfterCut=true;
UseSortHan=false;
Tetrode2Arena=true;
Cluster3D=true;
%% Add to path and load folders
FileDirNLX = dir([ComputerDir,'\',Date,'\',MouseName]); %Load NLX files
FileDirBioObserve = dir([BioObserveDir,'\',Date,'\',MouseName]);%Load Video files
addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\Field Trip Neuralynx'));
addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\NeuralynxMatlabImportExport_v6.0.0'));
end % CodeOptions
%% Collect files and create concatenated file
if ~SortedFile
for Collect_Files = 1:1
if CollectFiles
% Orgenize NLX folders
ConditionList= {};
for FileNumber = [3:length(FileDirNLX)] % creats ConditionList&ConditionFileOrder
for ExperimentalCondition=1:length(ExperimentalConditionList)% Run this for every possible condition
    if contains(char(FileDirNLX(FileNumber).name),char(ExperimentalConditionList(ExperimentalCondition)))
    ConditionList(FileNumber-2)={char(ExperimentalConditionList(ExperimentalCondition))};
    if contains(char(FileDirNLX(FileNumber).name),'Empty')
        ConditionFileOrder(1)=(FileNumber);
    elseif contains(char(FileDirNLX(FileNumber).name),'Chow')
     ConditionFileOrder(2)=(FileNumber);
    elseif contains(char(FileDirNLX(FileNumber).name),'Exposure')
     ConditionFileOrder(3)=(FileNumber);
     elseif contains(char(FileDirNLX(FileNumber).name),'OFF')
     ConditionFileOrder(4)=(FileNumber);
     end
    break
    end
end % for ExperimentalCondition
end % for FileNumber
% Orgenize Video files %%%%%
if Tetrode2Arena
    VideoConditionList={};
end % if Tetrode2Arena
%% import NLX files into Trial struct and concatenate into ConcatenateStruct
if ConcatenateFiles; ConcatenateStruct(1).name={'Concatenate'};ConcatenateStruct(1).Timestamps_spike=[];ConcatenateStruct(1).ScNumbers_spike=[];ConcatenateStruct(1).CellNumbers_spike=[];ConcatenateStruct(1).Features_spike=[];ConcatenateStruct(1).Samples_spike=[];end % for if ConcatenateFiles
for SelectedCondition=1:length(ConditionList) % import all the data into "Triel" struct
Trial(SelectedCondition).name={char(ConditionList(SelectedCondition))};
[Trial(SelectedCondition).Timestamps_spike, Trial(SelectedCondition).ScNumbers_spike, Trial(SelectedCondition).CellNumbers_spike, Trial(SelectedCondition).Features_spike, Trial(SelectedCondition).Samples_spike, Trial(SelectedCondition).Header_spike] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(SelectedCondition+2).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
Trial(SelectedCondition).Range=[Trial(SelectedCondition).Timestamps_spike(1) Trial(SelectedCondition).Timestamps_spike(end)];
Trial(SelectedCondition).Header_spike{33, 1}(17)='u';
if ConcatenateFiles
ConcatenateStruct(1).Timestamps_spike=[ConcatenateStruct(1).Timestamps_spike,Trial(SelectedCondition).Timestamps_spike];
ConcatenateStruct(1).ScNumbers_spike=[ConcatenateStruct(1).ScNumbers_spike,Trial(SelectedCondition).ScNumbers_spike];
ConcatenateStruct(1).CellNumbers_spike=[ConcatenateStruct(1).CellNumbers_spike,Trial(SelectedCondition).CellNumbers_spike];
ConcatenateStruct(1).Features_spike=[ConcatenateStruct(1).Features_spike,Trial(SelectedCondition).Features_spike];
ConcatenateStruct(1).Samples_spike=cat(3,ConcatenateStruct(1).Samples_spike,Trial(SelectedCondition).Samples_spike);
end % if ConcatenateFiles
end %for SelectedCondition
if ConcatenateFiles
ConcatenateStruct(1).Header_spike=[Trial(ConditionFileOrder(3)-2).Header_spike]; % take header from "Exposure" condition
% export concatenated file to .ntt file and save 
% if ~exist([ComputerDir,'\',Date,'\',MouseName,'\ConcatenatedFile'],'dir')
mkdir([ComputerDir,'\',Date,'\',MouseName],'ConcatenatedFile')
% end % if ~exist
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles_TT',num2str(TetrodeNumber),'.ntt'], 0, 1, [], [1 1 1 1 1 1], ...
    ConcatenateStruct(1).Timestamps_spike,ConcatenateStruct(1).ScNumbers_spike, ConcatenateStruct(1).CellNumbers_spike, ConcatenateStruct(1).Features_spike, ConcatenateStruct(1).Samples_spike, ConcatenateStruct(1).Header_spike);
end % if ConcatenateFiles
end % CollectFiles
end % Collect_files
%% create tetrode object
for Create_Tetrode_Obj=1:1
if CreateTetrodeObject
if SortConcatenatedFile
% SortFileFolder = [ComputerDir,'\',Date,'\',MouseName,'\ConcatenatedFile'];
% ConcatenatedFileName=['AllFiles',num2str(TetrodeNumber),'.ntt'];
% SortFileName=['AllFiles_TT',num2str(TetrodeNumber),'.ntt'];
Tetrode = tetrode([ComputerDir,'\',Date,'\',MouseName,'\ConcatenatedFile\AllFiles_TT',num2str(TetrodeNumber),'.ntt'],ClusteredFile);
% if working on not - concatenated
elseif SortEmpty||SortChow||SortJellyExposure||SortJellyOFF
if SortEmpty
    ConditionNumber=1
elseif SortChow
    ConditionNumber=2
elseif SortJellyExposure
    ConditionNumber=3
elseif SortJellyOFF
    ConditionNumber=4
end 
SortFileFolder = [ComputerDir,'\',Date,'\',MouseName];
Tetrode = tetrode([SortFileFolder,'\',char(FileNameList(ConditionFileOrder(ConditionNumber))),'\','TT',num2str(TetrodeNumber),'.ntt']);
end %if SortConcatenatedFile
end % CreateTetrodeObj
end % Create_Tetrode_Obj
%% Filter the concatenated file
%% Set an offline threshold
for Offline_Threshold=1:1
if OfflineThreshold
%% Remove noise (mostly lickometer)% pick 90% of the inputrange
Tetrode.remove_noise('clipping');close all;
%% Present the original figure and set interval
Tetrode.set_interval(interval);
[precluster_figure_1, precluster_figure_2] = Tetrode.present_figure(cluster_attribute, show_electrodes);
%% Set an offline threshold
if strcmp(threshold_method, '5sigma')
    Tetrode.set_5sigma_threshold();
    else
    % Keep increasing the threshold untill the user says stop.
    threshold = Tetrode.settings.ThresVal;
    keep_running = true;
    while(keep_running)
        threshold = threshold + 5;
        keep_running = Tetrode.set_offline_threshold(threshold);
    end
    disp(['Final threshold set to : ' num2str(threshold-5) ' nµV']);
end
% Save sorted cell numbers
%%%%% need this?
    SortHan=Tetrode.cells;
    save([ComputerDir,'\',Date,'\',MouseName,'\ConcatenatedFile\T',num2str(TetrodeNumber),' Sorted Cell Numbers Concatenated Han.mat'],'SortHan');
if UseSortHan
    try CellNumbersAll=SortHan; catch; disp('first run sort.m');end
    Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','ALLFILES_TT',num2str(TetrodeNumber),'_s_','Han','.ntt'], 0, 1, [], [1 1 1 1 1 1], TimestampsAll,ScNumbersAll, CellNumbersAll, FeaturesAll, SamplesAll, Header_spike_Jelly_Exposure);
end %if UseSortHan
%% the following code collects only the selected swips and makes a new .NTT file
if SaveAfterCut
    if SortedFile
ConcatenateStruct(2).name={'Clustered_Concatenate'};
[ConcatenateStruct(2).Timestamps_spike,ConcatenateStruct(2).ScNumbers_spike,ConcatenateStruct(2).CellNumbers_spike, ConcatenateStruct(2).Features_spike, ConcatenateStruct(2).Samples_spike, ConcatenateStruct(2).Header_spike] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles_TT',num2str(TetrodeNumber),'_s.ntt'], [1 1 1 1 1], 1, 1, [] ); 
UseData=2;
else
UseData=1;
end
% collect only the remaining cells:
IndexNan=~isnan(Tetrode.cells);
TimestampsThreshold=ConcatenateStruct(UseData).Timestamps_spike(IndexNan);
ScNumbersThreshold=ConcatenateStruct(UseData).ScNumbers_spike(IndexNan);
CellNumbersThreshold=ConcatenateStruct(UseData).CellNumbers_spike(IndexNan);
FeaturesThreshold=ConcatenateStruct(UseData).Features_spike(:,IndexNan);
SamplesThreshold=ConcatenateStruct(UseData).Samples_spike(:,:,IndexNan);
HeaderThreshold=ConcatenateStruct(UseData).Header_spike;
    if SortedFile
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles_TT',num2str(TetrodeNumber),'_s_Cut.ntt'], 0, 1, [], [1 1 1 1 1 1], TimestampsThreshold,ScNumbersThreshold, CellNumbersThreshold, FeaturesThreshold, SamplesThreshold, HeaderThreshold);
    else
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles_TT',num2str(TetrodeNumber),'_Cut.ntt'], 0, 1, [], [1 1 1 1 1 1], TimestampsThreshold,ScNumbersThreshold, CellNumbersThreshold, FeaturesThreshold, SamplesThreshold, HeaderThreshold);
    end
end
else 
threshold=[LastThreshold,LastThreshold,LastThreshold,LastThreshold];
end % if OfflineThreshold
end % for Offline_Threshold
%% Check if any of the cells are tagged
for Opto_Tagging=1:1
if opto_tagging
for Laser=3:length(FileDirNLX)
if contains(char(FileDirNLX(Laser).name),'Laser')
events = read_neuralynx_nev(([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(Laser).name),'\','Events.nev']));
break
end
end %for Laser
    for LineNumber=1:length(events)
        timestamps(LineNumber) = events(LineNumber).TimeStamp;
        TTL(LineNumber) = events(LineNumber).TTLValue;
    end
    stamps = timestamps(TTL==1);
    [OptoTaggingPlot_handle, OptoTaggingResults]=Tetrode.peri_event(stamps,'window',200*10^3);
end
end % if Opto_Tagging
end % if SortedFile
%% Sort concatenated file in 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster
for hanCluster=1:1
if ~Cluster3D
% if automated_clustering
%     automated_clustering = true; %for quick-select
%     new_cells = cluster_han(cluster_attribute, Tetrode, 'auto_gm', nr_clusters);
%     automated_clustering = false; % for next round
% else
%% open selected feature for manual clustering
    new_cells = cluster_han(cluster_attribute, Tetrode);
% end
% Figure out new NaN values
new_noise = isnan(new_cells);
new_noise(isnan(Tetrode.cells)) = false;
% Remove waveforms labeled as noise (but check with the user).
if sum(new_noise)>0
    Tetrode.remove_noise('indexer', new_noise);
end
% Update the cells property
Tetrode.cells(~new_noise) = new_cells(~new_noise);
Tetrode.nr_cells = max(Tetrode.cells); 
%    for SelectedCell=1:5
% Tetrode.Show_cell(SelectedCell)
% % ClusterSummary(SelectedCell).Quality=Tetrode.unit_quality
% end
%% Remove all but selected cluster
if RemoveAllButSelected 
Cluster0=Tetrode.cells==0;
Cluster1=Tetrode.cells==1;
Cluster2=Tetrode.cells==2;
Cluster3=Tetrode.cells==3;
Cluster4=Tetrode.cells==4;
for NotSelected=0:10
    if NotSelected~=SelectedCluster
Tetrode.remove_noise('remove cluster',NotSelected)
    end %if
end
end %if
%% Close the previous figures
close all
%% Replace values in Tetrode object
if Replace_Values
    Tetrode.raw_data.CellNumber=Tetrode.cells;
end
%% Save Tetrode object
if SaveResults
save([SortFileFolder,'\','T',num2str(TetrodeNumber),' Tetrode object Han.mat'],'Tetrode');
end
%% display analysis results
% Print the current clustogram
Tetrode.present_figure(cluster_attribute, show_electrodes);
if PlotHistogram    
Sec=1*10^6;
    ms=1*10^3;
    Time=Sec;
    Time2Inspect =2;%time in sec
    Time2Inspect=Time2Inspect*10^6;
    HistoFigure=figure;
    CatchCell=Tetrode.cells==SelectedCluster;
    CellToHistTimestamps=Tetrode.timestamps(CatchCell);
    CellToHist=CellToHistTimestamps-min(CellToHistTimestamps); % normalyse to minimal value
    NumberOfBins=max(CellToHist)/Time;% calculate the number of bins to get 1sec bins
    histogram(CellToHist,NumberOfBins)% plot histogram
    [NumberOfSpikes,BinEdges] = histcounts(CellToHist,NumberOfBins);% get values
    title([char(ConditionList(2)),' ',char(ConditionList(3)),' ',char(ConditionList(4)),' ',char(ConditionList(5))])
    saveas(HistoFigure,[ComputerDir,'\',Date,'\',MouseName,'\ConcatenatedFile\',['TimeHistogram_Cluster_',num2str(SelectedCluster)]],'jpg')    
   end % if PlotHistogram
end %Cluster3D
end % for hanCluster
%% Split files back to original location
 for Split_Files=1:1
% Import the sorted file 
if SplitFiles
try
ConcatenateStruct(3).name={'Clustered_Cut_Concatenate'};
[ConcatenateStruct(3).Timestamps_spike,ConcatenateStruct(3).ScNumbers_spike,ConcatenateStruct(3).CellNumbers_spike, ConcatenateStruct(3).Features_spike, ConcatenateStruct(3).Samples_spike, ConcatenateStruct(3).Header_spike] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles_TT',num2str(TetrodeNumber),'_Cut_s.ntt'], [1 1 1 1 1], 1, 1, [] ); 
try ConcatenateStruct(3).Threshold=threshold-5;catch;end
UseData=3;
FileSorted=true;
catch
    disp('you must sort the file first!!! save in format FileName_s')
FileSorted=false;
end
% Split the data according to the timestamp range of each trial
if FileSorted
clear TrialCut
for Condition=1:length(Trial)
clear IndexRange
IndexRange=ConcatenateStruct(3).Timestamps_spike<Trial(Condition).Timestamps_spike(end)&ConcatenateStruct(3).Timestamps_spike>Trial(Condition).Timestamps_spike(1);
TrialCut(Condition).name=Trial(Condition).name;
TrialCut(Condition).Timestamps_spike=ConcatenateStruct(3).Timestamps_spike(IndexRange);
TrialCut(Condition).ScNumbers_spike=ConcatenateStruct(3).ScNumbers_spike(IndexRange);
TrialCut(Condition).CellNumbers_spike=ConcatenateStruct(3).CellNumbers_spike(IndexRange);
TrialCut(Condition).Features_spike=ConcatenateStruct(3).Features_spike(:,IndexRange);
TrialCut(Condition).Samples_spike=ConcatenateStruct(3).Samples_spike(:,:,IndexRange);
TrialCut(Condition).Header_spike=Trial(Condition).Header_spike;
try TrialCut(Condition).threshold=threshold-5;catch ; TrialCut(Condition).threshold=50;end
% export the files to original locations
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(Condition+2).name),'\','TT',num2str(TetrodeNumber),'_Cut',num2str(threshold(1)-5),'.ntt'], 0, 1, [], [1 1 1 1 1 1],...
TrialCut(Condition).Timestamps_spike,TrialCut(Condition).ScNumbers_spike, TrialCut(Condition).CellNumbers_spike, TrialCut(Condition).Features_spike, TrialCut(Condition).Samples_spike, TrialCut(Condition).Header_spike);
if PlotFigures
TrialCut(Condition).Figure=figure;
% find the header line with ADBitVolts information and stop there
for LineNumber=1:length(TrialCut(Condition).Header_spike)
if contains([(char(TrialCut(Condition).Header_spike{LineNumber,:}))],'ADBitVolts')
break
   end %if
end %for
% split the string to individule values and store it in TrialCut.ADVoltsList
ADVoltsList=strsplit((char(TrialCut(Condition).Header_spike{LineNumber,:})),' ');
TrialCut(Condition).ADVoltsList =[str2double(char(ADVoltsList(2))),str2double(char(ADVoltsList(3))),str2double(char(ADVoltsList(4))),str2double(char(ADVoltsList(5)))];
% plot all the electrodes
for Electrode=2:5
TrialCut(Condition).Samples_MicroVolts(:,Electrode-1,:)=TrialCut(Condition).Samples_spike(:,Electrode-1,:)*(TrialCut(Condition).ADVoltsList(Electrode-1))*1000000;% 1000000 to show in microVolts
plot(mean(TrialCut(Condition).Samples_MicroVolts,3));
title(TrialCut(Condition).name);
end %for Electrode
end % if PlotFigures
end % for Condition
end % if FileSorted
end %if SplitFiles
end % for Split_Files
%% Analyze the data:
% Responsive unit detection
% Heat maps different food types + firing properties (run “Tetrode to arena.m”)
%Meta-analyse each condition - Merge each condition into one heatmap (run “average image.m”)
