
%% Variables
% addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\NeuralynxMatlabImportExport_v6.0.0'));
addpath(genpath('E:\analysis\Neuralynx Scripts\MatlabImportExport_v6.0.0'));
% ComputerDir='C:\Users\netas\Google Drive\Single Unit Data';
ComputerDir='E:\CheetahData\NG\Data';
Date='2019-03-06';
MouseName='SUBLAT1-5'; %%% specify mouse name
FileDir = dir([ComputerDir,'\',Date,'\',MouseName]);

%% Choose folders and tetrodes to run
ChooseFileNumber = [3]; % what files to run
ChooseTetrodeNumber=[1 2 3 4 5 6 7 8]; % what Tetrodes to run
% FileNameList=extractfield(FileDir,'name');char(FileNameList(FileNumber))

%% Function
for FileNumber = ChooseFileNumber+2
    FileDir(FileNumber).name % display file name
    FileFolder = [ComputerDir,'\',Date,'\',MouseName,'\',FileDir(FileNumber).name];
%now we split the name and date componantes:
%%%Excel format: [=D22&"_"&E22&" min_SU8-3_"&"z"&B22&"_"&F22&"_"&G22&"_"&H22&"_"&I22&"_"&round(J22,3)&" "&"gr"]
SplitName=strsplit(FileDir(FileNumber).name,'_');
TimeSplit = strsplit(char(SplitName(2)),'-');
%now we assign the variables with the componanets we got from the file name 
try DateStamp = char(SplitName(1)); catch; DateStamp = ''; end
try TimeStamp = [char(TimeSplit(1)),'-',char(TimeSplit(2))]; catch; TimeStamp = ''; end
try RecDuration = char(SplitName(3)); catch; RecDuration = ''; end
try MouseNumber = char(SplitName(4)); catch; MouseNumber = MouseName; end
try Zlocation = char(SplitName(5));catch; Zlocation = ''; end
try TTLLaser = char(SplitName(6));catch; TTLLaser = ''; end
try StimType=char(SplitName(7));catch; StimType = ''; end
try StimFreq=char(SplitName(8));catch; StimFreq = ''; end
try FoodType=char(SplitName(9));catch; FoodType = ''; end
try FoodConsumed=[char(SplitName(10)),' gr'];catch; FoodConsumed = ''; end

%% now we run the function for the selected tetrodes and clusters
%we build an arry to put in our data in:
DataSummary=zeros(20,8);
   count=0;
   
for TetrodeNumber = ChooseTetrodeNumber 
    Timer= tic;
    timeLimit = 5;
    try % first try and run a clusteres file
    [TimeStamps_spike, ScNumbers_spike, CellNumbers_spike, Features_spike, Samples_spike, Header_spike] = Nlx2MatSpike([FileFolder,'\',['TT',num2str(TetrodeNumber),'_s.ntt']], [1 1 1 1 1], 1, 1, [] ); %spike file
    ChooseClusters = [1 2 3 4 5 6];
    catch % If you can't find a clustered file, run the raw file
    [TimeStamps_spike, ScNumbers_spike, CellNumbers_spike, Features_spike, Samples_spike, Header_spike] = Nlx2MatSpike([FileFolder,'\',['TT',num2str(TetrodeNumber),'.ntt']], [1 1 1 1 1], 1, 1, [] ); %spike file
    ChooseClusters = [0];
    end
    try % the following try-catch loops are aimed to avoid a case of a disabled channle. each of the 4 tetrode channles can do the trick
    ChannleNumber = TetrodeNumber*4;
    [TimeStamps_csc, ChannelNumbers_csc, SampleFrequencies_csc, NumberOfValidSamples_csc, Samples_csc, Header_csc] = Nlx2MatCSC([FileFolder,'\',['CSC',num2str(ChannleNumber),'.ncs']], [1 1 1 1 1], 1, 1, [] ); %csc file
    catch 
        try
    ChannleNumber = TetrodeNumber*4-1
    [TimeStamps_csc, ChannelNumbers_csc, SampleFrequencies_csc, NumberOfValidSamples_csc, Samples_csc, Header_csc] = Nlx2MatCSC([FileFolder,'\',['CSC',num2str(ChannleNumber),'.ncs']], [1 1 1 1 1], 1, 1, [] ); %csc file
        catch
            try
    ChannleNumber = TetrodeNumber*4-2
    [TimeStamps_csc, ChannelNumbers_csc, SampleFrequencies_csc, NumberOfValidSamples_csc, Samples_csc, Header_csc] = Nlx2MatCSC([FileFolder,'\',['CSC',num2str(ChannleNumber),'.ncs']], [1 1 1 1 1], 1, 1, [] ); %csc file
            catch
    ChannleNumber = TetrodeNumber*4-3
    [TimeStamps_csc, ChannelNumbers_csc, SampleFrequencies_csc, NumberOfValidSamples_csc, Samples_csc, Header_csc] = Nlx2MatCSC([FileFolder,'\',['CSC',num2str(ChannleNumber),'.ncs']], [1 1 1 1 1], 1, 1, [] ); %csc file
            end
        end
    end

for clust = ChooseClusters; 
    
time = 60; %time (s) of interest

TimeStamps_csc_zeroed_s = (TimeStamps_csc-TimeStamps_csc(1))/1000000;

cell = cat(1,TimeStamps_spike,CellNumbers_spike); 
cell_index = find(cell(2,:) ~=clust); % "~=1" looks at cell#1, "~=2" looks at cell#2, etc.
cell(:,cell_index) = [];
if isempty(cell)
    continue
end
[c index] = min(abs(TimeStamps_csc_zeroed_s-time));

cell_index_time = find(cell(1,:) > TimeStamps_csc(index)); %finds the indices higher than 120s
cell(:,cell_index_time) = [];
cell = cell/1000000;
cell = cell(1,:) - cell(1,1);

cell(2,:) = 1;

%% runs burst detection algorithm

tn = cell(1,:);
limit = 0.1;
RSalpha = -log(0.05);

[archive_burst_RS, archive_burst_length, archive_burst_start] = burst(tn,limit,RSalpha);

%% creates an array ("burst_all_final") with burst indices

burst_end = zeros(1,length(archive_burst_start));

for i = 1:length(burst_end)
    burst_end(i) = archive_burst_start(i) + (archive_burst_length(i)-1);
end

burst_all_matrix = zeros(length(archive_burst_length),length(archive_burst_length));

for i = 1:length(archive_burst_start)
    burst_all_matrix(1:archive_burst_length(i),i) = archive_burst_start(i):burst_end(i);
end

burst_all_final = zeros(1,sum(archive_burst_length));

for i = 1:length(archive_burst_start)
    A = find(burst_all_final == 0);
    burst_all_final(A:A+(archive_burst_length(i)-1)) = burst_all_matrix(1:archive_burst_length(i),i);
end

%% plots indentified bursts

burst_time_s = zeros(1,length(burst_all_final));

for i = 1:length(burst_all_final)
    burst_time_s(i) = cell(1,burst_all_final(i));
end

burst_time_s(2,:) = 2;

% stem(cell(1,:),cell(2,:))
% hold on
% stem(burst_time_s(1,:),burst_time_s(2,:))    

%% firing frequency

% time = 300; %!!!!!!!!!!!!!!!!!!!!!!!delete this if you are looking at only one 2 min recording!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

frequency = length(cell)/time;
length(cell);

%% tonic firing frequency

burst_final_temp = burst_time_s(1,:);
cell_temp = cell(1,:);

[C,ia,ib] = intersect(cell_temp,burst_time_s);
cell_temp(ia) = [];
tonic_frequency = length(cell_temp)/time;

%% spikes_in_bursts

spikes_in_bursts = length(burst_time_s(1,:))/length(cell(1,:))*100;

%% spikes_per_burst

spikes_per_burst = mean(archive_burst_length);

%% interburst frequency
try
interburst_frequency = length(archive_burst_start) / (cell(1,archive_burst_start(end)) - cell(1,archive_burst_start(1)));
catch
    continue
end
%% intraburst frequency

intra_temp = zeros(1,length(archive_burst_start));

for i = 1:length(archive_burst_start)
    intra_temp(i) = archive_burst_length(i) / (cell(1,burst_end(i)) - cell(1,archive_burst_start(i)));
end

a = find (intra_temp >= 1000);
intra_temp(a) = [];

intraburst_frequency = mean(intra_temp);

%% display

U = ['Firing frequency ',num2str(frequency)];
V = ['Tonic frequency ',num2str(tonic_frequency)];
W = ['% of spikes in bursts ',num2str(spikes_in_bursts)];
X = ['Intraburst frequency ',num2str(intraburst_frequency)];
Y = ['Interburst frequency ',num2str(interburst_frequency)];
Z = ['Spikes per burst ',num2str(spikes_per_burst)];

% disp(U)
% disp(V)
% disp(W)
% disp(X)
% disp(Y)
% disp(Z)
count=count+1;
 DataSummary(count,1)=TetrodeNumber;
 DataSummary(count,2)=clust;
 DataSummary(count,3)=frequency;
 DataSummary(count,4)=tonic_frequency;
 DataSummary(count,5)=spikes_in_bursts;
 DataSummary(count,6)=intraburst_frequency;
 DataSummary(count,7)=interburst_frequency;
 DataSummary(count,8)=spikes_per_burst;
 
 
end
end

DataSummary=round(DataSummary,2)
%% Export and save all the data to an Excel file and save each condition in a specific folder
% try
SummaryFileName = ['DataSummary',' ',MouseName,' ',Date,' ',TTLLaser,' ',StimType,' ',StimFreq,' ',FoodType,' ',FoodConsumed];
Col_Title={'Tetrode #','cluster #','frequency','tonic frequency','spikes in bursts','intraburst frequency','interburst frequency','spikes per burst'};
ExperimentInfo={'Experiment info ',['mouse ',MouseNumber],['Date ',DateStamp],['Time Stamp ',TimeStamp],['recording Duration ',RecDuration],['Z=',Zlocation],[TTLLaser],['Stimulus: ',StimType,' ',StimFreq],['Food Type: ',FoodType],['Food Consumed: ',FoodConsumed]}';
        try SheetName= [num2str(FileNumber),'_',TimeStamp,'_',Zlocation,'_',TTLLaser(1:2),'_',FoodType];catch;SheetName=[num2str(FileNumber)];end
        try
        xlswrite(fullfile(FileFolder, SummaryFileName),DataSummary,SheetName,'B2');
        xlswrite(fullfile(FileFolder, SummaryFileName),Col_Title,SheetName,'B1');     %Write column title
        xlswrite(fullfile(FileFolder, SummaryFileName),ExperimentInfo,SheetName,'A1');     %Write column title
        catch
        SheetName= SheetName(1:31); 
        xlswrite(fullfile(FileFolder, SummaryFileName),DataSummary,SheetName,'B2');
        xlswrite(fullfile(FileFolder, SummaryFileName),Col_Title,SheetName,'B1');     %Write column title
        xlswrite(fullfile(FileFolder, SummaryFileName),ExperimentInfo,SheetName,'A1');     %Write column title
        end
% catch
%       'please close the file and try again'
% 
% end
%% save all to one excel file:
    FileFolder = [ComputerDir,'\',Date,'\',MouseName];
    SummaryFileName = ['mouse ',MouseNumber,'_',DateStamp];
try
Col_Title={'Tetrode #','cluster #','frequency','tonic frequency','spikes in bursts','intraburst frequency','interburst frequency','spikes per burst'};
ExperimentInfo={'Experiment info ',['mouse ',MouseNumber],['Date ',DateStamp],['Time Stamp ',TimeStamp],['recording Duration ',RecDuration],['Z=',Zlocation],[TTLLaser],['Stimulus: ',StimType,' ',StimFreq],['Food Type: ',FoodType],['Food Consumed: ',FoodConsumed]}';
        try SheetName= [num2str(FileNumber),'_',TimeStamp,'_',Zlocation,'_',TTLLaser(1:2),'_',FoodType];catch;SheetName=[num2str(FileNumber)];end
        try
        xlswrite(fullfile(FileFolder, SummaryFileName),DataSummary,SheetName,'B2');
        xlswrite(fullfile(FileFolder, SummaryFileName),Col_Title,SheetName,'B1');     %Write column title
        xlswrite(fullfile(FileFolder, SummaryFileName),ExperimentInfo,SheetName,'A1');     %Write column title
        catch try
        SheetName= SheetName(1:31); 
        xlswrite(fullfile(FileFolder, SummaryFileName),DataSummary,SheetName,'B2');
        xlswrite(fullfile(FileFolder, SummaryFileName),Col_Title,SheetName,'B1');     %Write column title
        xlswrite(fullfile(FileFolder, SummaryFileName),ExperimentInfo,SheetName,'A1');     %Write column title
            catch
               DataSummary 
            end
        end
catch
       'please close the file and try again'

end

end %For

   
   
   
   

%TableSummary = array2table(DataSummary,'VariableNames',{'TetrodeNumber','clust','frequency','tonic_frequency','spikes_in_bursts','intraburst_frequency','interburst_frequency','spikes_per_burst'})

%% spike autocorrelogram

% lag_time = 0.5; %in seconds
% 
% ISI_auto = zeros(length(tn),length(tn));
% 
% for i = 1:length(tn)
%     ISI_auto(i,:) = tn(i) - tn;
% end
% 
% ISI_auto = reshape(ISI_auto,1,[]);
% 
% k = find(ISI_auto >= -lag_time & ISI_auto <= lag_time);
% 
% ISI_auto = ISI_auto(k);
% 
% k = find(ISI_auto == 0);
% 
% ISI_auto(k) = [];
% 
% edges = [-lag_time:0.002:lag_time];
% histogram(ISI_auto,edges,'Normalization','probability')



