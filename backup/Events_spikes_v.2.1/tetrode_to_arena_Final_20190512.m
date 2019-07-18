%% Variables
close all
clear all
ChooseFileNumber = [2:5]; % what files to run %Exposure
ChooseTetrodeNumber=[8]; % what Tetrodes to run [1 2 3 4 5 6 7 8]
ChooseClustNumber=[3]; 
Date='2019-05-10';
MouseName='SUBLAT1-1'; %%% specify mouse name
ConditionList= {'Laser_2hz','Laser_20hz''Empty','Chow','Jelly_Exposure','Jelly_OFF'};
%% optional code changes
TTL=1; %indicate 'ON'(TTL=1) / 'OFF'(TTL=0)
SaveStatistics=0;%indicate 'YES'(1) / 'NO'(0)
SaveImages=1;%indicate 'YES'(1) / 'NO'(0)
Colorbar=0;%indicate 'YES'(1) / 'NO'(0)
SubplotFigures=1;%indicate 'YES'(1) / 'NO'(0)
%% other variables
% Gaussian kernel?
Gaussian_kernel = true;
sigma = 2; % Sigma of the kernel e^(-1*data/2*sigma^2)
% map binning
RangeMin = 0; % minimal value for heatmap
RangeMax = 80;% maximal value for heatmap
timestamp_interval = 5; %Sec
nr_bins = 20; %Actual the square root of the number of bins
% caxis([0 140])
% addpath(('C:\Users\admin\Documents\MATLAB\Field Trip Neuralynx'));
%arena = [200 -290 160 80];
%% add path to NLX import tools
addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\Field Trip Neuralynx'));
addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\NeuralynxMatlabImportExport_v6.0.0'));

%% Start code
SummaryList=zeros(9,9);
NameList={};
if SubplotFigures==1
SubplotFigure=figure('units','normalized','outerposition',[0 0 1 1]);
end
for FileNumber = ChooseFileNumber+2

%% File location
ComputerDir='E:\'; %location of NLX files my computer
BioObserveDir='G:\Single unit videos\Compressed';%location of Bioobserve video files my computer
HeatMapDir='E:\Maps';
%%

%% Choose folders and tetrodes to run
% FileNameList=extractfield(FileDir,'name');char(FileNameList(FileNumber))
%% write down which file numbers to run in the loop
FileDirNLX = dir([ComputerDir,'\',Date,'\',MouseName]);
FileDirBioObserve = dir([BioObserveDir,'\',Date,'\',MouseName]);
ConditionNumber=FileNumber-2;
VideoFileNumber=FileNumber-1;
FileName=FileDirNLX(FileNumber).name % display file name
FileFolder = [ComputerDir,'\',Date,'\',MouseName,'\',FileDirNLX(FileNumber).name];
%% now we split the name and date componantes:
%%%Excel format: [=D22&"_"&E22&" min_SU8-3_"&"z"&B22&"_"&F22&"_"&G22&"_"&H22&"_"&I22&"_"&round(J22,3)&" "&"gr"]
SplitName=strsplit(FileDirNLX(FileNumber).name,'_');
TimeSplit = strsplit(char(SplitName(2)),'-');
%now we assign the variables with the componanets we got from the file name 
try DateStamp = char(SplitName(1)); catch; DateStamp = ''; end
try TimeStamp = [char(TimeSplit(1)),'-',char(TimeSplit(2))]; catch; TimeStamp = ''; end
try RecDuration = char(SplitName(3)); catch; RecDuration = ''; end
try MouseNumber = char(SplitName(4)); catch; MouseNumber = MouseName; end
try Zlocation = char(SplitName(5));catch; Zlocation = ''; end
try TTLLaser = char(SplitName(6));catch; TTLLaser = ''; end
try FoodType=char(SplitName(7));catch; StimType = ''; end;
try FoodConsumed=char(SplitName(8));catch; StimFreq = ''; end

%% Import the files % write down which tetrodes to run in the loop
    SerialNumber=FileNumber-ChooseFileNumber(1)-1;
for TetrodeNumber = ChooseTetrodeNumber
    try
    % specifies filename: the correct filename format is "TT1_s.ntt","TT2_s.ntt" atc.
Tetrode = tetrode([FileFolder,'\','TT',num2str(TetrodeNumber),'_s.ntt']);
ClustNumber=ChooseClustNumber;  
    catch
        try
Tetrode = tetrode([FileFolder,'\','TT',num2str(TetrodeNumber),'.ntt']);
ClustNumber=0;  
disp('unclustered file analysed')
        catch
            disp('Cannot load file');
        continue
        end
    end
Events = read_neuralynx_nev([FileFolder,'\','Events.nev']);
Track = biobserve_import([BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',char(ConditionList(ConditionNumber)),'.csv']); %import video file
disp(['BioObserve file', [BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',char(ConditionList(ConditionNumber)),'.csv']])
%search for out of range values in the track and replace by middle of Arena (x,y) 
[HistTrack,HistBinCenter]=hist(Track(:,2),20);
CutOffBin=find(HistTrack==0,1);
OutOfRange=find(Track(:,2)>HistBinCenter(CutOffBin));
disp(['found ',num2str(length(OutOfRange)),' out-of-range values'])
MiddleX=(max(Track(:,1))+min(Track(:,1)))/2;
MiddleY=(max(Track(:,2))+min(Track(:,2)))/2;
for IndexOut=1:length(OutOfRange)
Track(OutOfRange(IndexOut),2)=MiddleY;
Track(OutOfRange(IndexOut),1)=MiddleX;
end
TrackOriginl=Track;
% %% TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_cells(1))/1000000; %for spikes
% %% TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_events(1))/1000000; %for laser pulses*
for Clust = ClustNumber 
%%
    if TTL==1
% Work on the timestamps from ephys
ev_indexer = false(length(Events),2);
for i=1:length(Events)
    pulses(i) = Events(i).TimeStamp;
    ev_indexer(i,1) = strcmp(Events(i).EventString(1:3),'TTL');
    ev_indexer(i,2) = logical(Events(i).TTLValue);
end

% Select only event related to TTL ONset
pulses = pulses(ev_indexer(:,1) & ev_indexer(:,2));
pulses = double(pulses)';
pulses = (pulses(1:length(pulses)))';

%% Work on biobserve timestamps
stamps = [timestamp_interval:timestamp_interval:Track(end,3)];
if length(stamps)<=length(pulses)
    pulses = pulses(1:length(stamps));
elseif length(stamps)>length(pulses)
    stamps = stamps(1:length(pulses));
end


%% Figure out how to allign the timestamps
disp(' ')
disp(' *** info about allignment of ephys and Viewer data *** ')
try
p = polyfit(pulses, stamps, 1);
catch
    pulses=pulses';
    p = polyfit(pulses, stamps, 1);
end
pulse_tester(:,1) = stamps;
pulse_tester(:,2) = pulses*p(1) + p(2);
pulse_tester(:,3) = pulse_tester(:,1) - pulse_tester(:,2);
disp(['Average error: ' num2str(mean(pulse_tester(:,3))), ' variation: ' num2str(std(pulse_tester(:,3)))])
disp(' ')


%% Change timestamps in the Tetrode file
Tetrode.timestamps = double(Tetrode.timestamps).*p(1) + p(2);

elseif  TTL==0 % this command is for recordings without TTL synchronisation
  
%find the recording time start from the BioObserve file:
[~,~,~,StartTimeBio] = biobserve_import([BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',char(ConditionList(ConditionNumber)),'.csv']); %import video file
SplitRecTime=strsplit(char(StartTimeBio),':');
SplitRecTime3=char(SplitRecTime(3));
SplitRecTime(3)={SplitRecTime3(1:2)};
StartTimeBioHH=str2double(char(SplitRecTime(1)));if StartTimeBioHH>12 ;StartTimeBioHH=StartTimeBioHH-12;end
StartTimeBioMM=str2double(char(SplitRecTime(2)));
StartTimeBioSS=strsplit(char(SplitRecTime(3)),' ');StartTimeBioSS=str2double(char(StartTimeBioSS(1)));
disp(['Bio started at: ',num2str(StartTimeBioHH),':',num2str(StartTimeBioMM),':',num2str(StartTimeBioSS)])
StartTimeBio=(StartTimeBioHH*3600+StartTimeBioMM*60+StartTimeBioSS);% time of recording start in seconds

%find the first timestamp from csc recording
ChannleNumber = TetrodeNumber*4;
[TimeStamps_csc, ~, ~, ~, ~, ~] = Nlx2MatCSC([FileFolder,'\',['CSC',num2str(ChannleNumber),'.ncs']], [1 1 1 1 1], 1, 1, [] ); %csc file

%find the recording time start from the NLX log file:
LogText= readtable([FileFolder,'\','CheetahLogFile.txt'],'HeaderLines',16);
LogText = table2cell(LogText);
count=1;
FinalArray={};
for i=1:size(LogText,1)
    TestSrt=strfind(char(LogText{i,1}),'StartRecording');
    if  ~isempty(TestSrt)
        FinalArray(count)=LogText(i);
        count=count+1;
    else
    end
end
FinalArray=FinalArray';
StartTimeLogFile=char(FinalArray(end));
StartTimeLogFileHH=str2double(StartTimeLogFile(14:15));if StartTimeLogFileHH>12 ;StartTimeLogFileHH=StartTimeLogFileHH-12;end
StartTimeLogFileMM=str2double(StartTimeLogFile(17:18));
StartTimeLogFileSS=str2double(StartTimeLogFile(20:21));
disp(['NLX started at: ',num2str(StartTimeLogFileHH),':',(StartTimeLogFile(17:18)),':',(StartTimeLogFile(20:21))])

StartTimeLogFile=(StartTimeLogFileHH*3600+StartTimeLogFileMM*60+StartTimeLogFileSS);% time of recording start in seconds

% % calculate the time difference between the Video and NLX files
TimeDiff=(double(StartTimeBio)-double(StartTimeLogFile));
TimeDiff=TimeDiff*10^6;% microseconds 

% Change timestamps new way - start here for video a
Tetrode.timestamps = (double(Tetrode.timestamps) - double(TimeStamps_csc(1))); %Make the first timestamp t0 by substracting the first timestamp in the csc file% double(StartTimeEvents)
Tetrode.timestamps = (double(Tetrode.timestamps) -TimeDiff); %reduce time difference between NLX and BioObserve 
Tetrode.timestamps=double(Tetrode.timestamps/1000000);
    end

%% Collect timestamps from one cell
cell_stamps = Tetrode.timestamps(Tetrode.cells==Clust)';


%% Find X and Y coordinates for all stamps
data = zeros(length(cell_stamps),2);

for i=1:length(cell_stamps)
    
    [~, index] = min(abs(Track(:,3)-cell_stamps(i)));
    
    data(i,:) = Track(index,1:2);
    
end


%% Density plot (firing rate)
bins = nr_bins;

bin_x = linspace(min(Track(:,1)), max(Track(:,1)), bins+1);
bin_y = linspace(min(Track(:,2)), max(Track(:,2)), bins+1);

map = zeros(bins,bins);

for i=1:bins
    for j=1:bins
        
        indexer = data(:,1)>=bin_x(i) & data(:,1)<bin_x(i+1) & data(:,2)>=bin_y(j) & data(:,2)<bin_y(j+1);
        APs = sum(indexer);
        
        indexer = Track(:,1)>=bin_x(i) & Track(:,1)<bin_x(i+1) & Track(:,2)>=bin_y(j) & Track(:,2)<bin_y(j+1);
        timespent = sum(indexer);
        
        map(j,i) = (APs/timespent)*25;
        
    end
end
% save z scored file for later 
MapZ=map;
MapZSigma=nanstd(MapZ(:));
MapZMean=nanmean(MapZ(:));
MapZ=(MapZ-MapZMean)./MapZSigma;
save([HeatMapDir,'\',char(ConditionList(ConditionNumber)),' ',Date,' T',num2str(TetrodeNumber),' C',num2str(Clust),' MapZ ',num2str(nr_bins^2),' bins']);

% mapSmooth = imgaussfilt(map);
%Calculate map parameters:
EmptySide=map(1:round(nr_bins/2),[1:round(nr_bins/4)]);
FoodSide=map([round((nr_bins/2)+1):nr_bins],[round(nr_bins*3/4+1):nr_bins]);
Statistics(1)=nanmean(map,'all'); % general average
Statistics(2)=nanmean(EmptySide,'all'); %empty average
Statistics(3)=nanmean(FoodSide,'all');  %food average
Statistics(4)=Statistics(3)/Statistics(2); %food/empty average
Statistics(5)=nanmedian(EmptySide,'all'); %empty median
Statistics(6)=nanmedian(FoodSide,'all'); %food median
Statistics(7)=Statistics(6)/Statistics(5); %food/empty median
Statistics(8)=nanstd(FoodSide,0,'all'); %stdev food
IndexFood=regexp(FoodConsumed,'\d'); FoodConsumed=FoodConsumed(IndexFood(1):IndexFood(end));
Statistics(9)=str2double(FoodConsumed); %Food Consumed
Statistics=Statistics';
Col_Title={'General Average','Empty zone Average','Food zone Average','Food/Empty Average', 'Empty zone median','Food zone median','Food/Empty median','standard diviation food side','Food consumed'}';
Identifier={[DateStamp,' ',MouseName,' T',num2str(TetrodeNumber),' C',num2str(Clust)]};

%% Make the figures
if SubplotFigures==1
    IsFigure=strfind(FileName,'Empty');
    if  ~isempty(IsFigure)
    SublpotLocation=1;
    SubplotTitle='Empty';
    else
        IsFigure=strfind(FileName,'Chow');
        if ~isempty(IsFigure)
        SublpotLocation=2;
        SubplotTitle='Chow';
        else
            IsFigure=strfind(FileName,'Exposure');
            if ~isempty(IsFigure)
            SublpotLocation=3;
            SubplotTitle='Jelly-Exposure';
            else
                IsFigure=strfind(FileName,'OFF');
                if ~isempty(IsFigure);
                SublpotLocation=4;
                SubplotTitle='Jelly-OFF';
                end
            end
        end
    end
end %if SubplotFigures

% Make the track figure and all APs
if SubplotFigures~=1
TrackFigure=figure;
SaveTrackFigureName = ['Track ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',char(ConditionList(ConditionNumber)),' ',' T',num2str(TetrodeNumber),' C',num2str(Clust)];
end
if SubplotFigures==1
   subplot(2,4,SublpotLocation)
end
plot(Track(:,1), Track(:,2))
axis equal
hold on
ax = gca;
ax.YDir = 'reverse';
if SubplotFigures==1
title(SubplotTitle);
else 
title(SaveTrackFigureName);

end

% scatter(data(:,1), data(:,2)); % plots firing over the physical distance

% % Make the heatplot figure 
 if SubplotFigures~=1
HeatMapFigure=figure;
    if TTL==1
    SaveHeatMapFigureName = ['Heatmap ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',char(ConditionList(ConditionNumber)),' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),' TTL'];
    elseif TTL==0
    SaveHeatMapFigureName = ['Heatmap ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',char(ConditionList(ConditionNumber)),' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),' NoTTL'];
    end
end
if SubplotFigures==1
subplot(2,4,(SublpotLocation+4))
end
%plot heatmap figure
imagesc(map);
gcf
caxis([RangeMin RangeMax])
if Colorbar==1
colorbar;
end
if SubplotFigures==1
title(SubplotTitle);
else
title(SaveHeatMapFigureName);
end
%%

%%% Now we save the image files in image(jpeg) and fig(MATLAB) format:
SaveFileFolder = [ComputerDir,'\',Date,'\',MouseName];
if SubplotFigures~=1
saveas(TrackFigure, fullfile(FileFolder, SaveTrackFigureName), 'jpeg'); % here you save the figure
saveas(TrackFigure, fullfile(FileFolder, SaveTrackFigureName), 'fig');% here you save the figure in MATLAB format
saveas(HeatMapFigure, fullfile(FileFolder, SaveHeatMapFigureName), 'jpeg'); % here you save the figure
saveas(HeatMapFigure, fullfile(FileFolder, SaveHeatMapFigureName), 'fig');% here you save the figure in MATLAB format

% saveas(HeatMapFigure, fullfile(SaveFileFolder, SaveHeatMapFigureName), 'jpeg'); % here you save the figure
% saveas(HeatMapFigure, fullfile(SaveFileFolder, SaveHeatMapFigureName), 'fig'); % here you save the figure
end%end SubplotFigures~=1
 
%%

% caxis([low high]); % click on the figure yet
% selection = map(30:50, 35:50); % will select row 30:50 and column 35:50
    if TTL==1
    SaveName= [num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp,' ',RecDuration,' ', Zlocation,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),' TTL2Video'];
    elseif TTL==0
    SaveName= [num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp,' ',RecDuration,' ', Zlocation,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),'_No_TTL2Video'];
    end
        if Colorbar==1
        SaveName=[SaveName,' colorbar'];
        end
% Col_Title={'Date','FileNumber','Tetrode #','cluster #','frequency','tonic frequency','spikes in bursts','intraburst frequency','interburst frequency','spikes per burst'};
% ExperimentInfo={'Experiment info ',['mouse ',MouseNumber],['Date ',DateStamp],['Time Stamp ',TimeStamp],['recording Duration ',RecDuration],['Z=',Zlocation],[TTLLaser],['Stimulus: ',StimType,' ',StimFreq],['Food Type: ',FoodType],['Food Consumed: ',FoodConsumed]}';
SheetName= char(ConditionList(ConditionNumber));
try
xlswrite(fullfile(SaveFileFolder, SaveName),map,SheetName,'A2');
xlswrite(fullfile(SaveFileFolder, SaveName),Col_Title,SheetName,'A23');     %Write column title
xlswrite(fullfile(SaveFileFolder, SaveName),Statistics,SheetName,'B23');
catch
SaveName=[SaveName,'_new'];
xlswrite(fullfile(SaveFileFolder, SaveName),map,SheetName,'A2');
xlswrite(fullfile(SaveFileFolder, SaveName),Col_Title,SheetName,'A23');     %Write column title
xlswrite(fullfile(SaveFileFolder, SaveName),Statistics,SheetName,'B23');
end %end try
end %ClusterNumber
end %TetrodeNumber
SummaryList(1,SerialNumber)=Statistics(1);
SummaryList(2,SerialNumber)=Statistics(2);
SummaryList(3,SerialNumber)=Statistics(3);
SummaryList(4,SerialNumber)=Statistics(4);
SummaryList(5,SerialNumber)=Statistics(5);
SummaryList(6,SerialNumber)=Statistics(6);
SummaryList(7,SerialNumber)=Statistics(7);
SummaryList(8,SerialNumber)=Statistics(8);
SummaryList(9,SerialNumber)=Statistics(9);
NameList{SerialNumber}=FoodType;

if SubplotFigures~=1
    close all
end

clearvars -except SubplotFigure TetrodeNumber Zlocation DateStamp SublpotLocation SubplotFigures FoodConsumed Colorbar SaveImages RangeMax RangeMin nr_bins timestamp_interval ConditionList MouseName Date FileFolder ChooseClustNumber Identifier Clust SaveStatistics TimeDiff ChooseFileNumber ChooseTetrodeNumber ClustNumber SerialNumber List NameList TTL SummaryList SaveFileFolder SaveName Col_Title;
end %FileNumber
% save subplot figure
if SubplotFigures==1
%     gcf
%     title([num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp,' T',num2str(TetrodeNumber),' C',num2str(Clust)])
    SaveName= [num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp, Zlocation,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),'_Subplotted'];
    saveas(SubplotFigure, fullfile(SaveFileFolder, SaveName), 'jpeg'); % here you save the figure
    saveas(SubplotFigure, fullfile(SaveFileFolder, SaveName), 'fig');% here you save the figure in MATLAB format

end
% save summary to excel file 
if SaveStatistics==1
    %%arrange the order of the table to be saved
    try
SummaryListTemp=zeros(size(SummaryList));
IndexEmpty = find(contains(NameList,'Empty'));
IndexChow = find(contains(NameList,'Chow'));
IndexJellyExposure = find(contains(NameList,'Exposure'));
IndexJellyOff = find(contains(NameList,'OFF'));
SummaryListTemp(:,1)=SummaryList(:,IndexEmpty);
SummaryListTemp(:,2)=SummaryList(:,IndexChow);
SummaryListTemp(:,3)=SummaryList(:,IndexJellyExposure);
SummaryListTemp(:,4)=SummaryList(:,IndexJellyOff);
for i=1:size(SummaryList,1)
SummaryListTemp(i,5)=SummaryListTemp(i,2)/SummaryListTemp(i,1);%Chow/Empty
SummaryListTemp(i,6)=SummaryListTemp(i,3)/SummaryListTemp(i,1);%Jelly Exposure/Empty
SummaryListTemp(i,7)=SummaryListTemp(i,3)/SummaryListTemp(i,2);%Jelly Exposure/Chow
SummaryListTemp(i,8)=SummaryListTemp(i,3)/SummaryListTemp(i,4);%Jelly Exposure/Jelly OFF
SummaryListTemp(i,9)=SummaryListTemp(i,4)/SummaryListTemp(i,1);%Jelly OFF/Empty end
end
NameList{1}='Empty';
NameList{2}='Chow';
NameList{3}='Jelly - Exposure';
NameList{4}='Jelly - OFF';
NameList{5}='Chow/Empty';
NameList{6}='Jelly Exposure/Empty';
NameList{7}='Jelly Exposure/Chow';
NameList{8}='Jelly Exposure/Jelly OFF';
NameList{9}='Jelly OFF/Empty';
SummaryList=SummaryListTemp;
    catch
        disp('cannot save statistics file without all conditions')
    end
xlswrite(fullfile(SaveFileFolder, SaveName),NameList,'Summary','B1');
xlswrite(fullfile(SaveFileFolder, SaveName),SummaryList,'Summary','B2');
xlswrite(fullfile(SaveFileFolder, SaveName),Col_Title,'Summary','A2');     %Write column title
xlswrite(fullfile(SaveFileFolder, SaveName),Identifier,'Summary','A1');     %Write identifier title

else 
disp('Summary file was not saved')
end
