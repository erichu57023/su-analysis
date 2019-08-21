close all; clearvars -except UseArena VideoFileOrder threshold VideoCount TrialCut DeafaultParams ClustersForFiringAnalysis threshold FinalNumber ConditionCount FileNumber TimeOfInterest RunPipeline CalculateFiringParameters CalculateFiringParameters RangeMin RangeMax Gaussian_kernel ClusteredFile HanSort FileNameList ComputerDir ConditionFileOrder ConditionList Date FileDirNLX MouseName TetrodeNumber
%% optional code changes
ClustNumber = [1];
RangeMax = 10;  % caxis([0 140])
% for Neta=[1 2 3 4 5 6];
% ClustNumber = Neta
Clust = ClustNumber; %Tagged 1 4 5 NonTagged 2 3 6 7
for CodeChanges=1:1
UseArena=true;
DeafaultParams=true;
TTL=true; 
SaveStatistics=true;
ConcatenateFiles=true; % indicate 1 for YES 0 for NO
SaveImages=true;%indicate 'YES'(1) / 'NO'(0)
Colorbar=false;%indicate 'YES'(1) / 'NO'(0)
SubplotFigures=true;%indicate 'YES'(1) / 'NO'(0)
CalculateFiringParameters=true;
PlotSmoothImages=false;
SplitFiles=true; % indicate 1 for YES 0 for NO
PlotFigures=true; % indicate if should plot average figures 0 for no 1 for Yes
ClusteredFile=true;
RunPipeline=true;
end %for CodeChanges
%% Variables
for Variables=1:1
% Left=138; Top=208;Width=165;Height=89;
if ~DeafaultParams
TetrodeNumber=[1]; 
Date='2019-07-16';
MouseName='SUBLAT1-1'; %%% specify mouse name
end %if ~DeafaultParams
ComputerDir='H:\'; %location of NLX files
RangeMin = 0; % minimal value for heatmap
Gaussian_kernel=true;
Sigma = 2; % Sigma of the kernel e^(-1*data/2*Sigma^2)
timestamp_interval = 5; %Sec
nr_bins = 30; %Actual the square root of the number of bins
TimeOfInterest=600; % time to calculate firing parameters
ClustersForFiringAnalysis=ClustNumber;
%arena = [200 -290 160 80];
% ConditionList= {'Laser','Laser','Empty','Chow','Jelly_Exposure','Jelly_OFF'};
end % for Variables
for Path=1:1
%% add path to NLX import tools
addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\Field Trip Neuralynx'));
addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\NeuralynxMatlabImportExport_v6.0.0'));
FileDirNLX = dir([ComputerDir,'\',Date,'\',MouseName]);
BioObserveDir='G:\Single unit videos\Compressed';%location of Bioobserve video files my computer
FileDirVideo = dir([BioObserveDir,'\',Date,'\',MouseName]);
HeatMapDir=[ComputerDir,'Maps'];
end %for Path
for Pipline=1:1
if RunPipeline
ConditionList= {};%'Laser','Empty','Chow','Jelly_Exposure','Jelly_OFF'};
FileNameList={};LaserSerial=0;
for FileNumber = [3:size(FileDirNLX,1)]
FileNameList(FileNumber-2) = {FileDirNLX(FileNumber).name};
if contains(char(FileNameList(FileNumber-2)),'Empty')
    ConditionFileOrder(1)=(FileNumber);
    ConditionList(FileNumber)={'Empty'};
end
if contains(char(FileNameList(FileNumber-2)),'Chow')
    ConditionFileOrder(2)=(FileNumber);
    ConditionList(FileNumber)={'Chow'};
end
if contains(char(FileNameList(FileNumber-2)),'Exposure')
    ConditionFileOrder(3)=(FileNumber);
    ConditionList(FileNumber)={'Jelly_Exposure'};
end
if contains(char(FileNameList(FileNumber-2)),'OFF')
    ConditionFileOrder(4)=(FileNumber);
    ConditionList(FileNumber)={'Jelly_OFF'};
end
if contains(char(FileNameList(FileNumber-2)),'Laser')
    LaserSerial=LaserSerial+1;
    ConditionFileOrder(4+LaserSerial)=(FileNumber);
    ConditionList(FileNumber)={'Laser'};
end
if contains(char(FileNameList(FileNumber-2)),'Concatenated')
    LaserSerial=LaserSerial+1;
    ConditionFileOrder(4+LaserSerial+1)=(FileNumber);
    ConditionList(FileNumber)={'Concatenated'};
end
end %for
ConditionList=ConditionList(3:end);
%% orgenize video recordings to the order: 'Empty','Chow',Jelly_Exposure','Jelly_OFF', 'Laser 1hz','Laser 20hz'..
% find file order in Video folder
VideoFileOrder={};
for csvFiles=1:length(FileDirVideo)
if contains(char(FileDirVideo(csvFiles).name),'csv');
if contains(char(FileDirVideo(csvFiles).name),'Empty');
VideoFileOrder(1)={char(FileDirVideo(csvFiles).name)};
elseif contains(char(FileDirVideo(csvFiles).name),'Chow');
VideoFileOrder(2)={char(FileDirVideo(csvFiles).name)};
elseif contains(char(FileDirVideo(csvFiles).name),'Exposure');
VideoFileOrder(3)={char(FileDirVideo(csvFiles).name)};
elseif contains(char(FileDirVideo(csvFiles).name),'OFF');
VideoFileOrder(4)={char(FileDirVideo(csvFiles).name)};
end
end
end % find video file names
end %if run pipeline
end % for Pipline
%% Start code
VideoCount=1;ConditionCount=0; FinalNumber=1;SummaryList=zeros(9,9);NameList={}; if SubplotFigures; SubplotFigure=figure('units','normalized','outerposition',[0 0 1 1]);end
for FileNumber = ConditionFileOrder(1:4)% Build heat maps 
DataFileName=FileDirNLX(FileNumber).name % display file name
DataFileFolder = [ComputerDir,'\',Date,'\',MouseName,'\',FileDirNLX(FileNumber).name];
BioObserveFile=[char(VideoFileOrder(VideoCount))];
%% Split the name and date componantes:
for SplitFileName=1:1
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
end % for SplitFileName
%% Import the files and create Obj with "Tetrode"
for CreateObj=1:1
    SerialNumber=FileNumber-(ConditionFileOrder(1)-3);
Tetrode = tetrode([DataFileFolder,'\',['TT',num2str(TetrodeNumber),'_Cut',num2str(threshold(1)-5),'.ntt']],ClusteredFile);
Events = read_neuralynx_nev([DataFileFolder,'\','Events.nev']);
Track = biobserve_import([BioObserveDir,'\',Date,'\',MouseName,'\',char(VideoFileOrder(VideoCount))]); %import video file
disp(['BioObserve file ', [BioObserveDir,'\',Date,'\',MouseName,'\',char(VideoFileOrder(VideoCount))]])
VideoCount=VideoCount+1;
% fix file names containing '-'
% % % if contains(char(ConditionList(FileNumber-2)),' - ')
% % % NewCondition = strrep(char(ConditionList(FileNumber-2)),' - ','_');
% % % Track = biobserve_import([BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',NewCondition,'.csv']); %import video file
% % % disp(['BioObserve file', [BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',NewCondition,'.csv']])
% % % BioObserveFile=[Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',NewCondition,'.csv'];
% % % else
% % % Track = biobserve_import([BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',char(ConditionList(FileNumber-2)),'.csv']); %import video file
% % % disp(['BioObserve file', [BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',char(ConditionList(FileNumber-2)),'.csv']])
% % % BioObserveFile=[Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',char(ConditionList(FileNumber-2)),'.csv'];
% % % end
% figure
% Tetrode.show_cell(Clust,'average')
end % for CreateObj
%% search for out of range values in the track and replace by middle of Arena (x,y) 
for OutValues=1:1
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
end % for OutValues
TrackOriginl=Track;
%% Sync with Video (2 methods)
for VideoSync=1:1 
%% Use this code if there was TTL synchronisation with video
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
% Work on biobserve timestamps
stamps = [timestamp_interval:timestamp_interval:Track(end,3)];
if length(stamps)<=length(pulses)
    pulses = pulses(1:length(stamps));
elseif length(stamps)>length(pulses)
    stamps = stamps(1:length(pulses));
end
% Figure out how to allign the timestamps
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
% Change timestamps in the Tetrode file
Tetrode.timestamps = double(Tetrode.timestamps).*p(1) + p(2);
%% Use this code if there was no TTL synchronisation with video
elseif  TTL==0 % this command is for recordings without TTL synchronisation
%find the recording time start from the BioObserve file:
[~,~,~,StartTimeBio] = biobserve_import([BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',char(ConditionList(FileNumber-2)),'.csv']); %import video file
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
[TimeStamps_csc, ~, ~, ~, ~, ~] = Nlx2MatCSC([DataFileFolder,'\',['CSC',num2str(ChannleNumber),'.ncs']], [1 1 1 1 1], 1, 1, [] ); %csc file
%find the recording time start from the NLX log file:
LogText= readtable([DataFileFolder,'\','CheetahLogFile.txt'],'HeaderLines',16);
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
% calculate the time difference between the Video and NLX files
TimeDiff=(double(StartTimeBio)-double(StartTimeLogFile));
TimeDiff=TimeDiff*10^6;% microseconds 
% Change timestamps new way - start here for video a
Tetrode.timestamps = (double(Tetrode.timestamps) - double(TimeStamps_csc(1))); %Make the first timestamp t0 by substracting the first timestamp in the csc file% double(StartTimeEvents)
Tetrode.timestamps = (double(Tetrode.timestamps) -TimeDiff); %reduce time difference between NLX and BioObserve 
Tetrode.timestamps=double(Tetrode.timestamps/1000000);
    end
end %for VideoSync 
%% Create and plot the Heatmaps
for CreateMap=1:1 
% Collect timestamps from one cell
cell_stamps = Tetrode.timestamps(Tetrode.cells==Clust)';
% Find X and Y coordinates for all stamps
data = zeros(length(cell_stamps),2);
for i=1:length(cell_stamps)
    
    [~, index] = min(abs(Track(:,3)-cell_stamps(i)));
    
    data(i,:) = Track(index,1:2);
end
% Density plot (firing rate)
if UseArena
try
ArenaInfo=table2struct(readtable([ComputerDir,Date,'\',MouseName,'\Arena.xlsx']));
catch
ArenaInfo=table2struct(readtable([ComputerDir,Date,'\',MouseName,'\',FileDirNLX(FileNumber).name,'\Arena.xlsx']));
end
Xmin = ArenaInfo(1).Val; % (Left)
Xmax = ArenaInfo(1).Val+ArenaInfo(3).Val; %(Left+Width)
Ymin = ArenaInfo(2).Val; %Top
Ymax = ArenaInfo(2).Val+ArenaInfo(4).Val; %Top+Height
% bin_x = linspace((Left), (Left+Width), bins+1);
% bin_y = linspace((Top), (Top+Height), bins+1);
disp('used Arena Info')
else
Xmin = min(Track(:,1));
Xmax = max(Track(:,1));
Ymin = min(Track(:,2));
Ymax = max(Track(:,2));
end 
bins = nr_bins;
bin_x = linspace(Xmin, Xmax, bins+1);
bin_y = linspace(Ymin, Ymax, bins+1);
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
save([ComputerDir,'\Maps','\',[FoodType,' ',MouseName,' ',Date,' T',num2str(TetrodeNumber),' C',num2str(Clust),' Q',num2str(Sigma),' Bin',num2str(nr_bins)]],'map')

% Add smooth image
if Gaussian_kernel
    nr_binsSmooth = 50;
    MapSmooth=zeros(nr_binsSmooth,nr_binsSmooth);
    bin_xSmooth = linspace(Xmin, Xmax, nr_binsSmooth+1);
    bin_ySmooth = linspace(Ymin, Ymax, nr_binsSmooth+1);
    % Change the bins to centerpoints instead
    bin_width_XSmooth = (Xmax - Xmin)/nr_binsSmooth;
    bin_width_YSmooth = (Ymax - Ymin)/nr_binsSmooth;
    bin_xSmooth = bin_xSmooth(1:end-1)+0.5*bin_width_XSmooth;
    bin_ySmooth = bin_ySmooth(1:end-1)+0.5*bin_width_YSmooth;
for i=1:nr_binsSmooth
    for j=1:nr_binsSmooth
            % First calculate the squared Euclidian distance to every
            % action potential
            AP_distanceSmooth = (bin_xSmooth(i)-data(:,1)).^2+(bin_ySmooth(j)-data(:,2)).^2;
            % Calculate the same for the time spend close to this pixel
            t_distanceSmooth = (bin_xSmooth(i)-Track(:,1)).^2+(bin_ySmooth(j)-Track(:,2)).^2;
            % Apply a Gaussian function
            AP_distanceSmooth = exp(-1.*AP_distanceSmooth/(2*Sigma^2));
            t_distanceSmooth = exp(-1.*t_distanceSmooth/(2*Sigma^2));
            % Divide the two to get a weighted firing frequency?
            MapSmooth(j,i) = sum(AP_distanceSmooth)/sum(t_distanceSmooth);
    end
end
mkdir([ComputerDir,'\Images to avearge'],'images')
save([ComputerDir,'\Images to avearge\images','\',[FoodType,' ',MouseName,' ',Date,' T',num2str(TetrodeNumber),' C',num2str(Clust),' Q',num2str(Sigma),' Bin',num2str(nr_binsSmooth)]],'MapSmooth')
end % if Gaussian

% Make the figures
if SubplotFigures
    IsFigure=strfind(DataFileName,'Empty');
    if  ~isempty(IsFigure)
    SublpotLocation=1;
    SubplotTitle='Empty';
    else
        IsFigure=strfind(DataFileName,'Chow');
        if ~isempty(IsFigure)
        SublpotLocation=2;
        SubplotTitle='Chow';
        else
            IsFigure=strfind(DataFileName,'Exposure');
            if ~isempty(IsFigure)
            SublpotLocation=3;
            SubplotTitle='Jelly-Exposure';
            else
                IsFigure=strfind(DataFileName,'OFF');
                if ~isempty(IsFigure)
                SublpotLocation=4;
                SubplotTitle='Jelly-OFF';
                end
            end
        end
    end
end %if SubplotFigures
% Make the track figure and all APs
if ~SubplotFigures
TrackFigure=figure;
SaveTrackFigureName = ['Track ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',char(ConditionList(FileNumber-2)),' ',' T',num2str(TetrodeNumber),' C',num2str(Clust)];
end
if SubplotFigures
   subplot(2,4,SublpotLocation)
end
plot(Track(:,1), Track(:,2))
% axis equal
xlim([Xmin Xmax]);
ylim([Ymin Ymax]);
hold on
ax = gca;
ax.YDir = 'reverse';
if SubplotFigures
title(SubplotTitle);
else 
title(SaveTrackFigureName);
end
% scatter(data(:,1), data(:,2)); % plots firing over the physical distance

% % Make the heatplot figure 
 if ~SubplotFigures
HeatMapFigure=figure;
    if TTL
    SaveHeatMapFigureName = ['Heatmap ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',char(ConditionList(FileNumber-2)),' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),' TTL'];
    elseif TTL==0
    SaveHeatMapFigureName = ['Heatmap ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',char(ConditionList(FileNumber-2)),' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),' NoTTL'];
    end
end
if SubplotFigures
subplot(2,4,(SublpotLocation+4))
end
%plot heatmap figure
if PlotSmoothImages
    imagesc(MapSmooth);
else
    imagesc(map);
end
gcf
caxis([RangeMin RangeMax])
if Colorbar
colorbar;
end
if SubplotFigures
title(SubplotTitle);
else
title(SaveHeatMapFigureName);
end
%%

%%% Now we save the image files in image(jpeg) and fig(MATLAB) format:
SaveFileFolder = [ComputerDir,'\',Date,'\',MouseName];
if ~SubplotFigures
saveas(TrackFigure, fullfile(DataFileFolder, SaveTrackFigureName), 'jpeg'); % here you save the figure
saveas(TrackFigure, fullfile(DataFileFolder, SaveTrackFigureName), 'fig');% here you save the figure in MATLAB format
saveas(HeatMapFigure, fullfile(DataFileFolder, SaveHeatMapFigureName), 'jpeg'); % here you save the figure
saveas(HeatMapFigure, fullfile(DataFileFolder, SaveHeatMapFigureName), 'fig');% here you save the figure in MATLAB format

% saveas(HeatMapFigure, fullfile(SaveFileFolder, SaveHeatMapFigureName), 'jpeg'); % here you save the figure
% saveas(HeatMapFigure, fullfile(SaveFileFolder, SaveHeatMapFigureName), 'fig'); % here you save the figure
end%end ~SubplotFigures
 
%%

% caxis([low high]); % click on the figure yet
% selection = map(30:50, 35:50); % will select row 30:50 and column 35:50
    if TTL
    SaveName= [num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp,' ',RecDuration,' ', Zlocation,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),' TTL2Video'];
    else
    SaveName= [num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp,' ',RecDuration,' ', Zlocation,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),'_No_TTL2Video'];
    end
        if Colorbar
        SaveName=[SaveName,' colorbar'];
        end
end % for CreateMap
 %% Calculate map parameters:
for MapParams=1:1
EmptySide=map([1:round(nr_bins/2)+1],[1:round(nr_bins*0.25)+1]); % take line 1:11, column 1:6
FoodSide=map([round(nr_bins*0.5):nr_bins],[round(nr_bins*0.75):nr_bins]); % take line 10:20, column 15:20
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
%% Save statistics
% Col_Title={'Date','FileNumber','Tetrode #','cluster #','frequency','tonic frequency','spikes in bursts','intraburst frequency','interburst frequency','spikes per burst'};
% ExperimentInfo={'Experiment info ',['mouse ',MouseNumber],['Date ',DateStamp],['Time Stamp ',TimeStamp],['recording Duration ',RecDuration],['Z=',Zlocation],[TTLLaser],['Stimulus: ',StimType,' ',StimFreq],['Food Type: ',FoodType],['Food Consumed: ',FoodConsumed]}';
SheetName= char(ConditionList(FileNumber-2));
try
xlswrite(fullfile(SaveFileFolder, SaveName),map,SheetName,'A2');
xlswrite(fullfile(SaveFileFolder, SaveName),Col_Title,SheetName,'A23');     %Write column title
xlswrite(fullfile(SaveFileFolder, SaveName),Statistics,SheetName,'B23');
xlswrite(fullfile(SaveFileFolder, SaveName),'something',SheetName,'A1');

catch
SaveName=[SaveName,'_new'];
xlswrite(fullfile(SaveFileFolder, SaveName),map,SheetName,'A2');
xlswrite(fullfile(SaveFileFolder, SaveName),Col_Title,SheetName,'A23');     %Write column title
xlswrite(fullfile(SaveFileFolder, SaveName),Statistics,SheetName,'B23');
xlswrite(fullfile(SaveFileFolder, SaveName),{char(BioObserveFile)},SheetName,'A1');
end %end try
XLSaveName=SaveName;
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
end %for MapParams
if ~SubplotFigures; close all;end
if CalculateFiringParameters
ConditionCount=ConditionCount+1;
ConditionList(ConditionCount+1)
% [TimeStamps_spike, ScNumbers_spike, ~, Features_spike, Samples_spike, Header_spike] = Nlx2MatSpike([DataFileFolder,'\',['TT',num2str(TetrodeNumber),'.ntt']], [1 1 1 1 1], 1, 1, [] ); %spike file
% [~, ~, CellNumbers_spike, ~, ~, ~] = Nlx2MatSpike([DataFileFolder,'\',['TT',num2str(TetrodeNumber),'_s.ntt']], [1 1 1 1 1], 1, 1, [] ); %spike file
[TimeStamps_spike, ScNumbers_spike, CellNumbers_spike, Features_spike, Samples_spike, Header_spike_Firing] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(FileNumber).name),'\','TT',num2str(TetrodeNumber),'_Cut',num2str(threshold(1)-5),'.ntt'], [1 1 1 1 1], 1, 1, [] ); %spike file
[TimeStamps_csc, ChannelNumbers_csc, SampleFrequencies_csc, NumberOfValidSamples_csc, Samples_csc, Header_csc] = Nlx2MatCSC([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(FileNumber).name),'\','CSC',num2str(TetrodeNumber*4),'.ncs'], [1 1 1 1 1], 1, 1, [] ); %csc file
for ClusterNumber=ClustersForFiringAnalysis;
try
TimeStamps_csc_zeroed_s = (TimeStamps_csc-TimeStamps_csc(1))/1000000;
cell = cat(1,TimeStamps_spike,CellNumbers_spike); 
cell_index = find(cell(2,:) ~=ClusterNumber); % "~=1" looks at cell#1, "~=2" looks at cell#2, etc.
cell(:,cell_index) = [];
[c,indexFiring] = min(abs(TimeStamps_csc_zeroed_s-TimeOfInterest));
cell_index_time = find(cell(1,:) > TimeStamps_csc(indexFiring)); %finds the indices higher than 120s
cell(:,cell_index_time) = [];
cell = cell/1000000;
cell = cell(1,:) - cell(1,1);
cell(2,:) = 1;
%% runs burst detection algorithm
tn = cell(1,:);
limit = 0.1;
RSalpha = -log(0.05);
[archive_burst_RS, archive_burst_length, archive_burst_start] = burst(tn,limit,RSalpha);
% creates an array ("burst_all_final") with burst indices
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
% plots indentified bursts
burst_time_s = zeros(1,length(burst_all_final));
for i = 1:length(burst_all_final)
    burst_time_s(i) = cell(1,burst_all_final(i));
end
burst_time_s(2,:) = 2;
% stem(cell(1,:),cell(2,:))
% hold on
% stem(burst_time_s(1,:),burst_time_s(2,:))    
% firing frequency
% TimeOfInterest = 300; %!!!!!!!!!!!!!!!!!!!!!!!delete this if you are looking at only one 2 min recording!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
frequency = length(cell)/TimeOfInterest;
length(cell);
% tonic firing frequency
burst_final_temp = burst_time_s(1,:);
cell_temp = cell(1,:);
[C,ia,ib] = intersect(cell_temp,burst_time_s);
cell_temp(ia) = [];
tonic_frequency = length(cell_temp)/TimeOfInterest;
% spikes_in_bursts
spikes_in_bursts = length(burst_time_s(1,:))/length(cell(1,:))*100;
if spikes_in_bursts~=0
% spikes_per_burst
spikes_per_burst = mean(archive_burst_length);
%% interburst frequency
interburst_frequency = length(archive_burst_start) / (cell(1,archive_burst_start(end)) - cell(1,archive_burst_start(1)));
%% intraburst frequency
intra_temp = zeros(1,length(archive_burst_start));
for i = 1:length(archive_burst_start)
    intra_temp(i) = archive_burst_length(i) / (cell(1,burst_end(i)) - cell(1,archive_burst_start(i)));
end
a = find (intra_temp >= 1000);
intra_temp(a) = [];
intraburst_frequency = mean(intra_temp);
else 
spikes_per_burst =0;
interburst_frequency=0;
intraburst_frequency=0;
end
%% display
U = ['Firing frequency ',num2str(frequency)];
V = ['Tonic frequency ',num2str(tonic_frequency)];
W = ['% of spikes in bursts ',num2str(spikes_in_bursts)];
X = ['Intraburst frequency ',num2str(intraburst_frequency)];
Y = ['Interburst frequency ',num2str(interburst_frequency)];
Z = ['Spikes per burst ',num2str(spikes_per_burst)];
disp(U);disp(V);disp(W);disp(X);disp(Y);
disp(Z)
Final(FinalNumber).Date=Date;
Final(FinalNumber).MouseName=MouseName;
Final(FinalNumber).TetrodeNumber=TetrodeNumber;
Final(FinalNumber).ClusterNumber=ClusterNumber;
Final(FinalNumber).frequency=frequency;
Final(FinalNumber).tonic_frequency=tonic_frequency;
Final(FinalNumber).spikes_in_bursts=spikes_in_bursts;
Final(FinalNumber).intraburst_frequency=intraburst_frequency;
Final(FinalNumber).spikes_per_burst=spikes_per_burst;
Final(FinalNumber).MedianEmpty=nanmedian(EmptySide,'all'); %empty median
Final(FinalNumber).Medianfood=nanmedian(FoodSide,'all'); %Food median
Final(FinalNumber).MeanEmptySide=nanmean(EmptySide,'all'); %empty median
Final(FinalNumber).MeanfoodSide=nanmean(FoodSide,'all'); %Food median
Final(FinalNumber).Food_Type=FoodType;
Final(FinalNumber).Food_Consumed=FoodConsumed;
Final(FinalNumber).Time_Of_Interest_Sec=TimeOfInterest;% save time analysed (sec)  
Final(FinalNumber).Zlocation=Zlocation;% save Z location 
Final(FinalNumber).FileName=char(FileDirNLX(FileNumber).name); % save file name
FinalNumber=FinalNumber+1;
catch
continue
end %try
end % for culsternumber
end %if CalculateFiringParameters
clearvars -except TrialCut UseArena DeafaultParams VideoCount VideoFileOrder ClustersForFiringAnalysis threshold CalculateFiringParameters Final Final FinalNumber ConditionCount TimeOfInterest CalculateFiringParameters RangeMin RangeMax PlotSmoothImages Sigma Gaussian_kernel ClusteredFile XLSaveName suptitle HeatMapDir FileNameList ComputerDir ConditionFileOrder ConditionList Date FileDirNLX MouseName TetrodeNumber BioObserveDir SubplotFigure TetrodeNumber Zlocation DateStamp SublpotLocation SubplotFigures FoodConsumed Colorbar SaveImages RangeMax RangeMin nr_bins timestamp_interval ConditionList MouseName Date FileFolder ChooseClustNumber Identifier Clust SaveStatistics TimeDiff ChooseFileNumber ChooseTetrodeNumber ClustNumber SerialNumber List NameList TTL SummaryList SaveFileFolder SaveName Col_Title;
end %for FileNumber
if CalculateFiringParameters
try writetable(struct2table(Final), [ComputerDir,Date,'\',MouseName,'\ConcatenatedFile\','Final_TT',num2str(TetrodeNumber),' Cluster ',num2str(ClustNumber),'.xlsx'])
catch writetable(struct2table(Final), [ComputerDir,Date,'\',MouseName,'\ConcatenatedFile\','Final_TT',num2str(TetrodeNumber),'_new.xlsx'])
end
end %if CalculateFiringParameters
for SaveSubplot=1:1
% save subplot figure
if SubplotFigures
%     suptitle('something');
%     gcf
%     title([num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp,' T',num2str(TetrodeNumber),' C',num2str(Clust)])
    SaveName= [num2str(nr_bins^2),'bins ',MouseName,' ',DateStamp, Zlocation,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust),'_Subplotted'];
    saveas(SubplotFigure, fullfile(SaveFileFolder, SaveName), 'jpeg'); % here you save the figure
    saveas(SubplotFigure, fullfile(SaveFileFolder, SaveName), 'fig');% here you save the figure in MATLAB format
end
end % for SaveSubplot
% save summary to excel file 
for SaveStat=1:1
if SaveStatistics
    %%arrange the order of the table to be saved
%     try
SummaryListTemp=zeros(size(SummaryList));
SummaryListTemp(:,1)=SummaryList(:,ConditionFileOrder(1)-1);
SummaryListTemp(:,2)=SummaryList(:,ConditionFileOrder(2)-1);
SummaryListTemp(:,3)=SummaryList(:,ConditionFileOrder(3)-1);
SummaryListTemp(:,4)=SummaryList(:,ConditionFileOrder(4)-1);
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
%     catch
%         disp('cannot save statistics file without all conditions')
%     end
xlswrite(fullfile(SaveFileFolder, XLSaveName),NameList,'Summary','B1');
xlswrite(fullfile(SaveFileFolder, XLSaveName),SummaryList,'Summary','B2');
xlswrite(fullfile(SaveFileFolder, XLSaveName),Col_Title,'Summary','A2');     %Write column title
xlswrite(fullfile(SaveFileFolder, XLSaveName),Identifier,'Summary','A1');     %Write identifier title
else 
disp('Summary file was not saved')
end
end % for SaveStat
clearvars -except UseArena threshold Tetrode VideoCount TrialCut VideoFileOrder DeafaultParams ClustersForFiringAnalysis FileNumber Final TimeOfInterest CalculateFiringParameters PlotSmoothImages Gaussian_kernel Sigma ClusteredFile SublpotLocation ClustNumber DeafaultParams FileNameList ComputerDir ConditionFileOrder ConditionList Date FileDirNLX MouseName TetrodeNumber
disp('if it comes out empty, change ClusteredFile=true in tetrode file')
% close all
% end % for Neta

