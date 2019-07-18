%% Start
close all
clear all
% addpath(('C:\Users\admin\Documents\MATLAB\Field Trip Neuralynx'));
addpath('C:\Users\admin\Documents\MATLAB\Events_spikes_v.2.1\Field Trip Neuralynx');
%% File location
ComputerDir='D:\CheetahData\NG\Data\3-7,4-5,1-1,4-9 Backup'; %location of NLX files
BioObserveDir='E:\Single unit videos\Compressed';%location of Bioobserve video files
%% Variables
Date='2019-04-30';
MouseName='SUBLAT1-1'; %%% specify mouse name
ChooseTetrodeNumber=[5]; % what Tetrodes to run [1 2 3 4 5 6 7 8]
% StartFile = 2;
% EndFile =5;
% for FileToRun=StartFile:EndFile
%     if FileToRun==2
% Condition='Empty';
% ChooseFileNumber = [2]; 
%  elseif FileToRun==3
% Condition='Chow';
% ChooseFileNumber = [3];
%  elseif FileToRun==4
% Condition='Jelly_Exposure';
% ChooseFileNumber = [4];
%  elseif FileToRun==5
% Condition='Jelly_OFF';
% ChooseFileNumber = [5];
% ChooseFileNumber=FileToRun;
% FileNameList=extractfield(FileDir,'name');char(FileNameList(FileNumber))

% other variables
timestamp_interval = 5; %Sec
%arena = [200 -290 160 80];
nr_bins = 20; %Actual the square root of the number of bins
RangeMin = 0; % minimal value for heatmap
RangeMax = 180;% maximal value for heatmap
%%
FileDirNLX = dir([ComputerDir,'\',Date,'\',MouseName]);
FileDirBioObserve = dir([BioObserveDir,'\',Date,'\',MouseName]);
ChooseFileNumber=4;
Condition='Jelly_Exposure';
for FileNumber = ChooseFileNumber+2
    VideoFileNumber=FileNumber-1;
        FileFolder = [ComputerDir,'\',Date,'\',MouseName,'\',FileDirNLX(FileNumber).name];
        FileName=FileDirNLX(FileNumber).name % display file name
        NameTest=strfind(FileName,'Empty');
    If strfind(FileName,'Empty')>0;
    Condition
%now we split the name and date componantes:
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
try StimType=char(SplitName(7));catch; StimType = ''; end
try StimFreq=char(SplitName(8));catch; StimFreq = ''; end
try FoodType=char(SplitName(9));catch; FoodType = ''; end
try FoodConsumed=[char(SplitName(10)),' gr'];catch; FoodConsumed = ''; end

%% Import the files % write down which tetrodes to run in the loop

for TetrodeNumber = ChooseTetrodeNumber
    try
    %%% specifies filename: the correct filename format is "TT1_s.ntt","TT2_s.ntt" atc.
Tetrode = tetrode([FileFolder,'\','TT',num2str(TetrodeNumber),'_s.ntt']);
Events = read_neuralynx_nev([FileFolder,'\','Events.nev']);
ClustNumber=[0 1 2 3 4 5 6 7 8];  
    catch
        try
Tetrode = tetrode([FileFolder,'\','TT',num2str(TetrodeNumber),'.ntt']);
Events = read_neuralynx_nev([FileFolder,'\','Events.nev']);
ClustNumber=0;          
        catch
            disp('Cannot load file');
        continue
        end
    end
Track = biobserve_import([BioObserveDir,'\',Date,'\',MouseName,'\',Date(1:4),Date(6:7),Date(9:10),'_',MouseName(1:7),'_',MouseName(9),'_',Condition,'.csv']); %import video file

%%
% %% TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_cells(1))/1000000; %for spikes
% %% TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_events(1))/1000000; %for laser pulses*
%% write down which clusters to run in the loop
for Clust = ClustNumber 
     




%% Work on the timestamps from ephys
ev_indexer = false(length(Events),2);
for i=1:length(Events)
    pulses(i) = Events(i).TimeStamp;
    ev_indexer(i,1) = strcmp(Events(i).EventString(1:3),'TTL');
    ev_indexer(i,2) = logical(Events(i).TTLValue);
end

% Select only event related to TTL ONset
pulses = pulses(ev_indexer(:,1) & ev_indexer(:,2));
pulses = double(pulses)';



%% Work on biobserve timestamps
stamps = [timestamp_interval:timestamp_interval:Track(end,3)];
stamps = stamps(1:length(pulses))';


%% Figure out how to allign the timestamps
disp(' ')
disp(' *** info about allignment of ephys and Viewer data *** ')
p = polyfit(pulses, stamps, 1)

pulse_tester(:,1) = stamps;
pulse_tester(:,2) = pulses*p(1) + p(2);
pulse_tester(:,3) = pulse_tester(:,1) - pulse_tester(:,2);
disp(['Average error: ' num2str(mean(pulse_tester(:,3))), ' variation: ' num2str(std(pulse_tester(:,3)))])
disp(' ')


%% Change timestamps in the Tetrode file
Tetrode.timestamps = double(Tetrode.timestamps).*p(1) + p(2);


%% Collect timestamps from one cell
cell_stamps = Tetrode.timestamps(Tetrode.cells==ClustNumber)';


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


%% Make the figures
% First the track and all APs
TrackFigure=figure;
% subplot(2,1,1)
plot(Track(:,1), Track(:,2))
axis equal
hold on
ax = gca;
ax.YDir = 'reverse';
% scatter(data(:,1), data(:,2)); % plots firing over the physical distance

% Then the heatplot
HeatMapFigure=figure;
% subplot(2,1,2)
imagesc(map);
gcf
caxis([RangeMin RangeMax])
%%% Now we save the image files in image(jpeg) and fig(MATLAB) format:
% Save Track image
SaveTrackFigureName = ['Track ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',Condition,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust)];
title(SaveTrackFigureName);
saveas(TrackFigure, fullfile([ComputerDir,'\',Date,'\',MouseName], SaveTrackFigureName), 'jpeg'); % here you save the figure
saveas(TrackFigure, fullfile([ComputerDir,'\',Date,'\',MouseName], SaveTrackFigureName), 'fig');% here you save the figure in MATLAB format

SaveHeatMapFigureName = ['Heatmap ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',Condition,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust)];
title(SaveHeatMapFigureName);
saveas(HeatMapFigure, fullfile([ComputerDir,'\',Date,'\',MouseName], SaveHeatMapFigureName), 'jpeg'); % here you save the figure
saveas(HeatMapFigure, fullfile([ComputerDir,'\',Date,'\',MouseName], SaveHeatMapFigureName), 'fig');% here you save the figure in MATLAB format



%% This code wil not run
% caxis([low high]); % click on the figure yet
% selection = map(30:50, 35:50); % will select row 30:50 and column 35:50



close all

end %ClusterNumber
end %TetrodeNumber
end %FileNumber
% end