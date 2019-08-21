%%%This code is designed to find sti`mulation responsive cells and produce figures for visual assessment of the cell response. 
%%% Please fill in the parameters:
close all
% clc
addpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\NeuralynxMatlabImportExport_v6.0.0');
close all; clearvars -except TrialCut DeafaultParams ClustersForFiringAnalysis threshold FinalNumber ConditionCount FileNumber TimeOfInterest RunPipeline CalculateFiringParameters CalculateFiringParameters RangeMin RangeMax Gaussian_kernel ClusteredFile HanSort FileNameList ComputerDir ConditionFileOrder ConditionList Date FileDirNLX MouseName TetrodeNumber
DeafaultParams=true;
Hz=2;
if ~DeafaultParams
clear all
ComputerDir='H:\'; %location of NLX files
Date='2019-05-02';
MouseName='SUKDLAT1-7'; %%% specify mouse name
FileDirNLX = dir([ComputerDir,'\',Date,'\',MouseName]);
ChooseTetrodeNumber=[3]; % what Tetrodes to run [1 2 3 4 5 6 7 8]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
ChooseTetrodeNumber=TetrodeNumber;
end
ClustNumber=[0:10];
%% Choose folders and tetrodes to run and graph features
ChooseFileNumber = [1]; % what files to run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Hz==20
ChooseMaxLimitForHistogram =0;% When this is set to 0 the code cooses the value. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChooseLimForHistogram = 0.02;%%%%%%%%%%%%%%%%%%%%%%%
ChooseEdgesForHistogram=[-ChooseLimForHistogram:0.0004:ChooseLimForHistogram]; %bin size is the middle value round((ChooseLimForHistogram/50),4)
ChooseTimeScaleForRaster=0.02;%%%%%%%%%%%%%%%%%%%%%%%
ChooseTickSizeForRaster=(ChooseTimeScaleForRaster/5); % deafault should beChooseTickSizeForRaster=(ChooseTimeScaleForRaster/5)
% FileNameList=extractfield(FileDir,'name');char(FileNameList(FileNumber))
else
ChooseMaxLimitForHistogram =50;% When this is set to 0 the code cooses the value. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChooseLimForHistogram = 0.2;%%%%%%%%%%%%%%%%%%%%%%%
ChooseEdgesForHistogram=[-ChooseLimForHistogram:0.004:ChooseLimForHistogram]; %bin size is the middle value round((ChooseLimForHistogram/50),4)
ChooseTimeScaleForRaster=0.2;%%%%%%%%%%%%%%%%%%%%%%%
ChooseTickSizeForRaster=(ChooseTimeScaleForRaster/5); % deafault should beChooseTickSizeForRaster=(ChooseTimeScaleForRaster/5)
% FileNameList=extractfield(FileDir,'name');char(FileNameList(FileNumber))
end
PlotRaster=1;
PlotHistogram=1;
%% write down which file numbers to run in the loop
for FileNumber = ChooseFileNumber+2
    disp(['File number ',num2str(FileNumber),'_',FileDirNLX(FileNumber).name]) % display file name
    disp(['bin size ',num2str(round((ChooseLimForHistogram/30),3))])
    FileFolder = [ComputerDir,'\',Date,'\',MouseName,'\',FileDirNLX(FileNumber).name];
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
%% write down which tetrodes to run in the loop
for TetrodeNumber = ChooseTetrodeNumber
    disp(['Tetrode number ', num2str(TetrodeNumber)])
 [TimeStamps_events, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV([FileFolder,'\','Events.nev'], [1 1 1 1 1], 1, 1, []); %event file
    try
    %% specifies filename: the correct filename format is "TT1_s.ntt","TT2_s.ntt" atc. 
[TimeStamps_cells, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike([FileFolder,'\',['TT',num2str(TetrodeNumber),'_Cut',num2str(threshold(1)-5),'.ntt']], [1 1 1 1 1], 1, 1, [] ); %clustering file
disp('Cut and Clustered file analysed')
    catch
try
[TimeStamps_cells, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike([FileFolder,'\',['TT',num2str(TetrodeNumber),'_s.ntt']], [1 1 1 1 1], 1, 1, [] ); %clustering file
disp('Custered file analysed')
catch
        try
[TimeStamps_cells, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike([FileFolder,'\',['TT',num2str(TetrodeNumber),'.ntt']], [1 1 1 1 1], 1, 1, [] ); %clustering file
ClustNumber=0; 
disp('uncustered file analysed')
        catch
        continue
        end
end
    end

%% write down which clusters to run in the loop
for Clust = ClustNumber 
    disp(['Cluster number ',num2str(Clust)])
     
if TimeStamps_cells(1) > TimeStamps_events(1)
    TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_events(1))/1000000;
    TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_events(1))/1000000;
else
    TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_cells(1))/1000000;
    TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_cells(1))/1000000; 
end

cell = cat(1,TimeStamps_cells_zeroed_s,CellNumbers); 
% select cluster to analyse
a = find(cell(2,:) ~=Clust); % "~=1" looks at cell#1, "~=2" looks at cell#2, etc.
cell(:,a) = [];
% ignor empty clusters:
if isempty(cell)
    continue
end
% open a figure and fit it to fullscrean display:
figure1 = figure('units','normalized','outerposition',[0 0 1 1]);
% use only TTL events that are not "0" 
laser = cat(1,TimeStamps_events_zeroed_s,TTLs);
laser(laser == 0) = NaN;
%count the number of laser pulses
eventcount = find(laser(2,:) > 0); 
eventcount = size(eventcount); 
% number the laser pulses from 1 to last
for i = 1:eventcount(:,2)-1 
    laser(2,2+2*i) = 1 + i;
end
for i = 1:eventcount(:,2)
    X(i,:) = [(laser(1,2*i)-ChooseTimeScaleForRaster):0.0001:(laser(1,2*i)+ChooseTimeScaleForRaster)]; %sets the range around stimulation -0.1 s to +0.5 s 
    sizeX = size(X(1,:)); %check the size of Xi array
    Y(i,(1:sizeX(:,2))) = i; %creates an Yi array for different trials Y1 - first trial Y2 - second trial etc.
        [v location_laser] = min(abs(X(i,:)-laser(1,2*i))); %finds the index of laser pulse time
        a = find(cell(1,:) >= laser(1,2*i)-ChooseTimeScaleForRaster & cell(1,:) <= laser(1,2*i)+ChooseTimeScaleForRaster); %finds the indices around the laser pulse (50ms before and 100 ms after)
    TF = isempty(a); %checks if it finds any values
    if TF == 0        
        numberofspikes = size(a);
        location_cell = zeros(1,numberofspikes(:,2));
        
        for j = 1:numberofspikes(:,2)
            [v location_cell(1,j)] = min(abs(X(i,:)-cell(1,a(j))));
        end
        
        X(i,:) = NaN; %changes all of the that are not event to NaN    
        X(i,location_laser) = 0; %makes the pulse location 0 
                for h = 1:numberofspikes(:,2); %inserts time values at the spike times
            X(i,location_cell(1,h)) = cell(1,a(h))-laser(1,2*i);
        end
    end
    if TF == 1
        X(i,:) = NaN; %changes all of the that are not event to NaN    
        X(i,location_laser) = 0; %makes the pulse location 0 
    end
      
    s(i) = scatter(X(i,:),Y(i,:));
        hold on
      
end
% Set the visual parameters and title of the plot:
try
set(s,'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k')
catch
    continue
end
ax = gca;
ax.XTick = [-ChooseTimeScaleForRaster:ChooseTickSizeForRaster:ChooseTimeScaleForRaster];
ax.XLim = [-ChooseTimeScaleForRaster ChooseTimeScaleForRaster];
% ax.YLim = [0 50];
% f1.InnerPosition = [680 558 800 350];
SaveFileName1 = [MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',TTLLaser,' ',FoodType,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust)];
title(SaveFileName1);

% Save the image files in image(jpeg) and fig(MATLAB) format:
saveas(figure1, fullfile(FileFolder, SaveFileName1), 'jpeg'); % here you save the figure
% saveas(figure1, fullfile(FileFolder, SaveFileName1), 'fig');% here you save the figure in MATLAB format

hold off

%% histogram
figure2 = figure;
k = find(isnan(X) == 1);
X(k) = 0;
k = find(X);
histo = X(k);
edges = ChooseEdgesForHistogram;
[HistoValues,HistoEdges]=histcounts(histo,edges);
[ValMax,LocationMax]=max(HistoValues);
LargestBinSizeToAverageBinRatio= 100*ValMax/mean(HistoValues);
disp(['Largest Bin ',num2str(ValMax)])
disp(['Largest Bin To Average Bin Ratio (%) ',num2str(LargestBinSizeToAverageBinRatio)])
LatancyMS=1000*HistoEdges(LocationMax);
disp(['Latency from light stimulation ',num2str(LatancyMS),' ms'])
histogram(histo,edges,'FaceColor','k','FaceAlpha',1);

ax = gca;
ax.XTick = [-ChooseLimForHistogram:(ChooseLimForHistogram/5):ChooseLimForHistogram];
ax.XLim = [-ChooseLimForHistogram ChooseLimForHistogram];
MaxLim=(round(1.2*(size(cell,2)/80)));
if MaxLim < 30;
    MaxLim =20;
end
if ChooseMaxLimitForHistogram ~=0;
    MaxLim =ChooseMaxLimitForHistogram;
end
ax.YLim = [0 MaxLim];
% figure2.InnerPosition = [680 558 800 350];
SaveFileName2 = ['Histogram',' ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',TTLLaser,' ',FoodType,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust)];
title(SaveFileName2);
saveas(figure2, fullfile(FileFolder, SaveFileName2), 'jpeg'); % here you save the figure
% saveas(figure2, fullfile(FileFolder, SaveFileName2), 'fig');% here you save the figure in MATLAB format

%% average response time

count_real = zeros(1,99);

for i = 1:99 %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
    a = find(histo(:,1) >= i/1000 - 0.001 & histo(:,1) <= i/1000 + 0.001);
    count_real(i) = size(a,1);
end

[Mr,Ir] = max(count_real); %M = highest number of spikes in interval, I = index of M
highestresponse = Mr;

bin_low = (Ir-1)/1000;
bin_up = (Ir+1)/1000;

a = find(histo(:,1) >= bin_low & histo(:,1) <= bin_up);
response_latency = mean(histo(a)); %gives average response latency

a = find(histo >= 0);
histo_pos = histo(a);

k = size(histo_pos,1);

bootstrap = zeros(1,10000);
for i = 1:10000;
    r = (0.1-0).*rand(size(histo_pos,1),1) + 0;
    count_shuffle = zeros(1,99);
        for j = 1:99 %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
            a = find(r(:,1) >= j/1000 - 0.001 & r(:,1) <= j/1000 + 0.001);
            count_shuffle(j) = size(a,1);
        end
    [Ms,Is] = max(count_shuffle);
    bootstrap(i) = Ms;
end

bootstrap = transpose(bootstrap);
pd = fitdist(bootstrap,'Normal');
x_values = 0:0.1:40;
y = pdf(pd,x_values);

f3 = figure
yyaxis left
histogram(bootstrap)
hold on
yyaxis right
plot(x_values,y,'LineWidth',2)

percentile = icdf(pd,0.999);
close all
end
end
end
close all
clearvars -except FileNameList ComputerDir ConditionFileOrder ConditionList Date FileDirNLX MouseName TetrodeNumber threshold
% clear 
%% Example for one laser pulse
    
% lasertime = laser(1,14);
% k = find(cell(1,:) >= laser(1,14)-0.05 & cell(1,:) <= laser(1,14)+0.95)
% 
% X = [(laser(1,14)-0.05):0.000001:(laser(1,14)+0.95)];
% sizeX = size(X(1,:));
% Y(1:sizeX(:,2)) = 1;
% [v location_l] = min(abs(X(1,:)-laser(1,14)));
% [v location_c1] = min(abs(X(1,:)-cell1(1,3)));
% [v location_c2] = min(abs(X(1,:)-cell1(1,4)));
% [v location_c3] = min(abs(X(1,:)-cell1(1,5)));
% 
% X(1,:) = NaN;
% X(1,location_l) = 0;
% X(1,location_c1) = cell1(1,3)-laser(1,14);
% X(1,location_c2) = cell1(1,4)-laser(1,14);
% X(1,location_c3) = cell1(1,5)-laser(1,14);
% scatter(X,Y)

