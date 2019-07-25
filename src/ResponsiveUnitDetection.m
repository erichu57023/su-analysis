%%% This code is designed to find stimulation-responsive cells and produce figures for visual assessment of the cell response. 
close all;
clear all;
clc;
%#ok<*UNRCH>

%% Parameters (please fill in)
% addpath(genpath('C:\Users\admin\Documents\MATLAB\NeuralynxMatlabImportExport_v6.0.0'));
ComputerDir = 'D:\CheetahData\NG\Data\Temp Storage';

Date = '2019-07-24';          % specify date
MouseName = 'SUBLAT1-7';      % specify mouse name

ChooseFileNumber = 1;         % choose files to run
ChooseTetrodeNumber = 2;    % choose tetrodes to run
ClustNumber = 0:10;

hz=20;

ViewEventOnlySpikesFlag = false;  % saves graphs of windowed spikes
window = [1 4];                      % define the window width (in ms)
SaveEventOnlyNTTFlag = true;    % saves windowed spikes into a seperate TT#_events.ntt file
TetrodesToExtract = ChooseTetrodeNumber;         % choose tetrodes to extract for new NTTs.

%% Setup
switch hz
    case 2
    ChooseMaxLimitForHistogram = 0; % When this is set to 0 the code chooses the value.
    ChooseLimForHistogram = 0.2;
    ChooseEdgesForHistogram = -ChooseLimForHistogram : 0.004 : ChooseLimForHistogram; % bin size is the middle value round((ChooseLimForHistogram/50),4)
    ChooseTimeScaleForRaster = 0.2;
    case 20
    ChooseMaxLimitForHistogram = 100; % When this is set to 0 the code chooses the value.
    ChooseLimForHistogram = 0.02;
    ChooseEdgesForHistogram = -ChooseLimForHistogram : 0.0004 : ChooseLimForHistogram; % bin size is the middle value round((ChooseLimForHistogram/50),4)
    ChooseTimeScaleForRaster = 0.02;
end

ChooseTickSizeForRaster = (ChooseTimeScaleForRaster / 5); % default ChooseTickSizeForRaster=(ChooseTimeScaleForRaster / 5)
FileDir = dir([ComputerDir,'\',Date,'\',MouseName]);

%% Loop through each file
for FileNumber = ChooseFileNumber + 2
    disp(['File number ', num2str(FileNumber-2), '_', FileDir(FileNumber).name]) % display file name
    disp(['bin size ', num2str(round((ChooseLimForHistogram / 30), 3))])
    FileFolder = [ComputerDir, '\', Date, '\', MouseName, '\', FileDir(FileNumber).name];
    %now we split the name and date componantes:
    %%%Excel format: [=D22&"_"&E22&" min_SU8-3_"&"z"&B22&"_"&F22&"_"&G22&"_"&H22&"_"&I22&"_"&round(J22,3)&" "&"gr"]
    SplitName = strsplit(FileDir(FileNumber).name, '_');
    TimeSplit = strsplit(char(SplitName(2)), '-');
    %now we assign the variables with the componanets we got from the file name 
    try DateStamp    = char(SplitName(1)); catch; DateStamp = ''; end
    try TimeStamp    = [char(TimeSplit(1)),'-',char(TimeSplit(2))]; catch; TimeStamp = ''; end
    try RecDuration  = char(SplitName(3)); catch; RecDuration = ''; end
    try MouseNumber  = char(SplitName(4)); catch; MouseNumber = MouseName; end
    try Zlocation    = char(SplitName(5)); catch; Zlocation = ''; end
    try TTLLaser     = char(SplitName(6)); catch; TTLLaser = ''; end
    try StimType     = char(SplitName(7)); catch; StimType = ''; end
    try StimFreq     = char(SplitName(8)); catch; StimFreq = ''; end
    try FoodType     = char(SplitName(9)); catch; FoodType = ''; end
    try FoodConsumed = [char(SplitName(10)),' gr']; catch; FoodConsumed = ''; end
    
    %% Loop through each tetrode
    for TetrodeNumber = ChooseTetrodeNumber
        disp(['Tetrode number ', num2str(TetrodeNumber)]);
        tetrodeEventsNotSaved = SaveEventOnlyNTTFlag; % used to make sure event-extracted files are only created once for each tetrode, regardless of cluster number.
        try
            % Specifies filename: the correct filename format is "TT1_s.ntt","TT2_s.ntt" etc. 
            [TimeStamps_events, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV([FileFolder,'\','Events.nev'], [1 1 1 1 1], 1, 1, []); %event file
            [TimeStamps_cells, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike([FileFolder,'\',['TT',num2str(TetrodeNumber),'_s.ntt']], [1 1 1 1 1], 1, 1, [] ); %clustering file
            if ClustNumber == 0 
                ClustNumber = 0:10;
            end
        catch
            try
                [TimeStamps_events, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV([FileFolder, '\', 'Events.nev'], [1 1 1 1 1], 1, 1, []); %event file
                [TimeStamps_cells, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike([FileFolder,'\', ['TT', num2str(TetrodeNumber), '.ntt']], [1 1 1 1 1], 1, 1, [] ); %clustering file
                ClustNumber=0;          
            catch
                continue
            end
        end
        % TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_cells(1)) / 1000000; %for spikes
        % TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_events(1)) / 1000000; %for laser pulses*
        
        %% Loop through each cluster
        for Clust = ClustNumber 
            disp(['  ', 'Cluster number ', num2str(Clust)]);

            if TimeStamps_cells(1) > TimeStamps_events(1)
                TimeStamps_events_zeroed_s = (TimeStamps_events - TimeStamps_events(1)) / 1000000;
                TimeStamps_cells_zeroed_s = (TimeStamps_cells - TimeStamps_events(1)) / 1000000;
            else
                TimeStamps_events_zeroed_s = (TimeStamps_events - TimeStamps_cells(1)) / 1000000;
                TimeStamps_cells_zeroed_s = (TimeStamps_cells - TimeStamps_cells(1)) / 1000000; 
            end

            cell = cat(1, TimeStamps_cells_zeroed_s, CellNumbers); 
            cell(:, CellNumbers ~= Clust) = []; % "~=1" looks at cell#1, "~=2" looks at cell#2, etc.
            
            %%% Now we test to ignore empty clusters:
            if isempty(cell)
                continue
            end
            
            %%% Now we open a figure and fit it to fullscrean display:
            figure1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
            hold on;
            
            laser = cat(1, TimeStamps_events_zeroed_s, TTLs);
            laser = laser(:, logical(TTLs));

            eventcount = sum(laser(2, :));
            laser(2, :) = 1:eventcount;
            
            Xall = zeros(eventcount, ChooseTimeScaleForRaster * 2e4 + 1);

            for i = 1 : eventcount
                minBound = laser(1,i) - ChooseTimeScaleForRaster;
                maxBound = laser(1,i) + ChooseTimeScaleForRaster;
                X = (minBound : 0.0001 : maxBound); % defines the range around stimulation -0.2 s to +0.2 s 
                sizeX = length(X); % check the size of Xi array
                Y = i * ones(1, sizeX); % aligns the plotted spikes of this event into one horizontal line
                [~, location_laser] = min(abs(X-laser(1,i))); % finds the index of laser pulse time

                a = find(cell(1,:) >= minBound & cell(1,:) <= maxBound); % finds the spikes around the laser pulse
                if ~isempty(a)  % checks if it finds any values
                    numberofspikes = size(a);
                    
                    location_cell = zeros(1, numberofspikes(:,2));
                    for j = 1:numberofspikes(:,2)
                        [~, location_cell(1,j)] = min(abs(X - cell(1,a(j))));
                    end

                    X(:) = NaN; % changes all of the x values that are not events to NaN    
                    X(location_laser) = 0; % cernters the pulse location on 0 

                    for h = 1:numberofspikes(:,2) % inserts time values at the spike times
                        X(location_cell(1,h)) = cell(1,a(h))-laser(1,i);
                    end
                else
                    X(:) = NaN; % changes all of the that are not event to NaN    
                    X(location_laser) = 0; % centers the pulse location on 0 
                end

                Xall(i, :) = X;
                line(X, Y, 'LineStyle', 'none', 'Marker', 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            end
            hold off;
            
            ax = gca;
            ax.XTick = -ChooseTimeScaleForRaster : ChooseTickSizeForRaster : ChooseTimeScaleForRaster;
            ax.XLim = [-ChooseTimeScaleForRaster ChooseTimeScaleForRaster];
            % ax.YLim = [0 50];
            % f1.InnerPosition = [680 558 800 350];
            SaveFileName1 = [MouseName, ' ', DateStamp, ' ',TimeStamp, ' ', RecDuration, ' ', Zlocation, ' ', TTLLaser, ' ', FoodType, ' ', ' T', num2str(TetrodeNumber), ' C', num2str(Clust)];
            title(SaveFileName1);           

            %% histogram
            figure2 = figure;
            Xall = nonzeros(Xall);
            edges = ChooseEdgesForHistogram;
            [HistoValues, HistoEdges] = histcounts(Xall, edges);
            [ValMax, LocationMax] = max(HistoValues);
            LargestBinSizeToAverageBinRatio = 100 * ValMax / mean(HistoValues);
            disp(['  ', 'Largest Bin ', num2str(ValMax)])
            disp(['  ', 'Largest Bin To Average Bin Ratio (%) ', num2str(LargestBinSizeToAverageBinRatio)])
            LatancyMS = 1000 * HistoEdges(LocationMax);
            disp(['  ', 'Latency from light stimulation ', num2str(LatancyMS), ' ms'])
            histogram(Xall, edges, 'FaceColor', 'k', 'FaceAlpha', 1);

            ax = gca;
            ax.XTick = -ChooseLimForHistogram : (ChooseLimForHistogram / 5) : ChooseLimForHistogram;
            ax.XLim = [-ChooseLimForHistogram ChooseLimForHistogram];
            MaxLim = round(1.2 * (size(cell, 2) / 80));
            if ChooseMaxLimitForHistogram ~= 0
                MaxLim = ChooseMaxLimitForHistogram;
            elseif MaxLim < 30
                MaxLim = 30;
            end
            ax.YLim = [0 MaxLim];
            % figure2.InnerPosition = [680 558 800 350];
            
            SaveFileName2 = ['Histogram', ' ', MouseName, ' ', DateStamp, ' ', TimeStamp, ' ', RecDuration, ' ', Zlocation, ' ', TTLLaser, ' ', FoodType, ' ', ' T', num2str(TetrodeNumber), ' C', num2str(Clust)];
            title(SaveFileName2);
            
            %% 
            if ViewEventOnlySpikesFlag && tetrodeEventsNotSaved
                if any(TetrodesToExtract == TetrodeNumber)
                    figure3 = ExtractEventSpikes(FileFolder, TetrodeNumber, Clust, window, true, true);
                else
                    figure3 = ExtractEventSpikes(FileFolder, TetrodeNumber, Clust, window, true, false);
                end
            elseif ViewEventOnlySpikesFlag
                figure3 = ExtractEventSpikes(FileFolder, TetrodeNumber, Clust, window, true, false);
            elseif tetrodeEventsNotSaved
                ExtractEventSpikes(FileFolder, TetrodeNumber, Clust, window, false, true);
            end
            tetrodeEventsNotSaved = false; % avoids saving _events.ntt files redundantly if there are multiple clusters in this tetrode.
            
            %%% Now we save the image files in image(jpeg) and fig(MATLAB) format:
            saveas(figure1, fullfile(FileFolder, SaveFileName1), 'jpeg');
            saveas(figure2, fullfile(FileFolder, SaveFileName2), 'jpeg');
            % saveas(figure1, fullfile(FileFolder, SaveFileName1), 'fig');% here you save the figure in MATLAB format
            % saveas(figure2, fullfile(FileFolder, SaveFileName2), 'fig');% here you save the figure in MATLAB format
            if ViewEventOnlySpikesFlag
                SaveFileName3 = ['EventOnlySpikes', ' ', MouseName, ' ', DateStamp, ' ', TimeStamp, ' ', RecDuration, ' ', Zlocation, ' ', TTLLaser, ' ', FoodType, ' ', ' T', num2str(TetrodeNumber), ' C', num2str(Clust)];
                sgtitle(SaveFileName3); 
                saveas(figure3, fullfile(FileFolder, SaveFileName3), 'jpeg');
                % saveas(figure3, fullfile(FileFolder, SaveFileName3), 'fig');% here you save the figure in MATLAB format
            end
                

            % % %% average response time
            % % 
            % % count_real = zeros(1,99);
            % % 
            % % for i = 1:99 %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
            % %     a = find(histo(:,1) >= (i/1000 - 0.001) & histo(:,1) <= (i/1000 + 0.001));
            % %     count_real(i) = size(a,1);
            % % end
            % % 
            % % [Mr,Ir] = max(count_real); %M = highest number of spikes in interval, I = index of M
            % % highestresponse = Mr;
            % % 
            % % bin_low = (Ir-1)/1000;
            % % bin_up = (Ir+1)/1000;
            % % 
            % % response_latency = mean(histo(histo(:,1) >= bin_low & histo(:,1) <= bin_up)); %gives average response latency
            % % 
            % % histo_pos = histo(histo >= 0);
            % % 
            % % k = size(histo_pos,1);
            % % 
            % % bootstrap = zeros(1,10000);
            % % for i = 1:10000
            % %     r = (0.1).*rand(size(histo_pos,1),1);
            % %     count_shuffle = zeros(1,99);
            % %         for j = 1:99 %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
            % %             a = find(r(:,1) >= j/1000 - 0.001 & r(:,1) <= j/1000 + 0.001);
            % %             count_shuffle(j) = size(a,1);
            % %         end
            % %     [Ms,Is] = max(count_shuffle);
            % %     bootstrap(i) = Ms;
            % % end
            % % 
            % %  
            % % bootstrap = transpose(bootstrap);
            % % pd = fitdist(bootstrap,'Normal');
            % % x_values = 0:0.1:40;
            % % y = pdf(pd,x_values);
            % % 
            % % figure3 = figure;
            % % yyaxis left
            % % histogram(bootstrap)
            % % hold on
            % % yyaxis right
            % % plot(x_values,y,'LineWidth',2)
            % % 
            % % percentile = icdf(pd,0.999);
            % % SaveFileName3 = ['Ave. response time',' ',MouseName,' ',DateStamp,' ',TimeStamp,' ',RecDuration,' ', Zlocation,' ',TTLLaser,' ',FoodType,' ',' T',num2str(TetrodeNumber),' C',num2str(Clust)];
            % % title(SaveFileName3);
            % % % saveas(figure3, fullfile(FileFolder, SaveFileName3), 'jpeg'); % here you save the figure
            close all
        end
    end
end

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

