%% Start
clear all
close all
%%
%Add experiment details for proper registration:
ConditionStim='Real Stim';%specify 'Fake Stim', Real Stim', 'No Stim' 
FoodType='NoFood';%specify 'Jelly', Chow', 'NoFood'
FoodTime=[];%specify '5', '10', '15'
FoodConsumed=[];%specify in gr.
%% Path
% addpath(genpath('C:\Users\netas\Documents\MATLAB\Events_spikes_v.2.1\Field Trip Neuralynx'));
addpath(genpath('E:\analysis\Events_spikes_v.2.1\Field Trip Neuralynx'));

%% Variables
% ComputerDir='C:\Users\netas\Google Drive\Single Unit Data';
ComputerDir='E:\CheetahData\NG\Data';
Date='2019-03-06';
MouseName='SUBLAT1-3'; %%% specify mouse name
FileDir = dir([ComputerDir,'\',Date,'\',MouseName]);
%% Choose folders and tetrodes to run
ChooseFileNumber = [1]; % what files to run
ChooseTetrodeNumber=[1]; % what Tetrodes to run [1 2 3 4 5 6 7 8]
min_pulse_interval = 50000; % (µs) minimum distance between pulses (there's some noise on the TTL)
attribute = 'pca_1'; % height, pca_1, etc....
%%

%%
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
for TetrodeNumber = ChooseTetrodeNumber
    try
    Tetrode = tetrode([FileFolder,'\',['TT',num2str(TetrodeNumber),'_s.ntt']]);
   catch
        try
    Tetrode = tetrode([FileFolder,'\',['TT',num2str(TetrodeNumber),'.ntt']]);
        catch 
              'Folder does not exist'
        continue
        end
    end
%% Import the laser pulses
events = read_neuralynx_nev([FileFolder,'\','Events.nev']);
for i=1:length(events)
    pulses(i) = double(events(i).TimeStamp)*double(events(i).TTLValue);
end
pulses = pulses(pulses>0);

% Kick out small intervals
intervals = pulses(2:end) - pulses(1:end-1);
indexer = [intervals>=min_pulse_interval, true];
pulses = pulses(indexer);


%% Plot the data
show_cell = []; % leave empty to plot all the cells
[waveform_figure, scatter_figure] = Tetrode.present_figure(attribute);
[handel_figure, results] = Tetrode.peri_event(pulses, 'cell', show_cell);
FigureOne=figure;
timeline = squeeze(results(1,:,1));
average = [];
for i = 1:size(results,3)
    temp_results = squeeze(results(2:end,:,i));
    average{i} = mean(temp_results);
    average{i} = average{i} - mean(average{i}(21:26));
    plot(timeline, average{i})
    hold on
end

SaveFileName = [char(FileNameList(FileNumber)),' TT',num2str(TetrodeNumber)];
% title(SaveFileName);
movegui(FigureOne,'south')
figure(FigureOne);
title([SaveFileName,' Norm freq'])
saveas(FigureOne, fullfile(FileFolder, [SaveFileName,' Norm freq']), 'jpeg');
movegui(handel_figure,'north')
figure(handel_figure);
saveas(handel_figure,fullfile(FileFolder, [SaveFileName,' Events']),'jpeg');
title([SaveFileName,' Events'])
movegui(waveform_figure,'northeast')
figure(waveform_figure);
title([SaveFileName,' Waveform'])
saveas(waveform_figure,fullfile(FileFolder, [SaveFileName,' Waveform']),'jpeg');
close
movegui(scatter_figure,'northwest')
figure(scatter_figure);
title([SaveFileName,' Scatter'])
saveas(scatter_figure,fullfile(FileFolder, [SaveFileName,' Scatter']),'jpeg');
close
% 
% 
%% Plot any cells
for i=0:Tetrode.nr_cells
  [figure_1, figure_2]= Tetrode.show_cell(i, 'average'); 
%     Tetrode.show_cell(i);
movegui(figure_1,'southeast')
figure(figure_1);
title([SaveFileName,' cell',num2str(i)])
saveas(figure_1, fullfile(FileFolder, [SaveFileName,' cell',num2str(i)]), 'jpeg');
movegui(figure_2,'southwest')
close
figure(figure_2);
title([SaveFileName,' ISI cell',num2str(i)])
saveas(figure_2, fullfile(FileFolder, [SaveFileName,' ISI cell',num2str(i)]), 'jpeg');
close 
end
end
end