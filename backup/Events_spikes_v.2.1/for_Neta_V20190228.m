%% Start
clear all
close all
addpath(genpath('E:\analysis\Events_spikes_v.2.1\Field Trip Neuralynx'));

%% Variables
ComputerDir='E:\CheetahData\NG\Data';
Date='2019-03-05';
MouseName='SU4-1'; %%% specify mouse name
%Add experiment details for proper registration:
ConditionStim='Real Stim';%specify 'Fake Stim', Real Stim', 'No Stim' 
FoodType='NoFood';%specify 'Jelly', Chow', 'NoFood'
FoodTime=[];%specify '5', '10', '15'
FoodConsumed=[];%specify in gr.
FileDir = dir([ComputerDir,'\',Date,'\',MouseName]);
FileNameList=extractfield(FileDir,'name');
%%
min_pulse_interval = 50000; % (µs) minimum distance between pulses (there's some noise on the TTL)
% tetrode_file = 'TT5_s.ntt'; % filename
% event_file = 'Events.nev'; % filename
attribute = 'pca_1'; % height, pca_1, etc....
%%
for FileNumber = ChooseFileNumber;
    char(FileNameList(FileNumber))
    FileFolder = [ComputerDir,'\',Date,'\',MouseName,'\',char(FileNameList(FileNumber))];
    %% Import tetrode
for TetrodeNumber = [1];
    try
    Tetrode = tetrode([FileFolder,'\',['TT',num2str(TetrodeNumber),'_s.ntt']]);
    catch 
    Tetrode = tetrode([FileFolder,'\',['TT',num2str(TetrodeNumber),'.ntt']]);
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