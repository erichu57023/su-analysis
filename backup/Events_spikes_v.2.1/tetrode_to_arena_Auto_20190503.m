%% Start
close all
clear all
% addpath(('C:\Users\admin\Documents\MATLAB\Field Trip Neuralynx'));

%% Variables
% filenames
filename_tetrode = 'TT5.ntt';
filename_biobserve = '20190431_SUBLAT1_1_Chow';
filename_eventfile = 'Events.nev';

% other variables
timestamp_interval = 5; %Sec
%arena = [200 -290 160 80];
nr_bins = 20; %Actual the square root of the number of bins


%% Import the files
Tetrode = tetrode(filename_tetrode);
Events = read_neuralynx_nev(filename_eventfile);
Track = biobserve_import(filename_biobserve);


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
ClusterNumber=3
cell_nr = ClusterNumber;
cell_stamps = Tetrode.timestamps(Tetrode.cells==ClusterNumber)';


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


%% Make a figure
% First the track and all APs
TrackFigure=figure;
% subplot(2,1,1)
plot(Track(:,1), Track(:,2))
axis equal
hold on
ax = gca;
ax.YDir = 'reverse';

scatter(data(:,1), data(:,2)); % plots firing over the physical distance

% Then the heatplot
HeatMapFigure=figure;
% subplot(2,1,2)
imagesc(map);


%% This code wil not run
% caxis([low high]); % click on the figure yet
% selection = map(30:50, 35:50); % will select row 30:50 and column 35:50







