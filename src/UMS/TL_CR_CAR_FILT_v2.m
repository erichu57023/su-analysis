function TL_CR_CAR_FILT_v2(pathname)
%% FUNCTION: Subtracts common average signal from all files within a folder 
% (assumes each file is a different channel recorded concurrently). Then,
% besselfilters each channel and saves the trace in a new folder.

% INPUTS:
% [pathname] : pathname to folder of experiment date
% VERSIONS:
% v2: chunks data so that CREF can be run on all good_channels and then
% bessel filtering does not need to occur twice
%% Set Up
if ~strcmp(pathname(end),filesep)
    pathname=[pathname filesep];
end

% Main folder holding .mat data, extract file names
pathname_mat = [pathname 'spikes\mat' filesep]; 
temp = dir([pathname_mat '*.mat']);
files = {temp(1:end).name}';
[files ~] = sort(files);
clear temp;

% load attributes

%% Extract attributes
load([pathname 'attributes\attributes.mat']);

badchs = attributes.experiment.Bad_Channels;

%% Calculate common signal from every other good channel
% (Computer crashes if run on 32 channels..)

% Make destination folder
if ~exist([pathname 'spikes\FILT'])
    mkdir([pathname 'spikes\FILT']);
end
pathname_FILT= [pathname 'spikes\FILT' filesep];

common_data = 0;
for ch = 1 : length(files)
    display([pathname_mat files{ch}]);
    load([pathname_mat files{ch}]);
    if ~isa(rawMat.data , 'double')
        rawMat.data = double(rawMat.data);
    end
    FILT.data = besselfilter(4 , 300 , 6000 , rawMat.sampling_rate , rawMat.data');
    FILT.data = single(FILT.data);
    FILT.sampling_rate = rawMat.sampling_rate;
    save([pathname_FILT files{ch}] , 'FILT');
    if ~ismember(ch , badchs)
        common_data = common_data + FILT.data;
    end
    clear FILT.data;
end
common_average = common_data / (length(files));
clear temp common_data;

% Make destination folder
if ~exist([pathname 'spikes\CREF'])
    mkdir([pathname 'spikes\CREF']);
end
pathname_CREF= [pathname 'spikes\CREF' filesep];

for ch = 1:length(files)
        display([files{ch}]);
    load([pathname_FILT files{ch}]);
    CREF.data = FILT.data - common_average;
    CREF.sampling_rate = FILT.sampling_rate;
    save([pathname_CREF files{ch}],'CREF');
end
