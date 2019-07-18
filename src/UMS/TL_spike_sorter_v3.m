function [spikes] = TL_spike_sorter_v3(pathname, threshold , cref_filt , tets , save_spikes)
%% This function spike sorts on single channels from tdt recordings processed via tdtCSV2mat -> CREF_tl

% INPUTS:
% [pathname] - path to recording date
% [threshold] - standard deviation threshold for spike cutting
% [cref] - use common average referenced signal (1), or just filtered
% signal (0)
% [tets] - index of which tetrodes to sort (in order of the tetrodes
% variable
% [save_spikes] = 1 to save data, 0 to not save

% NEED TO GET EVENT CHANNEL WORKING PROPERLY SO THAT DATA CAN BE ARRANGED
% FOR SORTING BY TRIAL (E.G. {TRIALS}(SAMPLES X CHANNELS)). OTHERWISE,
% 'PSEUDO TRIALS' WILL BE CREATED. CAN SET UP SO THAT IF THERE IS AN EVENTS
% FILE, IT CHECKS TRIAL START TIMES, OTHERWISE USES PSEUDO TRIALS

% v2: uses attributes file to extract trial start times
%% Set up parameters

if ~strcmp(pathname(end),filesep)
    pathname=[pathname filesep];
end
temp = strfind(pathname , '\');
animal_id = pathname(temp(end-2)+1:temp(end-1)-1);
rec_date = pathname(temp(end-1)+1:end-1); clear temp;

% Determine tetrode configuration (deep to superficial)
load([pathname 'attributes\attributes.mat']);
electrode = attributes.experiment.Electrode;
chmap = attributes.experiment.Channel_Map; % cell array, each cell is separate tetrode
badchs = attributes.experiment.Bad_Channels;

% tets = cell(length(chmap) , 1);
% [~ , indx , ~] = cellfun(@(z) setxor(z , badchs) , chmap , 'uniformoutput' , false);
% % Find the good channels on each tetrode to sort on
% for ch = 1 : length(indx)
%     tets{ch} = chmap{ch}(indx{ch});
% end

% Extract File names and set up destination folder
spiketype = {'CREF' , 'FILT'};
if strcmp(cref_filt , 'cref')
    data_pathname = [pathname 'spikes\CREF'];
    st = 1;
else if strcmp(cref_filt , 'filt')
        data_pathname = [pathname 'spikes\FILT'];
        st = 2;
    end
end

temp = dir(data_pathname);
files = {temp(3:end).name}; clear temp;

if save_spikes
    destination = [pathname 'spikes\SORTED']
    if exist(destination) ~= 7
        mkdir(destination);
    end
end

%% Sort individual channels and save
check_samprate = 0;
load([pathname 'attributes\attributes.mat']);
tri_starts = attributes.sweep.SweepOnsetSorting;
fshort = cellfun(@(x) x(end-5:end-4) , files , 'uniformoutput' , false);
for s = 1:length(tets)
    if ~isempty(chmap{tets(s)})
        name_str = [];
        for ch = 1:length(chmap{tets(s)})
            chstr = num2str(chmap{tets(s)}(ch));
            if length(chstr) == 1
                chstr = ['0' chstr];
            end
            chindx = strfind(fshort , chstr); chindx = ~cellfun(@isempty , chindx);
            chdata = load([data_pathname filesep files{chindx}]);
            %     DATA  -- matrix of data {trials}[ samples X channels ]
            %Only determine trial start and stop indices once (indices should
            %be the same on all channels)
            if s == 1 & ch == 1
                neuro_time = [1:length(chdata.(spiketype{st}).data)]/chdata.(spiketype{st}).sampling_rate;
                neuro_trials = zeros(length(tri_starts),2);
                for ts = 1 :length(tri_starts) - 1;
                    ind = find(neuro_time >= tri_starts(ts) & neuro_time < tri_starts(ts+1));
                    neuro_trials(ts,:) = [ind(1) ind(end)];
                    clear ind;
                end
                ind = find(neuro_time >= tri_starts(ts+1));
                neuro_trials(ts+1,:) = [ind(1) ind(end)];
                clear ind;
            end
            for ts = 1:length(tri_starts)
                DATA{ts}(:,ch) = chdata.(spiketype{st}).data(neuro_trials(ts,1):neuro_trials(ts,2));
                clear ind;
            end
        end
        
        spikes = ss_default_params_TL(chdata.(spiketype{st}).sampling_rate);
        spikes.params.display.dataset = [pathname(end-6:end-1) '_CH' num2str(s)];
        spikes.params.thresh = threshold;
        spikes  = ss_detect_BI(DATA,spikes);
        spikes = ss_align_BI(spikes);
        spikes = ss_kmeans_BI(spikes);
        spikes = ss_energy_parallel_BI(spikes);
        spikes = ss_aggregate_BI(spikes);
        spikes.info.channels = chmap{tets(s)};
        
        if save_spikes
            save([destination filesep animal_id '_' rec_date '_CH' num2str(tets(s)) '_SD' num2str(threshold) '_' cref_filt '_' datestr(now , 'yymmdd')] , ...
                'spikes');
        end
    end
end
end
