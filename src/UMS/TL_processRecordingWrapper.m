function TL_processRecordingWrapper(pathname , threshold , cref_filt , tets)
%% FUNCTION: Takes a folder 'YYMMDD' that contains csv files from igor (trials, header), 
% csv files from tdt containing neural data, and ddt files from tdt
% containing event information, and outputs automatically sorted units.
% Only variable that may need to be changed is the 'excel_path' in the Form
% Attributes section. This is the path to the user-inputted excel
% spreadsheet that contains experimental information from each recording

% INPUTS:
% [pathname] - pathname to YYMMDD folder. this folder must contain csv
% files from igor and tdt, and ddt files from tdt
% [threshold] - standard deviation threshold for spike cutting
% [cref_filt] - 'cref' to sort on common average referenced data, 'filt' to
% sort on just filtered data (no car)
% [tets] - vector of tetrode numbers to sort. In this wrapper, ch 1 is the
% deepest site

% Must remember to fill out excel spreadsheet containing experimental
% information, at least before running TL_form_attributes_v3
%% Arrange files and folders
% % Embedded functions: none
% TL_arrangefiles(pathname); 
% 
% % Convert recordings from csv to mat
% Embedded functions: none
% TL_tdtCSV2mat(pathname);
% 
% % Convert event channels from ddt to mat
% Embedded functions: ddtChanRead
% [Events] = TL_ddt2mat(pathname , chans , evt_names , save_events)
% TL_ddt2mat(pathname , [1 2 3] , {'Sweep_Start' , 'Sweep_Info' , 'Licks'} , 1);

%% Form attributes file from tdt events channels and igor tables
% Embedded functions: TL_DecodeRawMatEvents , TL_electrode_channels , replacevals 
excel_path = 'E:\Fmr1 Multi\Fmr1_Multi_Experiment_Info.xlsx'
% excel_path = 'E:\Autism Awake Data\Autism_Awake_Experiment_Info.xlsx'
TL_form_attributes_v3(pathname , excel_path , 1);

%% Filter and Common Average Reference spike recordings
% Embedded functions: none
TL_CR_CAR_FILT_v2(pathname);

%% Sort Recordings
% Embedded functions: ss_default_params_TL , ss_detect_BI , ss_align_BI ,
% ss_kmeans_BI , ss_energy_parallel_BI , ss_aggregate_BI , ums2k package
% functions
TL_spike_sorter_v3(pathname, threshold , cref_filt , tets , 1);







