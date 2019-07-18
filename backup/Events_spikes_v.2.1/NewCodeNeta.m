clear all; clc;
addpath(genpath('C:\Users\admin\Documents\MATLAB\NeuralynxMatlabImportExport_v6.0.0'));
% % % You can always type >>help Nlx2MatCSC or the appropriate function and get assistance.
% % % Filtering tool for use with .ncs files: http://neuralynx.com/software/NlxCSCFiltering.exe
% % % CSC Spike Extractor: http://neuralynx.com/software/CscSpikeExtractor_v301.zip

%% Import csc files to Matlab 

[TimestampsLaser1, ChannelNumbersLaser1, SampleFrequenciesLaser1, NumberOfValidSamplesLaser1, SamplesLaser1, HeaderLaser1] = Nlx2MatCSC('CSC17_Laser.ncs',[1 1 1 1 1], 1, 1, [] );
[TimestampsEmpty1, ChannelNumbersEmpty1, SampleFrequenciesEmpty1, NumberOfValidSamplesEmpty1, SamplesEmpty1, HeaderEmpty1] = Nlx2MatCSC('CSC17_Empty.ncs',[1 1 1 1 1], 1, 1, [] );
[TimestampsChow1, ChannelNumbersChow1, SampleFrequenciesChow1, NumberOfValidSamplesChow1, SamplesChow1, HeaderChow1] = Nlx2MatCSC('CSC17_Chow.ncs',[1 1 1 1 1], 1, 1, [] );
[TimestampsJelly_Exposure1, ChannelNumbersJelly_Exposure1, SampleFrequenciesJelly_Exposure1, NumberOfValidSamplesJelly_Exposure1, SamplesJelly_Exposure1, HeaderJelly_Exposure1] = Nlx2MatCSC('CSC17_Jelly_Exposure.ncs',[1 1 1 1 1], 1, 1, [] );
[TimestampsJelly_OFF1, ChannelNumbersJelly_OFF1, SampleFrequenciesJelly_OFF1, NumberOfValidSamplesJelly_OFF1, SamplesJelly_OFF1, HeaderJelly_OFF1] = Nlx2MatCSC('CSC17_Jelly_OFF.ncs',[1 1 1 1 1], 1, 1, [] );
%% Connect all files of one channel by timestamp to one "AllFiles" 
%find the timestamps (start/end)  for each condition 
RangeLaser=[TimestampsLaser1(1),TimestampsLaser1(end)];
RangeEmpty=[TimestampsEmpty1(1),TimestampsEmpty1(end)];
RangeChow=[TimestampsChow1(1),TimestampsChow1(end)];
RangeJelly_Exposure=[TimestampsJelly_Exposure1(1),TimestampsJelly_Exposure1(end)];
RangeJelly_OFF=[TimestampsJelly_OFF1(1),TimestampsJelly_OFF1(end)];
% Merge timestamps of all the files
TimestampsAllFiles1=[TimestampsLaser1,TimestampsEmpty1,TimestampsJelly_Exposure1,TimestampsJelly_OFF1,TimestampsChow1];
RangeAllFiles1=[TimestampsAllFiles1(1),TimestampsAllFiles1(end)];
SamplesAllFiles1=[SamplesLaser1,SamplesEmpty1];

ChannelNumbersAllFiles1=[ChannelNumbersLaser1,ChannelNumbersEmpty1];
SampleFrequenciesAllFiles1=[SampleFrequenciesLaser1,SampleFrequenciesEmpty1];
NumberOfValidSamplesAllFiles1=[NumberOfValidSamplesLaser1,NumberOfValidSamplesEmpty1];
HeaderAllFiles1=[HeaderLaser1,HeaderEmpty1];
%% 3. filter 300-6000 (the attached files are already filtered but for future recordings)

%% 4. collect spikes using threshold (numeric - 55 micro-volts or 4 STD from baseline)
%% 5. export large "sorting file" to "sorting file .ntt"  format (or other 3D compatible format)
%% 6. sort in 3D 
%% 7. import to Matlab 
%% 8. split  "sorting file .ntt"   to condition files and get 5 .ntt sorted files
