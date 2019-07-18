clear all; clc;
addpath(genpath('C:\Users\admin\Documents\MATLAB\NeuralynxMatlabImportExport_v6.0.0'));
%% Import NTT files to Matlab 


%% Variables
ComputerDir='D:\CheetahData\NG\Data\Temp Storage'; %location of NLX files
Date='2019-05-22';
MouseName='SUBLAT4-9'; %%% specify mouse name
TetrodeNumber=[7]; % what Tetrodes to run [1 2 3 4 5 6 7 8]
FileDirNLX = dir([ComputerDir,'\',Date,'\',MouseName]);
% orgenize recordings to the order: 'Empty','Chow',Jelly_Exposure','Jelly_OFF', 'Laser 2hz','Laser 20hz'..
FileNameList={};
for FileNumber = [3:size(FileDirNLX,1)]
FileNameList(FileNumber-2) = {FileDirNLX(FileNumber).name};
IsCondition=strfind(char(FileNameList(FileNumber-2)),'Empty');
if ~isempty(IsCondition)
    ConditionFileOrder(1)=(FileNumber);
end
IsCondition=strfind(char(FileNameList(FileNumber-2)),'Chow');
if ~isempty(IsCondition)
    ConditionFileOrder(2)=(FileNumber);
end
IsCondition=strfind(char(FileNameList(FileNumber-2)),'Exposure');
if ~isempty(IsCondition)
    ConditionFileOrder(3)=(FileNumber);
end
IsCondition=strfind(char(FileNameList(FileNumber-2)),'OFF');
if ~isempty(IsCondition)
    ConditionFileOrder(4)=(FileNumber);
end
IsCondition=strfind(char(FileNameList(FileNumber-2)),'2hz');
if ~isempty(IsCondition)
    ConditionFileOrder(5)=(FileNumber);
end
IsCondition=strfind(char(FileNameList(FileNumber-2)),'20hz');
if ~isempty(IsCondition)
    ConditionFileOrder(6)=(FileNumber);
end

end %for

%% Import files from NLX
[Timestamps_spike_Laser_2hz, ScNumbers_spike_Laser_2hz, CellNumbers_spike_Laser_2hz, Features_spike_Laser_2hz, Samples_spike_Laser_2hz, Header_spike_Laser_2hz] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(5)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Laser_20hz, ScNumbers_spike_Laser_20hz, CellNumbers_spike_Laser_20hz, Features_spike_Laser_20hz, Samples_spike_Laser_20hz, Header_spike_Laser_20hz] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(6)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Empty, ScNumbers_spike_Empty, CellNumbers_spike_Empty, Features_spike_Empty, Samples_spike_Empty, Header_spike_Empty] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(1)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Chow, ScNumbers_spike_Chow, CellNumbers_spike_Chow, Features_spike_Chow, Samples_spike_Chow, Header_spike_Chow] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(2)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Jelly_Exposure, ScNumbers_spike_Jelly_Exposure, CellNumbers_spike_Jelly_Exposure, Features_spike_Jelly_Exposure, Samples_spike_Jelly_Exposure, Header_spike_Jelly_Exposure] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(3)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Jelly_OFF, ScNumbers_spike_Jelly_OFF, CellNumbers_spike_Jelly_OFF, Features_spike_Jelly_OFF, Samples_spike_Jelly_OFF, Header_spike_Jelly_OFF] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(4)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
%% Concatenate all files 
TimestampsAll=[Timestamps_spike_Laser_2hz,Timestamps_spike_Laser_20hz,Timestamps_spike_Empty,Timestamps_spike_Chow,Timestamps_spike_Jelly_Exposure,Timestamps_spike_Jelly_OFF];
ScNumbersAll=[ScNumbers_spike_Laser_2hz,ScNumbers_spike_Laser_20hz,ScNumbers_spike_Empty,ScNumbers_spike_Chow,ScNumbers_spike_Jelly_Exposure,ScNumbers_spike_Jelly_OFF];
CellNumbersAll=[CellNumbers_spike_Laser_2hz,CellNumbers_spike_Laser_20hz,CellNumbers_spike_Empty,CellNumbers_spike_Chow,CellNumbers_spike_Jelly_Exposure,CellNumbers_spike_Jelly_OFF];
FeaturesAll=[Features_spike_Laser_2hz,Features_spike_Laser_20hz,Features_spike_Empty,Features_spike_Chow,Features_spike_Jelly_Exposure,Features_spike_Jelly_OFF];
SamplesAll=cat(3,Samples_spike_Laser_2hz,Samples_spike_Laser_20hz,Samples_spike_Empty,Samples_spike_Chow,Samples_spike_Jelly_Exposure,Samples_spike_Jelly_OFF);
HeaderAll=[Header_spike_Laser_2hz,Header_spike_Laser_20hz,Header_spike_Empty,Header_spike_Chow,Header_spike_Jelly_Exposure,Header_spike_Jelly_OFF];

%% Export the concatenated file to .ntt format and save
mkdir([ComputerDir,'\',Date,'\',MouseName],'ConcatenatedFile')
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles.ntt'], 0, 1, [], [1 1 1 1 1], TimestampsAll,ScNumbersAll, CellNumbersAll, FeaturesAll, SamplesAll, HeaderAll);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sort the large file in 3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import the sorted file 
[TimestampsAll, ScNumbersAll, CellNumbersAll, FeaturesAll, SamplesAll, HeaderAll] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles_s.ntt'], [1 1 1 1 1], 1, 1, [] ); 
%% Split the sorted file back to condition files
% save 2hz file
% Timestamps_spike_Laser_2hz=TimestampsAll(1:size(Timestamps_spike_Laser_2hz,2));
% ScNumbers_spike_Laser_2hz=ScNumbersAll(1:size(ScNumbers_spike_Laser_2hz,2));
CellNumbers_spike_Laser_2hz=CellNumbersAll(1:size(CellNumbers_spike_Laser_2hz,2));
% Features_spike_Laser_2hz=FeaturesAll(:,1:size(Features_spike_Laser_2hz,2));
% Samples_spike_Laser_2hz=SamplesAll(:,:,(1:size(Samples_spike_Laser_2hz,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(5)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Laser_2hz,ScNumbers_spike_Laser_2hz, CellNumbers_spike_Laser_2hz, Features_spike_Laser_2hz, Samples_spike_Laser_2hz, Header_spike_Laser_2hz);
% save Laser_20hz file
Timestamps_spike_Laser_20hz=TimestampsAll(1:size(Timestamps_spike_Laser_20hz,2));
ScNumbers_spike_Laser_20hz=ScNumbersAll(1:size(ScNumbers_spike_Laser_20hz,2));
CellNumbers_spike_Laser_20hz=CellNumbersAll(1:size(CellNumbers_spike_Laser_20hz,2));
Features_spike_Laser_20hz=FeaturesAll(:,1:size(Features_spike_Laser_20hz,2)); 
Samples_spike_Laser_20hz=SamplesAll(:,:,(1:size(Samples_spike_Laser_20hz,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(6)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Laser_20hz,ScNumbers_spike_Laser_20hz, CellNumbers_spike_Laser_20hz, Features_spike_Laser_20hz, Samples_spike_Laser_20hz, Header_spike_Laser_20hz);
% save Empty file
Timestamps_spike_Empty=TimestampsAll(1:size(Timestamps_spike_Empty,2));
ScNumbers_spike_Empty=ScNumbersAll(1:size(ScNumbers_spike_Empty,2));
CellNumbers_spike_Empty=CellNumbersAll(1:size(CellNumbers_spike_Empty,2));
Features_spike_Empty=FeaturesAll(:,1:size(Features_spike_Empty,2)); 
Samples_spike_Empty=SamplesAll(:,:,(1:size(Samples_spike_Empty,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(1)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Empty,ScNumbers_spike_Empty, CellNumbers_spike_Empty, Features_spike_Empty, Samples_spike_Empty, Header_spike_Empty);
% save Chow file
Timestamps_spike_Chow=TimestampsAll(1:size(Timestamps_spike_Chow,2));
ScNumbers_spike_Chow=ScNumbersAll(1:size(ScNumbers_spike_Chow,2));
CellNumbers_spike_Chow=CellNumbersAll(1:size(CellNumbers_spike_Chow,2));
Features_spike_Chow=FeaturesAll(:,1:size(Features_spike_Chow,2)); 
Samples_spike_Chow=SamplesAll(:,:,(1:size(Samples_spike_Chow,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(2)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Chow,ScNumbers_spike_Chow, CellNumbers_spike_Chow, Features_spike_Chow, Samples_spike_Chow, Header_spike_Chow);
% save Jelly_Exposure file
Timestamps_spike_Jelly_Exposure=TimestampsAll(1:size(Timestamps_spike_Jelly_Exposure,2));
ScNumbers_spike_Jelly_Exposure=ScNumbersAll(1:size(ScNumbers_spike_Jelly_Exposure,2));
CellNumbers_spike_Jelly_Exposure=CellNumbersAll(1:size(CellNumbers_spike_Jelly_Exposure,2));
Features_spike_Jelly_Exposure=FeaturesAll(:,1:size(Features_spike_Jelly_Exposure,2));
Samples_spike_Jelly_Exposure=SamplesAll(:,:,(1:size(Samples_spike_Jelly_Exposure,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(3)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Jelly_Exposure,ScNumbers_spike_Jelly_Exposure, CellNumbers_spike_Jelly_Exposure, Features_spike_Jelly_Exposure, Samples_spike_Jelly_Exposure, Header_spike_Jelly_Exposure);
% save Jelly_OFF file
Timestamps_spike_Jelly_OFF=TimestampsAll(1:size(Timestamps_spike_Jelly_OFF,2));
ScNumbers_spike_Jelly_OFF=ScNumbersAll(1:size(ScNumbers_spike_Jelly_OFF,2));
CellNumbers_spike_Jelly_OFF=CellNumbersAll(1:size(CellNumbers_spike_Jelly_OFF,2));
Features_spike_Jelly_OFF=FeaturesAll(:,1:size(Features_spike_Jelly_OFF,2)); 
Samples_spike_Jelly_OFF=SamplesAll(:,:,(1:size(Samples_spike_Jelly_OFF,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(4)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Jelly_OFF,ScNumbers_spike_Jelly_OFF, CellNumbers_spike_Jelly_OFF, Features_spike_Jelly_OFF, Samples_spike_Jelly_OFF, Header_spike_Jelly_OFF);
