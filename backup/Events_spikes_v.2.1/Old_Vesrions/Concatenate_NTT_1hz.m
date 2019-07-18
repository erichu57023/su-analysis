clear all; clc;close all
addpath(genpath('C:\Users\admin\Documents\MATLAB\NeuralynxMatlabImportExport_v6.0.0'));
%% Import NTT files to Matlab 

%% Variables
ComputerDir='G:\'; %location of NLX files
Date='2019-04-30';
MouseName='SUBLAT1-1'; %%% specify mouse name
TetrodeNumber=[5]; % what Tetrodes to run [1 2 3 4 5 6 7 8]
FileDirNLX = dir([ComputerDir,'\',Date,'\',MouseName]);
PlotFigures=0; % indicate if should plot average figures 0 for no 1 for Yes

%%orgenize recordings to the order: 'Empty','Chow',Jelly_Exposure','Jelly_OFF', 'Laser 2hz','Laser 20hz'..
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
IsCondition=strfind(char(FileNameList(FileNumber-2)),'Laser');
if ~isempty(IsCondition)
    ConditionFileOrder(5)=(FileNumber);
end
end %for

%% Import files from NLX
[Timestamps_spike_Laser_Laser, ScNumbers_spike_Laser_Laser, CellNumbers_spike_Laser_Laser, Features_spike_Laser_Laser, Samples_spike_Laser_Laser, Header_spike_Laser_Laser] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(5)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Empty, ScNumbers_spike_Empty, CellNumbers_spike_Empty, Features_spike_Empty, Samples_spike_Empty, Header_spike_Empty] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(1)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Chow, ScNumbers_spike_Chow, CellNumbers_spike_Chow, Features_spike_Chow, Samples_spike_Chow, Header_spike_Chow] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(2)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Jelly_Exposure, ScNumbers_spike_Jelly_Exposure, CellNumbers_spike_Jelly_Exposure, Features_spike_Jelly_Exposure, Samples_spike_Jelly_Exposure, Header_spike_Jelly_Exposure] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(3)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
[Timestamps_spike_Jelly_OFF, ScNumbers_spike_Jelly_OFF, CellNumbers_spike_Jelly_OFF, Features_spike_Jelly_OFF, Samples_spike_Jelly_OFF, Header_spike_Jelly_OFF] =  Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(4)).name),'\','TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); 
%% Concatenate all files 
TimestampsAll=[Timestamps_spike_Laser_Laser,Timestamps_spike_Empty,Timestamps_spike_Chow,Timestamps_spike_Jelly_Exposure,Timestamps_spike_Jelly_OFF];
ScNumbersAll=[ScNumbers_spike_Laser_Laser,ScNumbers_spike_Empty,ScNumbers_spike_Chow,ScNumbers_spike_Jelly_Exposure,ScNumbers_spike_Jelly_OFF];
CellNumbersAll=[CellNumbers_spike_Laser_Laser,CellNumbers_spike_Empty,CellNumbers_spike_Chow,CellNumbers_spike_Jelly_Exposure,CellNumbers_spike_Jelly_OFF];
FeaturesAll=[Features_spike_Laser_Laser,Features_spike_Empty,Features_spike_Chow,Features_spike_Jelly_Exposure,Features_spike_Jelly_OFF];
SamplesAll=cat(3,Samples_spike_Laser_Laser,Samples_spike_Empty,Samples_spike_Chow,Samples_spike_Jelly_Exposure,Samples_spike_Jelly_OFF);
HeaderAll={Header_spike_Laser_Laser,Header_spike_Empty,Header_spike_Chow,Header_spike_Jelly_Exposure,Header_spike_Jelly_OFF};

%% Export the concatenated file to .ntt format and save
mkdir([ComputerDir,'\',Date,'\',MouseName],'ConcatenatedFile')
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles.ntt'], 0, 1, [], [1 1 1 1 1], TimestampsAll,ScNumbersAll, CellNumbersAll, FeaturesAll, SamplesAll, HeaderAll);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sort the large file in 3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import the sorted file 
try
    [TimestampsAll, ScNumbersAll, CellNumbersAll, FeaturesAll, SamplesAll, HeaderAll] = Nlx2MatSpike([ComputerDir,'\',Date,'\',MouseName,'\','ConcatenatedFile','\','AllFiles_s.ntt'], [1 1 1 1 1], 1, 1, [] ); 
catch
    disp('you must sort the file first')
end
    %% Split the sorted file back to condition files
% save Laser file
Timestamps_spike_Laser_Laser=TimestampsAll(1:size(Timestamps_spike_Laser_Laser,2));
ScNumbers_spike_Laser_Laser=ScNumbersAll(1:size(ScNumbers_spike_Laser_Laser,2));
CellNumbers_spike_Laser_Laser=CellNumbersAll(1:size(CellNumbers_spike_Laser_Laser,2));
Features_spike_Laser_Laser=FeaturesAll(:,1:size(Features_spike_Laser_Laser,2));
Samples_spike_Laser_Laser=SamplesAll(:,:,(1:size(Samples_spike_Laser_Laser,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(5)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Laser_Laser,ScNumbers_spike_Laser_Laser, CellNumbers_spike_Laser_Laser, Features_spike_Laser_Laser, Samples_spike_Laser_Laser, Header_spike_Laser_Laser);
if PlotFigures==1
    Laser_Laser_Figure=figure;
for i=1:size(Header_spike_Laser_Laser,1)
IsADBitVolts=strfind((char(Header_spike_Laser_Laser{i,:})),'ADBitVolts');
if ~isempty(IsADBitVolts)
    ADVoltsList=Header_spike_Laser_Laser{i,:};
   end %if
end %for
 ADVoltsList=strsplit(ADVoltsList,' ');
         for Electrode=2:5
        ADBitVolts(Electrode-1) = str2double(ADVoltsList{1,Electrode});
        Samples_spike_Laser_Laser(:,Electrode-1,:)=Samples_spike_Laser_Laser(:,Electrode-1,:)*(ADBitVolts(Electrode-1));
        plot(mean(Samples_spike_Laser_Laser,3));
        title('Laser Laser');
        end %for
end %if
%% save Empty file
Timestamps_spike_Empty=TimestampsAll(1:size(Timestamps_spike_Empty,2));
ScNumbers_spike_Empty=ScNumbersAll(1:size(ScNumbers_spike_Empty,2));
CellNumbers_spike_Empty=CellNumbersAll(1:size(CellNumbers_spike_Empty,2));
Features_spike_Empty=FeaturesAll(:,1:size(Features_spike_Empty,2)); 
Samples_spike_Empty=SamplesAll(:,:,(1:size(Samples_spike_Empty,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(1)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Empty,ScNumbers_spike_Empty, CellNumbers_spike_Empty, Features_spike_Empty, Samples_spike_Empty, Header_spike_Empty);
if PlotFigures==1
Empty_Figure=figure;
for i=1:size(Header_spike_Empty,1)
IsADBitVolts=strfind((char(Header_spike_Empty{i,:})),'ADBitVolts');
if ~isempty(IsADBitVolts)
    ADVoltsList=Header_spike_Empty{i,:};
   end %if
end %for
 ADVoltsList=strsplit(ADVoltsList,' ');
        for Electrode=2:5
        ADBitVolts(Electrode-1) = str2double(ADVoltsList{1,Electrode});
        Samples_spike_Empty(:,Electrode-1,:)=Samples_spike_Empty(:,Electrode-1,:)*(ADBitVolts(Electrode-1));
        plot(mean(Samples_spike_Empty,3));
        title('Empty');
        end %for
        end %if
% save Chow file
Timestamps_spike_Chow=TimestampsAll(1:size(Timestamps_spike_Chow,2));
ScNumbers_spike_Chow=ScNumbersAll(1:size(ScNumbers_spike_Chow,2));
CellNumbers_spike_Chow=CellNumbersAll(1:size(CellNumbers_spike_Chow,2));
Features_spike_Chow=FeaturesAll(:,1:size(Features_spike_Chow,2)); 
Samples_spike_Chow=SamplesAll(:,:,(1:size(Samples_spike_Chow,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(2)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Chow,ScNumbers_spike_Chow, CellNumbers_spike_Chow, Features_spike_Chow, Samples_spike_Chow, Header_spike_Chow);
if PlotFigures==1
Chow_Figure=figure;
for i=1:size(Header_spike_Chow,1)
IsADBitVolts=strfind((char(Header_spike_Chow{i,:})),'ADBitVolts');
if ~isempty(IsADBitVolts)
    ADVoltsList=Header_spike_Chow{i,:};
   end %if
end %for
 ADVoltsList=strsplit(ADVoltsList,' ');
        for Electrode=2:5
        ADBitVolts(Electrode-1) = str2double(ADVoltsList{1,Electrode});
        Samples_spike_Chow(:,Electrode-1,:)=Samples_spike_Chow(:,Electrode-1,:)*(ADBitVolts(Electrode-1));
        plot(mean(Samples_spike_Chow,3));
        title('Chow');
        end %for
        end
% save Jelly_Exposure file
Timestamps_spike_Jelly_Exposure=TimestampsAll(1:size(Timestamps_spike_Jelly_Exposure,2));
ScNumbers_spike_Jelly_Exposure=ScNumbersAll(1:size(ScNumbers_spike_Jelly_Exposure,2));
CellNumbers_spike_Jelly_Exposure=CellNumbersAll(1:size(CellNumbers_spike_Jelly_Exposure,2));
Features_spike_Jelly_Exposure=FeaturesAll(:,1:size(Features_spike_Jelly_Exposure,2));
Samples_spike_Jelly_Exposure=SamplesAll(:,:,(1:size(Samples_spike_Jelly_Exposure,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(3)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Jelly_Exposure,ScNumbers_spike_Jelly_Exposure, CellNumbers_spike_Jelly_Exposure, Features_spike_Jelly_Exposure, Samples_spike_Jelly_Exposure, Header_spike_Jelly_Exposure);
if PlotFigures==1
Jelly_Exposure_Figure=figure;
for i=1:size(Header_spike_Jelly_Exposure,1)
IsADBitVolts=strfind((char(Header_spike_Jelly_Exposure{i,:})),'ADBitVolts');
if ~isempty(IsADBitVolts)
    ADVoltsList=Header_spike_Jelly_Exposure{i,:};
   end %if
end %for
 ADVoltsList=strsplit(ADVoltsList,' ');
        for Electrode=2:5
        ADBitVolts(Electrode-1) = str2double(ADVoltsList{1,Electrode});
        Samples_spike_Jelly_Exposure(:,Electrode-1,:)=Samples_spike_Jelly_Exposure(:,Electrode-1,:)*(ADBitVolts(Electrode-1));
        plot(mean(Samples_spike_Jelly_Exposure,3));
        title('Jelly Exposure');
        end %for
end %if
% save Jelly_OFF file
Timestamps_spike_Jelly_OFF=TimestampsAll(1:size(Timestamps_spike_Jelly_OFF,2));
ScNumbers_spike_Jelly_OFF=ScNumbersAll(1:size(ScNumbers_spike_Jelly_OFF,2));
CellNumbers_spike_Jelly_OFF=CellNumbersAll(1:size(CellNumbers_spike_Jelly_OFF,2));
Features_spike_Jelly_OFF=FeaturesAll(:,1:size(Features_spike_Jelly_OFF,2)); 
Samples_spike_Jelly_OFF=SamplesAll(:,:,(1:size(Samples_spike_Jelly_OFF,3)));
Mat2NlxSpike([ComputerDir,'\',Date,'\',MouseName,'\',char(FileDirNLX(ConditionFileOrder(4)).name),'\','TT',num2str(TetrodeNumber),'_s.ntt'], 0, 1, [], [1 1 1 1 1], Timestamps_spike_Jelly_OFF,ScNumbers_spike_Jelly_OFF, CellNumbers_spike_Jelly_OFF, Features_spike_Jelly_OFF, Samples_spike_Jelly_OFF, Header_spike_Jelly_OFF);
if PlotFigures==1
    Jelly_OFFFigure=figure;
for i=1:size(Header_spike_Jelly_OFF,1)
IsADBitVolts=strfind((char(Header_spike_Jelly_OFF{i,:})),'ADBitVolts');
if ~isempty(IsADBitVolts)
    ADVoltsList=Header_spike_Jelly_OFF{i,:};
   end %if
end %for
 ADVoltsList=strsplit(ADVoltsList,' ');
         for Electrode=2:5
        ADBitVolts(Electrode-1) = str2double(ADVoltsList{1,Electrode});
        Samples_spike_Jelly_OFF(:,Electrode-1,:)=Samples_spike_Jelly_OFF(:,Electrode-1,:)*(ADBitVolts(Electrode-1));
        plot(mean(Samples_spike_Jelly_OFF,3));
        title('Jelly OFF');
        end %for
end %if
        