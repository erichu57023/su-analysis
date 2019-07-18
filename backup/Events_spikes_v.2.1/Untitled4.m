clear all
clc
[Timestamps, ScNumbers, CellNumbers, Features, Samples_spike_Jelly_Exposure, Header] = Nlx2MatSpike(['TT6.ntt'], [1 1 1 1 1], 1, 1, [] ); 
% for i=1:size(Header,1)
% IsADBitVolts=strfind((char(Header{i,:})),'ADBitVolts');
% if ~isempty(IsADBitVolts)
%     ADVoltsList=Header{i,:};
%    end %if
% end %for
%  ADVoltsList=strsplit(ADVoltsList,' ');
%         for Electrode=2:5
%         ADBitVolts(Electrode-1) = str2double(ADVoltsList{1,Electrode});
%         Samples_spike_Jelly_Exposure(:,Electrode-1,:)=Samples_spike_Jelly_Exposure(:,Electrode-1,:)*(ADBitVolts(Electrode-1));
%        end %for
Mat2NlxSpike('TT5_export3.ntt', 0, 1, [], [1 1 1 1 1], Timestamps,ScNumbers, CellNumbers, Features, Samples_spike_Jelly_Exposure, Header);
