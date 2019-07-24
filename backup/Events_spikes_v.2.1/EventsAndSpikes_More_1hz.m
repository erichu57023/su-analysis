close all
for TetrodeNumber=[1 2 3 4]
[TimeStamps_events, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV('Events.nev', [1 1 1 1 1], 1, 1, []); %event file
[TimeStamps_cells, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(['TT',num2str(TetrodeNumber),'.ntt'], [1 1 1 1 1], 1, 1, [] ); %clustering file

%%
% TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_cells(1))/1000000; %for spikes
% TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_events(1))/1000000; %for laser pulses

if TimeStamps_cells(1) > TimeStamps_events(1)
    TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_events(1))/1000000;
    TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_events(1))/1000000;
else
    TimeStamps_evets_zeroed_s = (TimeStamps_events-TimeStamps_cells(1))/1000000;
    TimeStamps_cells_zeroed_s = (TimeStamps_cells-TimeStamps_cells(1))/1000000;
end

cell = cat(1,TimeStamps_cells_zeroed_s,CellNumbers); 
a = find(cell(2,:) ~=0); % "~=1" looks at cell#1, "~=2" looks at cell#2, etc.
cell(:,a) = [];
if isempty(cell)
    continue
end
laser = cat(1,TimeStamps_events_zeroed_s,TTLs);
laser(laser == 0) = NaN;

%%
eventcount = find(laser(2,:) == 1); %counts the number o laser pulses
eventcount = size(eventcount); 

for i = 1:eventcount(:,2)-1 %numbers the laser pulses from 1 to last
    laser(2,2+2*i) = 1 + i;
end

%% raster plot

f1 = figure;

for i = 1:eventcount(:,2)
    X(i,:) = [(laser(1,2*i)-0.2):0.0001:(laser(1,2*i)+0.2)]; %sets the range around stimulation -0.1 s to +0.5 s 
    sizeX = size(X(1,:)); %check the size of Xi array
    Y(i,(1:sizeX(:,2))) = i; %creates an Yi array for different trials Y1 - first trial Y2 - second trial etc.
    
    [v location_laser] = min(abs(X(i,:)-laser(1,2*i))); %finds the index of laser pulse time
    
    a = find(cell(1,:) >= laser(1,2*i)-0.05 & cell(1,:) <= laser(1,2*i)+0.1); %finds the indices around the laser pulse (50ms before and 100 ms after)
    TF = isempty(a); %checks if it finds any values
    if TF == 0        
        numberofspikes = size(a);
        location_cell = zeros(1,numberofspikes(:,2));
        
        for j = 1:numberofspikes(:,2)
            [v location_cell(1,j)] = min(abs(X(i,:)-cell(1,a(j))));
        end
        
        X(i,:) = NaN; %changes all of the that are not event to NaN    
        X(i,location_laser) = 0; %makes the pulse location 0 
        
        for h = 1:numberofspikes(:,2); %inserts time values at the spike times
            X(i,location_cell(1,h)) = cell(1,a(h))-laser(1,2*i);
        end
    end
    
    if TF == 1
        X(i,:) = NaN; %changes all of the that are not event to NaN    
        X(i,location_laser) = 0; %makes the pulse location 0 
    end
      
    s(i) = scatter(X(i,:),Y(i,:));
        hold on
end
try
set(s,'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k')
ax = gca;
ax.XTick = -0.05:0.01:0.1;
ax.YLim = [0 25];
f1.InnerPosition = [680 558 800 350];
catch
    continue
end

%% histogram

f2 = figure;

k = find(isnan(X) == 1);
X(k) = 0;
k = find(X);
histo = X(k);
edges = [-0.05:0.001:0.1]; 
histogram(histo,edges,'FaceColor','k','FaceAlpha',1)

ax = gca;
ax.XTick = -0.05:0.01:0.1;
ax.XLim = [-0.05 0.1];
ax.YLim = [0 30];
f2.InnerPosition = [680 558 800 350];

%% average response time

count_real = zeros(1,99);

for i = 1:99 %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
    a = find(histo(:,1) >= i/1000 - 0.001 & histo(:,1) <= i/1000 + 0.001);
    count_real(i) = size(a,1);
end

[Mr,Ir] = max(count_real); %M = highest number of spikes in interval, I = index of M
highestresponse = Mr;

bin_low = (Ir-1)/1000;
bin_up = (Ir+1)/1000;

a = find(histo(:,1) >= bin_low & histo(:,1) <= bin_up);
response_latency = mean(histo(a)); %gives average response latency

a = find(histo >= 0);
histo_pos = histo(a);

k = size(histo_pos,1);

bootstrap = zeros(1,10000);
for i = 1:10000;
    r = (0.1-0).*rand(size(histo_pos,1),1) + 0;
    count_shuffle = zeros(1,99);
        for j = 1:99 %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
            a = find(r(:,1) >= j/1000 - 0.001 & r(:,1) <= j/1000 + 0.001);
            count_shuffle(j) = size(a,1);
        end
    [Ms,Is] = max(count_shuffle);
    bootstrap(i) = Ms;
end

 
bootstrap = transpose(bootstrap);
pd = fitdist(bootstrap,'Normal');
x_values = 0:0.1:40;
y = pdf(pd,x_values);

f3 = figure
yyaxis left
histogram(bootstrap)
hold on
yyaxis right
plot(x_values,y,'LineWidth',2)

percentile = icdf(pd,0.999);
end

%% Example for one laser pulse
    
% lasertime = laser(1,14);
% k = find(cell(1,:) >= laser(1,14)-0.05 & cell(1,:) <= laser(1,14)+0.95)
% 
% X = [(laser(1,14)-0.05):0.000001:(laser(1,14)+0.95)];
% sizeX = size(X(1,:));
% Y(1:sizeX(:,2)) = 1;
% [v location_l] = min(abs(X(1,:)-laser(1,14)));
% [v location_c1] = min(abs(X(1,:)-cell1(1,3)));
% [v location_c2] = min(abs(X(1,:)-cell1(1,4)));
% [v location_c3] = min(abs(X(1,:)-cell1(1,5)));
% 
% X(1,:) = NaN;
% X(1,location_l) = 0;
% X(1,location_c1) = cell1(1,3)-laser(1,14);
% X(1,location_c2) = cell1(1,4)-laser(1,14);
% X(1,location_c3) = cell1(1,5)-laser(1,14);
% scatter(X,Y)
