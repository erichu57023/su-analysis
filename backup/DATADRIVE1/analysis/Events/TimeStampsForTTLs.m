[TimeStamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV('35ms1Hz50s.nev', [1 1 1 1 1], 1, 1, []);
TTLs(TTLs == 0) = NaN; %"deletes" 0 values so they wouldn't be plotted, 1 inidicate the TTL's start
TimeStamps_zeroed_s = (TimeStamps-TimeStamps(1))/1000000; %zeroes the timestamps and converts them to seconds
events = [TimeStamps_zeroed_s; TTLs];

eventcount = find(events(2,:) == 1); %counts the number of laser pulses
eventcount = size(eventcount); 
eventcount(:,2);fi([], 1, 0)

for i = 1:eventcount(:,2)-1 %numbers the laser pulses from 1 to last
    events(2,2+2*i) = 1 + i;
end

for i = 1:eventcount(:,2)
    X(i,:) = [(events(1,2*i)-0.05):0.0001:(events(1,2*i)+0.1)]; %sets the range around stimulation -0.05 s to +0.1 s 
    sizeX = size(X(1,:)); %check the size of Xi array
    Y(i,(1:sizeX(:,2))) = i; %creates an Yi array for different trials Y1 - first trial Y2 - second trial etc.
    [v location] = min(abs(X(i,:)-events(1,2*i))); %finds the index of laser pulse time
    X(i,(X(i,:) ~= X(i,location))) = NaN; %changes all of the that are not event to NaN
    X(i,(location)) = X(i,(location)) - X(i,(location)); %makes the pulse location 0
    scatter(X(i,:),Y(i,:))
    hold on
    axis([-50 100 0 i]) %sets the axis -50ms to +100ms
end

% X = [(events(1,2)-0.05):0.0001:(events(1,2)+0.1)]; %example for one event
% sizeX = size(X(1,:));
% Y(1:sizeX(:,2)) = 1;
% [v location] = min(abs(X(1,:)-events(1,2)));
% X(1,(X(1,:) ~= X(1,location))) = NaN;
% X(1,(location)) = X(1,(location)) - X(1,(location));
% scatter(X,Y)
% axis([-50 100 0 1])    

