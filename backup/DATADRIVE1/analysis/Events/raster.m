% Generate data
spikes=[1.44 4.36 9.94 15.32];
x=0:0.01:16;

% Set up XY coordinate matrix
% Row 1 = X coordinates
% Row 2 = Y coordinates
B=[x spikes;zeros(size(x)) ones(size(spikes))];

% Sort in ascending order of X coordinates
C=sortrows(B',1)';
% You can also SORT row 1 and reorder row 2 using the index you get from SORT

% Plot data
plot(C(1,:),C(2,:))