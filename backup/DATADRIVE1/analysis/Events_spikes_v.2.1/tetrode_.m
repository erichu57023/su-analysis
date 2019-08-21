classdef tetrode<handle
    %TETRODE Imports Neuralynx .NTT files and helps with clustering
    %   Current methods include:
    %
    %       - present_figure('atribute') e.g. obj.present_figure('height')
    %       will plot the height of all waveforms as a scatter plot on 3
    %       electrodes
    %       
    %       - show_cell(cell_number) e.g. obj.show_cell(1) will plot the
    %       waveforms on all 4 electrodes on cell 1 and obj.show_cell(nan)
    %       will simply plot all waveforms
    %
    %       - auto_cluster() e.g. obj.auto_cluster() will plot a 3D plot
    %       with the first 3 principle components of the entire dataset.
    %       The 3 axis represent the 3pca's, NOT 3 electrodes.
    %
    %       - peri_event() e.g. obj.peri_event(stamps); will plot a
    %       peri-event histogram of all units (and unclustered datapoints)
    %       around the timepoints in the array stamps. Note that stamps
    %       should be in 탎. You can also a window (in 탎) as a 3rd
    %       argument.
    %
    %   Tetrode uses part of Fieldtrip, (Copyright (C) 2005), Robert
    %   Oostenveld. For more information about Fieldtrip (including license
    %   information (GNU) see: http://www.fieldtriptoolbox.org.
    %
    %   Tetrode is part of Bearphys, Bearphys is made by Han de Jong,
    %   j.w.dejong@berkeley.edu
    
    
    properties
        % In principle all data properties are organized as follows:
        % (channels x data points x n)
        raw_data        % Stores imported tetrode (struct)
        data            % Data after manipulations (c,d,n)
        attributes      % Attributes of the sweeps in data (e.g. 'height')
        timestamps      % Timestamps of the waveforms in data (n)
        cells           % Clustered units(n);
        nr_cells        % The total number of clustered units
        handles         % Handles to figures and plots
        settings        % Struct with settings.
    end
    
    methods
        function obj = tetrode(filename)
            %TETRODE Construct an instance of this class
            %   Uses Fieldtrip m files to import a tetrode file
            
            obj.raw_data = read_neuralynx_ntt(filename, 1, inf);
            
            % Grab the data
            obj.data = obj.raw_data.dat;
            
            % Grab any pre-clustered cells
            obj.cells = obj.raw_data.CellNumber;
            
            % Find the total number of clusters
            obj.nr_cells = max(obj.cells);
            
            % Fill out the attributes
            obj.find_attributes;  
            
            % Fill out the timestamps
            obj.timestamps = obj.raw_data.TimeStamp;
            
            % Fill out the settings
            obj.settings.colorlist = ['k', 'g', 'r', 'b', 'm', 'y'];
            
        end
        
        function find_attributes(obj)
            % Find attributes updates all the sweep attributes
            
            for i=1:4
                obj.attributes(i).height = squeeze(max(obj.data(i,:,:)) - min(obj.data(i,:,:)));
                
                % run PCA
                [~ ,score,~ ,~ ,~ ] = pca(squeeze(obj.data(i,:,:))');
                obj.attributes(i).pca_1 = score(:,1);
                obj.attributes(i).pca_2 = score(:,2);
                obj.attributes(i).pca_3 = score(:,3);
                
            end
            
        end
        
        function outputArg = remove_noise(obj,inputArg)
            %REMOVE_NOISE Not written yet
            %   Detailed explanation goes here
            
            warning('Remove_noise method is not written yet.');
            
            outputArg = [];
             
        end
        
        function present_figure(obj, attribute, channel_list)
            %PRESENT_FIGURE presents a 3D scatterplot with the requested
            %atribute on 3 electrodes on 3 axis.
            % 
            %Examples:
            %
            %obj.present_figure('height') will plot the height on the first
            %3 electores.
            %
            %obj.present_figure('height', [1 3 4]) will plot the height on
            %electrode 1, 3 and 4.
            
            % Make a list of all attributes
            attribute_list = fieldnames(obj.attributes);
            
            % Colorlist for the markers
            colorlist = obj.settings.colorlist;
            
            % Check if the requested attribute is valid
            isvalid = false;
            for i=1:length(attribute_list)
                if strcmp(attribute_list{i},attribute)
                    isvalid=true;
                end
            end
            
            % If not valid, let the user know.
            if ~isvalid
                disp(' ')
                disp('Currently only the following attributes are supported:')
                for i=1:length(attribute_list)
                    disp(['    - ' attribute_list{i}])
                end
                return
            end
            
            % Grab the data
            for i=1:4
                XYZ(i,:) = getfield(obj.attributes(i),attribute);
            end
            
            % Figure out which three electrodes to show
            if nargin==3
                channel = channel_list;
            else
                channel = [1 2 3]; % show the first 3 electrodes
            end
            
            % Make the figure, plot unclustered datapoints
            figure
            
            % Plot identified clusters
            for i=0:obj.nr_cells
                indexer = obj.cells==i;
                
                % figure out the name of this dataset
                if i==0
                    dataname = 'Unclustered';
                else
                    dataname = ['Cell ' num2str(i)];
                end
                scatter3(XYZ(1,indexer), XYZ(2,indexer), XYZ(3,indexer),...
                    'DisplayName', dataname,...
                    'MarkerEdgeColor',colorlist(i+1));
                hold on
            end
            
            % Label figure and axis
            title(attribute)
            xlabel(['Electrode ' num2str(channel(1))])
            ylabel(['Electrode ' num2str(channel(2))])
            zlabel(['Electrode ' num2str(channel(3))])
            
        end
        
        function show_cell(obj, cell_number)
            %SHOW_CELL Presents the waveform of one cell on 4 electrodes
            
            % Does this cell exist?
            if cell_number>obj.nr_cells
                error(['This cell does not exist, there are/is only ' num2str(obj.nr_cells) ' cell(s) in this dataset.'])
            end
            
            % Make indexer for this cell
            if ~isnan(cell_number)
                indexer = obj.cells==cell_number;
            else % user put NaN, and want's all traces
                indexer = ones(size(obj.cells));
            end
            
            % If there are to many waveforms we will plot only 1% of them
            % to maintain performance
            cutoff = 10000;
            if sum(indexer)>cutoff
                new_indexer = zeros(size(indexer));
                new_indexer(1:100:end)=indexer(1:100:end);
                indexer = logical(new_indexer);
                disp(['Showing only 1% (n = ' num2str(sum(indexer==1)) ') of all waveforms.'])
            end
            
            % How to name the figure
            if isnan(cell_number)
                figure_name = 'All Waveforms';
            elseif cell_number==0
                figure_name = 'Unclustered';
            else
                figure_name = ['Cell ' num2str(cell_number)];
            end
          
            % Make a figure
            figure('Name',figure_name)
            for i=1:4
                subplot(2,2,i)
                plot(squeeze(obj.data(i,:,indexer)));
                title(['Electrode ' num2str(i)]);
            end
            
        end
        
        function [plot_handle, results] = peri_event(obj, stamps, window_size)
            % PERI_EVENT_PLOT_STAMPS present a peri-event plot of events
            % stamps. The stamps that will be used as the '0' timepoint are
            % in the variable 'stamps' (1D). All clusters as well as
            % non-clustered events will be ploted in different colors, they
            % are in the 3rd dimention in the results output.
            %
            % For now the results are also in 50 bins, the binsize is
            % dynamic
            %
            % Results is organized as follows: (trials, time bins, data
            % sets); There is a TIMELINE in the top row!
            %
            % Note that the stamps should the in the same unit of time
            % (e.g. us) as the data in the tetrode object.
            
            % Error handeling
            %       TO DO...  
            
            % CHECK IF DOUBLE AND UINT64 IS REALLY THE SAME
            
            % Colorlist for the markers
            colorlist = obj.settings.colorlist;
            
            % check the input arguments
            if strcmp(class(stamps),'uint64')
                stamps = double(stamps);
            end
            
            % Collect the data to be plotted
            units = cell(obj.nr_cells, 1);
            for i = 0:obj.nr_cells
                units{i+1} = double(obj.timestamps(obj.cells==i));      
            end
            
            % Error handeling on the intervals. Note that intervals 1000탎
            % are just probably some sort of error, but we should put out a
            % warning.
            intervals = stamps(2:end) - stamps(1:end-1);
            error_intervals = [intervals<1000, false];
            if sum(error_intervals)>0
                warning('Intervals <1000탎 between stamps not analyzed.')
                stamps = stamps(~error_intervals);
            end
            intervals = stamps(2:end) - stamps(1:end-1);
            
            % Figure out the appropriate window size
            if nargin==3 % user provided a window size
                window = window_size;
            else
                % Find the smalles inter-stamp interval
                window = round(0.5*(min(intervals)));
            end
            
            % Print the window size
            disp(['Window set to ' num2str(window) '탎'])
            
            % Figure out an appropriate bin size for the histogram
            bin_size = 2*window/50;
            timeline = round([-window:bin_size:window-bin_size]+0.5*bin_size,2);
            
            % Make the figure
            plot_handle = figure;
            subplot(2,1,1)
            
            % Make an empty results variable
            results = zeros(length(stamps) + 1, 50, length(units));
            
            % for every dataset
            for i = 1:length(units)
                
                % Add the timeline in the first row of the results
                results(1,:,i) = timeline;
                
                % for every trial
                % Start at row 2 because timeline in row 1
                for j = 2:length(stamps) + 1
                    
                    % Substract the trial from the input data to get the
                    % difference
                    temp_data = units{i} - stamps(j-1);

                    % for every bin
                    for k = 1:length(timeline)
                        
                        % Collect results
                        % note the biger-or-equal on one side, vs smaller
                        % on the other side
                        results(j,k,i) = sum(temp_data>=timeline(k)-0.5*bin_size & temp_data<timeline(k)+0.5*bin_size);
                    end
                    
                    % Plot the events
                    super_temp = temp_data(temp_data>= -window & temp_data< window);
                    plot(super_temp,ones(length(super_temp),1)*j,'.',...
                        'MarkerEdgeColor',colorlist(i));
                    hold on
                end
            end
            
            % Axis labeling etc.
            xlabel('Time (탎)')
            ylabel('Trial #')
            set(gca,'Ydir','reverse')
            xlim([-window, window])
            
            
            % Plot the histogram below
            subplot(2,1,2)
            
            % Collect the results
            results(2:end,:,:) = results(2:end,:,:)./bin_size;
            
            % Switch from 탎 to sec (we want the final results in Hz)
            results = results.*1000000;
            
            for i = 1:length(units)
                temp_results = results(2:end,:,i);
                sem_results=std(temp_results)./sqrt(length(stamps));
                fill([timeline';flipud(timeline')],[mean(temp_results)'-sem_results';flipud(mean(temp_results)'+sem_results')],[0 0 1],...
                    'linestyle','none',...
                    'FaceAlpha',0.1,...
                    'FaceColor',colorlist(i));
                hold on
                plot(timeline, mean(temp_results),'Color',colorlist(i)')
            end
            
            % Axis labeling etc.
            xlabel('Time (탎)')
            ylabel('Events (Hz)')
        end
        
        function output = auto_cluster(obj)
            %AUTO_CLUSTER does auto clustering on the basis of PCA
            %   Run a PCA on the data property and performs k-means
            %   clustering to possibly identify cells.
            
            %   The clustering does not work yet, but the PCA plot is
            %   informative.
            
            % Convert the data
            disp('Using only first 20 data points')
            new_data = [];
            for i=1:4
                new_data = [new_data; squeeze(obj.data(i,1:20,:))];
            end
            new_data = new_data';
            
            % Do the PCA
            [coeff,score,latent,tsquared,explained] = pca(new_data);
            
            % Print info
            disp(['First 3 coefficients explain ' num2str(sum(explained(1:3))) '% of the varation.'])
            
            % Make figure
            figure
            subplot(1,3,1)
            scatter(score(:,1), score(:,2));
            xlabel('1st Principal Component')
            ylabel('2nd Principal Component')
            subplot(1,3,2)
            scatter(score(:,1), score(:,3));
            xlabel('1st Principal Component')
            ylabel('3nd Principal Component')
            subplot(1,3,3)
            scatter(score(:,2), score(:,3));
            xlabel('2st Principal Component')
            ylabel('3nd Principal Component')
            
            % Plot scatter3
            figure
            scatter3(score(:,1),score(:,2),score(:,3))
            hold on
            
            % Plot identified clusters
            for i=1:obj.nr_cells
                indexer = obj.cells==i;
                scatter3(score(indexer,1), score(indexer,2), score(indexer,3))
            end
            
            % Print the output
            output = coeff;
            
        end
    end
end

