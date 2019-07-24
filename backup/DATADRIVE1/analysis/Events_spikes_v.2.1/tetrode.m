classdef tetrode<handle
    %TETRODE Imports Neuralynx .NTT files and helps with clustering
    %   Current methods include:
    %
    %       - remove_noise('noise type') e.g. obj.remove_noise('clipping');
    %       will remove all waveforms that clip the maximum value from the
    %       dataset. obj.remove_noise('max height', 100); will remove all
    %       waveforms with a height of >100µV.
    %
    %       - present_figure('atribute') e.g. obj.present_figure('height')
    %       will plot the height of all waveforms as a scatter plot on 3
    %       electrodes. You can also provide the electrodes you want to see
    %       like this: obj.present_figure('height',[1 2 4]); (For electrode
    %       1, 2 and 4.
    %       
    %       - show_cell(cell_number) e.g. obj.show_cell(1) will plot the
    %       waveforms on all 4 electrodes on cell 1 and obj.show_cell(nan)
    %       will simply plot all waveforms.
    %
    %       - full_pca() e.g. obj.full_pca() will plot a 3D plot
    %       with the first 3 principle components of the entire dataset.
    %       The 3 axis represent the 3pc's, NOT 3 electrodes.
    %
    %       - peri_event() e.g. obj.peri_event(stamps); will plot a
    %       peri-event histogram of all units (and unclustered datapoints)
    %       around the timepoints in the array stamps. Note that stamps
    %       should be in µs. You can also a window (in µs) as a 3rd
    %       argument.
    %
    %       - set_interval() e.g. obj.set_interval([0 20]); from now on the
    %       attributes (e.g. 'height', 'pca_1', etc..) are calulated using
    %       the first 20 datapoints of each waveform.
    %
    %       - save_data('filename') e.g. obj.save_data('filename') will
    %       save the tetrode object to the current folder.
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
            
            % Fill out the timestamps
            obj.timestamps = obj.raw_data.TimeStamp;
            
            % Find the recording interval
            obj.settings.recording_interval =...
                [obj.raw_data.hdr.FirstTimeStamp, obj.raw_data.hdr.LastTimeStamp];
            
            % Fill out the settings
            %obj.settings.interval = [1, size(obj.raw_data.dat,2)];
            obj.settings.interval = [1, 20];
            obj.settings.preferred_electrodes = [1 2 3]; % Because humans like three dimensions
            obj.settings.colorlist = ['kgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbm'];
            
            % Find the DBitVolts (from the header char...)
            header = obj.raw_data.hdr.Header;
            for i = 1:5000
                if strcmp('ADBitVolts', header(i:i+9))
                    break;
                end
            end
            obj.settings.ADBitVolts(1) = str2num(header(i+11:i+36));
            obj.settings.ADBitVolts(2) = str2num(header(i+38:i+63));
            obj.settings.ADBitVolts(3) = str2num(header(i+65:i+90));
            obj.settings.ADBitVolts(4) = str2num(header(i+92:i+117));
            
            % Figure out which electrodes are not recording anything at all
            % and label them 'false' in the working electrodes property
            for i=1:4
                if sum(sum(obj.data(i,:,:)))==0
                    obj.settings.working_electrodes(i) = false;
                else
                    obj.settings.working_electrodes(i) = true;
                end
            end
            
            % Find the Inputrange(from the header char...)
            for i = 1:5000
                if strcmp('InputRange', header(i:i+9)); break; end
            end
            new_char = header(i+11:i+30);
            for i=1:length(header)
                if strcmp(new_char(i),'-'); break; end
            end
            obj.settings.InputRange = str2num(new_char(1:i-1));
            
            % Find the Threshold Value(from the header char...)
            for i = 1:5000
                if strcmp('ThreshVal', header(i:i+8)); break; end
            end
            new_char = header(i+10:i+30);
            for i=1:length(header)
                if strcmp(new_char(i),'-'); break; end
            end
            obj.settings.ThresVal = str2num(new_char(1:i-1));
           
            % Fill out the attributes
            obj.find_attributes;
            
            % Set the cutoff for display performance issues
            obj.settings.disp_cutoff = 40000;
            
            % Figure out the prefered electrodes (which the most variation in height)
            for i=1:4
                devs(i) = std(obj.attributes(i).height);
            end
            [~, show_electrodes] = sort(devs, 'descend');
            show_electrodes = sort(show_electrodes(1:3));
            obj.settings.preferred_electrodes = show_electrodes;
            
            % That's it
        end
        
        function find_attributes(obj)
            % Find attributes updates all the sweep attributes
            
            % Analyse the following interval
            start = obj.settings.interval(1);
            stop = obj.settings.interval(2);
            
            % Grab the data from the 4 electrodes
            for i=1:4
                temp_data = squeeze(obj.data(i,start:stop,:)).*obj.settings.ADBitVolts(i)*10^6;
                
                % grab the peak
                obj.attributes(i).peak = max(temp_data);
                
                % grab the valley
                obj.attributes(i).valley = min(temp_data);
                
                % grab the height
                obj.attributes(i).height = obj.attributes(i).peak - obj.attributes(i).valley;
                
                % grab the energy
                obj.attributes(i).energy = sum(temp_data.^2);
                
                % run PCA, but only on 'working electrodes'
                if obj.settings.working_electrodes(i)
                    [~ ,score,~ ,~ ,~ ] = pca(squeeze(temp_data)');
                    obj.attributes(i).pca_1 = score(:,1)';
                    obj.attributes(i).pca_2 = score(:,2)';
                    obj.attributes(i).pca_3 = score(:,3)';
                else
                    obj.attributes(i).pca_1 = zeros(1,size(temp_data,2));
                    obj.attributes(i).pca_2 = zeros(1,size(temp_data,2));
                    obj.attributes(i).pca_3 = zeros(1,size(temp_data,2));
                end
                
            end
            
        end
        
        function removed = remove_noise(obj, varargin)
            %REMOVE_NOISE Removes noise on the basis of the arguments given
            %in varargin
            %
            %   Examples:
            %   >>remove_noise('max heigh', 100); Will remove all waveforms
            %   of which the height is bigger then 100µV on ANY electrode.
            %
            %   >>remove_noise('min height', 50); Will remove all waveforms
            %   of which the height is smaller then 100µV on EVERY
            %   electrode.
            %
            %   >>remove_noise('indexer',[true false false true]); will
            %   remove the 2th and 3rd waveform from the dataset. (In
            %   reality the logical will be much longer of course).
            
            % The output is boolean on wether or not any waveforms were
            % actually removed
            removed = false;
            
            % Deal with input arguments
            for i=1:2:length(varargin)
                
                switch varargin{i}
                    
                    % max height
                    case 'max height'
                            indexer = obj.attributes(1).height> varargin{i+1}...
                                | obj.attributes(2).height > varargin{i+1}...
                                | obj.attributes(3).height > varargin{i+1}...
                                | obj.attributes(4).height > varargin{i+1};
                            break;
                            
                    % max height
                    case 'min height'
                            indexer = obj.attributes(1).height< varargin{i+1}...
                                & obj.attributes(2).height < varargin{i+1}...
                                & obj.attributes(3).height < varargin{i+1}...
                                & obj.attributes(4).height < varargin{i+1};
                            break;
                    
                    % Remove isolated cluster        
                    case 'remove cluster'
                        indexer = obj.cells == varargin{i+1};
                        indexer2 = obj.cells>varargin{i+1};
                        obj.cells(indexer2) = obj.cells(indexer2)-1;
                        varargin{i+1} = 'skipp';
                        warning('Not a finished method, clicking no in the dialog will result in errors')
                    
                    % Remove all waveforms that clip the max value    
                    case 'clipping'
                        clip_value = obj.raw_data.hdr.ADMaxValue;
                        indexer = obj.data==clip_value | obj.data==-1*clip_value;
                        indexer = logical(squeeze(sum(sum(indexer))));
                    
                    % Remove all waveforms indexed in the following logical    
                    case 'indexer'
                        indexer = varargin{i+1};
                        varargin{i+1} = 'skipp';
                    
                    % Skipp this input argument    
                    case 'skipp'
                        X = 5; % Do nothing.
                        
                    otherwise
                        error('Unknown argument type, juse the help method to display example option.')
                end

            end
            
            % Check if the indexer even picked up on any waveforms
            if sum(indexer)==0
                disp('No waveforms met exclusion criteria.')
                return
            end
            
            % Show the removed noise
            [figure_1, figure_2] = obj.show_cell(indexer,'hypothetical');

            % Ask the user if they want to indeed remove these
            % waveforms?
            w_percentage = 100*sum(indexer)/length(obj.timestamps);
            input = questdlg(['Remove these ' num2str(sum(indexer)) '(' num2str(w_percentage) '% of total) waveforms?'],'Remove?','yes','no','no');
            if strcmp(input,'yes')
                % Invert the indexer to waveforms we want to keep
                indexer = ~indexer;
                
                % Just remove all these wavefroms (can always get them back
                % from raw data)
                obj.data = obj.data(:,:,indexer);
                obj.timestamps = obj.timestamps(indexer);
                obj.cells = obj.cells(indexer);
                obj.find_attributes;
                
                % Change the bool
                removed = true;
            end
            
            % Close the figures
            close(figure_1)
            close(figure_2)
            
        end
        
        function [waveform_figure, scatter_figure] = present_figure(obj, attribute, channel_list)
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
            
            
            % If the user did not give an attribute to plot, just plot
            % height
            if nargin==1
                attribute = 'height';
            end
            
            
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
                channel = obj.settings.preferred_electrodes;
            end
            
            
            % Make indexers for every cells
            for i=0:obj.nr_cells
                indexer{i+1} = obj.cells==i;
            end
            
            
            % If there are to many waveforms we will plot only 1% of them
            % to maintain performance
            cutoff = obj.settings.disp_cutoff;
            for i=1:length(indexer)
                if sum(indexer{i})>cutoff
                    new_indexer = zeros(size(indexer{i}));
                    new_indexer(1:20:end)=indexer{i}(1:20:end);
                    indexer{i} = logical(new_indexer);
                    
                    % Plot the waveforms
                    if i==1
                        disp(['Showing only 5% (n = ' num2str(sum(indexer{i}==1)) ') of unclustered waveforms/datapoints.'])
                    else
                        disp(['Showing only 5% (n = ' num2str(sum(indexer{i}==1)) ') of cell ' num2str(i-1) ' waveforms/datapoins.'])
                    end
                end
            end
            
            
            % Make the figure, plot unclustered datapoints
            scatter_figure = figure;
            scatter_figure.Position = [432, 300, 560, 420];
            
            % Plot identified clusters
            for i=0:obj.nr_cells
                
                % figure out the name of this dataset
                if i==0
                    dataname{i+1} = 'Unclustered';
                else
                    dataname{i+1} = ['Cell ' num2str(i)];
                end
                scatter3(XYZ(channel(1),indexer{i+1}), XYZ(channel(2),indexer{i+1}), XYZ(channel(3),indexer{i+1}),...
                    'DisplayName', dataname{i+1},...
                    'MarkerEdgeColor',colorlist(i+1));
                hold on
            end
            
            % Label figure and axis
            title(attribute)
            xlabel(['Electrode ' num2str(channel(1))])
            ylabel(['Electrode ' num2str(channel(2))])
            zlabel(['Electrode ' num2str(channel(3))])
            
            % ***** Plot the waveforms *****
            
            % Make a figure
            waveform_figure = figure('Name','Waveforms');
            waveform_figure.Position = [993 300, 560, 420];
            
            % Work on the xaxis
            interval = 10^6/obj.raw_data.hdr.SamplingFrequency; %µs
            timeline = [0:interval:1000-interval]';
            
            % 4 subplots for 4 electrodes
            for i=1:4
                subplot(2,2,i)
                
                % Plot all cells and unclustered
                for j=0:obj.nr_cells
                    % Plot ydata in µv
                    ydata = squeeze(obj.data(i,:,indexer{j+1}).*obj.settings.ADBitVolts(i))*10^6;
                    
                    % In the special case that there is only one waveform, make
                    % sure that Y is still oriented in the right way
                    if sum(indexer{j+1})==1
                        ydata = ydata';
                    end
                    
                    % Grab the Xdata in µs
                    xdata = repmat(timeline,1, size(ydata,2),1);
                    
                    % Convert the data for faster performance
                    xdata = [xdata; nan(1,size(xdata,2))];
                    xdata = xdata(:);
                    ydata = [ydata; nan(1,size(ydata,2))];
                    ydata = ydata(:);
                    
                    % Plot the data
                    plot(xdata, ydata, 'Color',colorlist(j+1),...
                        'DisplayName',dataname{j+1});
                    hold on
                end
                ylabel('Voltage (µv)')
                xlabel('Time (µs)')
                ylim([-obj.settings.InputRange(i) obj.settings.InputRange(i)])
                title(['Electrode ' num2str(i)]);
                
                % Plot the threshold line
                hold on
                line(get(gca,'XLim'), [obj.settings.ThresVal(i), obj.settings.ThresVal(i)],...
                    'Color',[0.2 0.2 0.2],...
                    'LineStyle','--',...
                    'DisplayName','Threshold');
            end
            
        end
        
        function [waveform_figure, histogram_figure] = show_cell(obj, cell_number, varargin)
            %SHOW_CELL Presents the waveform of one cell on 4 electrodes
            
            % Work on the flag
            plot_average = false;
            hypothetical_cell = false;
            disp_fraction = 1;
            for i = 1:length(varargin)
                switch varargin{i}
                    
                    % Skip this argument
                    case 'skipp'
                        hoihoi = 5;
                        
                    % We will plot the average and stdev instead of all
                    % waveforms   
                    case 'average' 
                        plot_average = true;
                    
                    % Plot a hypothetical cell, not an identified cluster
                    case 'hypothetical'
                        hypothetical_cell = true;
                        
                    % Plot only a fraction of the waveforms    
                    case 'fraction'
                        disp_fraction = round(100/varargin{i+1});
                        varargin{i+1} = 'skipp';
                        
                    otherwise
                        error('Unknown Argument')
                end
            end
            
            
            % Real unit or a hypothetical?
            if ~ hypothetical_cell
                % Does this cell exist?
                if cell_number>obj.nr_cells
                    error(['This cell does not exist, there are/is only ' num2str(obj.nr_cells) ' cell(s) in this dataset.'])
                end

                % Make indexer for this cell
                if ~isnan(cell_number)
                    indexer = obj.cells==cell_number;
                else % user put NaN, and want's all traces
                    indexer = ~isnan(obj.cells); %nan cells are noise
                end
            else
                indexer = cell_number;
                cell_number = 0;
            end
            
            
            % If there are to many waveforms we will plot only 5% of them
            % to maintain performance
            cutoff = obj.settings.disp_cutoff;
            original_indexer = indexer; % To be used for histogram for instance
            if sum(indexer)>cutoff && disp_fraction ==1
                disp_fraction = 20;
                disp(['Showing only 5% (n = ' num2str(sum(indexer==1)) ') of all waveforms.'])
            end
            
            
            % This is to plot only a certain fraction of all waveforms,
            % either because there are more waveforms than the cutoff or
            % because the user specifically requested it.
            new_indexer = zeros(size(indexer));
            new_indexer(1:disp_fraction:end)=indexer(1:disp_fraction:end);
            indexer = logical(new_indexer);
            
            
            % How to name the figure
            if isnan(cell_number)
                figure_name = 'All Waveforms';
            elseif cell_number==0
                figure_name = 'Unclustered';
            else
                figure_name = ['Cell ' num2str(cell_number)];
            end
            
            % Work on the xaxis
            interval = 10^6/obj.raw_data.hdr.SamplingFrequency; %µs
            timeline = [0:interval:1000-interval]';
            
            % Figure out the cell color
            cell_color = obj.settings.colorlist(cell_number+1);
            
            % Make a figure and plot the waveforms
            waveform_figure = figure('Name',figure_name,...
                'Position',[232, 300, 800, 600]);
            for i=1:4
                subplot(2,2,i)
                ydata = squeeze(obj.data(i,:,indexer)).*obj.settings.ADBitVolts(i)*10^6; % Data from bitvalues to µV
                
                % In the special case that there is only one waveform, make
                % sure that Y is still oriented in the right way
                if sum(indexer)==1
                    ydata = ydata';
                end
                
                % Plot average instead of all sweeps?
                if plot_average
                    stdev_y = std(ydata');
                    ydata = mean(ydata,2);
                    
                    % Plot the stdev shadow
                    fill([timeline', flipud(timeline)'], [ydata'-stdev_y, fliplr(ydata'+stdev_y)],[1 0 0],...
                        'FaceColor',cell_color,...
                        'linestyle','none',...
                        'FaceAlpha',0.1)
                    hold on   
                end
                
                % Plot the data
                xdata = repmat(timeline,1, size(ydata,2),1); % Timeline in µs
                
                % Convert the data for faster performance
                xdata = [xdata; nan(1,size(xdata,2))];
                xdata = xdata(:);
                ydata = [ydata; nan(1,size(ydata,2))];
                ydata = ydata(:);
                
                % Plot the data
                plot(xdata, ydata, 'Color', cell_color);
                hold on
                xlabel('Time (µs)')
                ylabel('Voltage (µv)')
                %ylim([-obj.settings.InputRange(i) obj.settings.InputRange(i)])
                title(['Electrode ' num2str(i)]);
                
                % Plot the threshold line
                line(get(gca,'XLim'), [obj.settings.ThresVal(i), obj.settings.ThresVal(i)],...
                    'Color',[0.2 0.2 0.2],...
                    'LineStyle','--');
            end
            
            % IMPORTANT NOTE: all analysis bellow this point should be done
            % using the original_indexer. The 'indexer' is only used for
            % display purpouses.
            
            % Work on the inter-spike histogram
            histogram_figure = figure('Name','cluster data',...
                'Position',[1033 300, 800, 600]);
            m_timestamps = obj.timestamps(original_indexer);
            intervals = double(m_timestamps(2:end)-m_timestamps(1:end-1));
            subplot(2,2,1);
            histogram(intervals,5000);
            set(gca,'XScale','log')
            xlabel('Interval (µs)')
            ylabel('Count #')
            title('ISI')
            grab_lim = xlim();
            xlim([10^2 grab_lim(2)]);
            line([1000 1000],ylim(),'Color',[1 0 0])
            
            
            % Plot the waveform count throughout the session (to check if
            % stable firing or appearing/disappearing througout the session.
            subplot(2,2,2);
            session_time = obj.settings.recording_interval;
            m_timestamps = double(m_timestamps - session_time(1)).*10^-6;
            session_time = double(session_time - session_time(1)).*10^-6;
            histogram(m_timestamps, 100)
            xlim(session_time)
            xlabel('Time (s)')
            ylabel('Count #')
            title('Firing Distribution')
     
            
            % Calculate the Mahalanobis distance between the center of this
            % cluster and every other spike and compare this to every other
            % spike in the tetrode object for good measure.
            start_index = obj.settings.interval(1);
            end_index = obj.settings.interval(2);
            ydata = [];
            for i=1:4
                % don't include broken electrodes
                if ~obj.settings.working_electrodes(i)
                    continue;
                end
                
                % grab all the ydata
                ydata_temp = squeeze(obj.data(i,start_index:end_index,:)).*obj.settings.ADBitVolts(i)*10^6; % Data from bitvalues to µV
                ydata = [ydata, ydata_temp']; % Matlab Mahalanobis function wants n x m.
            end
            
            % If there are enough sweeps we can plot the Mahalanobis
            % distance
            if size(ydata(original_indexer,:),1)>size(ydata(original_indexer,:),2)
                
                % Calculate Mahalanobis distance
                d2 = mahal(ydata, ydata(original_indexer,:));
                d2_cell = d2(original_indexer);
                d2_not_cell = d2(~original_indexer);

                % We have to trow out outliers because otherwise the histograms
                % are just un-interpetable. Note that we don't do this for the
                % d2_cell, because in that case the hole point is that we want
                % to see the outliers
                d2_not_cell = d2_not_cell(~isoutlier(d2_not_cell));

                % Plot the histograms
                subplot(2,2,3);
                h1 = histogram(d2_cell,...
                    'Normalization','probability',...
                    'DisplayName','Selected Unit');
                hold on
                h2 = histogram(d2_not_cell,...
                    'Normalization','probability',...
                    'DisplayName','Other Waveforms');

                xlim([0 max([mean(d2_cell)+3*std(d2_cell), mean(d2_not_cell)+3*std(d2_not_cell)])])
                xlabel('Mahalanobis distance')
                ylabel('probability')
                legend();
                title('Distance from unit')
                
                % Figure out which one has the smaller binsize and make
                % them similar
                binWidth = min([h1.BinWidth, h2.BinWidth]);
                h1.BinWidth = binWidth;
                h2.BinWidth = binWidth;
            end
            
            % And finally, we also plot peak hight througout the
            % session for every electrode.
            subplot(2,2,4)
            %smooth_factor = round(0.05 * sum(original_indexer));
            for i=1:4
                plot(m_timestamps, obj.attributes(i).peak(original_indexer),'.',...
                    'DisplayName',['Electrode ' num2str(i)])
                hold on
            end
            xlim(session_time)
            xlabel('Time (s)')
            ylabel('Peak (µV)')
            title(['Peak Stability'])
            legend();
        end
        
        function [plot_handle, results] = peri_event(obj, stamps, varargin)
            % PERI_EVENT_PLOT_STAMPS present a peri-event plot of events
            % stamps. The stamps that will be used as the '0' timepoint are
            % in the variable 'stamps' (1D). All clusters as well as
            % non-clustered events will be ploted in different colors, they
            % are in the 3rd dimention in the results output.
            %
            % For now the results are also in 50 bins, the binwith is
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
            
            % Deal with input arguments
            window = []; % Plot window
            cell_nr = []; % Cells to be analysed
            for i = 1:length(varargin)
                switch varargin{i}
                    
                    % Skipp this argument
                    case 'skipp'
                        skippy = 5;
                        
                    % Set the windowsize
                    case 'window'
                        window = varagin{i+1};
                        disp(['Window set to ' num2str(window) 'µs'])
                        varargin{i+1} = 'skipp';
                    
                    % Only a specific cell
                    case 'cell'
                        cell_nr = varargin{i+1};
                        if ~isempty(cell_nr)
                            disp(['Only plotting cell ' num2str(cell_nr)]);
                        else
                            disp(['Plotting all cells.'])
                        end
                        varargin{i+1} = 'skipp';
                        
                    otherwise
                        warning('Unknown argument,... ignored.')
                end
            end
            
            % Colorlist for the markers
            colorlist = obj.settings.colorlist;
            
            % check the input arguments
            if strcmp(class(stamps),'uint64')
                stamps = double(stamps);
            end
            
            % Collect the data to be plotted
            if isempty(cell_nr) % Do all cells
                cell_nr = [0:obj.nr_cells];
            end
            units = cell(length(cell_nr), 1);
            for i = 1:length(cell_nr)
                units{i} = double(obj.timestamps(obj.cells==cell_nr(i)));      
            end
            
            % Error handeling on the intervals. Note that intervals 1000µs
            % are just probably some sort of error, but we should put out a
            % warning.
            intervals = stamps(2:end) - stamps(1:end-1);
            error_intervals = [intervals<1000, false];
            if sum(error_intervals)>0
                warning('Intervals <1000µs between stamps not analyzed.')
                stamps = stamps(~error_intervals);
            end
            intervals = stamps(2:end) - stamps(1:end-1);
            
            % Figure out the appropriate window size
            if isempty(window)
                % Find the smalles inter-stamp interval
                window = round(0.5*(min(intervals)));
            end
            
            % Print the window size
            disp(['Window set to ' num2str(window) 'µs'])
            
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
            xlabel('Time (µs)')
            ylabel('Trial #')
            set(gca,'Ydir','reverse')
            xlim([-window, window])
            
            
            % Plot the histogram below
            subplot(2,1,2)
            
            % Collect the results
            results(2:end,:,:) = results(2:end,:,:)./bin_size;
            
            % Switch from µs to sec (we want the final results in Hz)
            results = results.*1000000;
            
            % Calculate mean and SEM and plot
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
            xlabel('Time (µs)')
            ylabel('Events (Hz)')
        end
        
        function set_interval(obj, interval)
            %SET_INTERVAL sets the interval that is used to calculate the
            %attributes (for instance only the first 20 datapoints: [0 20];
            
            %   Error handeling
            if length(interval)~=2
                error('Please provide an interval like this [start stop] (in datapoints)')
            end
            
            % Set the interval
            obj.settings.interval = interval;
            
            % Force recalculation of attributes
            obj.find_attributes;
            
        end
        
        function removed = set_offline_threshold(obj, threshold)
            % Applies an offline threshold to each electrode seperataly.
            % This threshold can not be lower than what was used to collect
            % the original tetrode file. The method will allow it, but
            % nothing will happen. The output 'removed' is a bool on wheter
            % or not any waveforms were actually removed from the data.
            %
            % Example: >> obj.set_offlinse_threshold([65 65 65 65]); Will
            % set a threshold of 65µV to each electrode and delete all
            % waveforms that don't make this threshold on any electrode.
            
            % Convert the data to µV
            m_data = obj.data;
            for i=1:4
                m_data(i,:,:) = m_data(i,:,:).*obj.settings.ADBitVolts(i)*10^6;
            end
            
            % See if any sweeps do not meet the new threshold
            indexer = true(size(m_data,3),1);
            for i = 1:4
                indexer(max(squeeze(m_data(i,:,:)))>=threshold(i)) = false;
            end
            
            % Indexer is now true for sweeps that do not meet the threshold
            % on any electrorde
            
            removed = obj.remove_noise('indexer', indexer);
            
            if removed
                obj.settings.ThresVal = threshold;
            end  
        end
        
        function [coeff, score, main_figure] = full_pca(obj)
            %FULL_PCA does a PCA over all data, independend of electrode
            %   First concatenates all electrodes into one long waveform
            %   for every sample, then performs pca over this data. This
            %   analylsis ignores the individual electrodes after and can
            %   even be perfomed on single-electrode recordings.
            
            % Find the interval
            start = obj.settings.interval(1);
            stop = obj.settings.interval(2);
            disp(['Using the datapoint between ' num2str(start) ' and ' num2str(stop) '.']);
            
            % Concatenate the data from the individual electrodes
            new_data = [];
            for i=1:4
                % don't include broken electrodes
                if ~obj.settings.working_electrodes(i)
                    continue;
                end
                
                new_data = [new_data, squeeze(obj.data(i,start:stop,:))'];
            end
            
            % Do the PCA
            [coeff, score, ~, ~, explained] = pca(new_data);
            
            % Print info
            disp(['First 3 coefficients explain ' num2str(sum(explained(1:3))) '% of the varation.'])
            
            % Make figure, print the waveforms in PCA space
            main_figure = figure();
            main_figure.Position(3:4) = [1200, 400];
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
            xlabel('1st Principal Component')
            ylabel('2nd Principal Component')
            zlabel('3nd Principal Component')
            
            % Plot identified clusters
            for i=1:obj.nr_cells
                indexer = obj.cells==i;
                scatter3(score(indexer,1), score(indexer,2), score(indexer,3))
            end 
        end
        
        function save_data(obj, filename)
            %SAVE will save the object to a .mat file with the name
            %'filename'.
            
            variable_name = genvarname(filename);
            eval([variable_name '=obj;']);
            
            
            save(filename, variable_name);
            
            
        end
        
        function remove_cell(obj, cell_number)
            %REMOVE_CELL labeles datapoints in that unit as unclustered
            
            % Grab the indexer
            indexer = obj.cells == cell_number;
            indexer2 = obj.cells>cell_number;
            
            % Uncluster the datapoints
            obj.cells(indexer) = 0;
            
            % Decrease all other cell numbers
            obj.cells(indexer2) = obj.cells(indexer2)-1;
            obj.nr_cells = obj.nr_cells-1;

            
        end
        
        function backwards_comp(obj)
            % BACKWARDS_COMP ensures backwards compatability to saved
            % objects that are made with an older version of this class.
            %
            % Mostly it runs code that is in the constructor in the current
            % version.
            
            % Find the recording interval
            obj.settings.recording_interval =...
                [obj.raw_data.hdr.FirstTimeStamp, obj.raw_data.hdr.LastTimeStamp];
            
            % Fill out the settings
            %obj.settings.interval = [1, size(obj.raw_data.dat,2)];
            obj.settings.interval = [1, 20];
            obj.settings.preferred_electrodes = [1 2 3]; % Because humans like three dimensions
            obj.settings.colorlist = ['kgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbm'];
            
            % Figure out which electrodes are not recording anything at all
            % and label them 'false' in the working electrodes property
            for i=1:4
                if sum(sum(obj.data(i,:,:)))==0
                    obj.settings.working_electrodes(i) = false;
                else
                    obj.settings.working_electrodes(i) = true;
                end
            end
            
        end
        
    end
end

