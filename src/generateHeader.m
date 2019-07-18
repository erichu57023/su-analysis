function header = generateHeader(range)
    % GENERATEHEADER(int[] Range)  Creates a generic Neuralynx Data File
    % Header containing only recording range data for voltage conversions.
    % For tetrode recordings, Range is a 1x4 int array.
    
    header = cell(6,1);
    header{1} = '######## Neuralynx Data File Header';
    header{2} = '-SamplingFrequency 32000';
    header{3} = '-ADMaxValue 32767';
    voltageString = [' ', num2str(range * 1e-6 / 32767, 16)];
    header{4} = ['-ADBitVolts', voltageString];
    rangeString = [' ', num2str(range)];
    header{5} = ['-InputRange', rangeString];
    header{6} = '-InputInverted True';
    
end
