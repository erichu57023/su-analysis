function modifyNTTFeatures(file, feature1, feature2)
    % MODIFYNTTFEATURES(file, feature1, feature2)
    % Adjusts the default features for the specified .ntt file.
    % file: char defining target file path
    % feature1: string representing first feature
    % feature2: string representing second feature

if (~contains(file, '.ntt'))
    disp([file, ' is not a .ntt file.']);
    return;
end

featureCount = 1;
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(file, [1 1 1 1 1], 1, 1, [] );
for ii = 1 : length(Header)
    splitString = strsplit(strtrim(Header{ii}));
    if strcmp(splitString{1}, '-Feature')
        if featureCount < 5
            splitString{2} = feature1;
        else
            splitString{2} = feature2;
        end
        featureCount = featureCount + 1;
    end
    Header{ii} = strjoin(splitString);
    if featureCount > 8
        break;
    end
end

Mat2NlxSpike(file, 0, 1, [], [1 1 1 1 1 1], Timestamps, ScNumbers, CellNumbers, Features, Samples, Header);
