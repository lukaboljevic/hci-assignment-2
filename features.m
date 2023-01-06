function [featVecFile, refFile] = features(subject, readWhich, dataDir, n, plot)
% Extract feature vectors for a given subject
% For a given subject, read 3 records, corresponding to one of 4 tasks, and
% extract feature vectors, using CSP, band pass filtering and log(var(.)).
    if readWhich < 3 || readWhich > 6
        error("Value for parameter `readWhich` has to be between 3 and 6!");
    end
    
    fprintf("\n");
    disp("====================================================");
    disp("Subject " + subject);
    disp("Reading records corresponding to task " + (readWhich - 2));
    disp("Band pass FIR filter order: " + n);
    fprintf("\n");

    eeg = 0;
    if nargin == 2
        eeg = 1;
        disp("Records will be downloaded directly from EEGMMI DB.");
    end

    subject = string(subject);


    % Which 3 subject records to read. These 3 records all correspond to
    % one of 4 tasks described on PhysioNet page for EEGMMI DB.
    records = [];
    for i=0:2
        records = [records string(num2str(readWhich + (i*4), "%02d"))];
    end
    fprintf("Record identifiers: ");
    disp(records);
    

    % Cell arrays to store T1 and T2 intervals for all records
    allT1 = {}; allT2 = {};
    

    % Read records and get their T0, T1 and T2 intervals
    if plot
        t1figure = figure();
    end
    for i=1:size(records, 2)
        if nargin == 2
            % Read record directly from EEGMMI DB. Slow!
            recordName = strcat("/eegmmidb/", subject, "/", subject, "R", records(i), ".edf");
        else
            % Record(s) for subject are already downloaded in directory
            % given with dataDir.
            cd(strcat(dataDir, "/", subject));
            recordName = strcat(subject, "R", records(i), ".edf");
        end
        disp("Record name: " + recordName);
        
        recordName = convertStringsToChars(recordName);
        [sig, fs, ~] = rdsamp(recordName, 1:64); 
        [~, T1, T2] = getIntervals(recordName, 'event', fs, size(sig, 1));
        sig = sig.';
        
        % If we didn't download from EEGMMI directly, but used local files
        if ~eeg
            cd ("../..");
        end

%         sec = 4.0;
%         fd = fs * sec;
%         [p, q] = size(sig);
        ignoreLast = 0.2 * fs;  % ignore last 0.2 seconds of any interval

        % Append all T1 intervals for this record
        for j=1:size(T1, 1)
            allT1{end+1} = sig(:, T1(j, 1):T1(j, 2)-ignoreLast);
        end
%         for j=1:size(T1, 1)
%             if (T1(j, 1) + fd) > q
%                 allT1{end+1} = sig(:, T1(j,1):q);
%             else
%                 allT1{end+1} = sig(:, T1(j,1):T1(j,1) + fd);
%             end
%         end
        
        % Append all T2 intervals for this record
        for j=1:size(T2, 1)
            allT2{end+1} = sig(:, T2(j, 1):T2(j, 2)-ignoreLast);
        end
%         for j=1:size(T2, 1)
%             if (T2(j, 1) + fd) > q
%                 allT2{end+1} = sig(:, T2(j,1):q);
%             else
%                 allT2{end+1} = sig(:, T2(j,1):T2(j,1) + fd);
%             end
%         end

        if plot
            displayName = strcat(subject, "R", records(i), ".edf");
            plot(allT1{end}(1, :), "DisplayName", displayName);
            hold on;
        end
    end
    fprintf("\n");
    
    if plot
        % Set figure position and width so we can see the T1s nicely
        figurePosition = t1figure.Position;  % current figure position and size
        y0 = figurePosition(2); figureHeight = figurePosition(4);
    
        screenSize = get(0, 'ScreenSize');  % screen size
        screenWidth = screenSize(3);
    
        title("One T1 interval for each record");
        legend;
        set(gcf, 'position', [0 y0-150 screenWidth figureHeight]);
    end


    % Calculate mean of first few T1 and T2 intervals, and use that for CSP
    N = 5;
    [allT1Mean] = calculateMean(allT1, N);
    [allT2Mean] = calculateMean(allT2, N);


    % CSP projection matrix
%     [W] = CSP(cell2mat(allT1(1)), cell2mat(allT2(1)));
%     allT1(1) = [];
%     allT2(1) = [];
    [W] = CSP(allT1Mean, allT2Mean);
    

    % Band pass filtering from 8 to 13Hz 
    % (mu rhythm, responsible for motor activities, lies in the frequency 
    % range 8-13Hz)
    f = [0 8 8 13 13 fs/2] / (fs/2);
    a = [0 0 1 1  0   0];
    b = firls(n, f, a);


    % Calculate feature vectors
    [features1] = featureVector(allT1, W, b); 
    [features2] = featureVector(allT2, W, b);


    % Scatter plot
%     figure;
%     scatter(features1(:, 1), features1(:, 2));
%     hold on;
%     scatter(features2(:, 1), features2(:, 2));
    

    % Save feature vectors and reference classes
    if ~isfolder("results")
        mkdir("results")
    end
    if ~isfolder(strcat("results/", subject))
        mkdir(strcat("results/", subject))
    end

    % readWhich-2 is the task number, as described on PhysioNet page for EEGMMI DB.
    featVecFile = strcat("results/", subject, "/", subject, "-task", num2str(readWhich-2), "-", num2str(n), "-featv.txt");
    refFile = strcat("results/", subject, "/", subject, "-task", num2str(readWhich-2), "-", num2str(n), "-ref.txt");

    ffv = fopen(featVecFile, "wt");
    frc = fopen(refFile, "wt");

    for i=1:size(features1, 1)
        fprintf(ffv, "%.8f %.8f\n", features1(i, 1), features1(i, 2));
        fprintf(frc, "T1\n");
    end

    for i=1:size(features2, 1)
        fprintf(ffv, "%.8f %.8f\n", features2(i, 1), features2(i, 2));
        fprintf(frc, "T2\n");
    end

    fclose(ffv);
    fclose(frc);
end


%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%

function [meanInterval] = calculateMean(intervals, N)
% Calculate the mean of first N T1 or T2 intervals
    meanInterval = cell2mat(intervals(1));
    for j = 2:N
        meanInterval = meanInterval + cell2mat(intervals(j));
    end
    meanInterval = meanInterval / N;
end


function [feats] = featureVector(intervals, W, b)
% Calculate feature vector from T1 or T2 intervals
    feats = [];
    for i=1:size(intervals, 2)
        % Filter using CSP projection matrix W
        tmp = W * cell2mat(intervals(i));

        % Take signals with the highest/lowest variance (entire point of CSP!).
        % We can also take e.g. the first 2/3, and last 2/3 signals, instead of
        % just the first and last. Taking the first/last 2 favours LDA, while
        % taking first/last 3 favours QDA (empirical). Taking just the first 
        % and last seems to be somewhere in between.
        tmp = [tmp(1,:).' tmp(size(tmp, 1), :).'].';
%         len = size(tmp, 1);
%         tmp = [tmp(1,:).' ...
%                tmp(2,:).' ...
% %                tmp(3,:).' ...
% %                tmp(len-2, :).' ...
%                tmp(len-1, :).' ...
%                tmp(len, :).'].';

        % Filter with band-pass filter (coefficients already calculated)
        tmp = filter(b, 1, tmp);

        % Apply log var operators to get features
        feats = [feats; log(var(tmp(1, :))) log(var(tmp(2, :)))];
    end
end


function [T0, T1, T2] = getIntervals(name, annotator, fs, samps)
% This function returns arrays of interval start and stop times in samples, 
% for the three different interval types - T0, T1, T2.
% Input:
% 
%   name: character vector containing record name, e.g 'S001R01.edf'
%   annotator: character vector containing annotator name e.g. 'event'
%   fs: sampling frequency of the record
%   samps: length of a data file in samples
%
% Output
%   
%   T0: array of start and end times in samples for interval T0
%   T1: array of start and end times in samples for interval T1
%   T2: array of start and end times in samples for interval T2

    T0 = [];
    T1 = [];
    T2 = [];

    [anot, ~, ~, ~, ~, cmt] = rdann(name, annotator);
    % anot: contains indexes where a certain interval (T0, T1, T2) starts
    % cmt: same length as anot, cmt(i) tells us which interval is at anot(i), and how long it lasts
    cmt = string(cmt);

    for i=1:length(anot)
        splitted = split(cmt(i), " ");
        duration = str2double(splitted(3)) * fs;  % duration of individual intervals (T0, T1, T2)

        % Way 1.
        intervalEnd = anot(i) + duration - 1;
        if (i < length(anot))
            intervalEnd = min(anot(i+1) - 1, intervalEnd);
        else
            intervalEnd = min(samps, intervalEnd);
        end

        % Way 2.
%         if (i<length(anot)) 
%             intervalEnd=anot(i+1) - 1;
%         else
%             intervalEnd=anot(i) + duration;
%         end
%         intervalEnd = min(intervalEnd, samps);

        if (splitted(1) == "T0")
            T0 = [T0; [anot(i) intervalEnd]];  % append start and end
        elseif (splitted(1) == "T1")
            T1 = [T1; [anot(i) intervalEnd]];
        elseif (splitted(1) == "T2")
            T2 = [T2; [anot(i) intervalEnd]];
        end
    end
end


function [result] = CSP(varargin)
% This code is for calulating the projection matrix for CSP 
% Haider Raza, Intelligent System Research Center, University of Ulster, Northern Ireland, UK.
%     Raza-H@email.ulster.ac.uk
%     Date: 03-Oct-2014
% Input:
%             
%       left:  left hand data
%       right: right hand data
% 
% Output:
%       left:  left hand data
%       right: right hand data    

    if (nargin ~= 2)
        disp('Must have 2 classes for CSP!')
    end
    
    Rsum=0;
    %finding the covariance of each class and composite covariance
    for i = 1:nargin 
        %mean here?
        R{i} = ((varargin{i}*varargin{i}')/trace(varargin{i}*varargin{i}'));%instantiate me before the loop!
        %Ramoser equation (2)
        Rsum=Rsum+R{i};
    end
   
    %   Find Eigenvalues and Eigenvectors of RC
    %   Sort eigenvalues in descending order
    [EVecsum,EValsum] = eig(Rsum);
    [EValsum,ind] = sort(diag(EValsum),'descend');
    EVecsum = EVecsum(:,ind);
    
    %   Find Whitening Transformation Matrix - Ramoser Equation (3)
    W = sqrt(inv(diag(EValsum))) * EVecsum';
    
    
    for k = 1:nargin
        S{k} = W * R{k} * W'; % Whiten Data Using Whiting Transform - Ramoser Equation (4)
    end
   
    %generalized eigenvectors/values
    [B,D] = eig(S{1},S{2});
    % Simultanous diagonalization
	% Should be equivalent to [B,D]=eig(S{1});
    
    [D,ind]=sort(diag(D));B=B(:,ind);
    
    %Resulting Projection Matrix-these are the spatial filter coefficients
    result = B'*W;

end