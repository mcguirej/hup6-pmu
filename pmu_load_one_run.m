function [ibi] = pmu_load_one_run(fname)
% load pulse-ox data for one run and return an ibi timeseries
%   (ibi = cardiac inter-beat interval)
% 
% Input:
%   fname (optional): string w/ path to the pmu data file for this subject
%       and run. If not provided, user is prompted to select a file.

try

% gui to select a file if not provided as input
if nargin<1
    [filename, pathname] = uigetfile('*','Please select a PMU data file');
    fname = fullfile(pathname,filename);
end    
    
% return if there is no datafile
if ~exist(fname,'file')
    fprintf('no pulse file exists.\n');
    ibi = [];
    return;
end

% load the file
% (4 lines of code adapted from siemens_PPGload.m)
fid = fopen(fname);
ignore = textscan(fid,'%s',4); % Ignore first 4 values.
data = textscan(fid,'%u16'); % Read data until end of u16 data.
footer = textscan(fid,'%s');   % Read in remaining data (time stamps and statistics).
fclose(fid);

% 'ignore', 'data', and 'footer' are 1-element cell arrays
%   'ignore' holds the first 4 data values
%       these aren't real data (unsure what they are)
%   'data' has a long vector of data
%   'footer' has a few potentially useful entries:
%       PULS Freq Per: 85 700
%       PULS Min Max Avg StdDiff: 693 700 696 2
%       LogStartMDHTime:  28866777
%       LogStopMDHTime:   31673655
%       LogStartMPCUTime: 28865552
%       LogStopMPCUTime:  31673615

% read timestamps from footer
% formatted as ms since midnight
fprintf('footer information:\n');
footStampLabels = {'LogStartMDHTime', 'LogStopMDHTime'};
for i = 1:length(footStampLabels)
    labelStr = footStampLabels{i};
    labelIdx = find(strcmp(footer{1},[labelStr,':']))+1; % position of this timestamp in 'footer'
    foot.(labelStr) = str2double(footer{1}{labelIdx}); % store the timestamp in a more convenient data structure
    fprintf('\t%s: %s\n',labelStr,footer{1}{labelIdx});
end % loop over labels

% marker values (5000, 6000, or 5003) are intermixed with real data (<5000)
% identify the real datapoints for the time grid
d = double(data{1}); % data vector
% index marker values of 5000 (only a subset of all the marker values)
markerIdx = (d==5000); % index elements of d that represent marker values
% identify the number of non-marker values
nSamps = sum(d<5000);

% evaluate timestamps relative to the number of values in the data record
% NOTE:
% this generally doesn't work out perfectly, but is
% usually off by less than 1 sec (50 samples) in either direction
% unclear why there is a mismatch here.
% we'll proceed on the premise that this LogStart is correct and that the
% inter-sample interval is exactly 20ms, but acknowledging a small amount
% of potential timing imprecision.
sampHz = 50;
sampPd = 1000/sampHz; % sampling period in msec (20ms)
mdhDurMsec = foot.LogStopMDHTime - foot.LogStartMDHTime;
mdhDurMins = mdhDurMsec/(60*1000);
predictedNSamps = floor(mdhDurMsec/sampPd);
fprintf('  MDH timestamps signal %1.2f mins of recording, implying %d datapoints.\n',...
    mdhDurMins,predictedNSamps);
fprintf('  Data record contains %d points (difference = %d).\n',...
    nSamps,nSamps-predictedNSamps);




keyboard
%%% resume here %%%





% load run onset times 
% (previously derived from dicom headers, using dicomContentTimes.m)
dicomInfo = load('runOnsetTimestamps_n22.mat');
subjIdx = strcmp(id,{dicomInfo.allOnsets.id});
%%% NOTE:
%%% there are 2 options here which give results differing by 2-3 sec
%%% it's not clear which is correct
%%% there are 2 dicom header fields that provide (different) timestamps
% runOnset = dicomInfo.allOnsets(subjIdx).runOnset; % to use CONTENTTIMES field
runOnset = dicomInfo.allOnsets(subjIdx).runOnset_AT; % to use ACQUISITIONTIMES field
nRuns = length(runOnset);

% evaluate the lag at the beginning and end of recording
runLengthSec = 615;
runLengthMsec = runLengthSec*1000;
logStartLagSec = (runOnset(1) - foot.LogStartMDHTime)/1000;
logStopLagSec = (foot.LogStopMDHTime - (runOnset(nRuns)+runLengthMsec))/1000;
fprintf('Estimated lag to start and stop recording:\n');
fprintf('  PMU log started %1.2f sec before start of first run.\n',logStartLagSec);
fprintf('  PMU log ended %1.2f sec after end of run %d.\n',logStopLagSec,nRuns);


%%% assign a timestamp to each marker value in the data file (i.e.,
%%% identify times at which a "beat" was logged)
%%% then convert to IBI

% 1. assign a timestamp to each data value (msec since midnight)
sampleTimes = (1:nSamps)*sampPd + foot.LogStartMDHTime; % non-markers only
dTimes = nan(size(d)); % has entries for both marker and non-marker values
dTimes(~markerIdx) = sampleTimes;

% 2. each marker is assigned the timestamp of the previous non-marker
% assume (but check) that there are never 2 consecutive markers
preMarkerPos = find(markerIdx) - 1;
preMarkerPos(preMarkerPos==0) = []; % in case 1st value is a marker
markerTimes = dTimes(preMarkerPos);
assert(~any(isnan(markerTimes)),'some marker times not assigned; possible consecutive marker values?');

% 3. assign each marker an IBI value corresponding to the distance to the
% PREVIOUS marker
ibiVals = markerTimes(2:end) - markerTimes(1:(end-1));
ibiTimes = markerTimes(2:end);

% 4. clean up outlier values in the ibi timeseries, if any.
% an undetected beat would create an interval ~2x the typical IBI
% also might have 2 beats "detected" close together
% censor out values 30% above or below the median
medianIBI = median(ibiVals);
outlierIdx = ibiVals > 1.3*medianIBI; % high outliers
outlierIdx = outlierIdx | (ibiVals < 0.7*medianIBI);
ibiVals(outlierIdx) = nan;
pctOutlier = 100*sum(outlierIdx)/length(outlierIdx);
fprintf('median IBI = %1.1fms; %d outlier values deleted (%1.2f%%).\n',...
    medianIBI,sum(outlierIdx),pctOutlier);
 

% 5. apply windowed-average smoothing to the IBI timeseries???
%%% NOTE: not being done for now, but could be added here.
%%% the downside is a (further) loss of temporal precision.
%%% also might consider assessing outliers relative to a running median.


%%% extract information specific to each run (with timestamps relative to
%%% run onsets).
ibiRuns = cell(nRuns,1);
for r = 1:nRuns
    
    % IBI output format:
    % col1 = time from run onset, in sec
    % col2 = ibi value that applies AFTER that time. 
    runIdx = (ibiTimes > runOnset(r)) & (ibiTimes < (runOnset(r) + runLengthMsec));
    if any(runIdx)
        ibiRuns{r}(:,1) = ibiTimes(runIdx) - runOnset(r);
        ibiRuns{r}(:,1) = ibiRuns{r}(:,1)/1000; % convert to sec
        ibiRuns{r}(:,2) = ibiVals(runIdx);
        ibiRuns{r}(end,2) = nan; % signal that there are no data after this point
    else
        ibiRuns{r} = []; % there is at least one run w/ no valid pmu data
    end
    
end % loop over runs


%%% convert to percent of median IBI 
%%% (calculating the median on data from during task runs)
% important because there are large differences in baseline
% IBI across individuals; we want to remove these differences to focus on
% changes over time, while retaining the ability to see differences between
% runs in the two conditions.
ibiAllRuns = cell2mat(ibiRuns);
medianIBI = nanmedian(ibiAllRuns(:,2));
for r = 1:nRuns
    if ~isempty(ibiRuns{r}) % one run has no data
        ibiRuns{r}(:,2) = 100*ibiRuns{r}(:,2)./medianIBI;
    end
end


catch ME
    
    disp(getReport(ME));
    keyboard;
    
end





