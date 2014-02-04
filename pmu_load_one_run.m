function [ibi] = pmu_load_one_run(fname)
% load pulse-ox data for one run and return an ibi timeseries
%   (ibi = cardiac inter-beat interval)
% 
% Input:
%   fname (optional): string w/ path to the pmu data file for this subject
%       and run. If not provided, user is prompted to select a file.
%
% Output:
%   ibi: 2-column matrix
%       col 1: heartbeat times, in s since run onset
%       col 2: preceding inter-beat interval in s

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
    fprintf('  %s: %s\n',labelStr,footer{1}{labelIdx});
end % loop over labels

% marker values (5000, 6000, or 5003) are intermixed with real data (<5000)
% identify the real datapoints for the time grid
d = double(data{1}); % data vector
% index marker values of 5000 (only a subset of all the marker values)
markerIdx = (d==5000); % index elements of d that represent marker values
% identify the number of non-marker values
sampIdx = (d<5000);
nSamps = sum(sampIdx);

% evaluate timestamps relative to the number of values in the data record
% NOTE:
% this generally doesn't work out perfectly, but is usually off by less 
% than 1 sec (50 samples) in either direction
% unclear why there is a mismatch here.
% we'll proceed on the premise that this LogStart is correct and that the
% inter-sample interval is exactly 20ms, but be aware of a small amount
% of potential timing imprecision.
sampHz = 50;
sampPd = 1000/sampHz; % sampling period in msec (20ms)
mdhDurMsec = foot.LogStopMDHTime - foot.LogStartMDHTime;
mdhDurMins = mdhDurMsec/(60*1000);
predictedNSamps = floor(mdhDurMsec/sampPd);
fprintf('  MDH timestamps signal %1.2f mins of recording.\n',mdhDurMins);
fprintf('  At %dHz this implies %d datapoints.\n',sampHz,predictedNSamps);
fprintf('  Data record contains %d points (difference = %d).\n',...
    nSamps,nSamps-predictedNSamps);

%%% align with dicom timing

% load results previously stored by dicom_times_one_run.m
dicom_fname = [fname,'.dicomOnset.mat'];
if ~exist(dicom_fname,'file')
    fprintf('Dicom onset info not found.\n');
    fprintf('Use dicom_times_one_run to create a file called %s\nExiting\n',dicom_fname);
    return;
end
dicomInfo = load(dicom_fname);

% use AcquisitionTimes dicom header field
runOnset = dicomInfo.runOnset_AT;

% evaluate the lag at the beginning and end of recording
logStartLagSec = (runOnset - foot.LogStartMDHTime)/1000;
logStopLagSec = (foot.LogStopMDHTime - runOnset)/1000;
fprintf('\nRun onset: %d\n',round(runOnset));
fprintf('  PMU log started %1.2f s before run onset.\n',logStartLagSec);
fprintf('  PMU log ended %1.2f s after run onsest.\n',logStopLagSec);

%%% assign a timestamp to each marker value in the data file (i.e.,
%%% identify times at which a "beat" was logged)
%%% then convert to IBI

% 1. assign a timestamp to each data value (ms since run onset)
sampleTimes = (1:nSamps)*sampPd + foot.LogStartMDHTime - runOnset; % non-markers only
sampleTimes = sampleTimes/1000; % convert to s since run onset
dTimes = nan(size(d)); % has entries for both marker and non-marker values
dTimes(sampIdx) = sampleTimes;

% 2. each marker is assigned the timestamp of the previous non-marker
% assume (but check) that there are never 2 consecutive markers
preMarkerPos = find(markerIdx) - 1;
preMarkerPos(preMarkerPos==0) = []; % in case 1st value is a marker
markerTimes = dTimes(preMarkerPos);
assert(~any(isnan(markerTimes)),'some marker times not assigned; possible consecutive marker values?');

% 3. assign each marker an IBI value corresponding to the distance to the
% previous marker
ibiVals = diff(markerTimes);
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
fprintf('\nmedian IBI = %1.3f s; %d outlier values deleted (%1.2f%%).\n',...
    medianIBI,sum(outlierIdx),pctOutlier);

% output matrix
ibi = [ibiTimes, ibiVals];

% may want to convert to percent of median IBI
% however, this should probably be done across all runs for a given subject
% (i.e., outside this function). 

% depending on the eventual analysis, may also want to apply smoothing to
% the ibi timecourse

catch ME
    
    disp(getReport(ME));
    keyboard;
    
end





