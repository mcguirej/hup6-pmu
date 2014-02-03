function [] = dicom_times_one_run(dicomDir,pmuFile,noPlot)
% reads ContentTime and AcquisitionTime field from dicom header files
% (CT and AT), and saves a run onset timestamp for future use.
%
% Inputs (all optional):
%   dicomDir (string): path to a dir with one run's dicom files
%       dicoms must be uncompressed
%       also must be numbered sequentially
%   pmuFile (string): path to the corresponding pmu file
%       run timestamp will be saved in a matching location
%   noPlot (logical): skip diagnostic plots
% 
% observations from examining an example run (TR = 2.5s):
% For each image, these 2 fields hold slightly different timestamps.
% For the first image in a run, AT is 2.5-3 sec earlier than CT. 
% Thereafter, the median spacing between adjacent images is 2510ms for AT
% and 2500 for CT. By the end, of a run, they nearly agree. 
% 
% Of note, it is not the case that the actual TR is 2510 ms; this would
% have made the experiment several seconds too long, but in fact its 
% duration was accurate within stopwatch precision. So AcquisitionTime
% entries are not perfectly accurate. However, the AT for the first image
% appears to be a good approximation of run onset time.
% 
% The CT is also not fully accurate; inter-image intervals vary quite a
% bit, but successive intervals tend to compensate so they basically agree
% on a common time grid. (CT might reflect the time the image is
% reconstructed).
%
% This function runs diagnostics on each set of timestamps, and saves
% estimated run start times (msec since midnight) based on CT and AT. 
% Work with PMU data should bear in mind that there may be a couple
% seconds of imprecision. 

% check dicomDir input: if absent, prompt to select a dir
if nargin<1 || isempty(dicomDir)
    fprintf('\nSelect a dicom directory...\n');
    dicomDir = uigetdir('','Select a dicom directory');
end
fprintf('Dicom dir:\n\t%s\n',dicomDir);
assert(exist(dicomDir,'dir')>0,'Dicom dir not found.');

% check pmuFile input: if absent, prompt to select a dir
if nargin<2 || isempty(pmuFile)
    fprintf('\nSelect the matching PMU file...\n');
    [pmuFile_name, pmuFile_path] = uigetfile('*','Select the matching PMU file');
    pmuFile = fullfile(pmuFile_path,pmuFile_name);
end
fprintf('PMU file:\n\t%s\n',pmuFile);
assert(exist(pmuFile,'file')>0,'PMU file not found.');

% check noPlot input: if absent, plots set to false (plots will be made)
if nargin<3, noPlot = false; end
assert(islogical(noPlot),'Input "noPlot" must be true or false');

% set the output filename on the basis of pmuFile
outFile = [pmuFile,'.dicomOnset.mat'];
fprintf('Output file:\n\t%s\n',outFile);
% do nothing if this output already exists
if exist(outFile,'file')
    fprintf('Output exists; skipping.\n');
    return;
end

% names of dicom header fields containing timestamps
tsFields = {'ContentTime', 'AcquisitionTime'};

try

% list files in the dicom directory
dicomList = dir(fullfile(dicomDir,'*'));
dicomNames = {dicomList(:).name}';

% remove hidden files (non-dicoms) from the list
isHidden = strncmp('.',dicomNames,1);
dicomList(isHidden) = [];
nTmpts = length(dicomList);
fprintf('reading headers for %d dicom images...',nTmpts);

% loop over timepoints
TR_all = nan(nTmpts,1); % will hold TR from each image header
for t = 1:nTmpts
    
    % check if the dicom is compressed
    % (should be uncompressed manually before calling this function)
    if strcmp('.gz',dicomList(t).name((end-2):end))
        fprintf('*** Dicoms appear to be gzipped ***\n*** Exiting ***\n');
        return;
    end
    
    % header info
    info = dicominfo(fullfile(dicomDir,dicomList(t).name));
    TR_all(t,1) = info.RepetitionTime;
    
    % process each of the 2 timestamps in the header
    % loop over 'ContentTime' and 'AcquisitionTime'
    for tsIdx = 1:2

        tsName = tsFields{tsIdx};
        timeStr = info.(tsName);

        % convert the timestamp from HHMMSS.SSSSSS to "msec since midnight"
        hrsInMsec = 60*60*1000*str2double(timeStr(1:2));
        minsInMsec = 60*1000*str2double(timeStr(3:4));
        secInMsec = 1000*str2double(timeStr(5:end));
        tstmp.(tsName)(t,1) = hrsInMsec + minsInMsec + secInMsec;

    end
    
end % loop over timepoints
fprintf('done.\n');

% check the TR
TR = unique(TR_all);
assert(numel(TR)==1,'Inconsistent TRs!');
fprintf('TR = %d\n',TR);

% check the TR implied by image timestamps
TR_CT = diff(tstmp.ContentTime);
TR_AT = diff(tstmp.AcquisitionTime);

% verify that timestamps are non-decreasing across images
% (i.e., that dicoms were loaded in the correct order)
assert(all(TR_CT>=0),'ContentTime values are out of order');
assert(all(TR_AT>=0),'AcquisitionTime values are out of order');

% run-onset estimate from ContentTime field (median across all images)
acqLags = (0:(nTmpts-1))'*TR; % each image's intended lag from run onset, in ms
onsetEstimates = tstmp.ContentTime - acqLags; % each image's estimate of run start time
runOnset_CT = median(onsetEstimates);

% run-onset estimate from AcquisitionTime field (first image)
runOnset_AT = tstmp.AcquisitionTime(1);

% save results
% output matfile includes both timestamps, and the name of the dicom
% directory on which they are based. 
save(outFile,'dicomDir','runOnset_CT','runOnset_AT');
fprintf('Saved results.\n\n');

% exit if diagnostic info was not requested
if noPlot, return; end



%%% optional diagnostic info and plots %%%

% difference between the two onset estimates
ct_minus_at_in_s = (runOnset_CT - runOnset_AT)/1000;
fprintf('Diagnostics:\n');
fprintf('  CT-based onset is %1.3f s later than AT-based onset.\n',ct_minus_at_in_s);

% median difference between successive image timestamps
fprintf('  ContentTime stamps imply median TR of %d ms\n',median(TR_CT));
fprintf('  AcquisitionTime stamps imply median TR of %d ms\n',median(TR_AT));

% plot the 2 timestamps for this run
figure(1); clf;
pred_ct = tstmp.ContentTime(1) + acqLags;
pred_at = tstmp.AcquisitionTime(1) + acqLags;
plot([tstmp.ContentTime, tstmp.AcquisitionTime]/1000,'o'); % convert from msec to sec
hold on;
plot([pred_ct, pred_at]/1000,'-');
set(gca,'Box','off','FontSize',20);
legend('ContentTime','AcquisitionTime','Predicted CT','Predicted AT');
ylabel('ContentTime minus AcquisitionTime (s)');
xlabel('Image number');

% plot the TR intervals implied by this run's dicom timestamps
figure(2); clf;
plotData = [TR_CT, TR_AT];
plot(plotData,'o-');
ylims = get(gca,'YLim');
set(gca,'Box','off','FontSize',20,'YLim',[0, ylims(2)]);
legend('ContentTime','AcquisitionTime');
ylabel('Inter-acquisition interval');
xlabel('Image number');

catch ME
    disp(getReport(ME));
    keyboard;
end

