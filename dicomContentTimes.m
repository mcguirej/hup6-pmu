function [] = dicomContentTimes()
% reads ContentTime and AcquisitionTime field from dicom header files
% (CT and AT)
% 
% observations from examining subject 1:
% For each image, these 2 fields hold slightly different timestamps.
% For the first image in a run, AT is 2.5-3 sec earlier than CT. 
% Thereafter, the median spacing between adjacent images is 2510ms for AT
% and 2500 for CT. By the end, of a run, they nearly agree. 
% 
% Of note, it is not the case that the actual TR is 2510 ms; this would
% have made the experiment several seconds too long, but in fact its 
% duration was accurate within stopwatch precision. So AcquisitionTime
% entries are not perfectly accurate. 
% 
% The CT is also not fully accurate; inter-image intervals vary quite a
% bit, but successive intervals tend to compensate so they basically agree
% on a common time grid. 
%
% This function runs diagnostics on each set of timestamps, and saves
% estimated run start times (msec since midnight) based on ContentTime. 
% Work with PMU data should bear in mind that there may be a couple
% seconds of imprecision. 



plotResults = false;
hdrFields = {'ContentTime', 'AcquisitionTime'};

addpath('../analysis'); % to access subject info
subjects = subInfo;
n = length(subjects);

% load existing file and append to it if possible
outputfname = sprintf('runOnsetTimestamps_n%d.mat',n);
if exist(outputfname,'file')
    d = load(outputfname);
    allOnsets = d.allOnsets;
    idsComplete = {allOnsets.id};
else
    allOnsets = struct([]);
    idsComplete = {};
end

try
for s = 1:n
% for s = 1
  id = subjects(s).id;
  fprintf('\n%s: \n',id);
  if any(strcmp(id,idsComplete)) % if this subject has been done already
      fprintf('  skipping.\n');
      continue;
  end
  
  nruns = subjects(s).nruns;
  onsetEstimates = cell(nruns,1);
  runOnset = nan(nruns,1);
  runOnset_AT = nan(nruns,1);
  
  for r = 1:nruns
    fprintf('  reading r%d... \n',r);
    [dicomPath, dicomDir] = fileparts(subjects(s).dicomEPI{r});
    dicomDir = fullfile('/Volumes','JTM01','jtm_cfn','qtask1','data','subjects',id,'raw',dicomDir);
    dicomList = dir(fullfile(dicomDir,'0*'));
    nTmpts = subjects(s).nTmpts;

    for t = 1:nTmpts(r)
      info = dicominfo(fullfile(dicomDir,dicomList(t).name));
      
      for hfIdx = 1:2 % for each if the 2 header fields containing timestamps
          
        tstampField = hdrFields{hfIdx};
        timeStr = info.(tstampField);

        % convert the timestamp to "msec since midnight" format
        hrsInMsec = 60*60*1000*str2num(timeStr(1:2));
        minsInMsec = 60*1000*str2num(timeStr(3:4));
        secInMsec = 1000*str2num(timeStr(5:end));
        tstmp.(tstampField){r}(t,1) = hrsInMsec + minsInMsec + secInMsec;
      
      end

    end % loop over timepoints
    
    % best estimate of run start time
    % use contenttime timestamp
    TR = 2.5*1000; % in msec
    acqLags = (0:(nTmpts(r)-1))'*TR; % in sec
    onsetEstimates{r} = tstmp.ContentTime{r} - acqLags;
    runOnset(r) = median(onsetEstimates{r});
    
    % diagnostic info to print:
    % 1. difference in between onset estimate and:
    %   - image 1 AT (expect AT will be ~2.5-3 sec earlier)
    %   - estimate based on final image AT (expect to be more similar)
    % 2. median TR implied by both CT and AT
    % run onset from 1st image acquisition time
    AT_first = tstmp.AcquisitionTime{r}(1); 
    runOnset_AT(r) = AT_first;
    fprintf('  1st-image AcquisitionTime disagrees by %1.2f sec.\n',abs(runOnset(r)-AT_first)/1000);
    % run onset from LAST image acquisition time
    AT_final = tstmp.AcquisitionTime{r}(end) - TR*(nTmpts(r)-1);
    fprintf('  last-image AcquisitionTime disagrees by %1.2f sec.\n',abs(runOnset(r)-AT_final)/1000);
    
    tr_ct = tstmp.ContentTime{r}(2:end) - tstmp.ContentTime{r}(1:(end-1));
    tr_at = tstmp.AcquisitionTime{r}(2:end) - tstmp.AcquisitionTime{r}(1:(end-1));
    fprintf('  ContentTime implies median TR of %dms\n',median(tr_ct));
    fprintf('  AcquisitionTime implies median TR of %dms\n',median(tr_at));
    
    if plotResults
    
        % plot difference between the 2 timestamps for this run
        figure(1);
        subplot(1,4,r);
        plotData = tstmp.ContentTime{r} - tstmp.AcquisitionTime{r};
        plotData = plotData/1000; % convert from msec to sec
        plot(plotData,'o');
        ylabel('content time vs. acquisition time (sec)');
        xlabel('acquisition number');
        title(sprintf('run %d',r));

        % plot the TR intervals implied by this run's timestamps
        figure(2);
        subplot(1,4,r);
        plotData = [tr_ct, tr_at];
        plot(plotData,'o-');
        legend('ContentTime','AcquisitionTime');
        ylabel('inter-acquisition interval');
        xlabel('acquisition number');
        
    end
    
  end % loop over runs
  
  % store the run onset estimates for this subject
  allOnsets(s).id = id;
  allOnsets(s).runOnset = runOnset;
  allOnsets(s).runOnset_AT = runOnset_AT;
  
  % display the results
  fprintf('  ** run onsets for %s (min since midnight) **\n  ',id);
  fprintf('%1.2f ',runOnset/60000); % (format as mins, not msec)
  fprintf('\n');

end % loop over subjects

% save out run onset timestamps for each subject. 
if length(allOnsets)==n % will not save for partial or test runs
    save(outputfname,'allOnsets');
end

catch ME
    
    disp(getReport(ME));
    keyboard;
    
end

