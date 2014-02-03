hup6-pmu
========

Code for processing physiological data acquired during MRI scans.

These functions load a PMU pulse ox data series and line it up with 
timestamps in the headers of the associated dicom files. 

dicom_times_one_run.m: This function prompts the user for a dicom directory
and the associated pulse file. It reads run onset information from the
dicoms (although there is potential ambiguity about this; see below) and
writes the results to a .mat file in the same location as the pulse file.
This function saves intermediate results so that it only needs to be run 
once for each scanning run (since it's somewhat slow, and dicom directories 
might not always be accessible).

pmu_load_one_run.m: Prompts the user for a pulse file. Reads timing 
information in the footer and lines up the data with dicom timing info
saved previously as described above. 

There are two reasons to believe the temporal alignment is imperfect.

(1) Dicom headers contain multiple timestamps, neither of which exactly
follows the expected pattern (increasing by the TR across successive
images). The "ContentTime" stamp is generally a couple seconds later than
the "AcquisitionTime" stamp. This might mean that ContentTime has more to 
do with the time the image was reconstructed. We are using the 
AcquisitionTime from the first image as the run onset timestamp. However,
of concern, AcquisitionTime values appear to increase across images by
slightly more than the TR. 

(We ruled out the possibility that the TR really is longer than intended
by stopwatch timing the run. The AcquisitionTime values imply a cumulative
lag of several seconds across a 10-min run, but the actual run duration 
as close to correct as we could measure.)

(2) Within the PMU data file, the number of data samples (after stripping
marker values) tends not to match exactly with the number implied by the
"start" and "stop" timestamps and the 50Hz sampling rate. The discrepancy
is reported by pmu_load_one_run.m. In test data it's usually <1s (<50
samples). 






