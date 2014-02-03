
dicomDir = '/Users/jtm/Documents/project_code/high-res-mri-protocol/consult_2013-12-05/dicoms/0007_BOLD_mb4_1.5x1.5x1.5_TR1000';
pmuFile = '~/Documents/repositories/hup6-pmu/example_pmu_data/BM971.puls';
noPlot = false;

% standard 3s TR dicoms
% dicomDir = '/Users/jtm/Documents/project_code/high-res-mri-protocol/consult_2013-12-05/dicoms/0003_BOLD_mb0_3x3x3_TR3000';

dicom_times_one_run(dicomDir,pmuFile,noPlot);

