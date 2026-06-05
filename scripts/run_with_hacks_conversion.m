% run_with_hacks_conversion
% Use this script after placing the raw HACKS table at:
%   input/raw/hacks_velocity_dispersions.txt

thisFile = mfilename('fullpath');
scriptDir = fileparts(thisFile);
rootDir = fileparts(scriptDir);
addpath(fullfile(rootDir,'src'));

convert_hacks_raw_to_csv(rootDir);
results = gaia_hst_validation_pipeline(rootDir);
assignin('base','results',results);
fprintf('
Saved results struct to workspace variable: results
');
