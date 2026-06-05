% run_gaia_hst_validation
% Execute the Gaia--HST validation workflow from the repository root.

thisFile = mfilename('fullpath');
scriptDir = fileparts(thisFile);
rootDir = fileparts(scriptDir);
addpath(fullfile(rootDir,'src'));

rawFile = fullfile(rootDir,'input','raw','hacks_velocity_dispersions.txt');
hacksCsv = fullfile(rootDir,'input','hacks_hst_dispersion_profiles.csv');
try
    H = readtable(hacksCsv,'VariableNamingRule','preserve');
catch
    H = table();
end
if height(H)==0 && isfile(rawFile)
    fprintf('Converting raw HACKS text to CSV...
');
    convert_hacks_raw_to_csv(rootDir);
end

results = gaia_hst_validation_pipeline(rootDir);
assignin('base','results',results);
fprintf('
Saved results struct to workspace variable: results
');
