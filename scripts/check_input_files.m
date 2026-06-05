% check_input_files
% Report whether the expected input files are present and readable.

thisFile = mfilename('fullpath');
scriptDir = fileparts(thisFile);
rootDir = fileparts(scriptDir);
files = {
 'input/gaia_edr3_profile_percentiles.csv'
 'input/gaia_edr3_profile_summary.csv'
 'input/gaia_edr3_catalogue_audit_170.csv'
 'input/hacks_hst_dispersion_profiles.csv'
 'input/apogee_dr17_gc_parameters.csv'
 'input/apogee_dr17_gc_member_summary.csv'
 'input/cluster_crosswalk_master.csv'
 'input/analysis_config.csv'
 'input/raw/hacks_velocity_dispersions.txt'
};
name = strings(numel(files),1); exists=false(numel(files),1); n_rows=nan(numel(files),1);
for i=1:numel(files)
    name(i)=string(files{i});
    f=fullfile(rootDir,files{i});
    exists(i)=isfile(f);
    if exists(i) && endsWith(f,'.csv')
        try
            T=readtable(f,'VariableNamingRule','preserve');
            n_rows(i)=height(T);
        catch
        end
    end
end
disp(table(name, exists, n_rows));
