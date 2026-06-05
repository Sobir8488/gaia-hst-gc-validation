function T = convert_hacks_raw_to_csv(rootDir)
% convert_hacks_raw_to_csv
% -------------------------------------------------------------------------
% Converts a whitespace-delimited HACKS velocity-dispersion text file to
% input/hacks_hst_dispersion_profiles.csv.
%
% Preferred raw input location:
%   input/raw/hacks_velocity_dispersions.txt
%
% Convenience fallback locations searched automatically:
%   hacks_velocity_dispersions.txt
%   Pasted text.txt
%   Pasted_text.txt
%   input/Pasted text.txt
%   input/raw/Pasted text.txt
% -------------------------------------------------------------------------
if nargin < 1 || isempty(rootDir); rootDir = pwd; end
rootDir = char(rootDir);
rawFile = findHacksRawFile(rootDir);
outFile = fullfile(rootDir,'input','hacks_hst_dispersion_profiles.csv');
if isempty(rawFile)
    error(['Raw HACKS file not found. Put the HACKS velocity-dispersion text at:\n  %s\n' ...
           'Then run:\n  convert_hacks_raw_to_csv(pwd)\n'], fullfile(rootDir,'input','raw','hacks_velocity_dispersions.txt'));
end
T = parseRawHacks(rawFile);
if height(T)==0
    error('The file was found but no HACKS rows were parsed: %s', rawFile);
end
if ~isfolder(fullfile(rootDir,'input')); mkdir(fullfile(rootDir,'input')); end
writetable(T,outFile);
fprintf('Parsed HACKS raw file:\n  %s\n', rawFile);
fprintf('Wrote %d HACKS rows for %d clusters to:\n  %s\n', height(T), numel(unique(T.cluster_id)), outFile);
end

function rawFile = findHacksRawFile(rootDir)
candidates = {
    fullfile(rootDir,'input','raw','hacks_velocity_dispersions.txt')
    fullfile(rootDir,'hacks_velocity_dispersions.txt')
    fullfile(rootDir,'Pasted text.txt')
    fullfile(rootDir,'Pasted_text.txt')
    fullfile(rootDir,'input','Pasted text.txt')
    fullfile(rootDir,'input','raw','Pasted text.txt')
};
rawFile = '';
for i=1:numel(candidates)
    if isfile(candidates{i})
        rawFile = candidates{i};
        return;
    end
end
end

function T = parseRawHacks(rawFile)
fid=fopen(rawFile,'r');
if fid<0; error('Cannot open raw HACKS file: %s', rawFile); end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
rows={};
while true
    line=fgetl(fid);
    if ~ischar(line); break; end
    line=strtrim(line);
    if isempty(line) || startsWith(line,'#'); continue; end
    parts=regexp(line,'\s+','split');
    if numel(parts)<10; continue; end
    cid=normalizeClusterIdLocal(parts{1});
    nums=str2double(parts(2:10));
    if any(~isfinite(nums)); continue; end
    rows(end+1,:)={string(cid),string(parts{1}),nums(1),nums(2),nums(3),nums(4),nums(5),nums(6),nums(7),nums(8),nums(9),"parsed_raw_hacks"}; %#ok<AGROW>
end
if isempty(rows)
    T = table(string.empty(0,1), string.empty(0,1), [], [], [], [], [], [], [], [], [], string.empty(0,1), ...
        'VariableNames',{'cluster_id','hacks_cluster_id','bin','n_stars','r_arcsec','sigma_pm_masyr','sigma_pm_err_masyr','sigma_radial_masyr','sigma_radial_err_masyr','sigma_tangential_masyr','sigma_tangential_err_masyr','source_note'});
else
    T=cell2table(rows,'VariableNames',{'cluster_id','hacks_cluster_id','bin','n_stars','r_arcsec','sigma_pm_masyr','sigma_pm_err_masyr','sigma_radial_masyr','sigma_radial_err_masyr','sigma_tangential_masyr','sigma_tangential_err_masyr','source_note'});
end
end

function cid=normalizeClusterIdLocal(s)
s=char(string(s)); u=upper(strtrim(s));
tok=regexp(u,'NGC[_\s-]*0*(\d+)','tokens','once'); if ~isempty(tok); cid=sprintf('NGC%04d',str2double(tok{1})); return; end
tok=regexp(u,'IC[_\s-]*0*(\d+)','tokens','once'); if ~isempty(tok); cid=sprintf('IC%04d',str2double(tok{1})); return; end
tok=regexp(u,'PAL[_\s-]*0*(\d+)','tokens','once'); if ~isempty(tok); cid=sprintf('PAL%02d',str2double(tok{1})); return; end
cid=regexprep(u,'[^A-Z0-9]+','_');
end
