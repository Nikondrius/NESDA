%% ==========================================================================
%  WAHN Phase 0: MATLAB Data Exploration
%  ==========================================================================
%  Author: Nikos Diederichs
%  Date: January 2025
%
%  Purpose: Explore WAHN data structure and NESDA reference before preprocessing
%
%  WICHTIG: Führe ZUERST phase0_exploration.sh aus!
%           Dieses Skript analysiert die Datenstrukturen im Detail.
%
%  Usage: Run in MATLAB on HPC
%  ==========================================================================

clear; clc;
fprintf('==========================================================================\n');
fprintf('  WAHN Phase 0: MATLAB Data Exploration\n');
fprintf('  Date: %s\n', datestr(now));
fprintf('==========================================================================\n\n');

%% Configuration Paths
WAHN_MRI_BASE    = '/volume/data/WAHN/MRI/19-Nov-2025/Data/Preprocessed/';
WAHN_EXCEL_DEMO  = '/volume/data/WAHN/DataDump/19-Nov-2025/WAHN_tabular_combined_data.xlsx';
WAHN_EXCEL_CAT12 = '/volume/data/WAHN/DataDump/19-Nov-2025/WAHN_baseline_data_cat12_r1207_tabular.xlsx';
NESDA_NM_PAT     = '/volume/projects/CV_NESDA/Data/NESDA_Waves/Wave_1/NM_Structures/NM_Pat.mat';
NESDA_NM_HC      = '/volume/projects/CV_NESDA/Data/NESDA_Waves/Wave_1/NM_Structures/NM_HC.mat';
BRAIN_MASK       = '/volume/projects/CV_NESDA/Analysis/bvFTD/Mask/brainmask_T1_2mm.nii';

%% ==========================================================================
%  SECTION 1: NESDA NM STRUCTURE REFERENCE
%  ==========================================================================
fprintf('\n========== SECTION 1: NESDA NM STRUCTURE REFERENCE ==========\n\n');

try
    fprintf('Loading NESDA NM_Pat.mat...\n');
    nesda_pat = load(NESDA_NM_PAT);

    fprintf('\nNM structure fields:\n');
    fields = fieldnames(nesda_pat.NM);
    for i = 1:length(fields)
        fprintf('  - %s\n', fields{i});
    end

    fprintf('\nY cell array contents:\n');
    if iscell(nesda_pat.NM.Y)
        for i = 1:length(nesda_pat.NM.Y)
            fprintf('  Y{%d}: size = [%d x %d]\n', i, size(nesda_pat.NM.Y{i}, 1), size(nesda_pat.NM.Y{i}, 2));
        end
    else
        fprintf('  Y: size = [%d x %d]\n', size(nesda_pat.NM.Y, 1), size(nesda_pat.NM.Y, 2));
    end

    % Check for other relevant fields
    if isfield(nesda_pat.NM, 'label')
        fprintf('\nLabel field:\n');
        fprintf('  Size: [%d x %d]\n', size(nesda_pat.NM.label, 1), size(nesda_pat.NM.label, 2));
        fprintf('  Unique values: %s\n', mat2str(unique(nesda_pat.NM.label)));
    end

    if isfield(nesda_pat.NM, 'id')
        fprintf('\nID field:\n');
        fprintf('  Size: [%d x %d]\n', size(nesda_pat.NM.id, 1), size(nesda_pat.NM.id, 2));
        if iscell(nesda_pat.NM.id)
            fprintf('  First 5 IDs: %s\n', strjoin(nesda_pat.NM.id(1:min(5,length(nesda_pat.NM.id))), ', '));
        end
    end

catch ME
    fprintf('ERROR loading NESDA NM_Pat: %s\n', ME.message);
end

% Also load HC for reference
try
    fprintf('\n--- Loading NESDA NM_HC.mat ---\n');
    nesda_hc = load(NESDA_NM_HC);
    if iscell(nesda_hc.NM.Y)
        fprintf('NM_HC Y{1} size: [%d x %d]\n', size(nesda_hc.NM.Y{1}, 1), size(nesda_hc.NM.Y{1}, 2));
    end
catch ME
    fprintf('ERROR loading NESDA NM_HC: %s\n', ME.message);
end

%% ==========================================================================
%  SECTION 2: BRAIN MASK ANALYSIS
%  ==========================================================================
fprintf('\n========== SECTION 2: BRAIN MASK ANALYSIS ==========\n\n');

try
    fprintf('Loading brain mask: %s\n', BRAIN_MASK);
    mask_hdr = spm_vol(BRAIN_MASK);
    mask_data = spm_read_vols(mask_hdr);

    fprintf('\nMask dimensions: [%d x %d x %d]\n', size(mask_data, 1), size(mask_data, 2), size(mask_data, 3));
    fprintf('Mask voxel size: [%.2f x %.2f x %.2f] mm\n', mask_hdr.mat(1,1), mask_hdr.mat(2,2), mask_hdr.mat(3,3));

    mask_idx = find(mask_data > 0);
    fprintf('Number of in-mask voxels: %d\n', length(mask_idx));
    fprintf('EXPECTED for Clara model: 71276\n');

    if length(mask_idx) == 71276
        fprintf('MATCH: Mask has correct number of features.\n');
    else
        fprintf('WARNING: Feature count mismatch!\n');
    end
catch ME
    fprintf('ERROR loading brain mask: %s\n', ME.message);
    fprintf('Make sure SPM is in MATLAB path.\n');
end

%% ==========================================================================
%  SECTION 3: WAHN EXCEL 1 - DEMOGRAPHICS
%  ==========================================================================
fprintf('\n========== SECTION 3: WAHN DEMOGRAPHICS (Excel 1) ==========\n\n');

try
    fprintf('Loading: %s\n', WAHN_EXCEL_DEMO);
    T1 = readtable(WAHN_EXCEL_DEMO);

    fprintf('\nTable size: %d rows x %d columns\n', height(T1), width(T1));

    fprintf('\nColumn names:\n');
    for i = 1:width(T1)
        fprintf('  %d. %s\n', i, T1.Properties.VariableNames{i});
    end

    fprintf('\nFirst 10 rows (selected columns):\n');
    % Try to display key columns
    if ismember('pident', T1.Properties.VariableNames)
        disp(T1(1:min(10, height(T1)), :));
    end

    % Analyze Age
    if ismember('Age', T1.Properties.VariableNames)
        fprintf('\nAge statistics:\n');
        fprintf('  N valid: %d\n', sum(~isnan(T1.Age)));
        fprintf('  Range: %.1f - %.1f years\n', min(T1.Age), max(T1.Age));
        fprintf('  Mean (SD): %.1f (%.1f)\n', nanmean(T1.Age), nanstd(T1.Age));
    end

    % Analyze Sex
    if ismember('Sexe', T1.Properties.VariableNames)
        fprintf('\nSex distribution:\n');
        sex_vals = T1.Sexe;
        fprintf('  Value 1: %d\n', sum(sex_vals == 1));
        fprintf('  Value 2: %d\n', sum(sex_vals == 2));
        fprintf('  Missing: %d\n', sum(isnan(sex_vals)));
    end

catch ME
    fprintf('ERROR loading Excel 1: %s\n', ME.message);
end

%% ==========================================================================
%  SECTION 4: WAHN EXCEL 2 - CAT12 METADATA
%  ==========================================================================
fprintf('\n========== SECTION 4: WAHN CAT12 METADATA (Excel 2) ==========\n\n');

try
    fprintf('Loading: %s\n', WAHN_EXCEL_CAT12);
    T2 = readtable(WAHN_EXCEL_CAT12);

    fprintf('\nTable size: %d rows x %d columns\n', height(T2), width(T2));

    fprintf('\nColumn names:\n');
    for i = 1:width(T2)
        fprintf('  %d. %s\n', i, T2.Properties.VariableNames{i});
    end

    fprintf('\nFirst 10 rows:\n');
    disp(T2(1:min(10, height(T2)), :));

    % Analyze TIV format
    if ismember('TIV', T2.Properties.VariableNames)
        fprintf('\n--- TIV FORMAT ANALYSIS ---\n');
        fprintf('TIV column type: %s\n', class(T2.TIV));

        fprintf('\nFirst 5 raw TIV values:\n');
        for i = 1:min(5, height(T2))
            if iscell(T2.TIV)
                fprintf('  %d: "%s"\n', i, T2.TIV{i});
            elseif isstring(T2.TIV)
                fprintf('  %d: "%s"\n', i, T2.TIV(i));
            else
                fprintf('  %d: %g\n', i, T2.TIV(i));
            end
        end

        % Try to convert if string/cell
        if iscell(T2.TIV) || isstring(T2.TIV)
            fprintf('\nAttempting TIV conversion (German decimal format)...\n');
            tiv_str = string(T2.TIV);
            tiv_str = strrep(tiv_str, ',', '.');
            tiv_numeric = str2double(tiv_str);

            fprintf('Converted TIV range: %.1f - %.1f\n', min(tiv_numeric), max(tiv_numeric));
            fprintf('Expected TIV range: 1200 - 1800 ml\n');

            if min(tiv_numeric) > 50000000
                fprintf('\nWARNING: TIV values seem to be in mm^3, not ml!\n');
                fprintf('To convert: TIV_ml = TIV_mm3 / 1000\n');
                fprintf('Converted range: %.1f - %.1f ml\n', min(tiv_numeric)/1000, max(tiv_numeric)/1000);
            end
        end
    end

    % Check mwp1 path format
    if ismember('CAT12_mwp1_name_unzipped', T2.Properties.VariableNames)
        fprintf('\n--- MWP1 PATH ANALYSIS ---\n');
        fprintf('First 5 mwp1 paths from Excel:\n');
        for i = 1:min(5, height(T2))
            if iscell(T2.CAT12_mwp1_name_unzipped)
                fprintf('  %d: %s\n', i, T2.CAT12_mwp1_name_unzipped{i});
            else
                fprintf('  %d: %s\n', i, string(T2.CAT12_mwp1_name_unzipped(i)));
            end
        end

        % Check if paths point to local Mac
        sample_path = string(T2.CAT12_mwp1_name_unzipped(1));
        if contains(sample_path, '/Users/') || contains(sample_path, 'C:\')
            fprintf('\nWARNING: Paths point to LOCAL machine, not HPC!\n');
            fprintf('These paths will not work. Need to find files on HPC.\n');
        end
    end

    % Check for duplicates
    if ismember('PSN', T2.Properties.VariableNames)
        fprintf('\n--- DUPLICATE CHECK ---\n');
        [unique_psn, ~, ic] = unique(T2.PSN);
        counts = accumarray(ic, 1);
        dups = unique_psn(counts > 1);

        fprintf('Total rows: %d\n', height(T2));
        fprintf('Unique PSN: %d\n', length(unique_psn));
        fprintf('Duplicates: %d\n', length(dups));

        if ~isempty(dups)
            fprintf('Duplicate PSNs:\n');
            for i = 1:length(dups)
                fprintf('  - %s (appears %dx)\n', string(dups(i)), counts(strcmp(unique_psn, dups(i))));
            end
        end
    end

catch ME
    fprintf('ERROR loading Excel 2: %s\n', ME.message);
end

%% ==========================================================================
%  SECTION 5: SEARCH FOR MWP1 FILES ON HPC
%  ==========================================================================
fprintf('\n========== SECTION 5: SEARCH FOR MWP1 FILES ==========\n\n');

fprintf('Searching for mwp1 files in: %s\n', WAHN_MRI_BASE);

% Method 1: Direct search with dir
mwp1_files = dir(fullfile(WAHN_MRI_BASE, '**', 'mwp1*.nii'));
fprintf('\nMethod 1 (dir mwp1*.nii): Found %d files\n', length(mwp1_files));

if ~isempty(mwp1_files)
    fprintf('First 10 mwp1 files:\n');
    for i = 1:min(10, length(mwp1_files))
        fprintf('  %d: %s\n', i, fullfile(mwp1_files(i).folder, mwp1_files(i).name));
    end
end

% Method 2: Broader search
mwp_files = dir(fullfile(WAHN_MRI_BASE, '**', 'mwp*.nii'));
fprintf('\nMethod 2 (dir mwp*.nii): Found %d files\n', length(mwp_files));

% Method 3: Check for mri subdirectories
mri_dirs = dir(fullfile(WAHN_MRI_BASE, '**', 'mri'));
mri_dirs = mri_dirs([mri_dirs.isdir]);
fprintf('\nMethod 3 (mri directories): Found %d directories\n', length(mri_dirs));

if ~isempty(mri_dirs)
    fprintf('Checking first mri directory for contents:\n');
    first_mri = fullfile(mri_dirs(1).folder, mri_dirs(1).name);
    mri_contents = dir(first_mri);
    for i = 1:length(mri_contents)
        if ~strcmp(mri_contents(i).name, '.') && ~strcmp(mri_contents(i).name, '..')
            fprintf('  - %s\n', mri_contents(i).name);
        end
    end
end

%% ==========================================================================
%  SECTION 6: TEST LOAD ONE MWP1 FILE (IF FOUND)
%  ==========================================================================
fprintf('\n========== SECTION 6: TEST LOAD MWP1 FILE ==========\n\n');

if ~isempty(mwp1_files)
    test_file = fullfile(mwp1_files(1).folder, mwp1_files(1).name);
    fprintf('Testing load of: %s\n', test_file);

    try
        vol_hdr = spm_vol(test_file);
        vol_data = spm_read_vols(vol_hdr);

        fprintf('\nFile dimensions: [%d x %d x %d]\n', size(vol_data, 1), size(vol_data, 2), size(vol_data, 3));
        fprintf('Voxel size: [%.2f x %.2f x %.2f] mm\n', abs(vol_hdr.mat(1,1)), abs(vol_hdr.mat(2,2)), abs(vol_hdr.mat(3,3)));
        fprintf('Data range: [%.4f, %.4f]\n', min(vol_data(:)), max(vol_data(:)));
        fprintf('Mean (non-zero): %.4f\n', mean(vol_data(vol_data > 0)));

        % Test masking
        if exist('mask_idx', 'var')
            masked_data = vol_data(mask_idx);
            fprintf('\nMasked data: %d values\n', length(masked_data));
            fprintf('Masked data range: [%.4f, %.4f]\n', min(masked_data), max(masked_data));
        end

    catch ME
        fprintf('ERROR loading test file: %s\n', ME.message);
    end
else
    fprintf('NO MWP1 FILES FOUND - Cannot test loading.\n');
    fprintf('\nACTION REQUIRED:\n');
    fprintf('  1. Contact Maja to upload mwp1 files to HPC\n');
    fprintf('  2. Or check if files are in a different location\n');
end

%% ==========================================================================
%  SECTION 7: SUMMARY AND NEXT STEPS
%  ==========================================================================
fprintf('\n========== SUMMARY AND NEXT STEPS ==========\n\n');

fprintf('CHECKLIST:\n');
fprintf('  [ ] mwp1 files found on HPC: ');
if ~isempty(mwp1_files)
    fprintf('YES (%d files)\n', length(mwp1_files));
else
    fprintf('NO - ACTION REQUIRED\n');
end

fprintf('  [ ] Brain mask features: ');
if exist('mask_idx', 'var')
    if length(mask_idx) == 71276
        fprintf('CORRECT (71276)\n');
    else
        fprintf('MISMATCH (%d vs 71276)\n', length(mask_idx));
    end
else
    fprintf('NOT CHECKED\n');
end

fprintf('  [ ] Demographics Excel: ');
if exist('T1', 'var')
    fprintf('LOADED (%d subjects)\n', height(T1));
else
    fprintf('ERROR\n');
end

fprintf('  [ ] CAT12 Excel: ');
if exist('T2', 'var')
    fprintf('LOADED (%d rows)\n', height(T2));
else
    fprintf('ERROR\n');
end

fprintf('\n');
if isempty(mwp1_files)
    fprintf('!!! CRITICAL: Cannot proceed without mwp1 files !!!\n');
    fprintf('Please contact Maja for file upload before running main script.\n');
else
    fprintf('Ready to proceed with: prepare_wahn_for_clara.m\n');
end

fprintf('\n==========================================================================\n');
fprintf('Exploration complete: %s\n', datestr(now));
fprintf('==========================================================================\n');
