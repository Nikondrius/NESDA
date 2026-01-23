%% ========================================================================
%  WAHN Data Preparation for Transition Model Application
%  ========================================================================
%  Author: Nikos Diederichs
%  Date: January 2025
%
%  Purpose: Create OOCV containers for Clara to apply Transition Model
%
%  WICHTIG: Führe ZUERST phase0_exploration.sh und phase0_exploration_matlab.m aus!
%
%  Outputs:
%    1. WAHN_Y_raw.mat           - Raw brain data matrix [N x 71276]
%    2. WAHN_Y_corrected.mat     - After Mean Offset Correction
%    3. WAHN_covariates.mat      - Age, Sex, TIV, Site (for ComBat)
%    4. WAHN_subject_table.csv   - Subject info for documentation
%    5. WAHN_validation_report.txt
%  ========================================================================

clear; clc; close all;

fprintf('========================================================================\n');
fprintf('  WAHN Data Preparation for Clara - Transition Model Application\n');
fprintf('  Date: %s\n', datestr(now));
fprintf('========================================================================\n\n');

%% ========================================================================
%  CONFIGURATION - PATHS
%  ========================================================================
fprintf('========== CONFIGURATION ==========\n\n');

% Input Paths - WAHN Data
WAHN_MRI_BASE    = '/volume/data/WAHN/MRI/19-Nov-2025/Data/Preprocessed/';
WAHN_EXCEL_DEMO  = '/volume/data/WAHN/DataDump/19-Nov-2025/WAHN_tabular_combined_data.xlsx';
WAHN_EXCEL_CAT12 = '/volume/data/WAHN/DataDump/19-Nov-2025/WAHN_baseline_data_cat12_r1207_tabular.xlsx';

% Reference - NESDA (for structure and Mean Offset)
NESDA_NM_PAT     = '/volume/projects/CV_NESDA/Data/NESDA_Waves/Wave_1/NM_Structures/NM_Pat.mat';
NESDA_NM_HC      = '/volume/projects/CV_NESDA/Data/NESDA_Waves/Wave_1/NM_Structures/NM_HC.mat';

% Brain Mask (SAME as NESDA - critical for feature alignment!)
BRAIN_MASK       = '/volume/projects/CV_NESDA/Analysis/bvFTD/Mask/brainmask_T1_2mm.nii';

% Output Directory
OUTPUT_DIR       = '/volume/projects/CV_NESDA/Data/WAHN/';

% Expected feature count (MUST match NESDA model!)
EXPECTED_FEATURES = 71276;

% Print configuration
fprintf('Input Paths:\n');
fprintf('  WAHN MRI Base:    %s\n', WAHN_MRI_BASE);
fprintf('  WAHN Demographics: %s\n', WAHN_EXCEL_DEMO);
fprintf('  WAHN CAT12 Excel:  %s\n', WAHN_EXCEL_CAT12);
fprintf('\nReference Paths:\n');
fprintf('  NESDA NM_Pat:     %s\n', NESDA_NM_PAT);
fprintf('  NESDA NM_HC:      %s\n', NESDA_NM_HC);
fprintf('  Brain Mask:       %s\n', BRAIN_MASK);
fprintf('\nOutput Directory:   %s\n', OUTPUT_DIR);
fprintf('Expected Features:  %d\n', EXPECTED_FEATURES);

% Create output directory if needed
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
    fprintf('\nCreated output directory: %s\n', OUTPUT_DIR);
end

%% ========================================================================
%  MODULE A: DATA DISCOVERY - FIND MWP1 FILES
%  ========================================================================
fprintf('\n========== MODULE A: DATA DISCOVERY ==========\n\n');

% A.1: Search for mwp1 files
fprintf('Searching for mwp1 files...\n');

% Primary search pattern
mwp1_files = dir(fullfile(WAHN_MRI_BASE, '**', 'mwp1*.nii'));

% Alternative pattern if primary fails
if isempty(mwp1_files)
    fprintf('No mwp1*.nii found, trying mwp*.nii...\n');
    mwp1_files = dir(fullfile(WAHN_MRI_BASE, '**', 'mwp*.nii'));
end

% Check in /mri/ subdirectories specifically
if isempty(mwp1_files)
    fprintf('Searching in /mri/ subdirectories...\n');
    mwp1_files = dir(fullfile(WAHN_MRI_BASE, '*', 'mri', 'mwp1*.nii'));
end

% Critical check
if isempty(mwp1_files)
    error(['\n!!! CRITICAL: No mwp1*.nii files found !!!\n\n' ...
           'Searched in: %s\n\n' ...
           'Action required:\n' ...
           '  1. Contact Maja - mwp1 files may be on her local Mac\n' ...
           '  2. Excel 2 shows paths like: /Users/maja/Desktop/...\n' ...
           '  3. Files need to be uploaded to HPC before proceeding\n\n' ...
           'STOP: Cannot proceed until mwp1 files are available.\n'], WAHN_MRI_BASE);
end

fprintf('Found %d mwp1 files\n', length(mwp1_files));

% A.2: Extract Subject IDs and paths
wahn_ids = cell(length(mwp1_files), 1);
wahn_paths = cell(length(mwp1_files), 1);

for i = 1:length(mwp1_files)
    wahn_paths{i} = fullfile(mwp1_files(i).folder, mwp1_files(i).name);

    % Extract ID: mwp1sub-XXXXX_run-01_T1w.nii -> XXXXX
    % Handle various naming patterns
    tokens = regexp(mwp1_files(i).name, 'sub-([A-Za-z0-9]+)', 'tokens');
    if ~isempty(tokens)
        wahn_ids{i} = tokens{1}{1};
    else
        % Fallback: use filename without extension
        [~, fname, ~] = fileparts(mwp1_files(i).name);
        wahn_ids{i} = fname;
    end
end

% A.3: Remove duplicates (keep first occurrence)
[unique_ids, unique_idx] = unique(wahn_ids, 'stable');
n_duplicates = length(wahn_ids) - length(unique_ids);

if n_duplicates > 0
    fprintf('\nWARNING: Found %d duplicate Subject IDs\n', n_duplicates);
    fprintf('Keeping first occurrence for each:\n');

    % Find and report duplicates
    [~, ~, ic] = unique(wahn_ids, 'stable');
    counts = accumarray(ic, 1);
    dup_idx = find(counts > 1);

    for i = 1:length(dup_idx)
        dup_id = unique_ids{dup_idx(i)};
        fprintf('  - %s (appeared %dx)\n', dup_id, counts(dup_idx(i)));
    end

    wahn_ids = wahn_ids(unique_idx);
    wahn_paths = wahn_paths(unique_idx);
end

fprintf('\nUnique subjects: %d\n', length(wahn_ids));

% Display first 10 IDs
fprintf('\nFirst 10 Subject IDs:\n');
for i = 1:min(10, length(wahn_ids))
    fprintf('  %d. %s\n', i, wahn_ids{i});
end

%% ========================================================================
%  MODULE B: LOAD DEMOGRAPHICS
%  ========================================================================
fprintf('\n========== MODULE B: LOAD DEMOGRAPHICS ==========\n\n');

% B.1: Load Excel 1 (main demographics)
fprintf('Loading demographics from Excel 1...\n');
try
    demo_table = readtable(WAHN_EXCEL_DEMO);
    fprintf('Loaded: %d rows, %d columns\n', height(demo_table), width(demo_table));
    fprintf('Columns: %s\n', strjoin(demo_table.Properties.VariableNames, ', '));
catch ME
    warning('Could not load demographics Excel: %s', ME.message);
    demo_table = table();
end

% B.2: Load Excel 2 (CAT12 metadata with TIV)
fprintf('\nLoading CAT12 metadata from Excel 2...\n');
try
    cat12_table = readtable(WAHN_EXCEL_CAT12);
    fprintf('Loaded: %d rows, %d columns\n', height(cat12_table), width(cat12_table));
    fprintf('Columns: %s\n', strjoin(cat12_table.Properties.VariableNames, ', '));

    % Fix TIV format (German decimal: "69,994" -> 69.994)
    if ismember('TIV', cat12_table.Properties.VariableNames)
        if iscell(cat12_table.TIV) || isstring(cat12_table.TIV)
            fprintf('\nConverting TIV from string format...\n');
            tiv_str = string(cat12_table.TIV);
            tiv_str = strrep(tiv_str, ',', '.');
            cat12_table.TIV = str2double(tiv_str);
        end

        % Check if values are in mm^3 instead of ml
        median_tiv = nanmedian(cat12_table.TIV);
        if median_tiv > 100000
            fprintf('TIV appears to be in mm^3, converting to ml...\n');
            cat12_table.TIV = cat12_table.TIV / 1000;
        end

        fprintf('TIV range: %.1f - %.1f ml\n', nanmin(cat12_table.TIV), nanmax(cat12_table.TIV));

        % Validate range
        if nanmin(cat12_table.TIV) < 500 || nanmax(cat12_table.TIV) > 2500
            warning('TIV values outside expected range (500-2500 ml)! Check format.');
        end
    end
catch ME
    warning('Could not load CAT12 Excel: %s', ME.message);
    cat12_table = table();
end

%% ========================================================================
%  MODULE C: LOAD REFERENCE DATA (NESDA & BRAIN MASK)
%  ========================================================================
fprintf('\n========== MODULE C: LOAD REFERENCE DATA ==========\n\n');

% C.1: Load NESDA NM structure for reference
fprintf('Loading NESDA NM_Pat for structure reference...\n');
try
    nesda_pat = load(NESDA_NM_PAT);
    fprintf('NESDA NM structure loaded\n');
    fprintf('  Fields: %s\n', strjoin(fieldnames(nesda_pat.NM), ', '));

    if iscell(nesda_pat.NM.Y)
        nesda_y_size = size(nesda_pat.NM.Y{1});
    else
        nesda_y_size = size(nesda_pat.NM.Y);
    end
    fprintf('  Y size: [%d x %d]\n', nesda_y_size(1), nesda_y_size(2));
catch ME
    warning('Could not load NESDA NM_Pat: %s', ME.message);
end

% C.2: Load NESDA HC for mean offset reference (optional)
fprintf('\nLoading NESDA NM_HC for mean offset reference...\n');
try
    nesda_hc = load(NESDA_NM_HC);
    if iscell(nesda_hc.NM.Y)
        nesda_hc_mean = mean(nesda_hc.NM.Y{1}, 1);
    else
        nesda_hc_mean = mean(nesda_hc.NM.Y, 1);
    end
    fprintf('NESDA HC mean computed: [1 x %d]\n', length(nesda_hc_mean));
    USE_NESDA_HC_MEAN = true;
catch ME
    warning('Could not load NESDA NM_HC: %s', ME.message);
    USE_NESDA_HC_MEAN = false;
end

% C.3: Load Brain Mask
fprintf('\nLoading brain mask...\n');
try
    mask_hdr = spm_vol(BRAIN_MASK);
    mask_data = spm_read_vols(mask_hdr);
    mask_idx = find(mask_data > 0);
    n_features = length(mask_idx);

    fprintf('Brain mask loaded\n');
    fprintf('  Dimensions: [%d x %d x %d]\n', size(mask_data, 1), size(mask_data, 2), size(mask_data, 3));
    fprintf('  In-mask voxels: %d\n', n_features);

    if n_features ~= EXPECTED_FEATURES
        error('Feature count mismatch! Expected %d, got %d. Check brain mask.', ...
              EXPECTED_FEATURES, n_features);
    end
    fprintf('  Feature count matches expected: %d\n', EXPECTED_FEATURES);
catch ME
    error('Failed to load brain mask: %s\nMake sure SPM is in MATLAB path.', ME.message);
end

%% ========================================================================
%  MODULE D: LOAD WAHN BRAIN DATA
%  ========================================================================
fprintf('\n========== MODULE D: LOAD BRAIN DATA ==========\n\n');

n_subjects = length(wahn_paths);
Y_raw = zeros(n_subjects, n_features);
load_success = true(n_subjects, 1);
load_errors = cell(n_subjects, 1);

fprintf('Loading %d brain volumes...\n', n_subjects);
fprintf('Progress: ');

tic;
for i = 1:n_subjects
    % Progress indicator
    if mod(i, 10) == 0
        fprintf('%d..', i);
    end

    try
        vol_hdr = spm_vol(wahn_paths{i});
        vol_data = spm_read_vols(vol_hdr);
        Y_raw(i, :) = vol_data(mask_idx)';

        % Basic QC: check for NaN or negative values
        if any(isnan(Y_raw(i, :)))
            warning('Subject %s has NaN values', wahn_ids{i});
        end

    catch ME
        load_success(i) = false;
        load_errors{i} = ME.message;
    end
end
load_time = toc;
fprintf('Done!\n');
fprintf('Loading time: %.1f seconds (%.2f s/subject)\n', load_time, load_time/n_subjects);

% Handle failed subjects
n_failed = sum(~load_success);
if n_failed > 0
    fprintf('\nWARNING: %d subjects failed to load:\n', n_failed);
    failed_idx = find(~load_success);
    for i = 1:length(failed_idx)
        idx = failed_idx(i);
        fprintf('  - %s: %s\n', wahn_ids{idx}, load_errors{idx});
    end

    fprintf('\nRemoving failed subjects from data matrix...\n');
    Y_raw = Y_raw(load_success, :);
    wahn_ids = wahn_ids(load_success);
    wahn_paths = wahn_paths(load_success);
    n_subjects = length(wahn_ids);
end

fprintf('\nBrain data loaded: [%d subjects x %d features]\n', size(Y_raw, 1), size(Y_raw, 2));
fprintf('Data range: [%.4f, %.4f]\n', min(Y_raw(:)), max(Y_raw(:)));
fprintf('Mean (all): %.4f\n', mean(Y_raw(:)));

%% ========================================================================
%  MODULE E: MEAN OFFSET CORRECTION
%  ========================================================================
fprintf('\n========== MODULE E: MEAN OFFSET CORRECTION ==========\n\n');

% Clara: "es ist einfach eine voxelweise Mittelwert Korrektur"
% Two options:
%   Option 1: Use WAHN mean as reference (within-cohort normalization)
%   Option 2: Use NESDA HC mean as reference (cross-cohort, more comparable)

fprintf('Computing Mean Offset Correction...\n');
fprintf('  Method: Voxel-wise mean subtraction\n');

% Compute WAHN mean (always computed for reporting)
mean_wahn = mean(Y_raw, 1);  % [1 x n_features]
fprintf('  WAHN mean computed: range [%.4f, %.4f]\n', min(mean_wahn), max(mean_wahn));

% Choose reference mean
if USE_NESDA_HC_MEAN
    fprintf('\n  Using NESDA HC mean as reference (cross-cohort correction)\n');
    mean_reference = nesda_hc_mean;
else
    fprintf('\n  Using WAHN mean as reference (within-cohort correction)\n');
    mean_reference = mean_wahn;
end

% Apply correction: Y_corrected = Y_raw - mean_reference
Y_corrected = Y_raw - mean_reference;

fprintf('\nCorrection applied:\n');
fprintf('  Raw data:       mean=%.6f, std=%.4f\n', mean(Y_raw(:)), std(Y_raw(:)));
fprintf('  Corrected data: mean=%.6f, std=%.4f\n', mean(Y_corrected(:)), std(Y_corrected(:)));

% Verify column means are approximately zero (for within-cohort)
if ~USE_NESDA_HC_MEAN
    col_means_corrected = mean(Y_corrected, 1);
    fprintf('  Column means after correction: mean=%.2e (should be ~0)\n', mean(col_means_corrected));
end

%% ========================================================================
%  MODULE F: MATCH AND PREPARE COVARIATES
%  ========================================================================
fprintf('\n========== MODULE F: PREPARE COVARIATES ==========\n\n');

% Initialize covariate arrays
age = nan(n_subjects, 1);
sex = nan(n_subjects, 1);
tiv = nan(n_subjects, 1);
site = ones(n_subjects, 1);  % WAHN is single-site, so all = 1

n_matched_demo = 0;
n_matched_tiv = 0;

fprintf('Matching subjects to demographics...\n');

for i = 1:n_subjects
    subj_id = wahn_ids{i};

    % Try matching in demo_table (Excel 1)
    if ~isempty(demo_table) && ismember('pident', demo_table.Properties.VariableNames)
        % Handle various ID column types
        if iscell(demo_table.pident)
            idx_demo = find(strcmp(demo_table.pident, subj_id));
        elseif isstring(demo_table.pident)
            idx_demo = find(strcmp(string(demo_table.pident), subj_id));
        else
            idx_demo = find(demo_table.pident == str2double(subj_id));
        end

        if ~isempty(idx_demo)
            idx_demo = idx_demo(1);  % Take first match

            if ismember('Age', demo_table.Properties.VariableNames)
                age(i) = demo_table.Age(idx_demo);
            end
            if ismember('Sexe', demo_table.Properties.VariableNames)
                sex(i) = demo_table.Sexe(idx_demo);
            end
            n_matched_demo = n_matched_demo + 1;
        end
    end

    % Try matching in cat12_table (Excel 2) for TIV
    if ~isempty(cat12_table) && ismember('PSN', cat12_table.Properties.VariableNames)
        if iscell(cat12_table.PSN)
            idx_cat12 = find(strcmp(cat12_table.PSN, subj_id));
        elseif isstring(cat12_table.PSN)
            idx_cat12 = find(strcmp(string(cat12_table.PSN), subj_id));
        else
            idx_cat12 = find(cat12_table.PSN == str2double(subj_id));
        end

        if ~isempty(idx_cat12) && ismember('TIV', cat12_table.Properties.VariableNames)
            idx_cat12 = idx_cat12(1);  % Take first match
            tiv(i) = cat12_table.TIV(idx_cat12);
            n_matched_tiv = n_matched_tiv + 1;
        end
    end
end

fprintf('Demographics matching:\n');
fprintf('  Age/Sex matched: %d/%d (%.1f%%)\n', n_matched_demo, n_subjects, 100*n_matched_demo/n_subjects);
fprintf('  TIV matched:     %d/%d (%.1f%%)\n', n_matched_tiv, n_subjects, 100*n_matched_tiv/n_subjects);

% Report covariate statistics
fprintf('\nCovariate statistics:\n');
fprintf('  Age:  N=%d, mean=%.1f, range=[%.0f, %.0f], missing=%d\n', ...
        sum(~isnan(age)), nanmean(age), nanmin(age), nanmax(age), sum(isnan(age)));
fprintf('  Sex:  male=%d, female=%d, missing=%d\n', ...
        sum(sex==1), sum(sex==2), sum(isnan(sex)));
fprintf('  TIV:  N=%d, mean=%.0f, range=[%.0f, %.0f], missing=%d\n', ...
        sum(~isnan(tiv)), nanmean(tiv), nanmin(tiv), nanmax(tiv), sum(isnan(tiv)));

% Create covariate matrix
covariates = [age, sex, tiv, site];
covariate_names = {'Age', 'Sex', 'TIV', 'Site'};

%% ========================================================================
%  MODULE G: SAVE OUTPUTS
%  ========================================================================
fprintf('\n========== MODULE G: SAVE OUTPUTS ==========\n\n');

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');

% G.1: Save Raw Y matrix
Y = Y_raw;  %#ok<NASGU>
outfile_raw = fullfile(OUTPUT_DIR, 'WAHN_Y_raw.mat');
save(outfile_raw, 'Y', '-v7.3');
fprintf('Saved: WAHN_Y_raw.mat [%d x %d]\n', size(Y_raw, 1), size(Y_raw, 2));

% G.2: Save Corrected Y matrix
Y = Y_corrected;  %#ok<NASGU>
outfile_corr = fullfile(OUTPUT_DIR, 'WAHN_Y_corrected.mat');
save(outfile_corr, 'Y', '-v7.3');
fprintf('Saved: WAHN_Y_corrected.mat [%d x %d]\n', size(Y_corrected, 1), size(Y_corrected, 2));

% G.3: Save Covariates
outfile_cov = fullfile(OUTPUT_DIR, 'WAHN_covariates.mat');
save(outfile_cov, 'covariates', 'covariate_names', 'age', 'sex', 'tiv', 'site');
fprintf('Saved: WAHN_covariates.mat\n');

% G.4: Save Subject Table
subject_table = table(wahn_ids, wahn_paths, age, sex, tiv, site, ...
    'VariableNames', {'SubjectID', 'MRI_Path', 'Age', 'Sex', 'TIV', 'Site'});
outfile_table = fullfile(OUTPUT_DIR, 'WAHN_subject_table.csv');
writetable(subject_table, outfile_table);
fprintf('Saved: WAHN_subject_table.csv\n');

% G.5: Save Mean Reference (for documentation)
outfile_mean = fullfile(OUTPUT_DIR, 'WAHN_mean_reference.mat');
mean_correction_method = 'NESDA_HC';
if ~USE_NESDA_HC_MEAN
    mean_correction_method = 'WAHN_within';
end
save(outfile_mean, 'mean_wahn', 'mean_reference', 'mean_correction_method');
fprintf('Saved: WAHN_mean_reference.mat (method: %s)\n', mean_correction_method);

% G.6: Save Validation Report
outfile_report = fullfile(OUTPUT_DIR, 'WAHN_validation_report.txt');
fid = fopen(outfile_report, 'w');

fprintf(fid, '========================================================================\n');
fprintf(fid, 'WAHN Data Preparation - Validation Report\n');
fprintf(fid, '========================================================================\n');
fprintf(fid, 'Date: %s\n', datestr(now));
fprintf(fid, 'Author: Nikos Diederichs\n');
fprintf(fid, '\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'DATA SUMMARY\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'Total subjects:     %d\n', n_subjects);
fprintf(fid, 'Features (voxels):  %d\n', n_features);
fprintf(fid, 'Expected features:  %d\n', EXPECTED_FEATURES);
fprintf(fid, 'Feature match:      %s\n', iif(n_features == EXPECTED_FEATURES, 'YES', 'NO'));
fprintf(fid, '\n');
fprintf(fid, 'Failed to load:     %d\n', n_failed);
fprintf(fid, 'Demographics match: %d/%d (%.1f%%)\n', n_matched_demo, n_subjects, 100*n_matched_demo/n_subjects);
fprintf(fid, 'TIV match:          %d/%d (%.1f%%)\n', n_matched_tiv, n_subjects, 100*n_matched_tiv/n_subjects);
fprintf(fid, '\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'DATA QUALITY\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'Y_raw:\n');
fprintf(fid, '  Dimensions: [%d x %d]\n', size(Y_raw, 1), size(Y_raw, 2));
fprintf(fid, '  Mean:       %.6f\n', mean(Y_raw(:)));
fprintf(fid, '  Std:        %.6f\n', std(Y_raw(:)));
fprintf(fid, '  Range:      [%.4f, %.4f]\n', min(Y_raw(:)), max(Y_raw(:)));
fprintf(fid, '  NaN count:  %d\n', sum(isnan(Y_raw(:))));
fprintf(fid, '\n');
fprintf(fid, 'Y_corrected:\n');
fprintf(fid, '  Dimensions: [%d x %d]\n', size(Y_corrected, 1), size(Y_corrected, 2));
fprintf(fid, '  Mean:       %.6f\n', mean(Y_corrected(:)));
fprintf(fid, '  Std:        %.6f\n', std(Y_corrected(:)));
fprintf(fid, '  Range:      [%.4f, %.4f]\n', min(Y_corrected(:)), max(Y_corrected(:)));
fprintf(fid, '  Correction: %s\n', mean_correction_method);
fprintf(fid, '\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'COVARIATES\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'Age:\n');
fprintf(fid, '  Valid:   %d\n', sum(~isnan(age)));
fprintf(fid, '  Missing: %d\n', sum(isnan(age)));
fprintf(fid, '  Mean:    %.1f years\n', nanmean(age));
fprintf(fid, '  Range:   [%.0f, %.0f]\n', nanmin(age), nanmax(age));
fprintf(fid, '\n');
fprintf(fid, 'Sex:\n');
fprintf(fid, '  Male (1):   %d\n', sum(sex==1));
fprintf(fid, '  Female (2): %d\n', sum(sex==2));
fprintf(fid, '  Missing:    %d\n', sum(isnan(sex)));
fprintf(fid, '\n');
fprintf(fid, 'TIV:\n');
fprintf(fid, '  Valid:   %d\n', sum(~isnan(tiv)));
fprintf(fid, '  Missing: %d\n', sum(isnan(tiv)));
fprintf(fid, '  Mean:    %.0f ml\n', nanmean(tiv));
fprintf(fid, '  Range:   [%.0f, %.0f]\n', nanmin(tiv), nanmax(tiv));
fprintf(fid, '\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'OUTPUT FILES\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'Location: %s\n', OUTPUT_DIR);
fprintf(fid, '\n');
fprintf(fid, '1. WAHN_Y_raw.mat         - Raw brain data [%d x %d]\n', size(Y_raw, 1), size(Y_raw, 2));
fprintf(fid, '2. WAHN_Y_corrected.mat   - Mean offset corrected [%d x %d]\n', size(Y_corrected, 1), size(Y_corrected, 2));
fprintf(fid, '3. WAHN_covariates.mat    - Age, Sex, TIV, Site\n');
fprintf(fid, '4. WAHN_subject_table.csv - Subject mapping\n');
fprintf(fid, '5. WAHN_mean_reference.mat- Mean vectors used for correction\n');
fprintf(fid, '6. WAHN_validation_report.txt - This report\n');
fprintf(fid, '\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, 'NEXT STEPS FOR CLARA\n');
fprintf(fid, '------------------------------------------------------------------------\n');
fprintf(fid, '1. Copy Y from WAHN_Y_raw.mat into NM OOCV structure\n');
fprintf(fid, '2. Apply Transition Model to raw data\n');
fprintf(fid, '3. Copy Y from WAHN_Y_corrected.mat into NM OOCV structure\n');
fprintf(fid, '4. Apply Transition Model to corrected data\n');
fprintf(fid, '5. (Optional) Use WAHN_covariates.mat for ComBat harmonization\n');
fprintf(fid, '\n');
fprintf(fid, '========================================================================\n');

fclose(fid);
fprintf('Saved: WAHN_validation_report.txt\n');

%% ========================================================================
%  FINAL SUMMARY
%  ========================================================================
fprintf('\n========================================================================\n');
fprintf('  WAHN DATA PREPARATION COMPLETE\n');
fprintf('========================================================================\n\n');

fprintf('Output location: %s\n\n', OUTPUT_DIR);

fprintf('Files for Clara:\n');
fprintf('  1. WAHN_Y_raw.mat         -> Apply model to raw data\n');
fprintf('  2. WAHN_Y_corrected.mat   -> Apply model to mean offset corrected\n');
fprintf('  3. WAHN_covariates.mat    -> For ComBat harmonization\n');
fprintf('  4. WAHN_subject_table.csv -> Subject documentation\n');

fprintf('\nData dimensions: [%d subjects x %d features]\n', n_subjects, n_features);
fprintf('Feature count matches NESDA: %s\n', iif(n_features == EXPECTED_FEATURES, 'YES', 'NO'));

fprintf('\nNext: Clara copies Y into NM OOCV structure and applies Transition Model.\n');

fprintf('\n========================================================================\n');
fprintf('Finished: %s\n', datestr(now));
fprintf('========================================================================\n');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
