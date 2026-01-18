%% ==========================================================================
%  NESDA DATA LOADER FOR NEUROMINER
%  ==========================================================================
%  Author: Claude AI Assistant (for Nikos Diederichs)
%  Date: January 2026
%
%  PURPOSE:
%  Loads the curated NESDA clinical data into MATLAB workspace in the format
%  expected by NeuroMiner for normative modeling with multi-label support.
%
%  USAGE:
%  1. Run this script on HPC/VNC
%  2. Start NeuroMiner: nm = nk_Initialize
%  3. Select "Load data for model discovery and cross-validation"
%  4. Use the pre-loaded workspace variables for data entry
%
%  OUTPUT VARIABLES (loaded into workspace):
%  - NM_features: [306 x 19] numeric matrix of clinical features
%  - NM_labels_Transition: [306 x 1] Transition-26 decision scores
%  - NM_labels_bvFTD: [306 x 1] bvFTD decision scores
%  - NM_subject_ids: {306 x 1} cell array of subject IDs (pident)
%  - NM_diagnosis_group: {306 x 1} cell array of diagnosis groups
%  - NM_feature_names: {1 x 19} cell array of feature names
%  - NM_covariates: struct with Age, Sex, Site for covariate correction
%  ==========================================================================

clear; clc;
fprintf('==========================================================\n');
fprintf('  NESDA DATA LOADER FOR NEUROMINER\n');
fprintf('==========================================================\n\n');

%% ==========================================================================
%  CONFIGURATION - UPDATE THESE PATHS IF NEEDED
%  ==========================================================================

% Most recent output folder from NESDA_Clinical_Associations.m
data_path = '/volume/projects/CV_NESDA/Analysis/Clinical_Associations_2026-01-12_17-00/Data/';

% Original NESDA data for covariates
nesda_data_path = '/volume/projects/CV_NESDA/Data/tabular_data/NESDA_tabular_combined_data.csv';

% Check if paths exist
if ~exist(data_path, 'dir')
    % Try to find the most recent Clinical_Associations folder
    base_path = '/volume/projects/CV_NESDA/Analysis/';
    folders = dir([base_path 'Clinical_Associations_*']);
    if ~isempty(folders)
        [~, idx] = max([folders.datenum]);
        data_path = [base_path folders(idx).name '/Data/'];
        fprintf('  Using most recent output folder:\n  %s\n\n', data_path);
    else
        error('Cannot find Clinical_Associations output folder!');
    end
end

%% ==========================================================================
%  SECTION 1: LOAD NM-READY FILES
%  ==========================================================================
fprintf('SECTION 1: LOADING NM-READY DATA FILES\n');
fprintf('----------------------------------------------------------\n\n');

% Load the all-subjects file (306 subjects: 69 HCs + 237 patients)
all_subjects_file = [data_path 'NM_Ready_ALL_Subjects_With_Labels.csv'];
fprintf('  Loading: %s\n', all_subjects_file);

if exist(all_subjects_file, 'file')
    data_all = readtable(all_subjects_file, 'VariableNamingRule', 'preserve');
    fprintf('  ✓ Loaded: [%d subjects × %d variables]\n\n', height(data_all), width(data_all));
else
    error('File not found: %s', all_subjects_file);
end

%% ==========================================================================
%  SECTION 2: EXTRACT COMPONENTS FOR NEUROMINER
%  ==========================================================================
fprintf('SECTION 2: EXTRACTING NEUROMINER COMPONENTS\n');
fprintf('----------------------------------------------------------\n\n');

% --- Subject IDs ---
NM_subject_ids = cellstr(num2str(data_all.pident));
fprintf('  Subject IDs: %d subjects extracted\n', length(NM_subject_ids));

% --- Diagnosis Group ---
NM_diagnosis_group = data_all.diagnosis_group;
n_hc = sum(strcmp(NM_diagnosis_group, 'HC'));
n_patients = height(data_all) - n_hc;
fprintf('  Diagnosis groups: %d HCs, %d Patients\n', n_hc, n_patients);

% --- Labels (Decision Scores) ---
NM_labels_Transition = data_all.Transition_26;
NM_labels_bvFTD = data_all.bvFTD;
fprintf('  Labels extracted:\n');
fprintf('    - Transition_26: n=%d valid, mean=%.3f, SD=%.3f\n', ...
    sum(~isnan(NM_labels_Transition)), nanmean(NM_labels_Transition), nanstd(NM_labels_Transition));
fprintf('    - bvFTD: n=%d valid, mean=%.3f, SD=%.3f\n', ...
    sum(~isnan(NM_labels_bvFTD)), nanmean(NM_labels_bvFTD), nanstd(NM_labels_bvFTD));

% --- Features (Clinical Variables) ---
% Exclude pident, diagnosis_group, and labels
exclude_vars = {'pident', 'diagnosis_group', 'Transition_26', 'bvFTD'};
all_vars = data_all.Properties.VariableNames;
feature_vars = setdiff(all_vars, exclude_vars, 'stable');

% Extract feature matrix
NM_features = zeros(height(data_all), length(feature_vars));
NM_feature_names = feature_vars;

for i = 1:length(feature_vars)
    var_data = data_all.(feature_vars{i});
    if isnumeric(var_data)
        NM_features(:, i) = var_data;
    else
        % Convert to numeric if needed
        NM_features(:, i) = double(var_data);
    end
end

fprintf('  Features extracted: [%d subjects × %d features]\n', size(NM_features, 1), size(NM_features, 2));
fprintf('    Feature names: %s\n', strjoin(NM_feature_names(1:5), ', '));
fprintf('    ... and %d more\n', length(NM_feature_names) - 5);

% Report NaN statistics
nan_per_feature = sum(isnan(NM_features), 1);
fprintf('  NaN per feature: min=%d, max=%d, mean=%.1f\n', ...
    min(nan_per_feature), max(nan_per_feature), mean(nan_per_feature));

%% ==========================================================================
%  SECTION 3: LOAD COVARIATES (Age, Sex, Site)
%  ==========================================================================
fprintf('\nSECTION 3: LOADING COVARIATES\n');
fprintf('----------------------------------------------------------\n\n');

if exist(nesda_data_path, 'file')
    fprintf('  Loading original NESDA data for covariates...\n');
    nesda_full = readtable(nesda_data_path, 'VariableNamingRule', 'preserve');

    % Match subjects
    [~, idx_nesda, idx_nm] = intersect(nesda_full.pident, data_all.pident);

    % Initialize covariates
    NM_covariates = struct();
    NM_covariates.Age = NaN(height(data_all), 1);
    NM_covariates.Sex = NaN(height(data_all), 1);
    NM_covariates.Site = NaN(height(data_all), 1);

    % Extract Age
    if ismember('Age', nesda_full.Properties.VariableNames)
        NM_covariates.Age(idx_nm) = nesda_full.Age(idx_nesda);
        fprintf('  ✓ Age: n=%d valid, mean=%.1f, SD=%.1f\n', ...
            sum(~isnan(NM_covariates.Age)), nanmean(NM_covariates.Age), nanstd(NM_covariates.Age));
    end

    % Extract Sex (coded as 1=Male, 2=Female in NESDA)
    if ismember('Sexe', nesda_full.Properties.VariableNames)
        NM_covariates.Sex(idx_nm) = nesda_full.Sexe(idx_nesda);
        n_male = sum(NM_covariates.Sex == 1);
        n_female = sum(NM_covariates.Sex == 2);
        fprintf('  ✓ Sex: %d Male (1), %d Female (2)\n', n_male, n_female);
    end

    % Extract Site (ascanloc = scan location, categorical variable)
    % NeuroMiner expects site as numeric codes for multi-site correction
    % Note: -2 in NESDA = missing scan location - reconstruct from alternatives
    if ismember('ascanloc', nesda_full.Properties.VariableNames)
        site_raw = nesda_full.ascanloc(idx_nesda);
        pidents = nesda_full.pident(idx_nesda);

        % Check if aarea (interview location) is available as backup
        has_aarea = ismember('aarea', nesda_full.Properties.VariableNames);
        if has_aarea
            aarea_raw = nesda_full.aarea(idx_nesda);
        end

        % Count and fix missing sites (-2 = missing in NESDA)
        n_missing_original = sum(site_raw == -2 | isnan(site_raw));
        n_fixed_aarea = 0;
        n_fixed_pident = 0;

        if n_missing_original > 0
            fprintf('  Fixing %d subjects with missing site (ascanloc=-2):\n', n_missing_original);

            for i = 1:length(site_raw)
                if isnan(site_raw(i)) || site_raw(i) < 0
                    % Option 1: Try aarea (interview location) as backup
                    if has_aarea && ~isnan(aarea_raw(i)) && aarea_raw(i) > 0 && aarea_raw(i) <= 3
                        site_raw(i) = aarea_raw(i);
                        n_fixed_aarea = n_fixed_aarea + 1;
                    else
                        % Option 2: Use first digit of pident (1, 2, 3 = valid NESDA sites)
                        pident_str = num2str(pidents(i));
                        first_digit = str2double(pident_str(1));
                        if first_digit >= 1 && first_digit <= 3
                            site_raw(i) = first_digit;
                            n_fixed_pident = n_fixed_pident + 1;
                        end
                        % If still invalid, remains NaN for exclusion
                    end
                end
            end

            fprintf('    - Fixed via aarea: %d\n', n_fixed_aarea);
            fprintf('    - Fixed via pident first digit: %d\n', n_fixed_pident);
            n_still_missing = sum(isnan(site_raw) | site_raw < 0);
            if n_still_missing > 0
                fprintf('    - Still missing (will be NaN): %d\n', n_still_missing);
                site_raw(site_raw < 0) = NaN;
            end
        end

        NM_covariates.Site(idx_nm) = site_raw;
        unique_sites = unique(NM_covariates.Site(~isnan(NM_covariates.Site)));
        fprintf('  ✓ Site (ascanloc): n=%d valid, %d unique sites\n', ...
            sum(~isnan(NM_covariates.Site)), length(unique_sites));
        fprintf('    Site codes: %s\n', num2str(unique_sites'));
    elseif ismember('site', nesda_full.Properties.VariableNames)
        site_raw = nesda_full.site(idx_nesda);

        % Handle missing site codes
        n_missing_site = sum(site_raw < 0);
        if n_missing_site > 0
            fprintf('  Note: %d subjects with negative site codes recoded to NaN\n', n_missing_site);
            site_raw(site_raw < 0) = NaN;
        end

        NM_covariates.Site(idx_nm) = site_raw;
        unique_sites = unique(NM_covariates.Site(~isnan(NM_covariates.Site)));
        fprintf('  ✓ Site: n=%d valid, %d unique sites\n', ...
            sum(~isnan(NM_covariates.Site)), length(unique_sites));
    else
        fprintf('  WARNING: Site variable (ascanloc/site) not found\n');
    end
else
    fprintf('  WARNING: Original NESDA data not found - covariates not loaded\n');
    NM_covariates = struct('Age', [], 'Sex', [], 'Site', []);
end

%% ==========================================================================
%  SECTION 4: CREATE BINARY GROUP LABELS FOR CLASSIFICATION
%  ==========================================================================
fprintf('\nSECTION 4: CREATING GROUP LABELS\n');
fprintf('----------------------------------------------------------\n\n');

% Binary HC vs Patient label (useful for some NM analyses)
NM_label_HC_vs_Patient = double(~strcmp(NM_diagnosis_group, 'HC'));
fprintf('  HC vs Patient: 0=HC (n=%d), 1=Patient (n=%d)\n', ...
    sum(NM_label_HC_vs_Patient == 0), sum(NM_label_HC_vs_Patient == 1));

% Create numeric diagnosis codes
NM_diagnosis_numeric = zeros(height(data_all), 1);
NM_diagnosis_numeric(strcmp(NM_diagnosis_group, 'HC')) = 0;
NM_diagnosis_numeric(strcmp(NM_diagnosis_group, 'Anxiety')) = 1;
NM_diagnosis_numeric(strcmp(NM_diagnosis_group, 'Depression')) = 2;
NM_diagnosis_numeric(strcmp(NM_diagnosis_group, 'Comorbid')) = 3;
fprintf('  Diagnosis codes: 0=HC, 1=Anxiety, 2=Depression, 3=Comorbid\n');

%% ==========================================================================
%  SECTION 5: SAVE WORKSPACE FOR QUICK RELOAD
%  ==========================================================================
fprintf('\nSECTION 5: SAVING WORKSPACE\n');
fprintf('----------------------------------------------------------\n\n');

workspace_file = [data_path 'NM_Workspace_Ready.mat'];
save(workspace_file, ...
    'NM_features', 'NM_feature_names', ...
    'NM_labels_Transition', 'NM_labels_bvFTD', ...
    'NM_subject_ids', 'NM_diagnosis_group', 'NM_diagnosis_numeric', ...
    'NM_label_HC_vs_Patient', 'NM_covariates', ...
    'data_path');
fprintf('  ✓ Saved: %s\n', workspace_file);
fprintf('    Quick reload: load(''%s'')\n', workspace_file);

%% ==========================================================================
%  SECTION 6: DISPLAY SUMMARY
%  ==========================================================================
fprintf('\n==========================================================\n');
fprintf('  WORKSPACE READY FOR NEUROMINER\n');
fprintf('==========================================================\n\n');

fprintf('VARIABLES LOADED:\n');
fprintf('----------------------------------------------------------\n');
fprintf('  NM_features           [%d × %d]  Clinical feature matrix\n', size(NM_features));
fprintf('  NM_feature_names      {1 × %d}   Feature names\n', length(NM_feature_names));
fprintf('  NM_labels_Transition  [%d × 1]   Transition-26 scores (regression)\n', length(NM_labels_Transition));
fprintf('  NM_labels_bvFTD       [%d × 1]   bvFTD scores (regression)\n', length(NM_labels_bvFTD));
fprintf('  NM_subject_ids        {%d × 1}   Subject IDs (pident)\n', length(NM_subject_ids));
fprintf('  NM_diagnosis_group    {%d × 1}   Diagnosis (HC/Anxiety/Depression/Comorbid)\n', length(NM_diagnosis_group));
fprintf('  NM_diagnosis_numeric  [%d × 1]   Numeric diagnosis (0-3)\n', length(NM_diagnosis_numeric));
fprintf('  NM_label_HC_vs_Patient[%d × 1]   Binary: 0=HC, 1=Patient\n', length(NM_label_HC_vs_Patient));
fprintf('  NM_covariates         struct     Age, Sex, Site for covariate correction\n\n');

fprintf('SAMPLE BREAKDOWN:\n');
fprintf('----------------------------------------------------------\n');
fprintf('  HC:         %3d (%.1f%%)\n', n_hc, 100*n_hc/height(data_all));
fprintf('  Anxiety:    %3d (%.1f%%)\n', sum(strcmp(NM_diagnosis_group, 'Anxiety')), ...
    100*sum(strcmp(NM_diagnosis_group, 'Anxiety'))/height(data_all));
fprintf('  Depression: %3d (%.1f%%)\n', sum(strcmp(NM_diagnosis_group, 'Depression')), ...
    100*sum(strcmp(NM_diagnosis_group, 'Depression'))/height(data_all));
fprintf('  Comorbid:   %3d (%.1f%%)\n', sum(strcmp(NM_diagnosis_group, 'Comorbid')), ...
    100*sum(strcmp(NM_diagnosis_group, 'Comorbid'))/height(data_all));
fprintf('  TOTAL:      %3d\n\n', height(data_all));

fprintf('NEUROMINER DATA ENTRY INSTRUCTIONS:\n');
fprintf('==========================================================\n\n');
fprintf('1. Start NeuroMiner:\n');
fprintf('   >> nm = nk_Initialize\n\n');
fprintf('2. Select: "Load data for model discovery and cross-validation"\n\n');
fprintf('3. For REGRESSION with Transition-26:\n');
fprintf('   - Input data: NM_features\n');
fprintf('   - Labels: NM_labels_Transition\n');
fprintf('   - ID variable: NM_subject_ids\n\n');
fprintf('4. For REGRESSION with bvFTD:\n');
fprintf('   - Input data: NM_features\n');
fprintf('   - Labels: NM_labels_bvFTD\n');
fprintf('   - ID variable: NM_subject_ids\n\n');
fprintf('5. For CLASSIFICATION (HC vs Patient):\n');
fprintf('   - Input data: NM_features\n');
fprintf('   - Labels: NM_label_HC_vs_Patient\n');
fprintf('   - Group 0 = HC, Group 1 = Patient\n\n');
fprintf('6. Add covariates (recommended for multi-site correction):\n');
fprintf('   - Age: NM_covariates.Age\n');
fprintf('   - Sex: NM_covariates.Sex\n');
fprintf('   - Site: NM_covariates.Site (categorical, for site harmonization)\n\n');
fprintf('==========================================================\n');
fprintf('  DATA LOADING COMPLETE - Ready for NeuroMiner!\n');
fprintf('==========================================================\n');
