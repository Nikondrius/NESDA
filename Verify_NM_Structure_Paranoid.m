%% ==========================================================================
%  PARANOID VERIFICATION OF NM STRUCTURE
%  ==========================================================================
%  Checks EVERYTHING - especially the 3 newly added demographic variables
%  ==========================================================================

fprintf('\n');
fprintf('##################################################################\n');
fprintf('#          PARANOID NM STRUCTURE VERIFICATION                   #\n');
fprintf('##################################################################\n\n');

errors_found = 0;
warnings_found = 0;

%% ==========================================================================
%  CHECK 1: NM STRUCTURE EXISTS
%% ==========================================================================
fprintf('CHECK 1: NM STRUCTURE EXISTS\n');
fprintf('----------------------------------------------------------\n');

if ~exist('NM', 'var')
    fprintf('  [FATAL] NM structure not found in workspace!\n');
    error('Cannot continue - load NM first!');
else
    fprintf('  [OK] NM structure exists\n\n');
end

%% ==========================================================================
%  CHECK 2: BASIC DIMENSIONS
%% ==========================================================================
fprintf('CHECK 2: BASIC DIMENSIONS\n');
fprintf('----------------------------------------------------------\n');

n_subjects = size(NM.Y{1}, 1);
n_features = size(NM.Y{1}, 2);
n_labels = size(NM.label, 2);

fprintf('  Subjects: %d\n', n_subjects);
fprintf('  Features: %d\n', n_features);
fprintf('  Labels: %d\n', n_labels);

if n_subjects ~= 306
    fprintf('  [ERROR] Expected 306 subjects, got %d!\n', n_subjects);
    errors_found = errors_found + 1;
else
    fprintf('  [OK] Subject count correct (306)\n');
end

if n_features ~= 21
    fprintf('  [ERROR] Expected 21 features, got %d!\n', n_features);
    errors_found = errors_found + 1;
else
    fprintf('  [OK] Feature count correct (21)\n');
end

if n_labels ~= 2
    fprintf('  [ERROR] Expected 2 labels (multi-label), got %d!\n', n_labels);
    errors_found = errors_found + 1;
else
    fprintf('  [OK] Label count correct (2 = Transition + bvFTD)\n');
end
fprintf('\n');

%% ==========================================================================
%  CHECK 3: FEATURE NAMES
%% ==========================================================================
fprintf('CHECK 3: FEATURE NAMES\n');
fprintf('----------------------------------------------------------\n');

if isfield(NM, 'featnames') && ~isempty(NM.featnames) && ~isempty(NM.featnames{1})
    feat_names = NM.featnames{1};
    n_names = length(feat_names);

    fprintf('  Feature names count: %d\n', n_names);

    if n_names ~= n_features
        fprintf('  [ERROR] Feature names (%d) != Feature count (%d)!\n', n_names, n_features);
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] Feature names match feature count\n');
    end

    % Check last 3 are the demographics
    expected_last_3 = {'abmi', 'aedu', 'aLCAsubtype'};
    actual_last_3 = feat_names(end-2:end);

    fprintf('\n  Last 3 feature names (should be demographics):\n');
    for i = 1:3
        idx = n_names - 3 + i;
        fprintf('    [%d] %s', idx, actual_last_3{i});
        if strcmp(actual_last_3{i}, expected_last_3{i})
            fprintf(' [OK]\n');
        else
            fprintf(' [ERROR] Expected: %s\n', expected_last_3{i});
            errors_found = errors_found + 1;
        end
    end
else
    fprintf('  [ERROR] NM.featnames not found or empty!\n');
    errors_found = errors_found + 1;
end
fprintf('\n');

%% ==========================================================================
%  CHECK 4: THE 3 NEW DEMOGRAPHIC VARIABLES (PARANOID MODE)
%% ==========================================================================
fprintf('CHECK 4: NEW DEMOGRAPHIC VARIABLES (PARANOID)\n');
fprintf('----------------------------------------------------------\n');

% Get the last 3 columns (should be abmi, aedu, aLCAsubtype)
new_vars = NM.Y{1}(:, end-2:end);
var_names = {'abmi', 'aedu', 'aLCAsubtype'};

% Expected ranges based on NESDA data
expected_ranges = struct();
expected_ranges.abmi = [15, 60];        % BMI realistic range
expected_ranges.aedu = [0, 25];         % Education years
expected_ranges.aLCAsubtype = [-3, 3];  % LCA subtype codes

% Expected statistics from the script output
expected_stats = struct();
expected_stats.abmi = struct('mean', 25.07, 'sd', 4.68, 'n', 306);
expected_stats.aedu = struct('mean', 12.76, 'sd', 3.20, 'n', 306);
expected_stats.aLCAsubtype = struct('mean', -0.40, 'sd', 1.97, 'n', 306);

for i = 1:3
    var_name = var_names{i};
    var_data = new_vars(:, i);

    fprintf('\n  --- %s (column %d) ---\n', var_name, n_features - 3 + i);

    % Basic stats
    n_valid = sum(~isnan(var_data));
    n_nan = sum(isnan(var_data));
    var_mean = nanmean(var_data);
    var_sd = nanstd(var_data);
    var_min = min(var_data);
    var_max = max(var_data);

    fprintf('  N valid: %d / %d\n', n_valid, n_subjects);
    fprintf('  N NaN: %d\n', n_nan);
    fprintf('  Mean: %.2f (expected: %.2f)\n', var_mean, expected_stats.(var_name).mean);
    fprintf('  SD: %.2f (expected: %.2f)\n', var_sd, expected_stats.(var_name).sd);
    fprintf('  Range: [%.2f, %.2f]\n', var_min, var_max);

    % Check for complete data
    if n_valid ~= expected_stats.(var_name).n
        fprintf('  [ERROR] Expected %d valid, got %d!\n', expected_stats.(var_name).n, n_valid);
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] All %d values present\n', n_valid);
    end

    % Check mean matches (tolerance 0.1)
    if abs(var_mean - expected_stats.(var_name).mean) > 0.1
        fprintf('  [ERROR] Mean mismatch! Expected %.2f, got %.2f\n', ...
            expected_stats.(var_name).mean, var_mean);
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] Mean matches expected value\n');
    end

    % Check SD matches (tolerance 0.1)
    if abs(var_sd - expected_stats.(var_name).sd) > 0.1
        fprintf('  [ERROR] SD mismatch! Expected %.2f, got %.2f\n', ...
            expected_stats.(var_name).sd, var_sd);
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] SD matches expected value\n');
    end

    % Check range is realistic
    exp_range = expected_ranges.(var_name);
    if var_min < exp_range(1) || var_max > exp_range(2)
        fprintf('  [WARNING] Values outside expected range [%d, %d]\n', exp_range(1), exp_range(2));
        warnings_found = warnings_found + 1;
    else
        fprintf('  [OK] Values within realistic range\n');
    end
end
fprintf('\n');

%% ==========================================================================
%  CHECK 5: LABELS (TRANSITION-26 + bvFTD)
%% ==========================================================================
fprintf('CHECK 5: LABELS (MULTI-LABEL)\n');
fprintf('----------------------------------------------------------\n');

label_names = {'Transition-26', 'bvFTD'};
for i = 1:2
    label_data = NM.label(:, i);
    n_valid = sum(~isnan(label_data));
    label_mean = nanmean(label_data);
    label_sd = nanstd(label_data);

    fprintf('  %s: n=%d valid, mean=%.3f, SD=%.3f\n', ...
        label_names{i}, n_valid, label_mean, label_sd);

    if n_valid ~= 306
        fprintf('  [WARNING] Not all subjects have %s scores\n', label_names{i});
        warnings_found = warnings_found + 1;
    end
end
fprintf('\n');

%% ==========================================================================
%  CHECK 6: COVARIATES
%% ==========================================================================
fprintf('CHECK 6: COVARIATES\n');
fprintf('----------------------------------------------------------\n');

if isfield(NM, 'covars') && ~isempty(NM.covars)
    n_covars = size(NM.covars, 2);
    fprintf('  Number of covariates: %d\n', n_covars);

    if n_covars ~= 3
        fprintf('  [ERROR] Expected 3 covariates (Age, Sex, Site), got %d!\n', n_covars);
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] Covariate count correct\n');
    end

    % Check covariate names if available
    if isfield(NM, 'covnames') && ~isempty(NM.covnames)
        fprintf('  Covariate names: %s\n', strjoin(NM.covnames, ', '));
    end

    % Check Site values (should be 1, 2, 3 only)
    site_col = NM.covars(:, 3);  % Assuming Site is 3rd covariate
    unique_sites = unique(site_col(~isnan(site_col)));
    fprintf('  Site values: %s\n', num2str(unique_sites'));

    if any(unique_sites < 1) || any(unique_sites > 3)
        fprintf('  [ERROR] Invalid site codes found!\n');
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] Site codes valid (1, 2, 3)\n');
    end
else
    fprintf('  [ERROR] NM.covars not found!\n');
    errors_found = errors_found + 1;
end
fprintf('\n');

%% ==========================================================================
%  CHECK 7: CASE IDs
%% ==========================================================================
fprintf('CHECK 7: CASE IDs (PIDENTS)\n');
fprintf('----------------------------------------------------------\n');

if isfield(NM, 'cases') && ~isempty(NM.cases)
    n_cases = length(NM.cases);
    fprintf('  Number of case IDs: %d\n', n_cases);

    if n_cases ~= n_subjects
        fprintf('  [ERROR] Case IDs (%d) != Subjects (%d)!\n', n_cases, n_subjects);
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] Case ID count matches subjects\n');
    end

    % Check for duplicates
    if iscell(NM.cases)
        unique_cases = unique(NM.cases);
    else
        unique_cases = unique(NM.cases);
    end

    if length(unique_cases) ~= n_cases
        fprintf('  [ERROR] Duplicate case IDs found!\n');
        errors_found = errors_found + 1;
    else
        fprintf('  [OK] All case IDs unique\n');
    end
else
    fprintf('  [ERROR] NM.cases not found!\n');
    errors_found = errors_found + 1;
end
fprintf('\n');

%% ==========================================================================
%  CHECK 8: CROSS-VERIFY DEMOGRAPHICS WITH ORIGINAL DATA
%% ==========================================================================
fprintf('CHECK 8: CROSS-VERIFY WITH ORIGINAL NESDA DATA\n');
fprintf('----------------------------------------------------------\n');

nesda_path = '/volume/projects/CV_NESDA/Data/tabular_data/NESDA_tabular_combined_data.csv';

if exist(nesda_path, 'file')
    fprintf('  Loading original NESDA data for verification...\n');
    nesda_orig = readtable(nesda_path, 'VariableNamingRule', 'preserve');

    % Get pidents from NM
    if iscell(NM.cases)
        nm_pidents = cellfun(@str2double, NM.cases);
    else
        nm_pidents = NM.cases;
    end

    % Spot check 5 random subjects
    rng(42);  % Reproducible
    check_idx = randperm(n_subjects, 5);

    fprintf('  Spot-checking 5 random subjects:\n\n');

    spot_check_ok = true;
    for i = 1:length(check_idx)
        subj_idx = check_idx(i);
        pident = nm_pidents(subj_idx);

        % Find in original data
        orig_idx = find(nesda_orig.pident == pident);

        if isempty(orig_idx)
            fprintf('  [ERROR] pident %d not found in original data!\n', pident);
            errors_found = errors_found + 1;
            spot_check_ok = false;
            continue;
        end

        fprintf('  Subject %d (pident=%d):\n', subj_idx, pident);

        % Check abmi
        nm_abmi = NM.Y{1}(subj_idx, end-2);
        orig_abmi = nesda_orig.abmi(orig_idx);
        match_abmi = abs(nm_abmi - orig_abmi) < 0.001 || (isnan(nm_abmi) && isnan(orig_abmi));
        fprintf('    abmi: NM=%.2f, Orig=%.2f %s\n', nm_abmi, orig_abmi, ternary(match_abmi, '[OK]', '[MISMATCH]'));
        if ~match_abmi, spot_check_ok = false; errors_found = errors_found + 1; end

        % Check aedu
        nm_aedu = NM.Y{1}(subj_idx, end-1);
        orig_aedu = nesda_orig.aedu(orig_idx);
        match_aedu = abs(nm_aedu - orig_aedu) < 0.001 || (isnan(nm_aedu) && isnan(orig_aedu));
        fprintf('    aedu: NM=%.2f, Orig=%.2f %s\n', nm_aedu, orig_aedu, ternary(match_aedu, '[OK]', '[MISMATCH]'));
        if ~match_aedu, spot_check_ok = false; errors_found = errors_found + 1; end

        % Check aLCAsubtype
        nm_lca = NM.Y{1}(subj_idx, end);
        orig_lca = nesda_orig.aLCAsubtype(orig_idx);
        match_lca = abs(nm_lca - orig_lca) < 0.001 || (isnan(nm_lca) && isnan(orig_lca));
        fprintf('    aLCAsubtype: NM=%.2f, Orig=%.2f %s\n', nm_lca, orig_lca, ternary(match_lca, '[OK]', '[MISMATCH]'));
        if ~match_lca, spot_check_ok = false; errors_found = errors_found + 1; end

        fprintf('\n');
    end

    if spot_check_ok
        fprintf('  [OK] All spot checks passed!\n');
    end
else
    fprintf('  [WARNING] Original NESDA data not accessible - skipping cross-verification\n');
    warnings_found = warnings_found + 1;
end
fprintf('\n');

%% ==========================================================================
%  FINAL SUMMARY
%% ==========================================================================
fprintf('##################################################################\n');
fprintf('#                    VERIFICATION SUMMARY                        #\n');
fprintf('##################################################################\n\n');

fprintf('  Errors found: %d\n', errors_found);
fprintf('  Warnings found: %d\n\n', warnings_found);

if errors_found == 0 && warnings_found == 0
    fprintf('  ============================================\n');
    fprintf('  =    ALL CHECKS PASSED - NM STRUCTURE OK   =\n');
    fprintf('  ============================================\n');
elseif errors_found == 0
    fprintf('  ============================================\n');
    fprintf('  =  PASSED WITH WARNINGS - Review above     =\n');
    fprintf('  ============================================\n');
else
    fprintf('  ============================================\n');
    fprintf('  =  ERRORS FOUND - DO NOT USE THIS FILE!    =\n');
    fprintf('  ============================================\n');
end

fprintf('\n');

%% ==========================================================================
%  HELPER FUNCTION
%% ==========================================================================
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
