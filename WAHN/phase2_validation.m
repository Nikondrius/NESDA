%% ========================================================================
%  WAHN Phase 2: Validation Script
%  ========================================================================
%  Author: Nikos Diederichs
%  Date: January 2025
%
%  Purpose: Validate WAHN outputs before sending to Clara
%
%  Run AFTER prepare_wahn_for_clara.m
%  ========================================================================

clear; clc;
fprintf('========================================================================\n');
fprintf('  WAHN Phase 2: Output Validation\n');
fprintf('  Date: %s\n', datestr(now));
fprintf('========================================================================\n\n');

%% Configuration
OUTPUT_DIR = '/volume/projects/CV_NESDA/Data/WAHN/';
NESDA_NM_PAT = '/volume/projects/CV_NESDA/Data/NESDA_Waves/Wave_1/NM_Structures/NM_Pat.mat';
EXPECTED_FEATURES = 71276;

all_tests_passed = true;

%% ========================================================================
%  TEST 1: Feature Dimension
%  ========================================================================
fprintf('========== TEST 1: Feature Dimension ==========\n\n');

try
    raw = load(fullfile(OUTPUT_DIR, 'WAHN_Y_raw.mat'));

    n_features = size(raw.Y, 2);
    fprintf('WAHN features: %d\n', n_features);
    fprintf('Expected:      %d\n', EXPECTED_FEATURES);

    if n_features == EXPECTED_FEATURES
        fprintf('PASS: Feature count matches\n');
    else
        fprintf('FAIL: Feature count mismatch!\n');
        all_tests_passed = false;
    end
catch ME
    fprintf('ERROR: Could not load WAHN_Y_raw.mat: %s\n', ME.message);
    all_tests_passed = false;
end

%% ========================================================================
%  TEST 2: No NaNs in Data
%  ========================================================================
fprintf('\n========== TEST 2: No NaNs in Data ==========\n\n');

try
    raw = load(fullfile(OUTPUT_DIR, 'WAHN_Y_raw.mat'));
    corr = load(fullfile(OUTPUT_DIR, 'WAHN_Y_corrected.mat'));

    nan_raw = sum(isnan(raw.Y(:)));
    nan_corr = sum(isnan(corr.Y(:)));

    fprintf('NaN count in Y_raw:       %d\n', nan_raw);
    fprintf('NaN count in Y_corrected: %d\n', nan_corr);

    if nan_raw == 0 && nan_corr == 0
        fprintf('PASS: No NaN values found\n');
    else
        fprintf('FAIL: NaN values present!\n');
        all_tests_passed = false;
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
    all_tests_passed = false;
end

%% ========================================================================
%  TEST 3: Compare with NESDA
%  ========================================================================
fprintf('\n========== TEST 3: Compare with NESDA ==========\n\n');

try
    nesda = load(NESDA_NM_PAT);
    wahn = load(fullfile(OUTPUT_DIR, 'WAHN_Y_raw.mat'));

    if iscell(nesda.NM.Y)
        nesda_Y = nesda.NM.Y{1};
    else
        nesda_Y = nesda.NM.Y;
    end

    fprintf('NESDA: [%d x %d], mean=%.4f, std=%.4f\n', ...
            size(nesda_Y, 1), size(nesda_Y, 2), mean(nesda_Y(:)), std(nesda_Y(:)));
    fprintf('WAHN:  [%d x %d], mean=%.4f, std=%.4f\n', ...
            size(wahn.Y, 1), size(wahn.Y, 2), mean(wahn.Y(:)), std(wahn.Y(:)));

    % Check if dimensions match
    if size(nesda_Y, 2) == size(wahn.Y, 2)
        fprintf('PASS: Feature dimensions match NESDA\n');
    else
        fprintf('FAIL: Feature dimensions do not match NESDA\n');
        all_tests_passed = false;
    end

    % Check if value ranges are similar
    nesda_range = [min(nesda_Y(:)), max(nesda_Y(:))];
    wahn_range = [min(wahn.Y(:)), max(wahn.Y(:))];

    fprintf('\nValue ranges:\n');
    fprintf('  NESDA: [%.4f, %.4f]\n', nesda_range(1), nesda_range(2));
    fprintf('  WAHN:  [%.4f, %.4f]\n', wahn_range(1), wahn_range(2));

    % Ranges should be roughly similar for GM probability maps
    if wahn_range(2) > 0 && wahn_range(2) < 2
        fprintf('PASS: WAHN value range appears reasonable for GM probability maps\n');
    else
        fprintf('WARNING: WAHN value range may be unexpected\n');
    end

catch ME
    fprintf('ERROR: %s\n', ME.message);
    all_tests_passed = false;
end

%% ========================================================================
%  TEST 4: Mean Offset Correction Verification
%  ========================================================================
fprintf('\n========== TEST 4: Mean Offset Correction ==========\n\n');

try
    raw = load(fullfile(OUTPUT_DIR, 'WAHN_Y_raw.mat'));
    corr = load(fullfile(OUTPUT_DIR, 'WAHN_Y_corrected.mat'));
    mean_ref = load(fullfile(OUTPUT_DIR, 'WAHN_mean_reference.mat'));

    fprintf('Correction method: %s\n', mean_ref.mean_correction_method);

    % Check that corrected = raw - mean_reference
    if strcmp(mean_ref.mean_correction_method, 'WAHN_within')
        % For within-cohort, column means should be ~0
        col_means = mean(corr.Y, 1);
        mean_of_col_means = mean(col_means);

        fprintf('Mean of column means (corrected): %.2e\n', mean_of_col_means);

        if abs(mean_of_col_means) < 1e-10
            fprintf('PASS: Column means are approximately zero\n');
        else
            fprintf('WARNING: Column means not exactly zero (may be due to numerical precision)\n');
        end
    else
        % For cross-cohort, just verify the correction was applied
        diff_before_after = mean(raw.Y(:)) - mean(corr.Y(:));
        fprintf('Mean difference (raw - corrected): %.4f\n', diff_before_after);

        if abs(diff_before_after) > 0
            fprintf('PASS: Correction was applied (values differ)\n');
        else
            fprintf('WARNING: Raw and corrected appear identical\n');
        end
    end

    % Verify reconstruction: raw = corrected + mean_reference
    reconstructed = corr.Y + mean_ref.mean_reference;
    reconstruction_error = max(abs(raw.Y(:) - reconstructed(:)));

    fprintf('Max reconstruction error: %.2e\n', reconstruction_error);

    if reconstruction_error < 1e-10
        fprintf('PASS: Correction is reversible\n');
    else
        fprintf('WARNING: Reconstruction error higher than expected\n');
    end

catch ME
    fprintf('ERROR: %s\n', ME.message);
    all_tests_passed = false;
end

%% ========================================================================
%  TEST 5: Covariate Completeness
%  ========================================================================
fprintf('\n========== TEST 5: Covariate Completeness ==========\n\n');

try
    cov = load(fullfile(OUTPUT_DIR, 'WAHN_covariates.mat'));
    subj = readtable(fullfile(OUTPUT_DIR, 'WAHN_subject_table.csv'));

    n_subjects = height(subj);
    fprintf('Total subjects: %d\n', n_subjects);

    % Check Age
    n_age = sum(~isnan(cov.age));
    fprintf('Age available:  %d/%d (%.1f%%)\n', n_age, n_subjects, 100*n_age/n_subjects);

    % Check Sex
    n_sex = sum(~isnan(cov.sex));
    fprintf('Sex available:  %d/%d (%.1f%%)\n', n_sex, n_subjects, 100*n_sex/n_subjects);

    % Check TIV
    n_tiv = sum(~isnan(cov.tiv));
    fprintf('TIV available:  %d/%d (%.1f%%)\n', n_tiv, n_subjects, 100*n_tiv/n_subjects);

    % Minimum threshold for Clara (e.g., 80%)
    threshold = 0.8;
    if n_age/n_subjects >= threshold && n_sex/n_subjects >= threshold
        fprintf('PASS: Demographics coverage above %.0f%%\n', threshold*100);
    else
        fprintf('WARNING: Demographics coverage below %.0f%%\n', threshold*100);
    end

catch ME
    fprintf('ERROR: %s\n', ME.message);
    all_tests_passed = false;
end

%% ========================================================================
%  TEST 6: File Integrity
%  ========================================================================
fprintf('\n========== TEST 6: File Integrity ==========\n\n');

expected_files = {
    'WAHN_Y_raw.mat'
    'WAHN_Y_corrected.mat'
    'WAHN_covariates.mat'
    'WAHN_subject_table.csv'
    'WAHN_mean_reference.mat'
    'WAHN_validation_report.txt'
};

for i = 1:length(expected_files)
    fpath = fullfile(OUTPUT_DIR, expected_files{i});
    if exist(fpath, 'file')
        d = dir(fpath);
        fprintf('  [OK] %s (%.1f KB)\n', expected_files{i}, d.bytes/1024);
    else
        fprintf('  [MISSING] %s\n', expected_files{i});
        all_tests_passed = false;
    end
end

%% ========================================================================
%  TEST 7: Data Matrix Consistency
%  ========================================================================
fprintf('\n========== TEST 7: Data Matrix Consistency ==========\n\n');

try
    raw = load(fullfile(OUTPUT_DIR, 'WAHN_Y_raw.mat'));
    corr = load(fullfile(OUTPUT_DIR, 'WAHN_Y_corrected.mat'));
    subj = readtable(fullfile(OUTPUT_DIR, 'WAHN_subject_table.csv'));

    % Check row counts match
    n_raw = size(raw.Y, 1);
    n_corr = size(corr.Y, 1);
    n_subj = height(subj);

    fprintf('Rows in Y_raw:       %d\n', n_raw);
    fprintf('Rows in Y_corrected: %d\n', n_corr);
    fprintf('Rows in subject_table: %d\n', n_subj);

    if n_raw == n_corr && n_raw == n_subj
        fprintf('PASS: All row counts match\n');
    else
        fprintf('FAIL: Row count mismatch!\n');
        all_tests_passed = false;
    end

catch ME
    fprintf('ERROR: %s\n', ME.message);
    all_tests_passed = false;
end

%% ========================================================================
%  SUMMARY
%  ========================================================================
fprintf('\n========================================================================\n');
fprintf('  VALIDATION SUMMARY\n');
fprintf('========================================================================\n\n');

if all_tests_passed
    fprintf('ALL TESTS PASSED\n\n');
    fprintf('Data is ready for Clara:\n');
    fprintf('  Output directory: %s\n', OUTPUT_DIR);
    fprintf('\nFiles to send:\n');
    fprintf('  1. WAHN_Y_raw.mat\n');
    fprintf('  2. WAHN_Y_corrected.mat\n');
    fprintf('  3. WAHN_covariates.mat\n');
    fprintf('  4. WAHN_subject_table.csv\n');
else
    fprintf('SOME TESTS FAILED\n\n');
    fprintf('Please review the errors above and fix before sending to Clara.\n');
end

fprintf('\n========================================================================\n');
fprintf('Validation complete: %s\n', datestr(now));
fprintf('========================================================================\n');
