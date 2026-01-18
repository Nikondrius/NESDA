%% ==========================================================================
%  ADD DEMOGRAPHIC VARIABLES TO EXISTING NM STRUCTURE
%  ==========================================================================
%  Author: Claude AI Assistant
%  Date: January 2026
%
%  PURPOSE:
%  Adds abmi, aedu, and aLCAsubtype to the existing NeuroMiner structure
%  by matching on pident (subject ID).
%
%  USAGE:
%  1. Load the existing NM structure in NeuroMiner (or have it in workspace)
%  2. Run this script
%  3. Save the updated NM structure
%  ==========================================================================

fprintf('==========================================================\n');
fprintf('  ADDING DEMOGRAPHIC VARIABLES TO NM STRUCTURE\n');
fprintf('==========================================================\n\n');

%% ==========================================================================
%  CONFIGURATION
%  ==========================================================================

% Path to original NESDA data with demographics
nesda_data_path = '/volume/projects/CV_NESDA/Data/tabular_data/NESDA_tabular_combined_data.csv';

% Variables to add
vars_to_add = {'abmi', 'aedu', 'aLCAsubtype'};

%% ==========================================================================
%  STEP 1: CHECK NM STRUCTURE EXISTS
%  ==========================================================================
fprintf('STEP 1: CHECKING NM STRUCTURE\n');
fprintf('----------------------------------------------------------\n\n');

if ~exist('NM', 'var')
    error(['NM structure not found in workspace!\n' ...
           'Please load it first:\n' ...
           '  1. Start NeuroMiner: nm\n' ...
           '  2. Select "Load NeuroMiner structure"\n' ...
           '  3. Navigate to your .mat file\n' ...
           '  4. Then run this script again']);
end

% Get current dimensions
n_subjects = size(NM.Y{1}, 1);
n_features_old = size(NM.Y{1}, 2);
fprintf('  Current NM structure: %d subjects × %d features\n\n', n_subjects, n_features_old);

%% ==========================================================================
%  STEP 2: GET SUBJECT IDs FROM NM
%  ==========================================================================
fprintf('STEP 2: EXTRACTING SUBJECT IDs FROM NM\n');
fprintf('----------------------------------------------------------\n\n');

% NM stores case IDs in NM.cases
if isfield(NM, 'cases') && ~isempty(NM.cases)
    nm_pidents = NM.cases;
    % Convert to numeric if stored as strings
    if iscell(nm_pidents)
        nm_pidents_num = cellfun(@(x) str2double(x), nm_pidents);
    else
        nm_pidents_num = nm_pidents;
    end
    fprintf('  Found %d subject IDs in NM.cases\n\n', length(nm_pidents_num));
else
    error('NM.cases not found - cannot match subjects!');
end

%% ==========================================================================
%  STEP 3: LOAD ORIGINAL NESDA DATA
%  ==========================================================================
fprintf('STEP 3: LOADING ORIGINAL NESDA DATA\n');
fprintf('----------------------------------------------------------\n\n');

if ~exist(nesda_data_path, 'file')
    error('NESDA data file not found: %s', nesda_data_path);
end

nesda_full = readtable(nesda_data_path, 'VariableNamingRule', 'preserve');
fprintf('  Loaded NESDA data: %d subjects × %d variables\n', height(nesda_full), width(nesda_full));

% Check which variables are available
available_vars = {};
for i = 1:length(vars_to_add)
    if ismember(vars_to_add{i}, nesda_full.Properties.VariableNames)
        available_vars{end+1} = vars_to_add{i};
        fprintf('  ✓ %s found\n', vars_to_add{i});
    else
        fprintf('  ✗ %s NOT found\n', vars_to_add{i});
    end
end
fprintf('\n');

if isempty(available_vars)
    error('None of the requested variables found in NESDA data!');
end

%% ==========================================================================
%  STEP 4: MATCH SUBJECTS AND EXTRACT NEW FEATURES
%  ==========================================================================
fprintf('STEP 4: MATCHING SUBJECTS BY PIDENT\n');
fprintf('----------------------------------------------------------\n\n');

% Initialize new feature matrix
n_new_features = length(available_vars);
new_features = NaN(n_subjects, n_new_features);

% Match each NM subject to NESDA data
n_matched = 0;
for i = 1:n_subjects
    pident = nm_pidents_num(i);
    idx_nesda = find(nesda_full.pident == pident);

    if ~isempty(idx_nesda)
        n_matched = n_matched + 1;
        for j = 1:n_new_features
            var_name = available_vars{j};
            value = nesda_full.(var_name)(idx_nesda(1));
            if isnumeric(value)
                new_features(i, j) = value;
            end
        end
    end
end

fprintf('  Matched: %d / %d subjects (%.1f%%)\n', n_matched, n_subjects, 100*n_matched/n_subjects);

% Report statistics for new features
fprintf('\n  New feature statistics:\n');
for j = 1:n_new_features
    var_name = available_vars{j};
    valid_n = sum(~isnan(new_features(:, j)));
    if valid_n > 0
        fprintf('    %s: n=%d valid, mean=%.2f, SD=%.2f\n', ...
            var_name, valid_n, nanmean(new_features(:, j)), nanstd(new_features(:, j)));
    else
        fprintf('    %s: NO valid values!\n', var_name);
    end
end
fprintf('\n');

%% ==========================================================================
%  STEP 5: UPDATE NM STRUCTURE
%  ==========================================================================
fprintf('STEP 5: UPDATING NM STRUCTURE\n');
fprintf('----------------------------------------------------------\n\n');

% Append new features to existing feature matrix
NM.Y{1} = [NM.Y{1}, new_features];
n_features_new = size(NM.Y{1}, 2);
fprintf('  Updated feature matrix: %d × %d → %d × %d\n', ...
    n_subjects, n_features_old, n_subjects, n_features_new);

% Update feature names
if isfield(NM, 'featnames') && ~isempty(NM.featnames) && ~isempty(NM.featnames{1})
    old_names = NM.featnames{1};
    if iscell(old_names)
        NM.featnames{1} = [old_names, available_vars];
    else
        % If stored differently, try to append
        NM.featnames{1} = [cellstr(old_names)', available_vars];
    end
    fprintf('  Updated feature names: added %s\n', strjoin(available_vars, ', '));
else
    fprintf('  WARNING: Could not update feature names (NM.featnames not found)\n');
    fprintf('           New features added at positions %d-%d\n', n_features_old+1, n_features_new);
end

%% ==========================================================================
%  STEP 6: VERIFICATION
%  ==========================================================================
fprintf('\nSTEP 6: VERIFICATION\n');
fprintf('----------------------------------------------------------\n\n');

fprintf('  Final NM structure:\n');
fprintf('    Subjects: %d\n', size(NM.Y{1}, 1));
fprintf('    Features: %d (was %d, added %d)\n', n_features_new, n_features_old, n_new_features);
fprintf('    Labels: %d targets\n', size(NM.label, 2));

if isfield(NM, 'covars') && ~isempty(NM.covars)
    fprintf('    Covariates: %d\n', size(NM.covars, 2));
end

%% ==========================================================================
%  DONE - REMIND USER TO SAVE
%  ==========================================================================
fprintf('\n==========================================================\n');
fprintf('  UPDATE COMPLETE!\n');
fprintf('==========================================================\n\n');

fprintf('  Added features: %s\n\n', strjoin(available_vars, ', '));

fprintf('  IMPORTANT: Save the updated NM structure!\n');
fprintf('  In NeuroMiner main menu:\n');
fprintf('    → Select "Save NeuroMiner structure"\n');
fprintf('    → Save as: NM_NESDA_Clinical_MultiLabel_Transition_bvFTD_with_Demographics.mat\n\n');

fprintf('==========================================================\n');
