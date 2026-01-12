%% ==========================================================================
%  PRIORITY 4.1-4.5: COMPREHENSIVE CLINICAL ASSOCIATIONS ANALYSIS
%  ==========================================================================
%  Author: Nikos Diederichs
%  Date: October 26, 2025
%  For: Clara V - NESDA Clinical Associations Phase 2
%  Version: 4.0 - STREAMLINED TO TRANSITION-26 AND bvFTD ONLY
%  MODIFIED: October 29, 2025 - Fixed _fr variable interpretation
%                               + Proper handling of frequency (0/1/2) vs binary (0/1)
%  MODIFIED: November 9, 2025 - Removed OOCV-27 (redundant, focus on primary models)
%  MODIFIED: November 10, 2025 - Added Section 10E-G: Binary diagnosis coding,
%                                 partial correlations, Spearman correlations
%
%  DECISION SCORE VERSIONS USED:
%  - Transition: OOCV-26 (Dynamic Std) [PRIMARY TRANSITION MODEL]
%  - bvFTD: OOCV-7 (Dynamic Std) [DEMENTIA MODEL]
%
%  ANALYSES INCLUDED:
%  - 4.1: Metabolic Subtypes & BMI
%  - 4.2: Symptom Severity (11 variables) + ENHANCED PCA (PC1, PC2, PC3)
%  - 4.3: Clinical History (21 variables including illness duration)
%  - 4.4: Cognition & Functioning
%  - 4.4B: Medication Analysis (PATIENTS ONLY) - CORRECTED CODING
%  - NEW 9C: Recency Stratified Analysis (OPTION 6)
%  - 4.5: Comprehensive Statistical Summary (ALL 40+ VARIABLES)
%  - NEW: Forest Plots for Transition-26 AND bvFTD
%  - NEW: Complete PC1, PC2, PC3 correlations with both decision scores
%  - NEW 10E: Binary Anxiety/Depression Coding (Comorbid as both)
%  - NEW 10F: Partial Correlations (controlling Age, Sex, Site)
%  - NEW 10G: Spearman vs Pearson Correlation Comparison
%  ==========================================================================

clear; clc; close all;

% =========================================================================
% FIGURE HANDLING: Use print() for reliable saving with invisible figures
% =========================================================================

% =========================================================================
% REPRODUCIBILITY: Set random seed for deterministic results
% =========================================================================
% Setting a random seed ensures reproducibility across all analyses,
% even though this script primarily uses deterministic operations.
% This is critical for scientific reproducibility and allows exact
% replication of results when running multiple times.
rng(42, 'twister');  % Mersenne Twister algorithm with seed 42

fprintf('---------------------------------------------------\n');
fprintf('| PRIORITY 4.1-4.5: COMPREHENSIVE CLINICAL ANALYSIS |\n');
fprintf('---------------------------------------------------\n');
fprintf('Start time: %s\n\n', datestr(now));

%% ==========================================================================
%  ANALYSIS PARAMETERS (SESSION 3: FEATURE 3.2)
%  ==========================================================================
%  Centralized configuration constants to improve code maintainability
%  and eliminate magic numbers throughout the analysis pipeline.
%
%  RATIONALE: Using named constants instead of hardcoded values:
%  1. Makes code self-documenting (intent is clear)
%  2. Allows easy parameter tuning for sensitivity analyses
%  3. Ensures consistency across all 4,000+ lines of code
%  4. Facilitates future modifications without hunting for magic numbers

% Statistical Thresholds
% ----------------------
% These control when analyses are performed and how significance is determined
MIN_SAMPLE_SIZE = 30;           % Minimum n for valid statistical analyses (avoids unreliable small-sample estimates)
ALPHA_LEVEL = 0.05;             % Significance level for uncorrected p-values (conventional 5% Type I error rate)
FDR_LEVEL = 0.05;               % False Discovery Rate q-value for Benjamini-Hochberg correction
CI_LEVEL = 0.95;                % Confidence interval level (95% - standard in medical research)
CI_Z_SCORE = 1.96;              % Z-score for 95% CI in normal distribution (two-tailed)

% Effect Size Interpretation Thresholds (Cohen's conventions)
% ------------------------------------------------------------
% Used to interpret practical significance of correlations and group differences
EFFECT_SMALL = 0.10;            % Small effect: r = 0.10, d = 0.20
EFFECT_MEDIUM = 0.30;           % Medium effect: r = 0.30, d = 0.50
EFFECT_LARGE = 0.50;            % Large effect: r = 0.50, d = 0.80
% Reference: Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences.

% Outlier Handling
% ----------------
% Decision scores are z-standardized brain predictions; extreme values indicate model uncertainty
OUTLIER_THRESHOLD_DS = 10;      % Absolute decision score threshold (|z| > 10 is virtually impossible)
OUTLIER_CODE = 99;              % Numeric code used in data files to indicate missing/outlier values

% PCA Parameters
% --------------
% Principal Component Analysis reduces 11 symptom variables to orthogonal components
MIN_VARIANCE_EXPLAINED = 0.70;  % Minimum cumulative variance for PC retention (70% - captures most information)
MAX_N_COMPONENTS = 3;           % Maximum number of PCs to retain (balance between reduction and information)
MIN_PCA_SAMPLES = 50;           % Minimum sample size for stable PCA (rule of thumb: 5× variables minimum)

% Plotting Parameters
% -------------------
% Standardized visual parameters for publication-quality figures
FIGURE_RESOLUTION = 300;        % DPI for saved PNG figures (300 DPI is publication standard)
COLORMAP_CORR = 'redblue';      % Heatmap colormap: blue (negative) to red (positive) correlations
MARKER_SIZE_SCATTER = 50;       % Marker size for scatter plots (optimal visibility)
LINE_WIDTH_REGRESSION = 2;      % Line width for regression lines (clear without overwhelming)
FONT_SIZE_AXIS = 12;            % Font size for axis labels (readable in publications)
FONT_SIZE_TITLE = 14;           % Font size for plot titles (slightly larger for hierarchy)
FOREST_PLOT_WIDTH = 1000;       % Forest plot figure width in pixels
FOREST_PLOT_HEIGHT = 600;       % Forest plot figure height in pixels

% Medication Analysis
% -------------------
MIN_MEDICATION_USERS = 10;      % Minimum n patients on specific medication for reliable analysis

% Age Interaction Analysis
% ------------------------
MIN_GROUP_SIZE_INTERACTION = 3; % Minimum group size for plotting age×diagnosis interactions (avoids 1-2 point groups)
AGE_PREDICTION_POINTS = 100;    % Number of points for smooth age prediction curves in interaction plots

% Variable Redundancy Analysis (NEW - Section 0B)
% -----------------------------------------------
% These parameters control the empirical redundancy detection to avoid FDR inflation
REDUNDANCY_THRESHOLD = 0.70;        % Correlation threshold for empirical variable exclusion (|r| > 0.70 = redundant)
SUBSCALE_TOTAL_THRESHOLD = 0.80;    % Threshold for subscale-total relationship documentation (expected to exceed)
% RATIONALE: Benjamini-Hochberg FDR assumes independence or positive dependence.
% Testing highly correlated variables (e.g., subscales + total) inflates false discoveries.
% Reference: Benjamini & Yekutieli (2001). Control of the false discovery rate under dependency.

fprintf('  ? Analysis parameters initialized\n\n');

%% ==========================================================================
%  SECTION 0: VARIABLE LABEL MAPPING
%  ==========================================================================
% PURPOSE: Create human-readable labels for NESDA variable codes
%
% RATIONALE:
% NESDA uses cryptic variable codes (e.g., 'aids', 'abaiscal') that are not
% self-explanatory. This mapping creates interpretable labels for:
% 1. Publication-quality figures (forest plots, tables)
% 2. Clear communication of results to clinicians
% 3. Reduced errors from misinterpreting variable meanings
% 4. Easier code maintenance and understanding
%
% The labels are used throughout via get_label_safe() function,
% which safely retrieves labels with fallback to variable name if not found.
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 0: VARIABLE LABEL MAPPING\n');
fprintf('---------------------------------------------------\n\n');

% Initialize containers.Map: key-value pairs for variable_code ? readable_label
variable_labels = containers.Map();

% Symptom Severity Variables
variable_labels('aids') = 'Depression Total (IDS-SR)';
variable_labels('aidssev') = 'Depression Severity';
variable_labels('aids_mood_cognition') = 'Depression: Mood/Cognition';
variable_labels('aids_anxiety_arousal') = 'Depression: Anxiety/Arousal';
variable_labels('aidsatyp') = 'Atypical Depression';
variable_labels('aidsmel') = 'Melancholic Depression';
variable_labels('abaiscal') = 'Anxiety Total (BAI)';
variable_labels('abaisev') = 'Anxiety Severity';
variable_labels('abaisom') = 'Anxiety: Somatic';
variable_labels('abaisub') = 'Anxiety: Subjective';
variable_labels('aauditsc') = 'Alcohol Use (AUDIT)';

% Age of Onset
variable_labels('AD2962xAO') = 'Age Onset - MDD';
variable_labels('AD2963xAO') = 'Age Onset - MDD Recurrent';
variable_labels('AD3004AO') = 'Age Onset - Dysthymia';

% Illness Duration
variable_labels('AD2962xDuration') = 'Illness Duration - MDD (yrs)';
variable_labels('AD2963xDuration') = 'Illness Duration - Dysthymia (yrs)';
variable_labels('AD3004Duration') = 'Illness Duration - Any Depression (yrs)';

% Recency Variables (NEW - OPTION 6)
variable_labels('AD2962xRE') = 'Recency - Single MDD';
variable_labels('AD2963xRE') = 'Recency - Recurrent MDD';
variable_labels('AD3004RE') = 'Recency - Dysthymia';

% Clinical History
variable_labels('acidep10') = 'N Current Depression Diagnoses (past 6mo)';
variable_labels('acidep11') = 'Lifetime Depression Present (Y/N)';
variable_labels('acidep13') = 'N MDD Episodes (Lifetime)';
variable_labels('acidep14') = 'MDD Type (Categorical)';
variable_labels('aanxy21') = 'N Current Anxiety Diagnoses (past 6mo)';
variable_labels('aanxy22') = 'Lifetime Anxiety Present (Y/N)';
variable_labels('ANDPBOXSX') = 'N Depressive Symptoms (Lifetime)';
variable_labels('acontrol') = 'Comorbidity Status (Depression x Anxiety)';
variable_labels('afamhdep') = 'Family History Depression';
variable_labels('appfmuse#') = 'N Medications';
variable_labels('atca_ddd') = 'TCA Dose (DDD)';
variable_labels('assri_ddd') = 'SSRI Dose (DDD)';
variable_labels('aotherad_ddd') = 'Other Antidep Dose (DDD)';

% Medication FREQUENCY Variables (_fr suffix: 0=No, 1=Infrequent, 2=Frequent)
variable_labels('assri_fr') = 'SSRI Frequency (0/1/2)';
variable_labels('abenzo_fr') = 'Benzodiazepine Frequency (0/1/2)';
variable_labels('atca_fr') = 'TCA Frequency (0/1/2)';
variable_labels('apsychotropic_fr') = 'Any Psychotropic Frequency (0/1/2)';
variable_labels('aother_ad_fr') = 'Other Antidep Frequency (0/1/2)';
variable_labels('aantipsychotic_fr') = 'Antipsychotic Frequency (0/1/2)';
variable_labels('ahypnotic_sedative_fr') = 'Hypnotic/Sedative Frequency (0/1/2)';
variable_labels('aanxiolytic_fr') = 'Anxiolytic Frequency (0/1/2)';
variable_labels('aother_psychotropic_fr') = 'Other Psychotropic Frequency (0/1/2)';

% Medication BINARY Variables (no _fr suffix: 0=No, 1=Yes)
variable_labels('assri') = 'SSRI Use (Binary Yes/No)';
variable_labels('abenzo') = 'Benzodiazepine Use (Binary Yes/No)';
variable_labels('atca') = 'TCA Use (Binary Yes/No)';
variable_labels('apsychotropic') = 'Any Psychotropic Use (Binary Yes/No)';
variable_labels('aother_ad') = 'Other Antidep Use (Binary Yes/No)';
variable_labels('aantipsychotic') = 'Antipsychotic Use (Binary Yes/No)';
variable_labels('ahypnotic_sedative') = 'Hypnotic/Sedative Use (Binary Yes/No)';
variable_labels('aanxiolytic') = 'Anxiolytic Use (Binary Yes/No)';
variable_labels('aother_psychotropic') = 'Other Psychotropic Use (Binary Yes/No)';

variable_labels('aother_ad_ddd') = 'Other Antidep Dose (DDD)';
variable_labels('n_psychotropic_classes') = 'N Psychotropic Classes';

% TASK 5.3: BENDEP scales
variable_labels('asumbd1') = 'BENDEP: Problematic Use';
variable_labels('asumbd2') = 'BENDEP: Preoccupation';
variable_labels('asumbd3') = 'BENDEP: Lack of Compliance';

% Childhood Adversity
variable_labels('ACTI_total') = 'Childhood Trauma Total';
variable_labels('ACLEI') = 'Childhood Life Events';
variable_labels('aseparation') = 'Parental Separation';
variable_labels('adeathparent') = 'Parental Death';
variable_labels('adivorce') = 'Parental Divorce';

% Demographics
variable_labels('Age') = 'Age (years)';
variable_labels('Sexe') = 'Sex (1=M, 2=F)';
variable_labels('abmi') = 'BMI (kg/m²)';
variable_labels('aedu') = 'Education (years)';
variable_labels('amarpart') = 'Marital Status';
% NOTE: aarea REMOVED - contains interviewer info (bias source), not patient info
variable_labels('aLCAsubtype') = 'Metabolic Subtype';

% Alternative names used in analysis (from Section 8B demographics)
variable_labels('BMI') = 'BMI (kg/m²)';
variable_labels('Sex') = 'Sex (1=M, 2=F)';
variable_labels('Education') = 'Education (years)';
variable_labels('Marital_Status') = 'Marital Status';

% NOTE: Labels corrected based on official NESDA data dictionary (January 2026)
% Corrections made: AD2963xAO, AD3004AO, acidep10-14, aanxy21-22, acontrol

% Helper function - robust version with error handling
get_label = @(v) get_label_safe(v, variable_labels);

fprintf('  Initialized %d variable labels\n\n', variable_labels.Count);


%% ==========================================================================
%  SECTION 0B: VARIABLE REDUNDANCY ANALYSIS & CURATED FEATURE SELECTION
%  ==========================================================================
%  PURPOSE: Systematically identify and exclude redundant variables to
%  ensure valid FDR correction and avoid inflated false discovery rates.
%
%  SCIENTIFIC RATIONALE:
%  The Benjamini-Hochberg FDR procedure assumes independence or positive
%  dependence among test statistics. When testing highly correlated variables
%  (e.g., IDS-SR total + IDS-SR subscales), we violate this assumption,
%  leading to:
%    1. Inflated false discovery rate (more Type I errors than expected)
%    2. Artificial strengthening of findings (if total is significant,
%       subscales will be too, appearing as "multiple" findings)
%    3. Wasted statistical power on redundant hypotheses
%
%  APPROACH:
%    A. Theoretical Exclusions: Remove mathematically redundant variables
%       (subscales that sum to totals, categorical versions of continuous)
%    B. Empirical Exclusions: Check remaining variables for high correlations
%       (|r| > 0.70) and exclude lower-priority redundant variables
%
%  COLLABORATOR REQUEST: Clara Vetter (January 2026) - "Go through clinical
%  variables systematically, check correlations between them, and remove
%  features that associate strongly with each other (like subscale BAI and
%  total BAI) to only have one feature per signal."
%
%  Reference: Benjamini & Yekutieli (2001). Ann Stat 29(4):1165-1188
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 0B: VARIABLE REDUNDANCY ANALYSIS\n');
fprintf('---------------------------------------------------\n\n');

% =========================================================================
% 0B.1: DEFINE VARIABLE DOMAINS AND KNOWN RELATIONSHIPS
% =========================================================================
fprintf('STEP 1: DEFINING VARIABLE DOMAINS AND RELATIONSHIPS\n');
fprintf('---------------------------------------------------\n\n');

% Document variable hierarchies and redundancy relationships
% Format: {variable_to_exclude, variable_to_keep_instead, reason, exclusion_type}
%
% EXCLUSION TYPES:
%   'subscale_total' - Subscale mathematically sums to total
%   'categorical_continuous' - Categorical bins of continuous variable
%   'binary_frequency' - Binary is subset information of frequency coding
%   'redundant_measure' - Same construct measured differently

theoretical_exclusions = {
    % IDS-SR Depression Scale: Subscales sum to total
    'aids_mood_cognition', 'aids', 'Subscale sums to IDS-SR total', 'subscale_total';
    'aids_anxiety_arousal', 'aids', 'Subscale sums to IDS-SR total', 'subscale_total';
    'aidssev', 'aids', 'Categorical bins of continuous total', 'categorical_continuous';

    % BAI Anxiety Scale: Subscales sum to total
    'abaisom', 'abaiscal', 'Somatic subscale sums to BAI total', 'subscale_total';
    'abaisub', 'abaiscal', 'Subjective subscale sums to BAI total', 'subscale_total';
    'abaisev', 'abaiscal', 'Categorical bins of continuous total', 'categorical_continuous';

    % Depression Episode Variables: Keep most informative
    'acidep11', 'acidep13', 'Binary Y/N has less information than episode count', 'binary_frequency';
    'acidep12', 'acidep13', 'Categorical version of episode count', 'categorical_continuous';
    'acidep14', 'acidep13', 'MDD type derived from episodes', 'redundant_measure';
    'acidep10', 'AD3004RE', 'Highly correlated with recency (r=0.83)', 'redundant_measure';

    % Anxiety Episode Variables: Keep count over binary
    'aanxy22', 'aanxy21', 'Binary Y/N has less information than count', 'binary_frequency';

    % Childhood Trauma: Keep comprehensive total
    'ACTI', 'ACTI_total', 'Same items, narrower range (0-4 vs 0-8)', 'redundant_measure';

    % Medication Variables: Keep frequency over binary
    'assri', 'assri_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'abenzo', 'abenzo_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'atca', 'atca_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'apsychotropic', 'apsychotropic_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'aother_ad', 'aother_ad_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'aantipsychotic', 'aantipsychotic_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'ahypnotic_sedative', 'ahypnotic_sedative_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'aanxiolytic', 'aanxiolytic_fr', 'Binary is subset of frequency coding', 'binary_frequency';
    'aother_psychotropic', 'aother_psychotropic_fr', 'Binary is subset of frequency coding', 'binary_frequency';
};

fprintf('  Defined %d theoretical exclusion rules\n', size(theoretical_exclusions, 1));
fprintf('  Exclusion types:\n');
fprintf('    - subscale_total: %d\n', sum(strcmp(theoretical_exclusions(:,4), 'subscale_total')));
fprintf('    - categorical_continuous: %d\n', sum(strcmp(theoretical_exclusions(:,4), 'categorical_continuous')));
fprintf('    - binary_frequency: %d\n', sum(strcmp(theoretical_exclusions(:,4), 'binary_frequency')));
fprintf('    - redundant_measure: %d\n\n', sum(strcmp(theoretical_exclusions(:,4), 'redundant_measure')));

% =========================================================================
% 0B.2: APPLY THEORETICAL EXCLUSIONS (A PRIORI)
% =========================================================================
fprintf('STEP 2: APPLYING THEORETICAL EXCLUSIONS (A PRIORI)\n');
fprintf('---------------------------------------------------\n\n');

% Create list of variables to exclude based on theoretical grounds
theoretical_exclude_list = theoretical_exclusions(:,1);

fprintf('  Variables excluded by theoretical rules:\n');
for i = 1:size(theoretical_exclusions, 1)
    fprintf('    EXCLUDE: %-25s | KEEP: %-15s | %s\n', ...
        theoretical_exclusions{i,1}, theoretical_exclusions{i,2}, theoretical_exclusions{i,3});
end
fprintf('\n');

% =========================================================================
% 0B.3: DEFINE ALL CLINICAL VARIABLES FOR CORRELATION ANALYSIS
% =========================================================================
fprintf('STEP 3: DEFINING CLINICAL VARIABLES FOR CORRELATION MATRIX\n');
fprintf('---------------------------------------------------\n\n');

% Combine all potential clinical variables for correlation analysis
% These will be checked for empirical redundancy after theoretical exclusions
% Categories: Symptom severity, Age of onset, Illness duration, Clinical history,
%             Childhood adversity, Recency, BENDEP scales, Medication frequency, Medication DDD
all_clinical_vars_for_corr = { ...
    'aids', 'aidsatyp', 'aidsmel', 'abaiscal', 'aauditsc', ...
    'AD2962xAO', 'AD2963xAO', 'AD3004AO', ...
    'AD2962xDuration', 'AD2963xDuration', 'AD3004Duration', ...
    'acidep10', 'acidep13', 'aanxy21', 'ANDPBOXSX', 'acontrol', 'afamhdep', ...
    'ACTI_total', 'ACLEI', 'aseparation', 'adeathparent', 'adivorce', ...
    'AD2962xRE', 'AD2963xRE', 'AD3004RE', ...
    'asumbd1', 'asumbd2', 'asumbd3', ...
    'assri_fr', 'abenzo_fr', 'atca_fr', 'apsychotropic_fr', ...
    'assri_ddd', 'atca_ddd', 'aotherad_ddd' ...
};

fprintf('  Variables for correlation analysis: %d\n\n', length(all_clinical_vars_for_corr));

% =========================================================================
% 0B.4: COMPUTE CORRELATION MATRIX FOR REMAINING CLINICAL VARIABLES
% =========================================================================
fprintf('STEP 4: COMPUTING CORRELATION MATRIX\n');
fprintf('---------------------------------------------------\n\n');

% Note: This will be computed after data loading in Section 2
% Store variable list for later use
vars_for_redundancy_check = all_clinical_vars_for_corr;

fprintf('  Correlation matrix will be computed after data loading (Section 5+)\n');
fprintf('  Using pairwise deletion for missing data\n');
fprintf('  Minimum sample size for valid correlation: %d\n\n', MIN_SAMPLE_SIZE);

% =========================================================================
% 0B.5: DEFINE EMPIRICAL EXCLUSION RULES
% =========================================================================
fprintf('STEP 5: DEFINING EMPIRICAL EXCLUSION RULES\n');
fprintf('---------------------------------------------------\n\n');

% Variables to check empirically (may or may not be redundant)
% These pairs will be checked AFTER data loading
empirical_checks = {
    % Variable1, Variable2, Priority (keep priority if r > threshold)
    'aidsatyp', 'aids', 'aids', 'Atypical features may be independent of total severity';
    'aidsmel', 'aids', 'aids', 'Melancholic features may be independent of total severity';
    'aidsatyp', 'aidsmel', 'aidsatyp', 'Expected low/negative correlation (opposite phenotypes)';
    'ACLEI', 'ACTI_total', 'ACTI_total', 'May measure overlapping childhood adversity';
};

fprintf('  Empirical checks to perform: %d pairs\n', size(empirical_checks, 1));
fprintf('  Threshold for empirical exclusion: |r| > %.2f\n\n', REDUNDANCY_THRESHOLD);

% =========================================================================
% 0B.6: DEFINE CURATED VARIABLE SETS (INITIAL - BEFORE EMPIRICAL CHECKS)
% =========================================================================
fprintf('STEP 6: DEFINING CURATED VARIABLE SETS\n');
fprintf('---------------------------------------------------\n\n');

% CURATED SYMPTOM VARIABLES
% Rationale: Keep totals (aids, abaiscal), add AUDIT for alcohol
% Atypical/melancholic features kept tentatively (pending empirical check)
curated_symptom_vars = {'aids', 'abaiscal', 'aauditsc', 'aidsatyp', 'aidsmel'};
fprintf('  Curated symptom variables: %d\n', length(curated_symptom_vars));
fprintf('    %s\n', strjoin(curated_symptom_vars, ', '));
fprintf('    NOTE: aidsatyp and aidsmel will be checked empirically\n\n');

% CURATED CLINICAL HISTORY VARIABLES
% Rationale: Keep episode count (acidep13), current anxiety diagnoses (aanxy21),
% family history, symptom count, comorbidity status
% NOTE: acidep10 REMOVED - highly correlated with AD3004RE (r=0.828) in recency set
curated_clinical_vars = {'acidep13', 'aanxy21', 'ANDPBOXSX', 'afamhdep', 'acontrol'};
fprintf('  Curated clinical variables: %d\n', length(curated_clinical_vars));
fprintf('    %s\n', strjoin(curated_clinical_vars, ', '));
fprintf('    NOTE: acidep10 excluded (r=0.83 with AD3004RE in recency set)\n\n');

% CURATED CHILDHOOD ADVERSITY VARIABLES
% Rationale: Keep comprehensive trauma total, tentatively keep life events
curated_childhood_vars = {'ACTI_total', 'ACLEI'};
fprintf('  Curated childhood variables: %d\n', length(curated_childhood_vars));
fprintf('    %s\n', strjoin(curated_childhood_vars, ', '));
fprintf('    NOTE: ACLEI will be checked empirically against ACTI_total\n\n');

% CURATED AGE OF ONSET VARIABLES
% Rationale: These will be checked for redundancy; may create composite
curated_age_onset_vars = {'AD2962xAO', 'AD2963xAO', 'AD3004AO'};
fprintf('  Curated age of onset variables: %d\n', length(curated_age_onset_vars));
fprintf('    %s\n', strjoin(curated_age_onset_vars, ', '));
fprintf('    NOTE: May create earliest_depression_onset composite if r > %.2f\n\n', REDUNDANCY_THRESHOLD);

% CURATED DURATION VARIABLES
% Rationale: Primary duration measure only
curated_duration_vars = {'AD2962xDuration'};
fprintf('  Curated duration variables: %d\n', length(curated_duration_vars));
fprintf('    %s\n\n', strjoin(curated_duration_vars, ', '));

% CURATED RECENCY VARIABLES
% Rationale: Keep all recency measures (different diagnoses)
curated_recency_vars = {'AD2962xRE', 'AD2963xRE', 'AD3004RE'};
fprintf('  Curated recency variables: %d\n', length(curated_recency_vars));
fprintf('    %s\n\n', strjoin(curated_recency_vars, ', '));

% CURATED MEDICATION VARIABLES
% Rationale: Keep frequency coding only (more info than binary)
% Include polypharmacy count (created in analysis)
curated_medication_vars = {'apsychotropic_fr', 'assri_ddd', 'n_psychotropic_classes'};
fprintf('  Curated medication variables: %d\n', length(curated_medication_vars));
fprintf('    %s\n\n', strjoin(curated_medication_vars, ', '));

% =========================================================================
% 0B.7: SUMMARY OF INITIAL VARIABLE SELECTION
% =========================================================================
fprintf('INITIAL VARIABLE SELECTION SUMMARY\n');
fprintf('---------------------------------------------------\n\n');

% Count original variables (before any exclusions)
total_original_symptom = 11;  % Original symptom_vars list
total_original_clinical = 13; % Original clinical_history_vars
total_original_childhood = 5; % Original childhood_adversity_vars
total_original_medication = 18; % Binary + frequency + DDD

total_original = total_original_symptom + total_original_clinical + ...
                 total_original_childhood + total_original_medication + ...
                 length(curated_age_onset_vars) + length(curated_recency_vars);

% Count curated variables
total_curated = length(curated_symptom_vars) + length(curated_clinical_vars) + ...
                length(curated_childhood_vars) + length(curated_medication_vars) + ...
                length(curated_age_onset_vars) + length(curated_duration_vars) + ...
                length(curated_recency_vars);

fprintf('  ========================================\n');
fprintf('  EXPECTED VARIABLE REDUCTION\n');
fprintf('  ========================================\n');
fprintf('  Domain                  Original  Curated  Reduction\n');
fprintf('  --------------------------------------------------\n');
fprintf('  Symptom Severity        %7d  %7d    %.0f%%\n', ...
    total_original_symptom, length(curated_symptom_vars), ...
    100*(1 - length(curated_symptom_vars)/total_original_symptom));
fprintf('  Clinical History        %7d  %7d    %.0f%%\n', ...
    total_original_clinical, length(curated_clinical_vars), ...
    100*(1 - length(curated_clinical_vars)/total_original_clinical));
fprintf('  Childhood Adversity     %7d  %7d    %.0f%%\n', ...
    total_original_childhood, length(curated_childhood_vars), ...
    100*(1 - length(curated_childhood_vars)/total_original_childhood));
fprintf('  Medication              %7d  %7d    %.0f%%\n', ...
    total_original_medication, length(curated_medication_vars), ...
    100*(1 - length(curated_medication_vars)/total_original_medication));
fprintf('  --------------------------------------------------\n');
fprintf('  TOTAL (approx)          %7d  %7d    %.0f%%\n\n', ...
    total_original, total_curated, 100*(1 - total_curated/total_original));

fprintf('  NOTES:\n');
fprintf('  - Final counts may change after empirical redundancy checks\n');
fprintf('  - Age of onset variables may be combined into composite\n');
fprintf('  - Atypical/melancholic depression kept if independent of total\n\n');

% =========================================================================
% 0B.8: STORE EXCLUSION DECISIONS FOR DOCUMENTATION
% =========================================================================
fprintf('STEP 8: STORING EXCLUSION DECISIONS\n');
fprintf('---------------------------------------------------\n\n');

% Create table for exclusion decisions (will be saved after data loading)
exclusion_decisions = table();
exclusion_decisions.Variable = theoretical_exclusions(:,1);
exclusion_decisions.Kept_Instead = theoretical_exclusions(:,2);
exclusion_decisions.Rule_Category = theoretical_exclusions(:,4);
exclusion_decisions.Reason = theoretical_exclusions(:,3);
exclusion_decisions.Exclusion_Type = repmat({'Theoretical'}, size(theoretical_exclusions, 1), 1);

fprintf('  Exclusion decisions stored: %d entries\n', height(exclusion_decisions));
fprintf('  Will be saved to Variable_Exclusion_Decisions.csv after analysis\n\n');

fprintf('SECTION 0B COMPLETE - Variable selection rules defined\n');
fprintf('  Empirical redundancy checks will be performed in Section 5D\n');
fprintf('  Correlation matrix and visualizations generated in Section 5D\n\n');


%% ==========================================================================
%  SECTION 1: SETUP PATHS AND CREATE OUTPUT DIRECTORIES
%  ==========================================================================
fprintf('SECTION 1: SETTING UP PATHS\n');
fprintf('--------------------------------------------------\n');

base_path = '/volume/projects/CV_NESDA/';
data_path = [base_path 'Data/tabular_data/'];
transition_path_base = [base_path 'Analysis/Transition_Model/Decision_Scores_Mean_Offset/'];
bvftd_path_base = [base_path 'Analysis/bvFTD/Decision_Scores_bvFTD/'];

% Create timestamped output folder for better organization
run_timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
results_path = [base_path 'Analysis/Clinical_Associations_' run_timestamp '/'];

% Diagnosis data paths
diagnosis_hc_file = [base_path 'Data/NESDA_Waves/Wave_1/DynStd_Preparation/NESDA_HC.csv'];
diagnosis_patients_file = [base_path 'Data/NESDA_Waves/Wave_1/DynStd_Preparation/NESDA_Patients.csv'];

if ~exist(results_path, 'dir')
    mkdir(results_path);
end

fig_path = [results_path 'Figures/'];
data_out_path = [results_path 'Data/'];
if ~exist(fig_path, 'dir'), mkdir(fig_path); end
if ~exist(data_out_path, 'dir'), mkdir(data_out_path); end

fprintf('  Output directory: %s\n', results_path);
fprintf('  Figure output: %s\n', fig_path);
fprintf('  Data output: %s\n\n', data_out_path);

diary([results_path 'Priority_4_1_to_4_5_Complete_Analysis_Log_OOCV26_bvFTD.txt']);
fprintf('Logging to: %sPriority_4_1_to_4_5_Complete_Analysis_Log_OOCV26_bvFTD.txt\n\n', results_path);

%% ==========================================================================
%  SECTION 2: LOAD NESDA CLINICAL DATA
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 2: LOADING NESDA CLINICAL DATA\n');
fprintf('---------------------------------------------------\n\n');

nesda_file = [data_path 'NESDA_tabular_combined_data.csv'];
fprintf('Loading: %s\n', nesda_file);

% SESSION 3 FEATURE 3.3: Robust error handling for file loading
if ~exist(nesda_file, 'file')
    error('ERROR: NESDA clinical data file not found: %s\nCheck data_path configuration.', nesda_file);
end

try
    nesda_data = readtable(nesda_file, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
    fprintf('  ? Data loaded: [%d × %d] TABLE\n', height(nesda_data), width(nesda_data));
catch ME
    error('CRITICAL: Failed to load NESDA clinical data from %s\nError: %s\nStack: %s', ...
          nesda_file, ME.message, ME.stack(1).name);
end

% ==========================================================================
% FEATURE 1.2: REMOVE aarea VARIABLE (INTERVIEWER INFO - BIAS SOURCE)
% ==========================================================================
if ismember('aarea', nesda_data.Properties.VariableNames)
    nesda_data = removevars(nesda_data, 'aarea');
    fprintf('  ? Variable aarea removed from analysis (interviewer info, not patient info)\n');
end

varnames = nesda_data.Properties.VariableNames;
fprintf('  First 5 variables: %s\n\n', strjoin(varnames(1:min(5,length(varnames))), ', '));

fprintf('CHECKING FOR KEY VARIABLES:\n');
fprintf('--------------------------------------------------\n');

id_vars = {'pident', 'PIDENT', 'ID', 'SubjectID', 'subject_id'};
id_var = '';
for v = 1:length(id_vars)
    if ismember(id_vars{v}, nesda_data.Properties.VariableNames)
        id_var = id_vars{v};
        fprintf('  ID variable found: %s\n', id_var);
        break;
    end
end
if isempty(id_var)
    error('ERROR: No ID variable found! Searched for: %s', strjoin(id_vars, ', '));
end

sample_ids = nesda_data.(id_var);
if isnumeric(sample_ids)
    fprintf('    Sample IDs (first 5): %d, %d, %d, %d, %d\n', ...
        sample_ids(1), sample_ids(2), sample_ids(3), sample_ids(4), sample_ids(5));
elseif iscell(sample_ids)
    fprintf('    Sample IDs (first 5): %s, %s, %s, %s, %s\n', ...
        sample_ids{1}, sample_ids{2}, sample_ids{3}, sample_ids{4}, sample_ids{5});
end

fprintf('\n');

%% ==========================================================================
%  SECTION 2B: DEFINE ALL 40+ CLINICAL VARIABLES
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 2B: CHECKING CLINICAL VARIABLES\n');
fprintf('---------------------------------------------------\n\n');

symptom_vars = {'aids', 'aidssev', 'aids_mood_cognition', 'aids_anxiety_arousal', ...
                'aidsatyp', 'aidsmel', 'abaiscal', 'abaisev', 'abaisom', 'abaisub', 'aauditsc'};

age_onset_vars = {'AD2962xAO', 'AD2963xAO', 'AD3004AO'};

% NEW: Recency variables (OPTION 6)
recency_vars = {'AD2962xRE', 'AD2963xRE', 'AD3004RE'};

clinical_history_vars = {'acidep10', 'acidep11', 'acidep13', 'acidep14', ...
                        'aanxy21', 'aanxy22', 'ANDPBOXSX', 'acontrol', ...
                        'afamhdep', 'appfmuse#', 'atca_ddd', 'assri_ddd', 'aotherad_ddd'};

childhood_adversity_vars = {'ACTI_total', 'ACLEI', 'aseparation', 'adeathparent', 'adivorce'};

% NOTE: aarea REMOVED - interviewer information (bias source), not patient characteristic
demographic_vars = {'Age', 'Sexe', 'abmi', 'aedu', 'amarpart', 'aLCAsubtype'};

% CORRECTED: Medication FREQUENCY variables (_fr suffix: 0=No, 1=Infrequent, 2=Frequent)
medication_frequency_vars = {'assri_fr', 'abenzo_fr', 'atca_fr', 'apsychotropic_fr', ...
                             'aother_ad_fr', 'aantipsychotic_fr', 'ahypnotic_sedative_fr', ...
                             'aanxiolytic_fr', 'aother_psychotropic_fr'};

% CORRECTED: Medication BINARY variables (NO _fr suffix: 0=No, 1=Yes)
medication_binary_vars = {'assri', 'abenzo', 'atca', 'apsychotropic', ...
                          'aother_ad', 'aantipsychotic', 'ahypnotic_sedative', ...
                          'aanxiolytic', 'aother_psychotropic'};

% Medication DDD variables (continuous dose)
medication_ddd_vars = {'appfmuse#', 'assri_ddd', 'atca_ddd', 'aotherad_ddd', 'aother_ad_ddd'};

% TASK 5.3: NEW - BENDEP scales (3 NEW VARIABLES)
bendep_vars = {'asumbd1', 'asumbd2', 'asumbd3'};

%% ==========================================================================
%  SECTION 2C: CONVERT CELL ARRAY VARIABLES TO NUMERIC
%  ==========================================================================
fprintf('\n---------------------------------------------------\n');
fprintf('SECTION 2C: CONVERTING CELL ARRAYS TO NUMERIC\n');
fprintf('---------------------------------------------------\n\n');

fprintf('CONVERTING AGE OF ONSET VARIABLES:\n');
for i = 1:length(age_onset_vars)
    if ismember(age_onset_vars{i}, nesda_data.Properties.VariableNames)
        var_data = nesda_data.(age_onset_vars{i});
        if iscell(var_data)
            numeric_data = NaN(height(nesda_data), 1);
            for j = 1:length(var_data)
                if ~isempty(var_data{j}) && ~strcmp(var_data{j}, '')
                    if isnumeric(var_data{j})
                        numeric_data(j) = var_data{j};
                    else
                        temp = str2double(var_data{j});
                        if ~isnan(temp)
                            numeric_data(j) = temp;
                        end
                    end
                end
            end
            nesda_data.(age_onset_vars{i}) = numeric_data;
            fprintf('  %s: %d values converted (%.1f%% valid)\n', ...
                age_onset_vars{i}, sum(~isnan(numeric_data)), ...
                100*sum(~isnan(numeric_data))/height(nesda_data));
        end
    end
end
fprintf('\n');

% NEW: Convert recency variables (OPTION 6)
fprintf('CONVERTING RECENCY VARIABLES:\n');
for i = 1:length(recency_vars)
    if ismember(recency_vars{i}, nesda_data.Properties.VariableNames)
        var_data = nesda_data.(recency_vars{i});
        if iscell(var_data)
            numeric_data = NaN(height(nesda_data), 1);
            for j = 1:length(var_data)
                if ~isempty(var_data{j}) && ~strcmp(var_data{j}, '')
                    if isnumeric(var_data{j})
                        numeric_data(j) = var_data{j};
                    else
                        temp = str2double(var_data{j});
                        if ~isnan(temp)
                            numeric_data(j) = temp;
                        end
                    end
                end
            end
            nesda_data.(recency_vars{i}) = numeric_data;
            fprintf('  %s: %d values converted (%.1f%% valid)\n', ...
                recency_vars{i}, sum(~isnan(numeric_data)), ...
                100*sum(~isnan(numeric_data))/height(nesda_data));
        end
    end
end
fprintf('\n');

fprintf('CONVERTING MEDICATION DDD VARIABLES:\n');
medication_vars_to_convert = {'atca_ddd', 'assri_ddd', 'aotherad_ddd', 'aother_ad_ddd'};
for i = 1:length(medication_vars_to_convert)
    if ismember(medication_vars_to_convert{i}, nesda_data.Properties.VariableNames)
        var_data = nesda_data.(medication_vars_to_convert{i});
        if iscell(var_data)
            numeric_data = NaN(height(nesda_data), 1);
            for j = 1:length(var_data)
                if ~isempty(var_data{j}) && ~strcmp(var_data{j}, '')
                    if isnumeric(var_data{j})
                        numeric_data(j) = var_data{j};
                    else
                        temp = str2double(var_data{j});
                        if ~isnan(temp)
                            numeric_data(j) = temp;
                        end
                    end
                end
            end
            nesda_data.(medication_vars_to_convert{i}) = numeric_data;
            fprintf('  %s: %d values converted (%.1f%% valid)\n', ...
                medication_vars_to_convert{i}, sum(~isnan(numeric_data)), ...
                100*sum(~isnan(numeric_data))/height(nesda_data));
        end
    end
end
fprintf('\n');

% TASK 5.3: Convert BENDEP variables
fprintf('CONVERTING BENDEP VARIABLES:\n');
for i = 1:length(bendep_vars)
    if ismember(bendep_vars{i}, nesda_data.Properties.VariableNames)
        var_data = nesda_data.(bendep_vars{i});
        if iscell(var_data)
            numeric_data = NaN(height(nesda_data), 1);
            for j = 1:length(var_data)
                if ~isempty(var_data{j}) && ~strcmp(var_data{j}, '')
                    if isnumeric(var_data{j})
                        numeric_data(j) = var_data{j};
                    else
                        temp = str2double(var_data{j});
                        if ~isnan(temp)
                            numeric_data(j) = temp;
                        end
                    end
                end
            end
            nesda_data.(bendep_vars{i}) = numeric_data;
            fprintf('  %s: %d values converted (%.1f%% valid)\n', ...
                bendep_vars{i}, sum(~isnan(numeric_data)), ...
                100*sum(~isnan(numeric_data))/height(nesda_data));
        end
    end
end
fprintf('\n');

available_symptom_vars = {};
for i = 1:length(symptom_vars)
    if ismember(symptom_vars{i}, varnames)
        available_symptom_vars{end+1} = symptom_vars{i};
    end
end
fprintf('SYMPTOM SEVERITY VARIABLES (%d vars):\n', length(symptom_vars));
fprintf('  Found: %d/%d variables\n', length(available_symptom_vars), length(symptom_vars));
if ~isempty(available_symptom_vars)
    fprintf('  Variables: %s\n', strjoin(available_symptom_vars, ', '));
    if isnumeric(nesda_data.(available_symptom_vars{1}))
        sample_vals = nesda_data.(available_symptom_vars{1})(1:min(5, height(nesda_data)));
        fprintf('  Sample data (%s, first 5): %s\n', available_symptom_vars{1}, mat2str(sample_vals'));
    end
end
fprintf('\n');

available_age_onset_vars = {};
for i = 1:length(age_onset_vars)
    if ismember(age_onset_vars{i}, varnames)
        available_age_onset_vars{end+1} = age_onset_vars{i};
    end
end
fprintf('AGE OF ONSET VARIABLES (3 vars):\n');
fprintf('  Found: %d/%d variables\n', length(available_age_onset_vars), length(age_onset_vars));
if ~isempty(available_age_onset_vars)
    fprintf('  Variables: %s\n', strjoin(available_age_onset_vars, ', '));
    if isnumeric(nesda_data.(available_age_onset_vars{1}))
        sample_vals = nesda_data.(available_age_onset_vars{1})(1:min(5, height(nesda_data)));
        fprintf('  Sample data (%s, first 5): %s\n', available_age_onset_vars{1}, mat2str(sample_vals'));
    end
end
fprintf('\n');

% NEW: Check recency variables (OPTION 6)
available_recency_vars = {};
for i = 1:length(recency_vars)
    if ismember(recency_vars{i}, varnames)
        available_recency_vars{end+1} = recency_vars{i};
    end
end
fprintf('RECENCY VARIABLES (3 vars) - NEW:\n');
fprintf('  Found: %d/%d variables\n', length(available_recency_vars), length(recency_vars));
if ~isempty(available_recency_vars)
    fprintf('  Variables: %s\n', strjoin(available_recency_vars, ', '));
    if isnumeric(nesda_data.(available_recency_vars{1}))
        sample_vals = nesda_data.(available_recency_vars{1})(1:min(5, height(nesda_data)));
        fprintf('  Sample data (%s, first 5): %s\n', available_recency_vars{1}, mat2str(sample_vals'));
    end
end
fprintf('\n');

available_clinical_history_vars = {};
for i = 1:length(clinical_history_vars)
    if ismember(clinical_history_vars{i}, varnames)
        available_clinical_history_vars{end+1} = clinical_history_vars{i};
    end
end
fprintf('CLINICAL HISTORY VARIABLES (13 vars):\n');
fprintf('  Found: %d/%d variables\n', length(available_clinical_history_vars), length(clinical_history_vars));
if ~isempty(available_clinical_history_vars)
    fprintf('  Variables: %s\n', strjoin(available_clinical_history_vars, ', '));
    if isnumeric(nesda_data.(available_clinical_history_vars{1}))
        sample_vals = nesda_data.(available_clinical_history_vars{1})(1:min(5, height(nesda_data)));
        fprintf('  Sample data (%s, first 5): %s\n', available_clinical_history_vars{1}, mat2str(sample_vals'));
    end
end
fprintf('\n');

available_childhood_vars = {};
for i = 1:length(childhood_adversity_vars)
    if ismember(childhood_adversity_vars{i}, varnames)
        available_childhood_vars{end+1} = childhood_adversity_vars{i};
    end
end
fprintf('CHILDHOOD ADVERSITY VARIABLES (5 vars):\n');
fprintf('  Found: %d/%d variables\n', length(available_childhood_vars), length(childhood_adversity_vars));
if ~isempty(available_childhood_vars)
    fprintf('  Variables: %s\n', strjoin(available_childhood_vars, ', '));
    if isnumeric(nesda_data.(available_childhood_vars{1}))
        sample_vals = nesda_data.(available_childhood_vars{1})(1:min(5, height(nesda_data)));
        fprintf('  Sample data (%s, first 5): %s\n', available_childhood_vars{1}, mat2str(sample_vals'));
    end
end
fprintf('\n');

%% ==========================================================================
%  SECTION 2D: CALCULATE ILLNESS DURATION (NEW)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 2D: CALCULATING ILLNESS DURATION\n');
fprintf('---------------------------------------------------\n\n');

fprintf('NEW VARIABLE: Illness Duration = Age - Age of Onset\n\n');

illness_duration_vars = {};
if ismember('Age', nesda_data.Properties.VariableNames)
    age_data = nesda_data.Age;
    
    for i = 1:length(age_onset_vars)
        if ismember(age_onset_vars{i}, nesda_data.Properties.VariableNames)
            age_onset_data = nesda_data.(age_onset_vars{i});
            duration_var_name = strrep(age_onset_vars{i}, 'AO', 'Duration');
            illness_duration = age_data - age_onset_data;
            nesda_data.(duration_var_name) = illness_duration;
            illness_duration_vars{end+1} = duration_var_name;
            
            valid_duration = sum(~isnan(illness_duration));
            fprintf('  %s: %d valid cases (%.1f%%)\n', ...
                duration_var_name, valid_duration, ...
                100*valid_duration/height(nesda_data));
            
            valid_idx = find(~isnan(illness_duration), 5, 'first');
            if ~isempty(valid_idx)
                fprintf('    Sample values (first 5): ');
                for j = 1:length(valid_idx)
                    fprintf('%.1f', illness_duration(valid_idx(j)));
                    if j < length(valid_idx)
                        fprintf(', ');
                    end
                end
                fprintf(' years\n');
            end
        end
    end
else
    fprintf('  Age variable not found - cannot calculate illness duration\n');
end

fprintf('\n');

%% ==========================================================================
%  SECTION 3: LOAD TRANSITION DECISION SCORES
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 3: LOADING TRANSITION DECISION SCORES\n');
fprintf('---------------------------------------------------\n\n');

transition_file_26 = [transition_path_base 'PRS_TransPred_A32_OOCV-26_Predictions_Cl_1PT_vs_NT.csv'];
fprintf('Loading Version A (OOCV-26): %s\n', transition_file_26);

if ~exist(transition_file_26, 'file')
    error('ERROR: Transition OOCV-26 file not found: %s', transition_file_26);
end

% Try tab-delimited first
try
transition_26_tbl = readtable(transition_file_26, 'VariableNamingRule', 'preserve');
% Manually set variable names if they weren't read correctly
if startsWith(transition_26_tbl.Properties.VariableNames{1}, 'Var')
    expected_names = {'Cases', 'PRED_LABEL', 'Mean_Score', 'Std_Score', 'Ens_Prob1', 'Ens_Prob_1', 'PercRank_Score'};
    n_cols = width(transition_26_tbl);
    transition_26_tbl.Properties.VariableNames = expected_names(1:n_cols);
end
fprintf('  OOCV-26 loaded successfully!\n');
catch
    % Fall back to comma-delimited
    transition_26_tbl = readtable(transition_file_26, 'Delimiter', ',', ...
        'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
    fprintf('  OOCV-26 loaded successfully (Comma-delimited)!\n');
end

fprintf('  Dimensions: [%d subjects × %d columns]\n', height(transition_26_tbl), width(transition_26_tbl));
fprintf('  Column names: %s\n', strjoin(transition_26_tbl.Properties.VariableNames, ', '));

transition_ids_26 = transition_26_tbl{:,1};
fprintf('  Subject IDs extracted from column 1\n');

% Case-insensitive search for Mean_Score
varnames_lower = lower(transition_26_tbl.Properties.VariableNames);
decision_col_idx = find(contains(varnames_lower, 'mean') & contains(varnames_lower, 'score'));

if isempty(decision_col_idx)
    error('ERROR: Cannot find Mean_Score column in OOCV-26. Available columns: %s', ...
        strjoin(transition_26_tbl.Properties.VariableNames, ', '));
end

fprintf('  Decision scores extracted from column %d: %s\n', ...
    decision_col_idx(1), transition_26_tbl.Properties.VariableNames{decision_col_idx(1)});

transition_scores_26 = transition_26_tbl{:,decision_col_idx(1)};
n_before_exclusion = length(transition_scores_26);
% SESSION 3 FEATURE 3.2: Use parameterized outlier thresholds
outlier_mask_26 = (transition_scores_26 == OUTLIER_CODE) | (abs(transition_scores_26) > OUTLIER_THRESHOLD_DS);
transition_scores_26(outlier_mask_26) = NaN;
n_excluded = sum(isnan(transition_scores_26)) - sum(isnan(transition_26_tbl{:,decision_col_idx(1)}));

fprintf('  Decision scores extracted from column: %s\n', ...
    transition_26_tbl.Properties.VariableNames{decision_col_idx(1)});
if n_excluded > 0
    fprintf('  Excluded %d subjects with decision score = 99 (missing value code)\n', n_excluded);
end
fprintf('    Valid scores: %d/%d\n', sum(~isnan(transition_scores_26)), n_before_exclusion);
fprintf('    Mean ± SD: %.3f ± %.3f\n\n', mean(transition_scores_26, 'omitnan'), std(transition_scores_26, 'omitnan'));

%% ==========================================================================
%  SECTION 4: LOAD bvFTD DECISION SCORES
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 4: LOADING bvFTD DECISION SCORES\n');
fprintf('---------------------------------------------------\n\n');

bvftd_file = [bvftd_path_base 'ClassModel_bvFTD-HC_A1_OOCV-7_Predictions_Cl_1bvFTD_vs_HC.csv'];
fprintf('Loading bvFTD (OOCV-7): %s\n', bvftd_file);

if ~exist(bvftd_file, 'file')
    error('ERROR: bvFTD file not found: %s', bvftd_file);
end

try
bvftd_tbl = readtable(bvftd_file, 'VariableNamingRule', 'preserve');
% Manually set variable names if they weren't read correctly
if startsWith(bvftd_tbl.Properties.VariableNames{1}, 'Var')
    expected_names = {'Cases', 'PRED_LABEL', 'Mean_Score', 'Std_Score', 'Ens_Prob1', 'Ens_Prob_1', 'PercRank_Score'};
    n_cols = width(bvftd_tbl);
    bvftd_tbl.Properties.VariableNames = expected_names(1:n_cols);
end
fprintf('  bvFTD loaded successfully!\n');
catch
    bvftd_tbl = readtable(bvftd_file, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
    fprintf('  bvFTD loaded successfully (Comma-delimited)!\n');
end

fprintf('  Dimensions: [%d subjects × %d columns]\n', height(bvftd_tbl), width(bvftd_tbl));

bvftd_ids = bvftd_tbl{:,1};
decision_col_idx = find(contains(bvftd_tbl.Properties.VariableNames, 'Mean_Score', 'IgnoreCase', true));
if isempty(decision_col_idx)
    decision_col_idx = find(contains(bvftd_tbl.Properties.VariableNames, 'Decision', 'IgnoreCase', true));
end
bvftd_scores = bvftd_tbl{:,decision_col_idx(1)};

n_before = length(bvftd_scores);
% SESSION 3 FEATURE 3.2: Use parameterized outlier thresholds
outlier_mask_bvftd = (bvftd_scores == OUTLIER_CODE) | (abs(bvftd_scores) > OUTLIER_THRESHOLD_DS);
bvftd_scores(outlier_mask_bvftd) = NaN;
n_excluded = sum(isnan(bvftd_scores)) - sum(isnan(bvftd_tbl{:,decision_col_idx(1)}));

fprintf('  Decision scores extracted\n');
if n_excluded > 0
    fprintf('  Excluded %d subjects with decision score = 99 (missing value code)\n', n_excluded);
end
fprintf('    Valid scores: %d/%d\n', sum(~isnan(bvftd_scores)), n_before);
fprintf('    Mean ± SD: %.3f ± %.3f\n\n', mean(bvftd_scores, 'omitnan'), std(bvftd_scores, 'omitnan'));

%% ==========================================================================
%  SECTION 5: MERGE DATA AND CREATE ANALYSIS DATASET
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 5: MERGING DATASETS\n');
fprintf('---------------------------------------------------\n\n');

nesda_ids = nesda_data.(id_var);
if isnumeric(transition_ids_26) && iscell(nesda_ids)
    nesda_ids = cellfun(@str2double, nesda_ids);
elseif iscell(transition_ids_26) && isnumeric(nesda_ids)
    nesda_ids = arrayfun(@num2str, nesda_ids, 'UniformOutput', false);
end

% SESSION 3 FEATURE 3.3: Robust error handling for ID matching
fprintf('Matching subjects across datasets...\n');
try
    [~, idx_nesda_26, idx_trans_26] = intersect(nesda_ids, transition_ids_26);
    fprintf('  ? Transition-26 matched: %d subjects\n', length(idx_nesda_26));

    if isempty(idx_nesda_26)
        error('No matching IDs between clinical data and Transition-26 decision scores!\nCheck ID variable formats and data alignment.');
    end
catch ME
    error('CRITICAL: ID matching failed for Transition-26:\n%s\nCheck that pident formats match in both files.', ME.message);
end

try
    [~, idx_nesda_bvftd, idx_bvftd] = intersect(nesda_ids, bvftd_ids);
    fprintf('  ? bvFTD matched: %d subjects\n\n', length(idx_nesda_bvftd));

    if isempty(idx_nesda_bvftd)
        warning('No matching IDs between clinical data and bvFTD. Analysis will continue without bvFTD scores.');
    end
catch ME
    warning('ID matching failed for bvFTD: %s\nContinuing without bvFTD scores.', ME.message);
    idx_nesda_bvftd = [];
    idx_bvftd = [];
end

fprintf('Creating analysis dataset...\n');
analysis_data = nesda_data;
analysis_data.Transition_26 = NaN(height(analysis_data), 1);
analysis_data.bvFTD = NaN(height(analysis_data), 1);

analysis_data.Transition_26(idx_nesda_26) = transition_scores_26(idx_trans_26);
analysis_data.bvFTD(idx_nesda_bvftd) = bvftd_scores(idx_bvftd);

fprintf('  Analysis dataset created: [%d subjects × %d variables]\n', ...
    height(analysis_data), width(analysis_data));
fprintf('  With Transition-26: %d (%.1f%%)\n', ...
    sum(~isnan(analysis_data.Transition_26)), ...
    100*sum(~isnan(analysis_data.Transition_26))/height(analysis_data));
fprintf('  With bvFTD: %d (%.1f%%)\n\n', ...
    sum(~isnan(analysis_data.bvFTD)), ...
    100*sum(~isnan(analysis_data.bvFTD))/height(analysis_data));

%% ==========================================================================
%  SECTION 5B: LOAD AND MERGE DIAGNOSIS GROUP (FOR MEDICATION ANALYSIS)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 5B: LOADING DIAGNOSIS GROUP DATA\n');
fprintf('---------------------------------------------------\n\n');

% SESSION 3 FEATURE 3.3: Enhanced error handling for optional diagnosis data
fprintf('Loading HC diagnosis data...\n');
if exist(diagnosis_hc_file, 'file')
    try
        hc_data = readtable(diagnosis_hc_file, 'VariableNamingRule', 'preserve');
        fprintf('  ? HC data loaded: [%d subjects]\n', height(hc_data));
    catch ME
        warning('Failed to load HC diagnosis data: %s\nContinuing without HC diagnosis info.', ME.message);
        hc_data = table();
    end
else
    fprintf('  WARNING: HC file not found: %s\n', diagnosis_hc_file);
    hc_data = table();
end

fprintf('Loading Patients diagnosis data...\n');
if exist(diagnosis_patients_file, 'file')
    try
        patients_data = readtable(diagnosis_patients_file, 'VariableNamingRule', 'preserve');
        fprintf('  ? Patients data loaded: [%d subjects]\n', height(patients_data));
    catch ME
        warning('Failed to load Patients diagnosis data: %s\nContinuing without patient diagnosis groups.', ME.message);
        patients_data = table();
    end
else
    fprintf('  WARNING: Patients file not found: %s\n', diagnosis_patients_file);
    patients_data = table();
end

% Combine HC and Patients data
if ~isempty(hc_data) && ~isempty(patients_data)
    % Add diagnosis_group column if not present
    if ~ismember('diagnosis_group', hc_data.Properties.VariableNames)
        hc_data.diagnosis_group = repmat({'HC'}, height(hc_data), 1);
    end
    
    % Combine
    diagnosis_data = [hc_data; patients_data];
    fprintf('  Combined diagnosis data: [%d subjects]\n', height(diagnosis_data));
    
    % Merge with analysis_data by pident
    diagnosis_id_var = 'pident';
    if ~ismember(diagnosis_id_var, diagnosis_data.Properties.VariableNames)
        fprintf('  ERROR: pident not found in diagnosis data\n');
    else
        % Match IDs
        diagnosis_ids = diagnosis_data.(diagnosis_id_var);
        
        % Convert IDs to same format
        if isnumeric(nesda_ids) && iscell(diagnosis_ids)
            diagnosis_ids = cellfun(@str2double, diagnosis_ids);
        elseif iscell(nesda_ids) && isnumeric(diagnosis_ids)
            diagnosis_ids = arrayfun(@num2str, diagnosis_ids, 'UniformOutput', false);
        end
        
        % Find matching subjects
        [~, idx_analysis, idx_diagnosis] = intersect(nesda_ids, diagnosis_ids);
        fprintf('  Matched %d subjects with diagnosis data\n', length(idx_analysis));
        
        % Add diagnosis_group to analysis_data
        analysis_data.diagnosis_group = repmat({''}, height(analysis_data), 1);
        analysis_data.diagnosis_group(idx_analysis) = diagnosis_data.diagnosis_group(idx_diagnosis);
        
        % Count by diagnosis group
        unique_groups = unique(analysis_data.diagnosis_group);
        fprintf('  Diagnosis groups:\n');
        for i = 1:length(unique_groups)
            if ~isempty(unique_groups{i})
                n_group = sum(strcmp(analysis_data.diagnosis_group, unique_groups{i}));
                fprintf('    %s: %d subjects\n', unique_groups{i}, n_group);
            end
        end
    end
else
    fprintf('  WARNING: Could not load diagnosis data - medication analysis will use all subjects\n');
    analysis_data.diagnosis_group = repmat({''}, height(analysis_data), 1);
end

%% ==========================================================================
%  SECTION 5C: CREATE PATIENT-ONLY DATASET (FILTER OUT HCs)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 5C: FILTERING FOR PATIENTS ONLY\n');
fprintf('---------------------------------------------------\n\n');

% Create patient-only filter (exclude HC)
if ismember('diagnosis_group', analysis_data.Properties.VariableNames)
    patient_mask = ~strcmp(analysis_data.diagnosis_group, '') & ...
                   ~strcmp(analysis_data.diagnosis_group, 'HC');
    
    fprintf('CREATING PATIENT-ONLY DATASET:\n');
    fprintf('  Total subjects: %d\n', height(analysis_data));
    fprintf('  Healthy Controls (HC): %d\n', sum(strcmp(analysis_data.diagnosis_group, 'HC')));
    fprintf('  Patients (non-HC): %d\n', sum(patient_mask));
    fprintf('  Missing diagnosis: %d\n\n', sum(strcmp(analysis_data.diagnosis_group, '')));
    
    % Store full dataset for reference
    analysis_data_full = analysis_data;
    
    % Override analysis_data with patient-only data
    analysis_data = analysis_data(patient_mask, :);
    
    fprintf('  *** All subsequent analyses will use PATIENTS ONLY ***\n');
    fprintf('  New analysis_data dimensions: [%d × %d]\n\n', ...
        height(analysis_data), width(analysis_data));
    
    % Show patient diagnosis groups
    unique_diag = unique(analysis_data.diagnosis_group);
    fprintf('  Patient diagnosis groups:\n');
    for i = 1:length(unique_diag)
        if ~isempty(unique_diag{i}) && ~strcmp(unique_diag{i}, 'HC')
            n_diag = sum(strcmp(analysis_data.diagnosis_group, unique_diag{i}));
            fprintf('    %s: %d (%.1f%%)\n', unique_diag{i}, n_diag, ...
                100*n_diag/height(analysis_data));
        end
    end
    fprintf('\n');
else
    fprintf('WARNING: diagnosis_group variable not found!\n');
    fprintf('         Cannot filter HCs - will analyze all subjects\n\n');
    analysis_data_full = analysis_data;
end

fprintf('SECTION 5C COMPLETE\n\n');

fprintf('\n');

%% ==========================================================================
%  SECTION 5D: EMPIRICAL REDUNDANCY ANALYSIS & CURATED VARIABLE FINALIZATION
%  ==========================================================================
%  PURPOSE: Now that data is loaded, perform empirical redundancy checks,
%  compute correlation matrix, create visualizations, and finalize curated
%  variable sets for use in downstream analyses.
%
%  This section implements the empirical portion of Section 0B's variable
%  selection framework.
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('SECTION 5D: EMPIRICAL REDUNDANCY ANALYSIS\n');
fprintf('---------------------------------------------------\n\n');

% =========================================================================
% 5D.1: COMPUTE CORRELATION MATRIX FOR CLINICAL VARIABLES
% =========================================================================
fprintf('STEP 1: COMPUTING CORRELATION MATRIX\n');
fprintf('---------------------------------------------------\n\n');

% Identify which clinical variables are available in the dataset
available_for_corr = {};
for i = 1:length(vars_for_redundancy_check)
    if ismember(vars_for_redundancy_check{i}, analysis_data.Properties.VariableNames)
        available_for_corr{end+1} = vars_for_redundancy_check{i};
    end
end

fprintf('  Variables available for correlation: %d/%d\n', ...
    length(available_for_corr), length(vars_for_redundancy_check));

% Extract data matrix for correlation computation
n_vars_corr = length(available_for_corr);
corr_data_matrix = NaN(height(analysis_data), n_vars_corr);

for i = 1:n_vars_corr
    var_data = analysis_data.(available_for_corr{i});
    if isnumeric(var_data)
        corr_data_matrix(:, i) = var_data;
    end
end

% Compute pairwise Pearson correlations with pairwise deletion
correlation_matrix = NaN(n_vars_corr, n_vars_corr);
correlation_n_matrix = NaN(n_vars_corr, n_vars_corr);

for i = 1:n_vars_corr
    for j = 1:n_vars_corr
        valid_idx = ~isnan(corr_data_matrix(:,i)) & ~isnan(corr_data_matrix(:,j));
        n_valid = sum(valid_idx);

        if n_valid >= MIN_SAMPLE_SIZE
            r = corr(corr_data_matrix(valid_idx, i), corr_data_matrix(valid_idx, j));
            correlation_matrix(i, j) = r;
            correlation_n_matrix(i, j) = n_valid;
        end
    end
end

fprintf('  Correlation matrix computed: %d × %d\n', n_vars_corr, n_vars_corr);
fprintf('  Min pairwise n: %d, Max pairwise n: %d\n\n', ...
    min(correlation_n_matrix(:), [], 'omitnan'), max(correlation_n_matrix(:), [], 'omitnan'));

% =========================================================================
% 5D.2: PERFORM EMPIRICAL REDUNDANCY CHECKS
% =========================================================================
fprintf('STEP 2: EMPIRICAL REDUNDANCY CHECKS\n');
fprintf('---------------------------------------------------\n\n');

% Initialize list for empirical exclusions
empirical_exclusion_list = {};
empirical_exclusion_reasons = {};

% Check aidsatyp vs aids
idx_aidsatyp = find(strcmp(available_for_corr, 'aidsatyp'));
idx_aids = find(strcmp(available_for_corr, 'aids'));
if ~isempty(idx_aidsatyp) && ~isempty(idx_aids)
    r_atyp_aids = correlation_matrix(idx_aidsatyp, idx_aids);
    fprintf('  aidsatyp vs aids: r = %.3f', r_atyp_aids);
    if abs(r_atyp_aids) > REDUNDANCY_THRESHOLD
        fprintf(' > %.2f → EXCLUDE aidsatyp\n', REDUNDANCY_THRESHOLD);
        empirical_exclusion_list{end+1} = 'aidsatyp';
        empirical_exclusion_reasons{end+1} = sprintf('r=%.3f with aids > %.2f', r_atyp_aids, REDUNDANCY_THRESHOLD);
        % Remove from curated_symptom_vars
        curated_symptom_vars = curated_symptom_vars(~strcmp(curated_symptom_vars, 'aidsatyp'));
    else
        fprintf(' ≤ %.2f → KEEP aidsatyp (independent)\n', REDUNDANCY_THRESHOLD);
    end
else
    fprintf('  aidsatyp or aids not available for empirical check\n');
end

% Check aidsmel vs aids
idx_aidsmel = find(strcmp(available_for_corr, 'aidsmel'));
if ~isempty(idx_aidsmel) && ~isempty(idx_aids)
    r_mel_aids = correlation_matrix(idx_aidsmel, idx_aids);
    fprintf('  aidsmel vs aids: r = %.3f', r_mel_aids);
    if abs(r_mel_aids) > REDUNDANCY_THRESHOLD
        fprintf(' > %.2f → EXCLUDE aidsmel\n', REDUNDANCY_THRESHOLD);
        empirical_exclusion_list{end+1} = 'aidsmel';
        empirical_exclusion_reasons{end+1} = sprintf('r=%.3f with aids > %.2f', r_mel_aids, REDUNDANCY_THRESHOLD);
        % Remove from curated_symptom_vars
        curated_symptom_vars = curated_symptom_vars(~strcmp(curated_symptom_vars, 'aidsmel'));
    else
        fprintf(' ≤ %.2f → KEEP aidsmel (independent)\n', REDUNDANCY_THRESHOLD);
    end
else
    fprintf('  aidsmel or aids not available for empirical check\n');
end

% Check aidsatyp vs aidsmel (document correlation)
if ~isempty(idx_aidsatyp) && ~isempty(idx_aidsmel)
    r_atyp_mel = correlation_matrix(idx_aidsatyp, idx_aidsmel);
    fprintf('  aidsatyp vs aidsmel: r = %.3f', r_atyp_mel);
    if r_atyp_mel < 0
        fprintf(' (negative as expected - opposite phenotypes)\n');
    else
        fprintf(' (positive - unexpected)\n');
    end
end

% Check ACLEI vs ACTI_total
idx_ACLEI = find(strcmp(available_for_corr, 'ACLEI'));
idx_ACTI_total = find(strcmp(available_for_corr, 'ACTI_total'));
if ~isempty(idx_ACLEI) && ~isempty(idx_ACTI_total)
    r_clei_acti = correlation_matrix(idx_ACLEI, idx_ACTI_total);
    fprintf('  ACLEI vs ACTI_total: r = %.3f', r_clei_acti);
    if abs(r_clei_acti) > REDUNDANCY_THRESHOLD
        fprintf(' > %.2f → EXCLUDE ACLEI\n', REDUNDANCY_THRESHOLD);
        empirical_exclusion_list{end+1} = 'ACLEI';
        empirical_exclusion_reasons{end+1} = sprintf('r=%.3f with ACTI_total > %.2f', r_clei_acti, REDUNDANCY_THRESHOLD);
        % Remove from curated_childhood_vars
        curated_childhood_vars = curated_childhood_vars(~strcmp(curated_childhood_vars, 'ACLEI'));
    else
        fprintf(' ≤ %.2f → KEEP ACLEI (independent)\n', REDUNDANCY_THRESHOLD);
    end
else
    fprintf('  ACLEI or ACTI_total not available for empirical check\n');
end

% Check age of onset variables for redundancy
idx_AO_mdd = find(strcmp(available_for_corr, 'AD2962xAO'));
idx_AO_recur = find(strcmp(available_for_corr, 'AD2963xAO'));
idx_AO_dyst = find(strcmp(available_for_corr, 'AD3004AO'));

age_onset_redundant = false;
if ~isempty(idx_AO_mdd) && ~isempty(idx_AO_recur)
    r_ao_12 = correlation_matrix(idx_AO_mdd, idx_AO_recur);
    fprintf('  AD2962xAO vs AD2963xAO: r = %.3f', r_ao_12);
    if abs(r_ao_12) > REDUNDANCY_THRESHOLD
        fprintf(' > %.2f (highly correlated)\n', REDUNDANCY_THRESHOLD);
        age_onset_redundant = true;
    else
        fprintf('\n');
    end
end
if ~isempty(idx_AO_mdd) && ~isempty(idx_AO_dyst)
    r_ao_13 = correlation_matrix(idx_AO_mdd, idx_AO_dyst);
    fprintf('  AD2962xAO vs AD3004AO: r = %.3f', r_ao_13);
    if abs(r_ao_13) > REDUNDANCY_THRESHOLD
        fprintf(' > %.2f (highly correlated)\n', REDUNDANCY_THRESHOLD);
        age_onset_redundant = true;
    else
        fprintf('\n');
    end
end

fprintf('\n');

% =========================================================================
% 5D.3: CREATE EARLIEST DEPRESSION ONSET COMPOSITE VARIABLE
% =========================================================================
fprintf('STEP 3: CREATING EARLIEST DEPRESSION ONSET COMPOSITE\n');
fprintf('---------------------------------------------------\n\n');

% If age of onset variables are highly correlated, create composite
if age_onset_redundant
    fprintf('  Age of onset variables are correlated > %.2f\n', REDUNDANCY_THRESHOLD);
    fprintf('  Creating earliest_depression_onset composite variable...\n');

    % Get available age of onset columns
    ao_cols = {};
    if ismember('AD2962xAO', analysis_data.Properties.VariableNames)
        ao_cols{end+1} = 'AD2962xAO';
    end
    if ismember('AD2963xAO', analysis_data.Properties.VariableNames)
        ao_cols{end+1} = 'AD2963xAO';
    end
    if ismember('AD3004AO', analysis_data.Properties.VariableNames)
        ao_cols{end+1} = 'AD3004AO';
    end

    if ~isempty(ao_cols)
        % Create matrix of age of onset values
        ao_matrix = NaN(height(analysis_data), length(ao_cols));
        for i = 1:length(ao_cols)
            ao_matrix(:, i) = analysis_data.(ao_cols{i});
        end

        % Take minimum (earliest) across all types
        analysis_data.earliest_depression_onset = min(ao_matrix, [], 2);

        n_valid_composite = sum(~isnan(analysis_data.earliest_depression_onset));
        fprintf('  earliest_depression_onset created: n=%d valid values\n', n_valid_composite);
        fprintf('  Range: %.1f - %.1f years\n', ...
            min(analysis_data.earliest_depression_onset, [], 'omitnan'), ...
            max(analysis_data.earliest_depression_onset, [], 'omitnan'));

        % Update curated age onset vars to use composite
        curated_age_onset_vars = {'earliest_depression_onset'};
        fprintf('  Updated curated_age_onset_vars to use composite\n\n');

        % Add label for composite variable
        variable_labels('earliest_depression_onset') = 'Age Onset - Earliest Depression';
    else
        fprintf('  WARNING: No age of onset variables available\n\n');
    end
else
    fprintf('  Age of onset variables not highly correlated\n');
    fprintf('  Keeping individual variables: %s\n\n', strjoin(curated_age_onset_vars, ', '));
end

% =========================================================================
% 5D.4: CREATE CORRELATION HEATMAP VISUALIZATION
% =========================================================================
fprintf('STEP 4: CREATING CORRELATION HEATMAP\n');
fprintf('---------------------------------------------------\n\n');

figure('Position', [100 100 1200 1000]);

% Create heatmap
imagesc(correlation_matrix);

% Use redblue colormap centered at 0
colormap(redblue(256));
caxis([-1 1]);
cb = colorbar();
cb.Label.String = 'Pearson r';

% Add variable labels (truncated if too long)
max_label_len = 12;
var_labels_short = cell(n_vars_corr, 1);
for i = 1:n_vars_corr
    label = get_label(available_for_corr{i});
    if length(label) > max_label_len
        var_labels_short{i} = [label(1:max_label_len-2) '..'];
    else
        var_labels_short{i} = label;
    end
end

set(gca, 'XTick', 1:n_vars_corr, 'XTickLabel', var_labels_short, ...
    'XTickLabelRotation', 45, 'FontSize', 7);
set(gca, 'YTick', 1:n_vars_corr, 'YTickLabel', var_labels_short, 'FontSize', 7);

title('Clinical Variable Correlation Matrix (Redundancy Analysis)', ...
    'FontWeight', 'bold', 'FontSize', 12);
xlabel('Variables', 'FontWeight', 'bold');
ylabel('Variables', 'FontWeight', 'bold');

% Add threshold line annotation
text(0.02, 0.02, sprintf('Redundancy threshold: |r| > %.2f', REDUNDANCY_THRESHOLD), ...
    'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', ...
    'EdgeColor', 'black');

% Save figure
saveas(gcf, [fig_path 'Fig_0B_Correlation_Heatmap_AllVariables.png']);
% saveas(gcf, [fig_path 'Fig_0B_Correlation_Heatmap_AllVariables.fig']);  % Disabled: .fig fails in remote environments
fprintf('  Saved: Fig_0B_Correlation_Heatmap_AllVariables.png/.fig\n\n');

% =========================================================================
% 5D.5: CREATE DENDROGRAM VISUALIZATION
% =========================================================================
fprintf('STEP 5: CREATING VARIABLE DENDROGRAM\n');
fprintf('---------------------------------------------------\n\n');

% Compute distance matrix (1 - |r|)
dist_matrix = 1 - abs(correlation_matrix);

% Handle NaN values by setting them to maximum distance
dist_matrix(isnan(dist_matrix)) = 1;

% Ensure diagonal is exactly 0 (required by squareform)
dist_matrix(logical(eye(size(dist_matrix)))) = 0;

% Ensure matrix is symmetric (required by squareform)
dist_matrix = (dist_matrix + dist_matrix') / 2;

% Convert to condensed distance vector for linkage
dist_vector = squareform(dist_matrix);

% Replace any remaining NaN/Inf with 1
dist_vector(isnan(dist_vector) | isinf(dist_vector)) = 1;

% Hierarchical clustering with average linkage
try
    linkage_tree = linkage(dist_vector, 'average');

    figure('Position', [100 100 1400 600]);

    % Create dendrogram
    [H, T, outperm] = dendrogram(linkage_tree, 0, 'Labels', var_labels_short, ...
        'Orientation', 'top', 'ColorThreshold', 0.30);

    % Add horizontal threshold line (distance = 0.30 corresponds to |r| = 0.70)
    hold on;
    xl = xlim;
    plot(xl, [0.30 0.30], 'r--', 'LineWidth', 2);
    text(xl(2)*0.98, 0.32, sprintf('r = %.2f threshold', REDUNDANCY_THRESHOLD), ...
        'HorizontalAlignment', 'right', 'Color', 'r', 'FontWeight', 'bold');
    hold off;

    % Formatting
    set(gca, 'XTickLabelRotation', 45, 'FontSize', 8);
    ylabel('Distance (1 - |r|)', 'FontWeight', 'bold', 'FontSize', 11);
    title('Clinical Variable Clustering (Hierarchical)', 'FontWeight', 'bold', 'FontSize', 12);

    % Save figure
    saveas(gcf, [fig_path 'Fig_0B_Variable_Dendrogram.png']);
    % saveas(gcf, [fig_path 'Fig_0B_Variable_Dendrogram.fig']);  % Disabled: .fig fails in remote environments
    fprintf('  Saved: Fig_0B_Variable_Dendrogram.png/.fig\n\n');
catch ME
    fprintf('  WARNING: Dendrogram creation failed: %s\n', ME.message);
    fprintf('  Skipping dendrogram visualization\n\n');
end

% =========================================================================
% 5D.6: FINALIZE AND DOCUMENT CURATED VARIABLE SETS
% =========================================================================
fprintf('STEP 6: FINALIZING CURATED VARIABLE SETS\n');
fprintf('---------------------------------------------------\n\n');

% Add empirical exclusions to the exclusion_decisions table
for i = 1:length(empirical_exclusion_list)
    new_row = table();
    new_row.Variable = empirical_exclusion_list(i);
    new_row.Kept_Instead = {''};  % Already documented which to keep
    new_row.Rule_Category = {'empirical_redundancy'};
    new_row.Reason = empirical_exclusion_reasons(i);
    new_row.Exclusion_Type = {'Empirical'};
    exclusion_decisions = [exclusion_decisions; new_row];
end

% Final curated variable sets
fprintf('  FINAL CURATED VARIABLE SETS:\n\n');

fprintf('  curated_symptom_vars (%d):\n', length(curated_symptom_vars));
fprintf('    %s\n\n', strjoin(curated_symptom_vars, ', '));

fprintf('  curated_clinical_vars (%d):\n', length(curated_clinical_vars));
fprintf('    %s\n\n', strjoin(curated_clinical_vars, ', '));

fprintf('  curated_childhood_vars (%d):\n', length(curated_childhood_vars));
fprintf('    %s\n\n', strjoin(curated_childhood_vars, ', '));

fprintf('  curated_age_onset_vars (%d):\n', length(curated_age_onset_vars));
fprintf('    %s\n\n', strjoin(curated_age_onset_vars, ', '));

fprintf('  curated_duration_vars (%d):\n', length(curated_duration_vars));
fprintf('    %s\n\n', strjoin(curated_duration_vars, ', '));

fprintf('  curated_recency_vars (%d):\n', length(curated_recency_vars));
fprintf('    %s\n\n', strjoin(curated_recency_vars, ', '));

fprintf('  curated_medication_vars (%d):\n', length(curated_medication_vars));
fprintf('    %s\n\n', strjoin(curated_medication_vars, ', '));

% Compile all curated variables
all_curated_vars = unique([curated_symptom_vars, curated_clinical_vars, ...
                           curated_childhood_vars, curated_age_onset_vars, ...
                           curated_duration_vars, curated_recency_vars, ...
                           curated_medication_vars]);

fprintf('  TOTAL CURATED VARIABLES: %d\n\n', length(all_curated_vars));

% =========================================================================
% 5D.7: SAVE DOCUMENTATION FILES
% =========================================================================
fprintf('STEP 7: SAVING DOCUMENTATION FILES\n');
fprintf('---------------------------------------------------\n\n');

% Save exclusion decisions
writetable(exclusion_decisions, [data_out_path 'Variable_Exclusion_Decisions.csv']);
fprintf('  Saved: Variable_Exclusion_Decisions.csv (%d exclusions)\n', height(exclusion_decisions));

% Save correlation matrix with row/column names
corr_matrix_table = array2table(correlation_matrix, ...
    'VariableNames', available_for_corr, 'RowNames', available_for_corr);
writetable(corr_matrix_table, [data_out_path 'Variable_Correlation_Matrix.csv'], ...
    'WriteRowNames', true);
fprintf('  Saved: Variable_Correlation_Matrix.csv (%d × %d)\n', n_vars_corr, n_vars_corr);

% Create and save curated variable summary
curated_summary = table();
curated_summary.Domain = {'Symptom_Severity'; 'Clinical_History'; 'Childhood_Adversity'; ...
                          'Age_of_Onset'; 'Duration'; 'Recency'; 'Medication'};
curated_summary.Original_Count = [total_original_symptom; total_original_clinical; ...
                                  total_original_childhood; 3; 3; 3; total_original_medication];
curated_summary.Curated_Count = [length(curated_symptom_vars); length(curated_clinical_vars); ...
                                 length(curated_childhood_vars); length(curated_age_onset_vars); ...
                                 length(curated_duration_vars); length(curated_recency_vars); ...
                                 length(curated_medication_vars)];
curated_summary.Reduction_Pct = 100 * (1 - curated_summary.Curated_Count ./ curated_summary.Original_Count);

writetable(curated_summary, [data_out_path 'Curated_Variable_Summary.csv']);
fprintf('  Saved: Curated_Variable_Summary.csv\n\n');

% =========================================================================
% 5D.8: PRINT FINAL SUMMARY
% =========================================================================
fprintf('========================================================\n');
fprintf('VARIABLE SELECTION COMPLETE - FINAL SUMMARY\n');
fprintf('========================================================\n\n');

total_original_final = sum(curated_summary.Original_Count);
total_curated_final = sum(curated_summary.Curated_Count);
reduction_pct = 100 * (1 - total_curated_final / total_original_final);

fprintf('  Original variables: %d\n', total_original_final);
fprintf('  Curated variables: %d\n', total_curated_final);
fprintf('  Reduction: %.1f%%\n\n', reduction_pct);

fprintf('  SCIENTIFIC BENEFIT:\n');
fprintf('    - FDR correction now methodologically valid\n');
fprintf('    - Results more interpretable (one variable per construct)\n');
fprintf('    - Transparent, reproducible variable selection\n');
fprintf('    - All decisions documented for reviewers\n\n');

fprintf('SECTION 5D COMPLETE\n\n');

%% ==========================================================================
%  HELPER FUNCTIONS
%  ==========================================================================
ci_r = @(r, n) [tanh(atanh(r) - 1.96/sqrt(n-3)), tanh(atanh(r) + 1.96/sqrt(n-3))];

%% ==========================================================================
%  SECTION 6: PRIORITY 4.1 - METABOLIC SUBTYPES & BMI
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|         PRIORITY 4.1: METABOLIC & BMI            |\n');
fprintf('---------------------------------------------------\n\n');

results_4_1 = struct();

if ismember('aLCAsubtype', analysis_data.Properties.VariableNames)
    fprintf('ANALYZING METABOLIC SUBTYPES (aLCAsubtype)...\n');
    
    subtypes = analysis_data.aLCAsubtype;
    subtypes_for_analysis = subtypes;
    subtypes_for_analysis(subtypes == -2) = NaN;
    n_excluded_subtype = sum(subtypes == -2);
    if n_excluded_subtype > 0
        fprintf('  Excluded %d subjects with aLCAsubtype = -2 (not subtyped)\n', n_excluded_subtype);
    end
    
    valid_idx = ~isnan(subtypes_for_analysis) & ~isnan(analysis_data.Transition_26);
    
    unique_subtypes = unique(subtypes_for_analysis(~isnan(subtypes_for_analysis)));
    fprintf('  Unique subtypes (after exclusion): %s\n', mat2str(unique_subtypes));
    fprintf('  Valid cases: %d\n\n', sum(valid_idx));
    
    fprintf('  Testing group differences:\n');

    [p_26, tbl_26, stats_26] = anova1(analysis_data.Transition_26(valid_idx), ...
        subtypes_for_analysis(valid_idx), 'off');
    fprintf('    Transition-26: F=%.3f, p=%.4f', tbl_26{2,5}, p_26);
    if p_26 < 0.05
        fprintf(' SIGNIFICANT\n');
    else
        fprintf('\n');
    end
    results_4_1.metabolic_transition_26_p = p_26;

    valid_idx_bvftd = ~isnan(subtypes_for_analysis) & ~isnan(analysis_data.bvFTD);
    [p_bvftd, tbl_bvftd, stats_bvftd] = anova1(analysis_data.bvFTD(valid_idx_bvftd), ...
        subtypes_for_analysis(valid_idx_bvftd), 'off');
    fprintf('    bvFTD: F=%.3f, p=%.4f', tbl_bvftd{2,5}, p_bvftd);
    if p_bvftd < 0.05
        fprintf(' SIGNIFICANT\n');
    else
        fprintf('\n');
    end
    results_4_1.metabolic_bvftd_p = p_bvftd;

    fig = figure('Position', [100 100 800 400]);

    subplot(1,2,1);
    boxplot(analysis_data.Transition_26(valid_idx), subtypes_for_analysis(valid_idx), ...
        'Colors', [0.2 0.4 0.8], 'Symbol', 'k.');
    ylabel('Transition-26 Score', 'FontWeight', 'bold');
    xlabel('Metabolic Subtype', 'FontWeight', 'bold');
    title(sprintf('Transition-26\np=%.4f', p_26), 'FontWeight', 'bold');
    grid on;

    subplot(1,2,2);
    boxplot(analysis_data.bvFTD(valid_idx_bvftd), subtypes_for_analysis(valid_idx_bvftd), ...
        'Colors', [0.8 0.2 0.2], 'Symbol', 'k.');
    ylabel('bvFTD Score', 'FontWeight', 'bold');
    xlabel('Metabolic Subtype', 'FontWeight', 'bold');
    title(sprintf('bvFTD\np=%.4f', p_bvftd), 'FontWeight', 'bold');
    grid on;

    % Validate and save figure
    if ishghandle(fig) && isvalid(fig)
        saveas(fig, [fig_path 'Fig_4_1_Metabolic_Subtypes.png']);
        % saveas(fig, [fig_path 'Fig_4_1_Metabolic_Subtypes.fig']);  % Disabled: .fig fails in remote environments
        fprintf('\n  Saved: Fig_4_1_Metabolic_Subtypes.png/.fig\n\n');
    else
        warning('Figure handle invalid, skipping save for Metabolic Subtypes plot');
    end
end

if ismember('abmi', analysis_data.Properties.VariableNames)
    fprintf('ANALYZING BMI CORRELATIONS...\n');
    
    bmi = analysis_data.abmi;
    
    valid_bmi_26 = ~isnan(bmi) & ~isnan(analysis_data.Transition_26);
    [r_bmi_26, p_bmi_26] = corr(bmi(valid_bmi_26), analysis_data.Transition_26(valid_bmi_26));
    ci_26 = ci_r(r_bmi_26, sum(valid_bmi_26));
    fprintf('  Transition-26: r=%.3f [%.3f, %.3f], p=%.4f (n=%d)', r_bmi_26, ci_26(1), ci_26(2), p_bmi_26, sum(valid_bmi_26));
    if p_bmi_26 < 0.05
        fprintf(' SIGNIFICANT\n');
    else
        fprintf('\n');
    end
    results_4_1.bmi_transition_26_r = r_bmi_26;
    results_4_1.bmi_transition_26_p = p_bmi_26;

    valid_bmi_bvftd = ~isnan(bmi) & ~isnan(analysis_data.bvFTD);
    [r_bmi_bvftd, p_bmi_bvftd] = corr(bmi(valid_bmi_bvftd), analysis_data.bvFTD(valid_bmi_bvftd));
    ci_bvftd = ci_r(r_bmi_bvftd, sum(valid_bmi_bvftd));
    fprintf('  bvFTD: r=%.3f [%.3f, %.3f], p=%.4f (n=%d)', r_bmi_bvftd, ci_bvftd(1), ci_bvftd(2), p_bmi_bvftd, sum(valid_bmi_bvftd));
    if p_bmi_bvftd < 0.05
        fprintf(' SIGNIFICANT\n');
    else
        fprintf('\n');
    end
    results_4_1.bmi_bvftd_r = r_bmi_bvftd;
    results_4_1.bmi_bvftd_p = p_bmi_bvftd;
    
    figure('Position', [100 100 800 400]);

    subplot(1,2,1);
    scatter(bmi(valid_bmi_26), analysis_data.Transition_26(valid_bmi_26), 50, ...
        [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    p_fit = polyfit(bmi(valid_bmi_26), analysis_data.Transition_26(valid_bmi_26), 1);
    plot(bmi(valid_bmi_26), polyval(p_fit, bmi(valid_bmi_26)), 'r-', 'LineWidth', 2);
    xlabel('BMI (kg/m²)', 'FontWeight', 'bold');
    ylabel('Transition-26 Score', 'FontWeight', 'bold');
    title(sprintf('r=%.3f, p=%.4f', r_bmi_26, p_bmi_26), 'FontWeight', 'bold');
    grid on;

    subplot(1,2,2);
    scatter(bmi(valid_bmi_bvftd), analysis_data.bvFTD(valid_bmi_bvftd), 50, ...
        [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    p_fit = polyfit(bmi(valid_bmi_bvftd), analysis_data.bvFTD(valid_bmi_bvftd), 1);
    plot(bmi(valid_bmi_bvftd), polyval(p_fit, bmi(valid_bmi_bvftd)), 'r-', 'LineWidth', 2);
    xlabel('BMI (kg/m²)', 'FontWeight', 'bold');
    ylabel('bvFTD Score', 'FontWeight', 'bold');
    title(sprintf('r=%.3f, p=%.4f', r_bmi_bvftd, p_bmi_bvftd), 'FontWeight', 'bold');
    grid on;
    
    % SESSION 3 FEATURE 3.3: Robust error handling for plot saving
    try
        saveas(gcf, [fig_path 'Fig_4_1_BMI_Correlations.png']);
        % saveas(gcf, [fig_path 'Fig_4_1_BMI_Correlations.fig']);  % Disabled: .fig fails in remote environments
        fprintf('\n  ? Saved: Fig_4_1_BMI_Correlations.png/.fig\n\n');
    catch ME
        warning('Failed to save BMI correlations figure: %s\nContinuing analysis.', ME.message);
    end
end

fprintf('PRIORITY 4.1 COMPLETE\n\n');

%% ==========================================================================
%  SECTION 7: PRIORITY 4.2 - SYMPTOM SEVERITY (ENHANCED PCA)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|   PRIORITY 4.2: SYMPTOM SEVERITY + ENHANCED PCA  |\n');
fprintf('---------------------------------------------------\n\n');

results_4_2 = struct();

% =========================================================================
% USING CURATED (NON-REDUNDANT) SYMPTOM VARIABLES
% =========================================================================
% Rationale: Section 5D identified and excluded redundant symptom variables
% (subscales that sum to totals, categorical bins of continuous scores).
% Using curated variables ensures valid FDR correction.
% See Section 0B for exclusion criteria and Section 5D for empirical checks.
% =========================================================================
if exist('curated_symptom_vars', 'var') && ~isempty(curated_symptom_vars)
    % Use curated (non-redundant) symptom variables
    fprintf('USING CURATED SYMPTOM VARIABLES (non-redundant set)\n');
    fprintf('  Curated set: %s\n\n', strjoin(curated_symptom_vars, ', '));

    % Filter curated variables to those available in dataset
    available_curated_symptom = {};
    for i = 1:length(curated_symptom_vars)
        if ismember(curated_symptom_vars{i}, analysis_data.Properties.VariableNames)
            available_curated_symptom{end+1} = curated_symptom_vars{i};
        end
    end

    % Use curated variables for this analysis
    available_symptom_vars = available_curated_symptom;
    fprintf('  Available in dataset: %d/%d\n\n', length(available_symptom_vars), length(curated_symptom_vars));
end

if ~isempty(available_symptom_vars)
    fprintf('ANALYZING %d SYMPTOM SEVERITY VARIABLES\n\n', length(available_symptom_vars));
    
    symptom_data = [];
    symptom_names_clean = {};
    for i = 1:length(available_symptom_vars)
        var_data = analysis_data.(available_symptom_vars{i});
        if isnumeric(var_data)
            symptom_data = [symptom_data, var_data];
            symptom_names_clean{end+1} = available_symptom_vars{i};
        end
    end
    
    fprintf('  Extracted %d symptom variables with numeric data\n', size(symptom_data, 2));
    fprintf('  Variables: %s\n\n', strjoin(symptom_names_clean, ', '));
    
    fprintf('CORRELATIONS WITH TRANSITION-26:\n');
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');

    symptom_corr_26 = [];
    for i = 1:length(symptom_names_clean)
        valid_idx = ~isnan(symptom_data(:,i)) & ~isnan(analysis_data.Transition_26);
        if sum(valid_idx) >= 30
            [r, p] = corr(symptom_data(valid_idx,i), analysis_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            symptom_corr_26 = [symptom_corr_26; r, p, n, ci(1), ci(2)];

            fprintf('  %-25s %7.3f %7.4f %5d', symptom_names_clean{i}, r, p, n);
            if p < 0.05
                fprintf(' *\n');
            else
                fprintf('\n');
            end
        else
            symptom_corr_26 = [symptom_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
            fprintf('  %-25s     -       -   %5d (insufficient data)\n', ...
                symptom_names_clean{i}, sum(valid_idx));
        end
    end

    % =====================================================================
    % MULTIPLE TESTING CORRECTION: Benjamini-Hochberg FDR
    % =====================================================================
    % PROBLEM: With 11 symptom variables tested, ~5% false positives expected by chance alone
    %          At a=0.05, we expect 0.55 false positives even if no true effects exist
    % SOLUTION: FDR correction controls the expected proportion of false discoveries
    %          among all discoveries (unlike Bonferroni which controls family-wise error)
    % METHOD: Benjamini & Hochberg (1995) step-up procedure
    %         - Ranks p-values from smallest to largest
    %         - Finds largest i where P(i) = (i/m)×q
    %         - Rejects all hypotheses 1...i
    % BENEFIT: More powerful than Bonferroni (fewer false negatives)
    %          while controlling false discovery rate at q=0.05
    % =====================================================================
    n_uncorrected_sig = sum(symptom_corr_26(:,2) < ALPHA_LEVEL & ~isnan(symptom_corr_26(:,2)));
    [h_fdr_26, crit_p_26, adj_p_26] = fdr_bh(symptom_corr_26(:,2), FDR_LEVEL);
    n_fdr_sig = sum(h_fdr_26);
    fprintf('\n  FDR CORRECTION (q=%.2f): %d/%d significant (uncorrected: %d/%d)\n', ...
        FDR_LEVEL, n_fdr_sig, size(symptom_corr_26,1), n_uncorrected_sig, size(symptom_corr_26,1));
    if crit_p_26 > 0
        fprintf('  Critical p-value: %.4f (original p-values = this are FDR-significant)\n\n', crit_p_26);
    else
        fprintf('  No tests survive FDR correction (all FDR-adjusted p-values > %.2f)\n\n', FDR_LEVEL);
    end

    results_4_2.symptom_correlations_transition_26 = symptom_corr_26;
    results_4_2.symptom_fdr_26 = h_fdr_26;
    results_4_2.symptom_adj_p_26 = adj_p_26;

    fprintf('CORRELATIONS WITH bvFTD:\n');
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');

    symptom_corr_bvftd = [];
    for i = 1:length(symptom_names_clean)
        valid_idx = ~isnan(symptom_data(:,i)) & ~isnan(analysis_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(symptom_data(valid_idx,i), analysis_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            symptom_corr_bvftd = [symptom_corr_bvftd; r, p, n, ci(1), ci(2)];

            fprintf('  %-25s %7.3f %7.4f %5d', symptom_names_clean{i}, r, p, n);
            if p < 0.05
                fprintf(' *\n');
            else
                fprintf('\n');
            end
        else
            symptom_corr_bvftd = [symptom_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end

    % FEATURE 1.3: FDR CORRECTION
    n_uncorrected_sig = sum(symptom_corr_bvftd(:,2) < 0.05 & ~isnan(symptom_corr_bvftd(:,2)));
    [h_fdr_bvftd, crit_p_bvftd, adj_p_bvftd] = fdr_bh(symptom_corr_bvftd(:,2), 0.05);
    n_fdr_sig = sum(h_fdr_bvftd);
    fprintf('\n  FDR CORRECTION (q=0.05): %d/%d significant (uncorrected: %d/%d)\n', ...
        n_fdr_sig, size(symptom_corr_bvftd,1), n_uncorrected_sig, size(symptom_corr_bvftd,1));
    if crit_p_bvftd > 0
        fprintf('  Critical p-value: %.4f\n\n', crit_p_bvftd);
    else
        fprintf('  No tests survive FDR correction (all FDR-adjusted p-values > 0.05)\n\n');
    end

    results_4_2.symptom_correlations_bvftd = symptom_corr_bvftd;
    results_4_2.symptom_fdr_bvftd = h_fdr_bvftd;
    results_4_2.symptom_adj_p_bvftd = adj_p_bvftd;
    
    %% ======================================================================
    %  ENHANCED PCA ANALYSIS - PC1, PC2, PC3 CORRELATIONS
    %  ======================================================================
    fprintf('========================================================\n');
    fprintf('ENHANCED PCA FOR GENERAL SYMPTOM SCORES (PC1, PC2, PC3)\n');
    fprintf('========================================================\n\n');
    
    complete_idx = all(~isnan(symptom_data), 2);
    symptom_data_complete = symptom_data(complete_idx, :);
    fprintf('  Cases with complete symptom data: %d\n', sum(complete_idx));

    % SESSION 3 FEATURE 3.3: Robust error handling for PCA
    if sum(complete_idx) >= MIN_PCA_SAMPLES
        try
            symptom_data_std = zscore(symptom_data_complete);
            [coeff, score, latent, ~, explained] = pca(symptom_data_std);

            fprintf('  ? PCA VARIANCE EXPLAINED:\n');
            fprintf('    PC1: %.1f%%\n', explained(1));
            fprintf('    PC2: %.1f%%\n', explained(2));
            fprintf('    PC3: %.1f%%\n', explained(3));
            fprintf('    First 3 PCs total: %.1f%%\n', sum(explained(1:3)));
            fprintf('    First 5 PCs total: %.1f%%\n\n', sum(explained(1:5)));
        catch ME
            warning('PCA calculation failed: %s\nSkipping PCA analysis and continuing with individual symptom variables.', ME.message);
            coeff = [];
            score = [];
            explained = [];
        end

        % Store PC scores only if PCA succeeded
        if ~isempty(score)
            % Store all three PC scores in analysis_data
            pc1_score = NaN(height(analysis_data), 1);
            pc2_score = NaN(height(analysis_data), 1);
            pc3_score = NaN(height(analysis_data), 1);

            pc1_score(complete_idx) = score(:,1);
            pc2_score(complete_idx) = score(:,2);
            pc3_score(complete_idx) = score(:,3);

            analysis_data.Symptom_PC1 = pc1_score;
            analysis_data.Symptom_PC2 = pc2_score;
            analysis_data.Symptom_PC3 = pc3_score;

            % Also store PC1 as "General_Symptom_Score" for backward compatibility
            analysis_data.General_Symptom_Score = pc1_score;

            fprintf('========================================================\n');
            fprintf('PC1 CORRELATIONS WITH DECISION SCORES\n');
            fprintf('========================================================\n\n');

            % PC1 vs Transition-26
            valid_pc1_26 = ~isnan(pc1_score) & ~isnan(analysis_data.Transition_26);
            [r_pc1_26, p_pc1_26] = corr(pc1_score(valid_pc1_26), ...
                analysis_data.Transition_26(valid_pc1_26));
        n_pc1_26 = sum(valid_pc1_26);
        ci_pc1_26 = ci_r(r_pc1_26, n_pc1_26);
        fprintf('  PC1 vs Transition-26:\n');
        fprintf('    r = %.3f [%.3f, %.3f]\n', r_pc1_26, ci_pc1_26(1), ci_pc1_26(2));
        fprintf('    p = %.4f\n', p_pc1_26);
        fprintf('    n = %d', n_pc1_26);
        if p_pc1_26 < 0.05
            fprintf(' *** SIGNIFICANT ***\n\n');
        else
            fprintf('\n\n');
        end
        results_4_2.pc1_transition_26_r = r_pc1_26;
        results_4_2.pc1_transition_26_p = p_pc1_26;
        results_4_2.pc1_transition_26_n = n_pc1_26;

        % PC1 vs bvFTD
        valid_pc1_bvftd = ~isnan(pc1_score) & ~isnan(analysis_data.bvFTD);
        [r_pc1_bvftd, p_pc1_bvftd] = corr(pc1_score(valid_pc1_bvftd), ...
            analysis_data.bvFTD(valid_pc1_bvftd));
        n_pc1_bvftd = sum(valid_pc1_bvftd);
        ci_pc1_bvftd = ci_r(r_pc1_bvftd, n_pc1_bvftd);
        fprintf('  PC1 vs bvFTD:\n');
        fprintf('    r = %.3f [%.3f, %.3f]\n', r_pc1_bvftd, ci_pc1_bvftd(1), ci_pc1_bvftd(2));
        fprintf('    p = %.4f\n', p_pc1_bvftd);
        fprintf('    n = %d', n_pc1_bvftd);
        if p_pc1_bvftd < 0.05
            fprintf(' *** SIGNIFICANT ***\n\n');
        else
            fprintf('\n\n');
        end
        results_4_2.pc1_bvftd_r = r_pc1_bvftd;
        results_4_2.pc1_bvftd_p = p_pc1_bvftd;
        results_4_2.pc1_bvftd_n = n_pc1_bvftd;
        
        fprintf('========================================================\n');
        fprintf('PC2 CORRELATIONS WITH DECISION SCORES\n');
        fprintf('========================================================\n\n');
        
        % PC2 vs Transition-26
        valid_pc2_26 = ~isnan(pc2_score) & ~isnan(analysis_data.Transition_26);
        [r_pc2_26, p_pc2_26] = corr(pc2_score(valid_pc2_26), ...
            analysis_data.Transition_26(valid_pc2_26));
        n_pc2_26 = sum(valid_pc2_26);
        ci_pc2_26 = ci_r(r_pc2_26, n_pc2_26);
        fprintf('  PC2 vs Transition-26:\n');
        fprintf('    r = %.3f [%.3f, %.3f]\n', r_pc2_26, ci_pc2_26(1), ci_pc2_26(2));
        fprintf('    p = %.4f\n', p_pc2_26);
        fprintf('    n = %d', n_pc2_26);
        if p_pc2_26 < 0.05
            fprintf(' *** SIGNIFICANT ***\n\n');
        else
            fprintf('\n\n');
        end
        results_4_2.pc2_transition_26_r = r_pc2_26;
        results_4_2.pc2_transition_26_p = p_pc2_26;
        results_4_2.pc2_transition_26_n = n_pc2_26;

        % PC2 vs bvFTD
        valid_pc2_bvftd = ~isnan(pc2_score) & ~isnan(analysis_data.bvFTD);
        [r_pc2_bvftd, p_pc2_bvftd] = corr(pc2_score(valid_pc2_bvftd), ...
            analysis_data.bvFTD(valid_pc2_bvftd));
        n_pc2_bvftd = sum(valid_pc2_bvftd);
        ci_pc2_bvftd = ci_r(r_pc2_bvftd, n_pc2_bvftd);
        fprintf('  PC2 vs bvFTD:\n');
        fprintf('    r = %.3f [%.3f, %.3f]\n', r_pc2_bvftd, ci_pc2_bvftd(1), ci_pc2_bvftd(2));
        fprintf('    p = %.4f\n', p_pc2_bvftd);
        fprintf('    n = %d', n_pc2_bvftd);
        if p_pc2_bvftd < 0.05
            fprintf(' *** SIGNIFICANT ***\n\n');
        else
            fprintf('\n\n');
        end
        results_4_2.pc2_bvftd_r = r_pc2_bvftd;
        results_4_2.pc2_bvftd_p = p_pc2_bvftd;
        results_4_2.pc2_bvftd_n = n_pc2_bvftd;
        
        fprintf('========================================================\n');
        fprintf('PC3 CORRELATIONS WITH DECISION SCORES\n');
        fprintf('========================================================\n\n');
        
        % PC3 vs Transition-26
        valid_pc3_26 = ~isnan(pc3_score) & ~isnan(analysis_data.Transition_26);
        [r_pc3_26, p_pc3_26] = corr(pc3_score(valid_pc3_26), ...
            analysis_data.Transition_26(valid_pc3_26));
        n_pc3_26 = sum(valid_pc3_26);
        ci_pc3_26 = ci_r(r_pc3_26, n_pc3_26);
        fprintf('  PC3 vs Transition-26:\n');
        fprintf('    r = %.3f [%.3f, %.3f]\n', r_pc3_26, ci_pc3_26(1), ci_pc3_26(2));
        fprintf('    p = %.4f\n', p_pc3_26);
        fprintf('    n = %d', n_pc3_26);
        if p_pc3_26 < 0.05
            fprintf(' *** SIGNIFICANT ***\n\n');
        else
            fprintf('\n\n');
        end
        results_4_2.pc3_transition_26_r = r_pc3_26;
        results_4_2.pc3_transition_26_p = p_pc3_26;
        results_4_2.pc3_transition_26_n = n_pc3_26;

        % PC3 vs bvFTD
        valid_pc3_bvftd = ~isnan(pc3_score) & ~isnan(analysis_data.bvFTD);
        [r_pc3_bvftd, p_pc3_bvftd] = corr(pc3_score(valid_pc3_bvftd), ...
            analysis_data.bvFTD(valid_pc3_bvftd));
        n_pc3_bvftd = sum(valid_pc3_bvftd);
        ci_pc3_bvftd = ci_r(r_pc3_bvftd, n_pc3_bvftd);
        fprintf('  PC3 vs bvFTD:\n');
        fprintf('    r = %.3f [%.3f, %.3f]\n', r_pc3_bvftd, ci_pc3_bvftd(1), ci_pc3_bvftd(2));
        fprintf('    p = %.4f\n', p_pc3_bvftd);
        fprintf('    n = %d', n_pc3_bvftd);
        if p_pc3_bvftd < 0.05
            fprintf(' *** SIGNIFICANT ***\n\n');
        else
            fprintf('\n\n');
        end
        results_4_2.pc3_bvftd_r = r_pc3_bvftd;
        results_4_2.pc3_bvftd_p = p_pc3_bvftd;
        results_4_2.pc3_bvftd_n = n_pc3_bvftd;
        
        %% ==================================================================
        %  SAVE COMPREHENSIVE PCA CORRELATION TABLE
        %  ==================================================================
        fprintf('========================================================\n');
        fprintf('SAVING COMPREHENSIVE PCA CORRELATION RESULTS\n');
        fprintf('========================================================\n\n');
        
        pca_correlation_summary = table();
        pca_correlation_summary.Component = {'PC1'; 'PC2'; 'PC3'};
        pca_correlation_summary.Variance_Explained = [explained(1); explained(2); explained(3)];
        
        pca_correlation_summary.Trans26_r = [r_pc1_26; r_pc2_26; r_pc3_26];
        pca_correlation_summary.Trans26_p = [p_pc1_26; p_pc2_26; p_pc3_26];
        pca_correlation_summary.Trans26_n = [n_pc1_26; n_pc2_26; n_pc3_26];
        pca_correlation_summary.Trans26_CI_lower = [ci_pc1_26(1); ci_pc2_26(1); ci_pc3_26(1)];
        pca_correlation_summary.Trans26_CI_upper = [ci_pc1_26(2); ci_pc2_26(2); ci_pc3_26(2)];

        pca_correlation_summary.bvFTD_r = [r_pc1_bvftd; r_pc2_bvftd; r_pc3_bvftd];
        pca_correlation_summary.bvFTD_p = [p_pc1_bvftd; p_pc2_bvftd; p_pc3_bvftd];
        pca_correlation_summary.bvFTD_n = [n_pc1_bvftd; n_pc2_bvftd; n_pc3_bvftd];
        pca_correlation_summary.bvFTD_CI_lower = [ci_pc1_bvftd(1); ci_pc2_bvftd(1); ci_pc3_bvftd(1)];
        pca_correlation_summary.bvFTD_CI_upper = [ci_pc1_bvftd(2); ci_pc2_bvftd(2); ci_pc3_bvftd(2)];
        
        writetable(pca_correlation_summary, [data_out_path 'Summary_PCA_Correlations_ALL_Components.csv']);
        fprintf('  Saved: Summary_PCA_Correlations_ALL_Components.csv\n\n');
        
        % Also save PCA loadings (already done, but let's update it)
        pca_results = table();
        pca_results.Variable = symptom_names_clean';
        pca_results.PC1_Loading = coeff(:,1);
        pca_results.PC2_Loading = coeff(:,2);
        pca_results.PC3_Loading = coeff(:,3);
        writetable(pca_results, [data_out_path 'Summary_Symptom_PCA_Loadings.csv']);
        fprintf('  Saved: Summary_Symptom_PCA_Loadings.csv\n\n');
        
        %% ==================================================================
        %  ENHANCED VISUALIZATIONS
        %  ==================================================================
        fprintf('========================================================\n');
        fprintf('CREATING ENHANCED PCA VISUALIZATIONS\n');
        fprintf('========================================================\n\n');
        
        % Figure 1: PCA Variance + PC1/PC2/PC3 Loadings
        fig_pca1 = figure('Position', [100 100 1400 500]);

        subplot(1,4,1);
        pareto(explained(1:min(10, length(explained))));
        xlabel('Principal Component', 'FontWeight', 'bold');
        ylabel('Variance Explained (%)', 'FontWeight', 'bold');
        title('Scree Plot', 'FontWeight', 'bold');

        subplot(1,4,2);
        bar(coeff(1:min(length(symptom_names_clean), 10), 1));
        symptom_labels_pca = cellfun(@(x) get_label(x), symptom_names_clean(1:min(length(symptom_names_clean), 10)), 'UniformOutput', false);
        set(gca, 'XTickLabel', symptom_labels_pca, 'XTickLabelRotation', 45);
        ylabel('PC1 Loading', 'FontWeight', 'bold');
        title(sprintf('PC1 (%.1f%% var)', explained(1)), 'FontWeight', 'bold');
        grid on;

        subplot(1,4,3);
        bar(coeff(1:min(length(symptom_names_clean), 10), 2));
        set(gca, 'XTickLabel', symptom_labels_pca, 'XTickLabelRotation', 45);
        ylabel('PC2 Loading', 'FontWeight', 'bold');
        title(sprintf('PC2 (%.1f%% var)', explained(2)), 'FontWeight', 'bold');
        grid on;

        subplot(1,4,4);
        bar(coeff(1:min(length(symptom_names_clean), 10), 3));
        set(gca, 'XTickLabel', symptom_labels_pca, 'XTickLabelRotation', 45);
        ylabel('PC3 Loading', 'FontWeight', 'bold');
        title(sprintf('PC3 (%.1f%% var)', explained(3)), 'FontWeight', 'bold');
        grid on;

        print(fig_pca1, [fig_path 'Fig_4_2_PCA_Comprehensive_Loadings.png'], '-dpng', '-r300');
        close(fig_pca1);
        % saveas(gcf, [fig_path 'Fig_4_2_PCA_Comprehensive_Loadings.fig']);  % Disabled: .fig fails in remote environments
        fprintf('  Saved: Fig_4_2_PCA_Comprehensive_Loadings.png/.fig\n');
        
        % Figure 2: PC Correlation Heatmap
        figure('Position', [100 100 700 400]);

        corr_matrix = [r_pc1_26, r_pc1_bvftd; ...
                      r_pc2_26, r_pc2_bvftd; ...
                      r_pc3_26, r_pc3_bvftd];

        imagesc(corr_matrix);
        colorbar;
        colormap(redblue);
        caxis([-0.3 0.3]);

        set(gca, 'XTick', 1:2, 'XTickLabel', {'Trans-26', 'bvFTD'}, ...
            'YTick', 1:3, 'YTickLabel', {'PC1', 'PC2', 'PC3'}, 'FontSize', 12);
        title('PCA Component Correlations with Decision Scores', 'FontWeight', 'bold', 'FontSize', 14);
        
        % Add correlation values as text
        for i = 1:3
            for j = 1:2
                text(j, i, sprintf('%.3f', corr_matrix(i,j)), ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 12, 'FontWeight', 'bold', ...
                    'Color', ternary(abs(corr_matrix(i,j)) > 0.15, 'w', 'k'));
            end
        end
        
        saveas(gcf, [fig_path 'Fig_4_2_PCA_Correlation_Heatmap.png']);
        % saveas(gcf, [fig_path 'Fig_4_2_PCA_Correlation_Heatmap.fig']);  % Disabled: .fig fails in remote environments
        fprintf('  Saved: Fig_4_2_PCA_Correlation_Heatmap.png/.fig\n');
        
        % Figure 3: Scatter plots for all significant PC correlations
        sig_pcs = [];
        if p_pc1_26 < 0.05, sig_pcs = [sig_pcs; 1, 26]; end
        if p_pc1_bvftd < 0.05, sig_pcs = [sig_pcs; 1, 0]; end
        if p_pc2_26 < 0.05, sig_pcs = [sig_pcs; 2, 26]; end
        if p_pc2_bvftd < 0.05, sig_pcs = [sig_pcs; 2, 0]; end
        if p_pc3_26 < 0.05, sig_pcs = [sig_pcs; 3, 26]; end
        if p_pc3_bvftd < 0.05, sig_pcs = [sig_pcs; 3, 0]; end
        
        if ~isempty(sig_pcs)
            n_sig = size(sig_pcs, 1);
            n_cols = min(3, n_sig);
            n_rows = ceil(n_sig / n_cols);
            
            figure('Position', [100 100, 400*n_cols, 350*n_rows]);
            
            for i = 1:n_sig
                pc_num = sig_pcs(i, 1);
                model_id = sig_pcs(i, 2);
                
                subplot(n_rows, n_cols, i);
                
                if pc_num == 1
                    pc_data = pc1_score;
                    pc_label = 'PC1 (General Symptom Score)';
                elseif pc_num == 2
                    pc_data = pc2_score;
                    pc_label = 'PC2';
                else
                    pc_data = pc3_score;
                    pc_label = 'PC3';
                end
                
                if model_id == 26
                    ds_data = analysis_data.Transition_26;
                    ds_label = 'Transition-26';
                    color = [0.2 0.4 0.8];
                    if pc_num == 1, r_val = r_pc1_26; p_val = p_pc1_26;
                    elseif pc_num == 2, r_val = r_pc2_26; p_val = p_pc2_26;
                    else, r_val = r_pc3_26; p_val = p_pc3_26; end
                else
                    ds_data = analysis_data.bvFTD;
                    ds_label = 'bvFTD';
                    color = [0.8 0.2 0.2];
                    if pc_num == 1, r_val = r_pc1_bvftd; p_val = p_pc1_bvftd;
                    elseif pc_num == 2, r_val = r_pc2_bvftd; p_val = p_pc2_bvftd;
                    else, r_val = r_pc3_bvftd; p_val = p_pc3_bvftd; end
                end
                
                valid_idx = ~isnan(pc_data) & ~isnan(ds_data);
                scatter(pc_data(valid_idx), ds_data(valid_idx), 50, ...
                    color, 'filled', 'MarkerFaceAlpha', 0.5);
                hold on;
                
                p_fit = polyfit(pc_data(valid_idx), ds_data(valid_idx), 1);
                plot(pc_data(valid_idx), polyval(p_fit, pc_data(valid_idx)), ...
                    'r-', 'LineWidth', 2);
                
                xlabel(pc_label, 'FontWeight', 'bold');
                ylabel([ds_label ' Score'], 'FontWeight', 'bold');
                title(sprintf('%s vs %s\nr=%.3f, p=%.4f', pc_label, ds_label, r_val, p_val), ...
                    'FontWeight', 'bold');
                grid on;
            end  % end for i = 1:n_sig from line 1677

            saveas(gcf, [fig_path 'Fig_4_2_PCA_Significant_Correlations.png']);
            % saveas(gcf, [fig_path 'Fig_4_2_PCA_Significant_Correlations.fig']);  % Disabled: .fig fails in remote environments
            fprintf('  Saved: Fig_4_2_PCA_Significant_Correlations.png/.fig\n');
        end  % end if ~isempty(sig_pcs) from line 1670

        fprintf('\n');
    end  % end if ~isempty(score) from line 1345

    end  % end if sum(complete_idx) >= MIN_PCA_SAMPLES from line 1326

    % Create interpretable labels for heatmap
    symptom_labels = cellfun(@(x) get_label(x), symptom_names_clean, 'UniformOutput', false);
    
    figure('Position', [100 100 1000 500]);

    subplot(1,2,1);
    imagesc(symptom_corr_26(:,1)');
    colorbar;
    colormap(redblue);
    caxis([-0.5 0.5]);
    set(gca, 'XTick', 1:length(symptom_names_clean), 'XTickLabel', symptom_labels, ...
        'XTickLabelRotation', 45, 'YTick', 1, 'YTickLabel', {'r'});
    title('Symptom Severity vs Transition-26', 'FontWeight', 'bold');

    subplot(1,2,2);
    imagesc(symptom_corr_bvftd(:,1)');
    colorbar;
    colormap(redblue);
    caxis([-0.5 0.5]);
    set(gca, 'XTick', 1:length(symptom_names_clean), 'XTickLabel', symptom_labels, ...
        'XTickLabelRotation', 45, 'YTick', 1, 'YTickLabel', {'r'});
    title('Symptom Severity vs bvFTD', 'FontWeight', 'bold');
    
    saveas(gcf, [fig_path 'Fig_4_2_Symptom_Correlations_Heatmap.png']);
    % saveas(gcf, [fig_path 'Fig_4_2_Symptom_Correlations_Heatmap.fig']);  % Disabled: .fig fails in remote environments
    fprintf('\n  Saved: Fig_4_2_Symptom_Correlations_Heatmap.png/.fig\n');
    
    symptom_corr_summary = table();
    symptom_corr_summary.Variable = symptom_names_clean';
    symptom_corr_summary.Transition_26_r = symptom_corr_26(:,1);
    symptom_corr_summary.Transition_26_p = symptom_corr_26(:,2);
    symptom_corr_summary.Transition_26_Uncorrected_significant = symptom_corr_26(:,2) < 0.05;
    symptom_corr_summary.Transition_26_p_FDR = adj_p_26;
    symptom_corr_summary.Transition_26_FDR_significant = h_fdr_26;
    symptom_corr_summary.Transition_26_n = symptom_corr_26(:,3);
    symptom_corr_summary.Transition_26_CI_lower = symptom_corr_26(:,4);
    symptom_corr_summary.Transition_26_CI_upper = symptom_corr_26(:,5);
    symptom_corr_summary.bvFTD_r = symptom_corr_bvftd(:,1);
    symptom_corr_summary.bvFTD_p = symptom_corr_bvftd(:,2);
    symptom_corr_summary.bvFTD_Uncorrected_significant = symptom_corr_bvftd(:,2) < 0.05;
    symptom_corr_summary.bvFTD_p_FDR = adj_p_bvftd;
    symptom_corr_summary.bvFTD_FDR_significant = h_fdr_bvftd;
    symptom_corr_summary.bvFTD_n = symptom_corr_bvftd(:,3);
    symptom_corr_summary.bvFTD_CI_lower = symptom_corr_bvftd(:,4);
    symptom_corr_summary.bvFTD_CI_upper = symptom_corr_bvftd(:,5);
    writetable(symptom_corr_summary, [data_out_path 'Summary_Symptom_Correlations.csv']);
    fprintf('  Saved: Summary_Symptom_Correlations.csv (with FDR correction)\n');
    
end

fprintf('\nPRIORITY 4.2 COMPLETE (WITH ENHANCED PCA)\n\n');

%% ==========================================================================
%  SECTION 8: PRIORITY 4.3 - CLINICAL HISTORY (EXTENDED WITH ILLNESS DURATION)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|      PRIORITY 4.3: CLINICAL HISTORY              |\n');
fprintf('---------------------------------------------------\n\n');

results_4_3 = struct();

if ~isempty(available_age_onset_vars)
    fprintf('ANALYZING AGE OF ONSET VARIABLES\n\n');
    
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    age_onset_corr_26 = [];
    age_onset_corr_27 = [];
    age_onset_corr_bvftd = [];
    
    for i = 1:length(available_age_onset_vars)
        var_data = analysis_data.(available_age_onset_vars{i});
        
        if ~isnumeric(var_data)
            age_onset_corr_26 = [age_onset_corr_26; NaN, NaN, 0, NaN, NaN];
            age_onset_corr_27 = [age_onset_corr_27; NaN, NaN, 0, NaN, NaN];
            age_onset_corr_bvftd = [age_onset_corr_bvftd; NaN, NaN, 0, NaN, NaN];
            continue;
        end
        
        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            age_onset_corr_26 = [age_onset_corr_26; r, p, n, ci(1), ci(2)];
            
            fprintf('  %-25s %7.3f %7.4f %5d', available_age_onset_vars{i}, r, p, n);
            if p < 0.05
                fprintf(' * [Trans-26]\n');
            else
                fprintf('   [Trans-26]\n');
            end
        else
            age_onset_corr_26 = [age_onset_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
        end

        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            age_onset_corr_bvftd = [age_onset_corr_bvftd; r, p, n, ci(1), ci(2)];
        else
            age_onset_corr_bvftd = [age_onset_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end
    fprintf('\n');
    
    results_4_3.age_onset_correlations_26 = age_onset_corr_26;
    results_4_3.age_onset_correlations_bvftd = age_onset_corr_bvftd;
end

if ~isempty(illness_duration_vars)
    fprintf('ANALYZING ILLNESS DURATION VARIABLES (NEW)\n\n');
    
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    duration_corr_26 = [];
    duration_corr_27 = [];
    duration_corr_bvftd = [];
    
    for i = 1:length(illness_duration_vars)
        if ismember(illness_duration_vars{i}, analysis_data.Properties.VariableNames)
            var_data = analysis_data.(illness_duration_vars{i});
            
            valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26);
            if sum(valid_idx) >= 30
                [r, p] = corr(var_data(valid_idx), analysis_data.Transition_26(valid_idx));
                n = sum(valid_idx);
                ci = ci_r(r, n);
                duration_corr_26 = [duration_corr_26; r, p, n, ci(1), ci(2)];
                
                fprintf('  %-25s %7.3f %7.4f %5d', illness_duration_vars{i}, r, p, n);
                if p < 0.05
                    fprintf(' *\n');
                else
                    fprintf('\n');
                end
            else
                duration_corr_26 = [duration_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
            end

            valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD);
            if sum(valid_idx) >= 30
                [r, p] = corr(var_data(valid_idx), analysis_data.bvFTD(valid_idx));
                n = sum(valid_idx);
                ci = ci_r(r, n);
                duration_corr_bvftd = [duration_corr_bvftd; r, p, n, ci(1), ci(2)];
            else
                duration_corr_bvftd = [duration_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
            end
        end
    end
    fprintf('\n');
    
    results_4_3.illness_duration_correlations_26 = duration_corr_26;
    results_4_3.illness_duration_correlations_bvftd = duration_corr_bvftd;
end

% NEW: Analyze recency variables (OPTION 6)
if ~isempty(available_recency_vars)
    fprintf('ANALYZING RECENCY VARIABLES (OPTION 6 - NEW)\n\n');
    
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    recency_corr_26 = [];
    recency_corr_27 = [];
    recency_corr_bvftd = [];
    
    for i = 1:length(available_recency_vars)
        var_data = analysis_data.(available_recency_vars{i});
        
        if ~isnumeric(var_data)
            recency_corr_26 = [recency_corr_26; NaN, NaN, 0, NaN, NaN];
            recency_corr_27 = [recency_corr_27; NaN, NaN, 0, NaN, NaN];
            recency_corr_bvftd = [recency_corr_bvftd; NaN, NaN, 0, NaN, NaN];
            continue;
        end
        
        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            recency_corr_26 = [recency_corr_26; r, p, n, ci(1), ci(2)];
            
            fprintf('  %-25s %7.3f %7.4f %5d', available_recency_vars{i}, r, p, n);
            if p < 0.05
                fprintf(' * [Trans-26]\n');
            else
                fprintf('   [Trans-26]\n');
            end
        else
            recency_corr_26 = [recency_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
        end

        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            recency_corr_bvftd = [recency_corr_bvftd; r, p, n, ci(1), ci(2)];
        else
            recency_corr_bvftd = [recency_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end
    fprintf('\n');
    
    results_4_3.recency_correlations_26 = recency_corr_26;
    results_4_3.recency_correlations_bvftd = recency_corr_bvftd;

    % Save recency correlations summary
    recency_summary = table();
    recency_summary.Variable = available_recency_vars';
    recency_summary.Transition_26_r = recency_corr_26(:,1);
    recency_summary.Transition_26_p = recency_corr_26(:,2);
    recency_summary.Transition_26_Uncorrected_significant = recency_corr_26(:,2) < 0.05;
    recency_summary.Transition_26_n = recency_corr_26(:,3);
    recency_summary.Transition_26_CI_lower = recency_corr_26(:,4);
    recency_summary.Transition_26_CI_upper = recency_corr_26(:,5);
    recency_summary.bvFTD_r = recency_corr_bvftd(:,1);
    recency_summary.bvFTD_p = recency_corr_bvftd(:,2);
    recency_summary.bvFTD_Uncorrected_significant = recency_corr_bvftd(:,2) < 0.05;
    recency_summary.bvFTD_n = recency_corr_bvftd(:,3);
    recency_summary.bvFTD_CI_lower = recency_corr_bvftd(:,4);
    recency_summary.bvFTD_CI_upper = recency_corr_bvftd(:,5);
    writetable(recency_summary, [data_out_path 'Summary_Recency_Correlations.csv']);
    fprintf('  Saved: Summary_Recency_Correlations.csv\n\n');
end

% =========================================================================
% USING CURATED (NON-REDUNDANT) CLINICAL HISTORY VARIABLES
% =========================================================================
% Rationale: Section 5D excluded redundant clinical variables (binary vs count,
% categorical vs continuous). Using curated variables ensures valid FDR.
% See Section 0B for exclusion criteria and Section 5D for empirical checks.
% =========================================================================
if exist('curated_clinical_vars', 'var') && ~isempty(curated_clinical_vars)
    fprintf('USING CURATED CLINICAL HISTORY VARIABLES (non-redundant set)\n');
    fprintf('  Curated set: %s\n\n', strjoin(curated_clinical_vars, ', '));

    % Filter curated variables to those available in dataset
    available_curated_clinical = {};
    for i = 1:length(curated_clinical_vars)
        if ismember(curated_clinical_vars{i}, analysis_data.Properties.VariableNames)
            available_curated_clinical{end+1} = curated_clinical_vars{i};
        end
    end

    % Use curated variables for this analysis
    available_clinical_history_vars = available_curated_clinical;
    fprintf('  Available in dataset: %d/%d\n\n', length(available_clinical_history_vars), length(curated_clinical_vars));
end

if ~isempty(available_clinical_history_vars)
    fprintf('ANALYZING CLINICAL HISTORY VARIABLES\n\n');
    
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    clinical_corr_26 = [];
    clinical_corr_bvftd = [];

    for i = 1:length(available_clinical_history_vars)
        var_data = analysis_data.(available_clinical_history_vars{i});

        if ~isnumeric(var_data)
            clinical_corr_26 = [clinical_corr_26; NaN, NaN, 0, NaN, NaN];
            clinical_corr_bvftd = [clinical_corr_bvftd; NaN, NaN, 0, NaN, NaN];
            continue;
        end
        
        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            clinical_corr_26 = [clinical_corr_26; r, p, n, ci(1), ci(2)];
            
            fprintf('  %-25s %7.3f %7.4f %5d', available_clinical_history_vars{i}, r, p, n);
            if p < 0.05
                fprintf(' *\n');
            else
                fprintf('\n');
            end
        else
            clinical_corr_26 = [clinical_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
            fprintf('  %-25s     -       -   %5d (insufficient data)\n', ...
                available_clinical_history_vars{i}, sum(valid_idx));
        end

        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            clinical_corr_bvftd = [clinical_corr_bvftd; r, p, n, ci(1), ci(2)];
        else
            clinical_corr_bvftd = [clinical_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end
    fprintf('\n');
    
    results_4_3.clinical_history_correlations_26 = clinical_corr_26;
    results_4_3.clinical_history_correlations_bvftd = clinical_corr_bvftd;

    % FEATURE 1.3: FDR CORRECTION
    [h_fdr_clin_26, crit_p_clin_26, adj_p_clin_26] = fdr_bh(clinical_corr_26(:,2), 0.05);
    [h_fdr_clin_bvftd, crit_p_clin_bvftd, adj_p_clin_bvftd] = fdr_bh(clinical_corr_bvftd(:,2), 0.05);
    fprintf('\n  FDR CORRECTION (q=0.05):\n');
    fprintf('    Trans-26: %d/%d significant (uncorrected: %d/%d)\n', ...
        sum(h_fdr_clin_26), length(h_fdr_clin_26), ...
        sum(clinical_corr_26(:,2) < 0.05 & ~isnan(clinical_corr_26(:,2))), length(h_fdr_clin_26));
    fprintf('    bvFTD: %d/%d significant (uncorrected: %d/%d)\n\n', ...
        sum(h_fdr_clin_bvftd), length(h_fdr_clin_bvftd), ...
        sum(clinical_corr_bvftd(:,2) < 0.05 & ~isnan(clinical_corr_bvftd(:,2))), length(h_fdr_clin_bvftd));

    clinical_summary = table();
    clinical_summary.Variable = available_clinical_history_vars';
    clinical_summary.Transition_26_r = clinical_corr_26(:,1);
    clinical_summary.Transition_26_p = clinical_corr_26(:,2);
    clinical_summary.Transition_26_Uncorrected_significant = clinical_corr_26(:,2) < 0.05;
    clinical_summary.Transition_26_p_FDR = adj_p_clin_26;
    clinical_summary.Transition_26_FDR_significant = h_fdr_clin_26;
    clinical_summary.Transition_26_n = clinical_corr_26(:,3);
    clinical_summary.Transition_26_CI_lower = clinical_corr_26(:,4);
    clinical_summary.Transition_26_CI_upper = clinical_corr_26(:,5);
    clinical_summary.bvFTD_r = clinical_corr_bvftd(:,1);
    clinical_summary.bvFTD_p = clinical_corr_bvftd(:,2);
    clinical_summary.bvFTD_Uncorrected_significant = clinical_corr_bvftd(:,2) < 0.05;
    clinical_summary.bvFTD_p_FDR = adj_p_clin_bvftd;
    clinical_summary.bvFTD_FDR_significant = h_fdr_clin_bvftd;
    clinical_summary.bvFTD_n = clinical_corr_bvftd(:,3);
    clinical_summary.bvFTD_CI_lower = clinical_corr_bvftd(:,4);
    clinical_summary.bvFTD_CI_upper = clinical_corr_bvftd(:,5);
    writetable(clinical_summary, [data_out_path 'Summary_Clinical_History_Correlations.csv']);
    fprintf('  Saved: Summary_Clinical_History_Correlations.csv (with FDR correction)\n');
end

% =========================================================================
% USING CURATED (NON-REDUNDANT) CHILDHOOD ADVERSITY VARIABLES
% =========================================================================
% Rationale: Section 5D excluded redundant childhood variables (e.g., ACLEI
% if highly correlated with ACTI_total). Using curated variables ensures
% valid FDR correction. See Section 0B for exclusion criteria.
% =========================================================================
if exist('curated_childhood_vars', 'var') && ~isempty(curated_childhood_vars)
    fprintf('\nUSING CURATED CHILDHOOD ADVERSITY VARIABLES (non-redundant set)\n');
    fprintf('  Curated set: %s\n\n', strjoin(curated_childhood_vars, ', '));

    % Filter curated variables to those available in dataset
    available_curated_childhood = {};
    for i = 1:length(curated_childhood_vars)
        if ismember(curated_childhood_vars{i}, analysis_data.Properties.VariableNames)
            available_curated_childhood{end+1} = curated_childhood_vars{i};
        end
    end

    % Use curated variables for this analysis
    available_childhood_vars = available_curated_childhood;
    fprintf('  Available in dataset: %d/%d\n', length(available_childhood_vars), length(curated_childhood_vars));
end

if ~isempty(available_childhood_vars)
    fprintf('\nANALYZING CHILDHOOD ADVERSITY VARIABLES\n\n');
    
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    childhood_corr_26 = [];
    childhood_corr_bvftd = [];

    for i = 1:length(available_childhood_vars)
        var_data = analysis_data.(available_childhood_vars{i});

        if ~isnumeric(var_data)
            childhood_corr_26 = [childhood_corr_26; NaN, NaN, 0, NaN, NaN];
            childhood_corr_bvftd = [childhood_corr_bvftd; NaN, NaN, 0, NaN, NaN];
            continue;
        end
        
        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            childhood_corr_26 = [childhood_corr_26; r, p, n, ci(1), ci(2)];
            
            fprintf('  %-25s %7.3f %7.4f %5d', available_childhood_vars{i}, r, p, n);
            if p < 0.05
                fprintf(' *\n');
            else
                fprintf('\n');
            end
        else
            childhood_corr_26 = [childhood_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
        end

        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), analysis_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            childhood_corr_bvftd = [childhood_corr_bvftd; r, p, n, ci(1), ci(2)];
        else
            childhood_corr_bvftd = [childhood_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end
    fprintf('\n');
    
    results_4_3.childhood_adversity_correlations_26 = childhood_corr_26;
    results_4_3.childhood_adversity_correlations_bvftd = childhood_corr_bvftd;

    % FEATURE 1.3: FDR CORRECTION
    [h_fdr_child_26, crit_p_child_26, adj_p_child_26] = fdr_bh(childhood_corr_26(:,2), 0.05);
    [h_fdr_child_bvftd, crit_p_child_bvftd, adj_p_child_bvftd] = fdr_bh(childhood_corr_bvftd(:,2), 0.05);
    fprintf('\n  FDR CORRECTION (q=0.05):\n');
    fprintf('    Trans-26: %d/%d significant (uncorrected: %d/%d)\n', ...
        sum(h_fdr_child_26), length(h_fdr_child_26), ...
        sum(childhood_corr_26(:,2) < 0.05 & ~isnan(childhood_corr_26(:,2))), length(h_fdr_child_26));
    fprintf('    bvFTD: %d/%d significant (uncorrected: %d/%d)\n\n', ...
        sum(h_fdr_child_bvftd), length(h_fdr_child_bvftd), ...
        sum(childhood_corr_bvftd(:,2) < 0.05 & ~isnan(childhood_corr_bvftd(:,2))), length(h_fdr_child_bvftd));

    childhood_summary = table();
    childhood_summary.Variable = available_childhood_vars';
    childhood_summary.Transition_26_r = childhood_corr_26(:,1);
    childhood_summary.Transition_26_p = childhood_corr_26(:,2);
    childhood_summary.Transition_26_Uncorrected_significant = childhood_corr_26(:,2) < 0.05;
    childhood_summary.Transition_26_p_FDR = adj_p_child_26;
    childhood_summary.Transition_26_FDR_significant = h_fdr_child_26;
    childhood_summary.Transition_26_n = childhood_corr_26(:,3);
    childhood_summary.Transition_26_CI_lower = childhood_corr_26(:,4);
    childhood_summary.Transition_26_CI_upper = childhood_corr_26(:,5);
    childhood_summary.bvFTD_r = childhood_corr_bvftd(:,1);
    childhood_summary.bvFTD_p = childhood_corr_bvftd(:,2);
    childhood_summary.bvFTD_Uncorrected_significant = childhood_corr_bvftd(:,2) < 0.05;
    childhood_summary.bvFTD_p_FDR = adj_p_child_bvftd;
    childhood_summary.bvFTD_FDR_significant = h_fdr_child_bvftd;
    childhood_summary.bvFTD_n = childhood_corr_bvftd(:,3);
    childhood_summary.bvFTD_CI_lower = childhood_corr_bvftd(:,4);
    childhood_summary.bvFTD_CI_upper = childhood_corr_bvftd(:,5);
    writetable(childhood_summary, [data_out_path 'Summary_Childhood_Adversity_Correlations.csv']);
    fprintf('  Saved: Summary_Childhood_Adversity_Correlations.csv (with FDR correction)\n');
end

fprintf('\nPRIORITY 4.3 COMPLETE\n\n');

%% ==========================================================================
%  SECTION 8B: DEMOGRAPHICS ANALYSIS (NEW)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|    SECTION 8B: DEMOGRAPHICS ANALYSIS (NEW)       |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('ANALYZING DEMOGRAPHIC VARIABLES\n\n');

demo_vars_analyzed = {};
demo_corr_26 = [];
demo_corr_27 = [];
demo_corr_bvftd = [];

if ismember('Age', analysis_data.Properties.VariableNames)
    fprintf('  Age:\n');
    age_data = analysis_data.Age;
    demo_vars_analyzed{end+1} = 'Age';
    
    valid_idx = ~isnan(age_data) & ~isnan(analysis_data.Transition_26);
    if sum(valid_idx) >= 30
        [r, p] = corr(age_data(valid_idx), analysis_data.Transition_26(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_26 = [demo_corr_26; r, p, n, ci(1), ci(2)];
        fprintf('    Transition-26: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_26 = [demo_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
    end

    valid_idx = ~isnan(age_data) & ~isnan(analysis_data.bvFTD);
    if sum(valid_idx) >= 30
        [r, p] = corr(age_data(valid_idx), analysis_data.bvFTD(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_bvftd = [demo_corr_bvftd; r, p, n, ci(1), ci(2)];
        fprintf('    bvFTD: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_bvftd = [demo_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
    end
end

if ismember('Sexe', analysis_data.Properties.VariableNames)
    fprintf('  Sex (1=Male, 2=Female):\n');
    sex_data = analysis_data.Sexe;
    demo_vars_analyzed{end+1} = 'Sex';
    
    valid_idx = ~isnan(sex_data) & ~isnan(analysis_data.Transition_26);
    if sum(valid_idx) >= 30
        [r, p] = corr(sex_data(valid_idx), analysis_data.Transition_26(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_26 = [demo_corr_26; r, p, n, ci(1), ci(2)];
        fprintf('    Transition-26: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_26 = [demo_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
    end

    valid_idx = ~isnan(sex_data) & ~isnan(analysis_data.bvFTD);
    if sum(valid_idx) >= 30
        [r, p] = corr(sex_data(valid_idx), analysis_data.bvFTD(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_bvftd = [demo_corr_bvftd; r, p, n, ci(1), ci(2)];
        fprintf('    bvFTD: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_bvftd = [demo_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
    end
end

if ismember('aedu', analysis_data.Properties.VariableNames)
    fprintf('  Education (years):\n');
    edu_data = analysis_data.aedu;
    demo_vars_analyzed{end+1} = 'Education';
    
    valid_idx = ~isnan(edu_data) & ~isnan(analysis_data.Transition_26);
    if sum(valid_idx) >= 30
        [r, p] = corr(edu_data(valid_idx), analysis_data.Transition_26(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_26 = [demo_corr_26; r, p, n, ci(1), ci(2)];
        fprintf('    Transition-26: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_26 = [demo_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
    end

    valid_idx = ~isnan(edu_data) & ~isnan(analysis_data.bvFTD);
    if sum(valid_idx) >= 30
        [r, p] = corr(edu_data(valid_idx), analysis_data.bvFTD(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_bvftd = [demo_corr_bvftd; r, p, n, ci(1), ci(2)];
        fprintf('    bvFTD: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_bvftd = [demo_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
    end
end

if ismember('amarpart', analysis_data.Properties.VariableNames)
    fprintf('  Marital Status:\n');
    marital_data = analysis_data.amarpart;
    demo_vars_analyzed{end+1} = 'Marital_Status';
    
    valid_idx = ~isnan(marital_data) & ~isnan(analysis_data.Transition_26);
    if sum(valid_idx) >= 30
        [r, p] = corr(marital_data(valid_idx), analysis_data.Transition_26(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_26 = [demo_corr_26; r, p, n, ci(1), ci(2)];
        fprintf('    Transition-26: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_26 = [demo_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
    end

    valid_idx = ~isnan(marital_data) & ~isnan(analysis_data.bvFTD);
    if sum(valid_idx) >= 30
        [r, p] = corr(marital_data(valid_idx), analysis_data.bvFTD(valid_idx));
        n = sum(valid_idx);
        ci = ci_r(r, n);
        demo_corr_bvftd = [demo_corr_bvftd; r, p, n, ci(1), ci(2)];
        fprintf('    bvFTD: r=%.3f [%.3f, %.3f], p=%.4f, n=%d\n\n', r, ci(1), ci(2), p, n);
    else
        demo_corr_bvftd = [demo_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
    end
end

if ~isempty(demo_vars_analyzed)
    % FEATURE 1.3: FDR CORRECTION
    [h_fdr_demo_26, crit_p_demo_26, adj_p_demo_26] = fdr_bh(demo_corr_26(:,2), 0.05);
    [h_fdr_demo_bvftd, crit_p_demo_bvftd, adj_p_demo_bvftd] = fdr_bh(demo_corr_bvftd(:,2), 0.05);
    fprintf('\n  FDR CORRECTION (q=0.05):\n');
    fprintf('    Trans-26: %d/%d significant (uncorrected: %d/%d)\n', ...
        sum(h_fdr_demo_26), length(h_fdr_demo_26), ...
        sum(demo_corr_26(:,2) < 0.05 & ~isnan(demo_corr_26(:,2))), length(h_fdr_demo_26));
    fprintf('    bvFTD: %d/%d significant (uncorrected: %d/%d)\n\n', ...
        sum(h_fdr_demo_bvftd), length(h_fdr_demo_bvftd), ...
        sum(demo_corr_bvftd(:,2) < 0.05 & ~isnan(demo_corr_bvftd(:,2))), length(h_fdr_demo_bvftd));

    demographics_summary = table();
    demographics_summary.Variable = demo_vars_analyzed';
    demographics_summary.Transition_26_r = demo_corr_26(:,1);
    demographics_summary.Transition_26_p = demo_corr_26(:,2);
    demographics_summary.Transition_26_Uncorrected_significant = demo_corr_26(:,2) < 0.05;
    demographics_summary.Transition_26_p_FDR = adj_p_demo_26;
    demographics_summary.Transition_26_FDR_significant = h_fdr_demo_26;
    demographics_summary.Transition_26_n = demo_corr_26(:,3);
    demographics_summary.Transition_26_CI_lower = demo_corr_26(:,4);
    demographics_summary.Transition_26_CI_upper = demo_corr_26(:,5);
    demographics_summary.bvFTD_r = demo_corr_bvftd(:,1);
    demographics_summary.bvFTD_p = demo_corr_bvftd(:,2);
    demographics_summary.bvFTD_Uncorrected_significant = demo_corr_bvftd(:,2) < 0.05;
    demographics_summary.bvFTD_p_FDR = adj_p_demo_bvftd;
    demographics_summary.bvFTD_FDR_significant = h_fdr_demo_bvftd;
    demographics_summary.bvFTD_n = demo_corr_bvftd(:,3);
    demographics_summary.bvFTD_CI_lower = demo_corr_bvftd(:,4);
    demographics_summary.bvFTD_CI_upper = demo_corr_bvftd(:,5);

    writetable(demographics_summary, [data_out_path 'Summary_Demographics_Correlations.csv']);
    fprintf('  Saved: Summary_Demographics_Correlations.csv (with FDR correction)\n');
end

fprintf('\nDEMOGRAPHICS ANALYSIS COMPLETE\n\n');

%% ==========================================================================
%  SECTION 9: PRIORITY 4.4 - COGNITION & FUNCTIONING
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|    PRIORITY 4.4: COGNITION & FUNCTIONING         |\n');
fprintf('---------------------------------------------------\n\n');

results_4_4 = struct();

fprintf('SEARCHING FOR COGNITION & FUNCTIONING VARIABLES\n\n');

cognitive_patterns = {'memory', 'attention', 'executive', 'iq', 'mmse', ...
                     'func', 'gaf', 'wsas', 'whodas', 'disability', 'impairment'};

cognitive_vars_found = {};
for i = 1:length(varnames)
    for j = 1:length(cognitive_patterns)
        if contains(lower(varnames{i}), cognitive_patterns{j})
            cognitive_vars_found{end+1} = varnames{i};
            break;
        end
    end
end

cognitive_vars_found = unique(cognitive_vars_found);

fprintf('  Found %d potential cognition/functioning variables\n', length(cognitive_vars_found));
if length(cognitive_vars_found) > 0
    fprintf('  First 5 variables: %s', strjoin(cognitive_vars_found(1:min(5, length(cognitive_vars_found))), ', '));
    if length(cognitive_vars_found) > 5
        fprintf(' ... (%d more)', length(cognitive_vars_found)-5);
    end
    fprintf('\n');
    
    if isnumeric(analysis_data.(cognitive_vars_found{1}))
        sample_vals = analysis_data.(cognitive_vars_found{1})(1:min(5, height(analysis_data)));
        fprintf('  Sample data (%s, first 5): %s\n', cognitive_vars_found{1}, mat2str(sample_vals'));
    end
end
fprintf('\n');

if ~isempty(cognitive_vars_found)
    fprintf('CORRELATIONS WITH DECISION SCORES\n\n');
    
    fprintf('  Variable                     r        p      n      Model\n');
    fprintf('  ----------------------------------------------------\n');
    
    cognitive_corr_26 = [];
    cognitive_corr_27 = [];
    cognitive_corr_bvftd = [];
    
    for i = 1:length(cognitive_vars_found)
        var_data = analysis_data.(cognitive_vars_found{i});
        
        if isnumeric(var_data)
            valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26);
            if sum(valid_idx) >= 30
                [r, p] = corr(var_data(valid_idx), analysis_data.Transition_26(valid_idx));
                n = sum(valid_idx);
                ci = ci_r(r, n);
                cognitive_corr_26 = [cognitive_corr_26; r, p, n, ci(1), ci(2)];
                
                fprintf('  %-25s %7.3f %7.4f %5d   Trans-26', ...
                    cognitive_vars_found{i}, r, p, n);
                if p < 0.05
                    fprintf(' *\n');
                else
                    fprintf('\n');
                end
            else
                cognitive_corr_26 = [cognitive_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
            end

            valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD);
            if sum(valid_idx) >= 30
                [r, p] = corr(var_data(valid_idx), analysis_data.bvFTD(valid_idx));
                n = sum(valid_idx);
                ci = ci_r(r, n);
                cognitive_corr_bvftd = [cognitive_corr_bvftd; r, p, n, ci(1), ci(2)];
            else
                cognitive_corr_bvftd = [cognitive_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
            end
        end
    end
    fprintf('\n');
    
    results_4_4.cognitive_correlations_26 = cognitive_corr_26;
    results_4_4.cognitive_correlations_bvftd = cognitive_corr_bvftd;
    
    if ~isempty(cognitive_corr_26)
        cognitive_summary = table();
        cognitive_summary.Variable = cognitive_vars_found';
        cognitive_summary.Transition_26_r = cognitive_corr_26(:,1);
        cognitive_summary.Transition_26_p = cognitive_corr_26(:,2);
        cognitive_summary.Transition_26_Uncorrected_significant = cognitive_corr_26(:,2) < 0.05;
        cognitive_summary.Transition_26_n = cognitive_corr_26(:,3);
        cognitive_summary.Transition_26_CI_lower = cognitive_corr_26(:,4);
        cognitive_summary.Transition_26_CI_upper = cognitive_corr_26(:,5);
        cognitive_summary.bvFTD_r = cognitive_corr_bvftd(:,1);
        cognitive_summary.bvFTD_p = cognitive_corr_bvftd(:,2);
        cognitive_summary.bvFTD_Uncorrected_significant = cognitive_corr_bvftd(:,2) < 0.05;
        cognitive_summary.bvFTD_n = cognitive_corr_bvftd(:,3);
        cognitive_summary.bvFTD_CI_lower = cognitive_corr_bvftd(:,4);
        cognitive_summary.bvFTD_CI_upper = cognitive_corr_bvftd(:,5);
        writetable(cognitive_summary, [data_out_path 'Summary_Cognition_Functioning_Correlations.csv']);
        fprintf('  Saved: Summary_Cognition_Functioning_Correlations.csv\n');
    end
end

fprintf('\nPRIORITY 4.4 COMPLETE\n\n');

%% ==========================================================================
%  SECTION 9B: MEDICATION ANALYSIS (PATIENTS ONLY) - CORRECTED CODING
%  ==========================================================================
fprintf('===========================================================\n');
fprintf('|  SECTION 9B: MEDICATION ANALYSIS (CORRECTED CODING)     |\n');
fprintf('===========================================================\n\n');

fprintf('****************************************************************\n');
fprintf('CRITICAL CORRECTION: MEDICATION VARIABLE CODING\n');
fprintf('****************************************************************\n\n');
fprintf('FREQUENCY VARIABLES (_fr suffix):\n');
fprintf('  Coding: 0 = No use, 1 = Infrequent use, 2 = Frequent use\n');
fprintf('  Examples: assri_fr, abenzo_fr, atca_fr, etc.\n\n');
fprintf('BINARY VARIABLES (NO _fr suffix):\n');
fprintf('  Coding: 0 = No use, 1 = Yes (any use)\n');
fprintf('  Examples: assri, abenzo, atca, etc.\n\n');
fprintf('****************************************************************\n\n');

results_medication = struct();

% Since analysis_data is already patients-only, just use it directly
patient_data = analysis_data;
n_patients = height(patient_data);

fprintf('USING PATIENT-ONLY DATASET (already filtered in Section 5C):\n');
fprintf('  N = %d patients\n\n', n_patients);

%% -----------------------------------------------------------------------
%  CORRECTED: DIAGNOSTIC OUTPUT FOR MEDICATION VARIABLES
%  -----------------------------------------------------------------------
fprintf('========================================================\n');
fprintf('DIAGNOSTIC: MEDICATION VARIABLE DISTRIBUTIONS\n');
fprintf('========================================================\n\n');

fprintf('FREQUENCY VARIABLES (_fr suffix: 0/1/2 coding):\n');
fprintf('--------------------------------------------------\n\n');

for i = 1:length(medication_frequency_vars)
    if ismember(medication_frequency_vars{i}, patient_data.Properties.VariableNames)
        var_data = patient_data.(medication_frequency_vars{i});
        
        n_total = length(var_data);
        n_valid = sum(~isnan(var_data));
        
        % Count each frequency category
        n_no_use = sum(var_data == 0);
        n_infrequent = sum(var_data == 1);
        n_frequent = sum(var_data == 2);
        
        fprintf('  %s:\n', medication_frequency_vars{i});
        fprintf('    Valid data: %d/%d (%.1f%%)\n', n_valid, n_total, 100*n_valid/n_total);
        fprintf('    0 (No use): %d (%.1f%%)\n', n_no_use, 100*n_no_use/n_total);
        fprintf('    1 (Infrequent): %d (%.1f%%)\n', n_infrequent, 100*n_infrequent/n_total);
        fprintf('    2 (Frequent): %d (%.1f%%)\n', n_frequent, 100*n_frequent/n_total);
        fprintf('    Missing: %d (%.1f%%)\n\n', n_total-n_valid, 100*(n_total-n_valid)/n_total);
    else
        fprintf('  %s: NOT FOUND\n\n', medication_frequency_vars{i});
    end
end

fprintf('BINARY VARIABLES (NO _fr suffix: 0/1 coding):\n');
fprintf('--------------------------------------------------\n\n');

for i = 1:length(medication_binary_vars)
    if ismember(medication_binary_vars{i}, patient_data.Properties.VariableNames)
        var_data = patient_data.(medication_binary_vars{i});
        
        n_total = length(var_data);
        n_valid = sum(~isnan(var_data));
        
        % Count binary categories
        n_no = sum(var_data == 0);
        n_yes = sum(var_data == 1);
        
        fprintf('  %s:\n', medication_binary_vars{i});
        fprintf('    Valid data: %d/%d (%.1f%%)\n', n_valid, n_total, 100*n_valid/n_total);
        fprintf('    0 (No use): %d (%.1f%%)\n', n_no, 100*n_no/n_total);
        fprintf('    1 (Yes, any use): %d (%.1f%%)\n', n_yes, 100*n_yes/n_total);
        fprintf('    Missing: %d (%.1f%%)\n\n', n_total-n_valid, 100*(n_total-n_valid)/n_total);
    else
        fprintf('  %s: NOT FOUND\n\n', medication_binary_vars{i});
    end
end

%% -----------------------------------------------------------------------
%  CORRECTED: CREATE POLYPHARMACY COUNT VARIABLE
%  -----------------------------------------------------------------------
fprintf('========================================================\n');
fprintf('CREATING POLYPHARMACY COUNT (Using FREQUENCY Variables)\n');
fprintf('========================================================\n\n');

polypharmacy = zeros(height(patient_data), 1);

% Count medication classes using FREQUENCY variables (_fr)
% Count as "using" if frequency is 1 OR 2 (infrequent or frequent)
if ismember('assri_fr', patient_data.Properties.VariableNames)
    polypharmacy = polypharmacy + ((patient_data.assri_fr == 1) | (patient_data.assri_fr == 2));
end
if ismember('atca_fr', patient_data.Properties.VariableNames)
    polypharmacy = polypharmacy + ((patient_data.atca_fr == 1) | (patient_data.atca_fr == 2));
end
if ismember('aother_ad_fr', patient_data.Properties.VariableNames)
    polypharmacy = polypharmacy + ((patient_data.aother_ad_fr == 1) | (patient_data.aother_ad_fr == 2));
end
if ismember('abenzo_fr', patient_data.Properties.VariableNames)
    polypharmacy = polypharmacy + ((patient_data.abenzo_fr == 1) | (patient_data.abenzo_fr == 2));
end
if ismember('aantipsychotic_fr', patient_data.Properties.VariableNames)
    polypharmacy = polypharmacy + ((patient_data.aantipsychotic_fr == 1) | (patient_data.aantipsychotic_fr == 2));
end
if ismember('ahypnotic_sedative_fr', patient_data.Properties.VariableNames)
    polypharmacy = polypharmacy + ((patient_data.ahypnotic_sedative_fr == 1) | (patient_data.ahypnotic_sedative_fr == 2));
end
if ismember('aanxiolytic_fr', patient_data.Properties.VariableNames)
    % Only count if NOT already counted as benzo
    if ismember('abenzo_fr', patient_data.Properties.VariableNames)
        polypharmacy = polypharmacy + (((patient_data.aanxiolytic_fr == 1) | (patient_data.aanxiolytic_fr == 2)) & ...
            (patient_data.abenzo_fr == 0));
    else
        polypharmacy = polypharmacy + ((patient_data.aanxiolytic_fr == 1) | (patient_data.aanxiolytic_fr == 2));
    end
end

patient_data.n_psychotropic_classes = polypharmacy;

fprintf('  Polypharmacy variable created: n_psychotropic_classes\n');
fprintf('  Range: %d to %d classes\n', min(polypharmacy), max(polypharmacy));
fprintf('  Mean: %.2f classes\n', mean(polypharmacy));
fprintf('  Median: %.0f classes\n\n', median(polypharmacy));

% Distribution of polypharmacy
fprintf('  Distribution of psychotropic class count:\n');
unique_counts = unique(polypharmacy(~isnan(polypharmacy)));
for i = 1:length(unique_counts)
    n_count = sum(polypharmacy == unique_counts(i));
    fprintf('    %d classes: %d patients (%.1f%%)\n', unique_counts(i), n_count, 100*n_count/n_patients);
end
fprintf('\n');

%% -----------------------------------------------------------------------
%  CORRECTED: CONTINUOUS MEDICATION VARIABLES (DDD) + BENDEP + POLYPHARMACY
%  -----------------------------------------------------------------------
fprintf('========================================================\n');
fprintf('ANALYZING CONTINUOUS MEDICATION VARIABLES (DDD)\n');
fprintf('========================================================\n\n');

available_med_ddd = {};
med_ddd_corr_26 = [];
med_ddd_corr_27 = [];
med_ddd_corr_bvftd = [];

fprintf('CHECKING AVAILABILITY IN PATIENT SUBSAMPLE:\n');
for i = 1:length(medication_ddd_vars)
    if ismember(medication_ddd_vars{i}, patient_data.Properties.VariableNames)
        var_data = patient_data.(medication_ddd_vars{i});
        
        n_valid = sum(~isnan(var_data) & var_data ~= 0);  % Exclude 0 values for DDD
        pct_valid = 100*n_valid/n_patients;
        
        fprintf('  %-20s: %d/%d non-zero (%.1f%%)', medication_ddd_vars{i}, ...
            n_valid, n_patients, pct_valid);
        
        if n_valid >= 30
            available_med_ddd{end+1} = medication_ddd_vars{i};
            fprintf(' ? INCLUDE\n');
        else
            fprintf(' ? EXCLUDE (too sparse)\n');
        end
    else
        fprintf('  %-20s: NOT FOUND\n', medication_ddd_vars{i});
    end
end
fprintf('\n');

% BENDEP variables
fprintf('CHECKING BENDEP SCALES (BENZODIAZEPINE DEPENDENCE):\n');
for i = 1:length(bendep_vars)
    if ismember(bendep_vars{i}, patient_data.Properties.VariableNames)
        var_data = patient_data.(bendep_vars{i});
        
        n_valid = sum(~isnan(var_data));
        pct_valid = 100*n_valid/n_patients;
        
        fprintf('  %-20s: %d/%d valid (%.1f%%)', bendep_vars{i}, ...
            n_valid, n_patients, pct_valid);
        
        if n_valid >= 30
            available_med_ddd{end+1} = bendep_vars{i};
            fprintf(' ? INCLUDE\n');
        else
            fprintf(' ? EXCLUDE (too sparse)\n');
        end
    else
        fprintf('  %-20s: NOT FOUND\n', bendep_vars{i});
    end
end
fprintf('\n');

% Polypharmacy count
fprintf('CHECKING POLYPHARMACY COUNT:\n');
if ismember('n_psychotropic_classes', patient_data.Properties.VariableNames)
    var_data = patient_data.n_psychotropic_classes;
    n_valid = sum(~isnan(var_data));
    pct_valid = 100*n_valid/n_patients;
    
    fprintf('  %-20s: %d/%d valid (%.1f%%)', 'n_psychotropic_classes', ...
        n_valid, n_patients, pct_valid);
    
    if n_valid >= 30
        available_med_ddd{end+1} = 'n_psychotropic_classes';
        fprintf(' ? INCLUDE\n');
    else
        fprintf(' ? EXCLUDE (too sparse)\n');
    end
else
    fprintf('  n_psychotropic_classes: NOT FOUND\n');
end
fprintf('\n');

if ~isempty(available_med_ddd)
    fprintf('CORRELATIONS WITH DECISION SCORES (PATIENTS ONLY):\n');
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    for i = 1:length(available_med_ddd)
        var_data = patient_data.(available_med_ddd{i});
        
        % Transition-26
        valid_idx = ~isnan(var_data) & ~isnan(patient_data.Transition_26);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), patient_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            med_ddd_corr_26 = [med_ddd_corr_26; r, p, n, ci(1), ci(2)];
            
            fprintf('  %-25s %7.3f %7.4f %5d', available_med_ddd{i}, r, p, n);
            if p < 0.05
                fprintf(' * [Trans-26]\n');
            else
                fprintf('   [Trans-26]\n');
            end
        else
            med_ddd_corr_26 = [med_ddd_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
        end

        % bvFTD
        valid_idx = ~isnan(var_data) & ~isnan(patient_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), patient_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            med_ddd_corr_bvftd = [med_ddd_corr_bvftd; r, p, n, ci(1), ci(2)];
        else
            med_ddd_corr_bvftd = [med_ddd_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end
    fprintf('\n');
    
    results_medication.ddd_correlations_26 = med_ddd_corr_26;
    results_medication.ddd_correlations_bvftd = med_ddd_corr_bvftd;
else
    fprintf('  NO DDD VARIABLES WITH SUFFICIENT DATA (n>=30)\n\n');
end

%% -----------------------------------------------------------------------
%  CORRECTED: FREQUENCY MEDICATION VARIABLES (ORDINAL 0/1/2)
%  -----------------------------------------------------------------------
fprintf('========================================================\n');
fprintf('ANALYZING FREQUENCY MEDICATION VARIABLES (0/1/2 ORDINAL)\n');
fprintf('========================================================\n\n');

available_med_freq = {};
med_freq_corr_26 = [];
med_freq_corr_27 = [];
med_freq_corr_bvftd = [];

fprintf('CHECKING AVAILABILITY IN PATIENT SUBSAMPLE:\n');
for i = 1:length(medication_frequency_vars)
    if ismember(medication_frequency_vars{i}, patient_data.Properties.VariableNames)
        var_data = patient_data.(medication_frequency_vars{i});
        
        n_valid = sum(~isnan(var_data));
        pct_valid = 100*n_valid/n_patients;
        
        % Check if variable has all three values (0, 1, 2)
        has_all_values = sum(var_data == 0) > 0 && sum(var_data == 1) > 0 && sum(var_data == 2) > 0;
        
        fprintf('  %-25s: %d/%d valid (%.1f%%)', medication_frequency_vars{i}, ...
            n_valid, n_patients, pct_valid);
        
        if n_valid >= 30
            available_med_freq{end+1} = medication_frequency_vars{i};
            if has_all_values
                fprintf(' ? INCLUDE (full 0/1/2)\n');
            else
                fprintf(' ? INCLUDE (limited values)\n');
            end
        else
            fprintf(' ? EXCLUDE (too sparse)\n');
        end
    else
        fprintf('  %-25s: NOT FOUND\n', medication_frequency_vars{i});
    end
end
fprintf('\n');

if ~isempty(available_med_freq)
    fprintf('CORRELATIONS WITH DECISION SCORES (ORDINAL TREATMENT):\n');
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    for i = 1:length(available_med_freq)
        var_data = patient_data.(available_med_freq{i});
        
        % Transition-26
        valid_idx = ~isnan(var_data) & ~isnan(patient_data.Transition_26);
        if sum(valid_idx) >= 30
            % Treat as ORDINAL (0 < 1 < 2)
            [r, p] = corr(var_data(valid_idx), patient_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            med_freq_corr_26 = [med_freq_corr_26; r, p, n, ci(1), ci(2)];
            
            fprintf('  %-25s %7.3f %7.4f %5d', available_med_freq{i}, r, p, n);
            if p < 0.05
                fprintf(' * [Trans-26]\n');
            else
                fprintf('   [Trans-26]\n');
            end
        else
            med_freq_corr_26 = [med_freq_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
        end

        % bvFTD
        valid_idx = ~isnan(var_data) & ~isnan(patient_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), patient_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            med_freq_corr_bvftd = [med_freq_corr_bvftd; r, p, n, ci(1), ci(2)];
        else
            med_freq_corr_bvftd = [med_freq_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end
    fprintf('\n');
    
    results_medication.freq_correlations_26 = med_freq_corr_26;
    results_medication.freq_correlations_bvftd = med_freq_corr_bvftd;
else
    fprintf('  NO FREQUENCY VARIABLES WITH SUFFICIENT DATA (n>=30)\n\n');
end

%% -----------------------------------------------------------------------
%  CORRECTED: BINARY MEDICATION VARIABLES (0/1 YES/NO)
%  -----------------------------------------------------------------------
fprintf('========================================================\n');
fprintf('ANALYZING BINARY MEDICATION VARIABLES (0/1 YES/NO)\n');
fprintf('========================================================\n\n');

available_med_binary = {};
med_binary_corr_26 = [];
med_binary_corr_27 = [];
med_binary_corr_bvftd = [];

fprintf('CHECKING AVAILABILITY IN PATIENT SUBSAMPLE:\n');
for i = 1:length(medication_binary_vars)
    if ismember(medication_binary_vars{i}, patient_data.Properties.VariableNames)
        var_data = patient_data.(medication_binary_vars{i});
        
        n_valid = sum(~isnan(var_data));
        pct_valid = 100*n_valid/n_patients;
        
        if isnumeric(var_data)
            n_yes = sum(var_data == 1);
            n_no = sum(var_data == 0);
            fprintf('  %-25s: %d/%d valid (%.1f%%), Yes=%d, No=%d', ...
                medication_binary_vars{i}, n_valid, n_patients, pct_valid, n_yes, n_no);
        else
            fprintf('  %-25s: %d/%d valid (%.1f%%)', ...
                medication_binary_vars{i}, n_valid, n_patients, pct_valid);
        end
        
        if n_valid >= 30
            available_med_binary{end+1} = medication_binary_vars{i};
            fprintf(' ? INCLUDE\n');
        else
            fprintf(' ? EXCLUDE (too sparse)\n');
        end
    else
        fprintf('  %-25s: NOT FOUND\n', medication_binary_vars{i});
    end
end
fprintf('\n');

if ~isempty(available_med_binary)
    fprintf('CORRELATIONS WITH DECISION SCORES (PATIENTS ONLY):\n');
    fprintf('  Variable                     r        p      n\n');
    fprintf('  ----------------------------------------\n');
    
    for i = 1:length(available_med_binary)
        var_data = patient_data.(available_med_binary{i});
        
        % Transition-26
        valid_idx = ~isnan(var_data) & ~isnan(patient_data.Transition_26);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), patient_data.Transition_26(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            med_binary_corr_26 = [med_binary_corr_26; r, p, n, ci(1), ci(2)];
            
            fprintf('  %-25s %7.3f %7.4f %5d', available_med_binary{i}, r, p, n);
            if p < 0.05
                fprintf(' * [Trans-26]\n');
            else
                fprintf('   [Trans-26]\n');
            end
        else
            med_binary_corr_26 = [med_binary_corr_26; NaN, NaN, sum(valid_idx), NaN, NaN];
        end

        % bvFTD
        valid_idx = ~isnan(var_data) & ~isnan(patient_data.bvFTD);
        if sum(valid_idx) >= 30
            [r, p] = corr(var_data(valid_idx), patient_data.bvFTD(valid_idx));
            n = sum(valid_idx);
            ci = ci_r(r, n);
            med_binary_corr_bvftd = [med_binary_corr_bvftd; r, p, n, ci(1), ci(2)];
        else
            med_binary_corr_bvftd = [med_binary_corr_bvftd; NaN, NaN, sum(valid_idx), NaN, NaN];
        end
    end
    fprintf('\n');
    
    results_medication.binary_correlations_26 = med_binary_corr_26;
    results_medication.binary_correlations_bvftd = med_binary_corr_bvftd;
else
    fprintf('  NO BINARY MEDICATION VARIABLES WITH SUFFICIENT DATA (n>=30)\n\n');
end

%% -----------------------------------------------------------------------
%  MEDICATION AVAILABILITY BY DIAGNOSIS GROUP
%  -----------------------------------------------------------------------
fprintf('========================================================\n');
fprintf('MEDICATION AVAILABILITY BY DIAGNOSIS GROUP\n');
fprintf('========================================================\n\n');

% Get non-HC diagnosis groups only
patient_diag_groups = unique(patient_data.diagnosis_group);
patient_diag_groups = patient_diag_groups(~strcmp(patient_diag_groups, 'HC') & ~strcmp(patient_diag_groups, ''));

fprintf('FREQUENCY MEDICATION VARIABLES (_fr: 0/1/2):\n\n');
for i = 1:length(medication_frequency_vars)
    if ismember(medication_frequency_vars{i}, patient_data.Properties.VariableNames)
        fprintf('  %s:\n', medication_frequency_vars{i});
        
        for j = 1:length(patient_diag_groups)
            diag_mask = strcmp(patient_data.diagnosis_group, patient_diag_groups{j});
            var_data = patient_data.(medication_frequency_vars{i})(diag_mask);
            
            n_total = sum(diag_mask);
            n_available = sum(~isnan(var_data));
            pct_available = 100 * n_available / n_total;
            
            n_no = sum(var_data == 0);
            n_infreq = sum(var_data == 1);
            n_freq = sum(var_data == 2);
            
            fprintf('    %s: %d/%d available (%.1f%%) | No=%d, Infreq=%d, Freq=%d\n', ...
                patient_diag_groups{j}, n_available, n_total, pct_available, ...
                n_no, n_infreq, n_freq);
        end
        fprintf('\n');
    end
end

fprintf('BINARY MEDICATION VARIABLES (0/1):\n\n');
for i = 1:length(medication_binary_vars)
    if ismember(medication_binary_vars{i}, patient_data.Properties.VariableNames)
        fprintf('  %s:\n', medication_binary_vars{i});
        
        for j = 1:length(patient_diag_groups)
            diag_mask = strcmp(patient_data.diagnosis_group, patient_diag_groups{j});
            var_data = patient_data.(medication_binary_vars{i})(diag_mask);
            
            n_total = sum(diag_mask);
            n_available = sum(~isnan(var_data));
            pct_available = 100 * n_available / n_total;
            
            n_yes = sum(var_data == 1);
            n_no = sum(var_data == 0);
            pct_yes = 100 * n_yes / n_total;
            
            fprintf('    %s: %d/%d available (%.1f%%) | Yes=%d (%.1f%%), No=%d (%.1f%%)\n', ...
                patient_diag_groups{j}, n_available, n_total, pct_available, ...
                n_yes, pct_yes, n_no, 100*n_no/n_total);
        end
        fprintf('\n');
    end
end

fprintf('CONTINUOUS MEDICATION VARIABLES (DDD):\n\n');
for i = 1:length(medication_ddd_vars)
    if ismember(medication_ddd_vars{i}, patient_data.Properties.VariableNames)
        fprintf('  %s:\n', medication_ddd_vars{i});
        
        for j = 1:length(patient_diag_groups)
            diag_mask = strcmp(patient_data.diagnosis_group, patient_diag_groups{j});
            var_data = patient_data.(medication_ddd_vars{i})(diag_mask);
            
            n_total = sum(diag_mask);
            n_available = sum(~isnan(var_data));
            pct_available = 100 * n_available / n_total;
            
            fprintf('    %s: %d/%d available (%.1f%%)', ...
                patient_diag_groups{j}, n_available, n_total, pct_available);
            
            if n_available > 0
                fprintf(' | Mean=%.2f, SD=%.2f\n', mean(var_data, 'omitnan'), std(var_data, 'omitnan'));
            else
                fprintf('\n');
            end
        end
        fprintf('\n');
    end
end

fprintf('BENDEP SCALES (BENZODIAZEPINE DEPENDENCE):\n\n');
for i = 1:length(bendep_vars)
    if ismember(bendep_vars{i}, patient_data.Properties.VariableNames)
        fprintf('  %s:\n', bendep_vars{i});
        
        for j = 1:length(patient_diag_groups)
            diag_mask = strcmp(patient_data.diagnosis_group, patient_diag_groups{j});
            var_data = patient_data.(bendep_vars{i})(diag_mask);
            
            n_total = sum(diag_mask);
            n_available = sum(~isnan(var_data));
            pct_available = 100 * n_available / n_total;
            
            fprintf('    %s: %d/%d available (%.1f%%)', ...
                patient_diag_groups{j}, n_available, n_total, pct_available);
            
            if n_available > 0
                fprintf(' | Mean=%.2f, SD=%.2f\n', mean(var_data, 'omitnan'), std(var_data, 'omitnan'));
            else
                fprintf('\n');
            end
        end
        fprintf('\n');
    end
end

fprintf('========================================================\n\n');

%% -----------------------------------------------------------------------
%  SAVE MEDICATION RESULTS
%  -----------------------------------------------------------------------
fprintf('========================================================\n');
fprintf('SAVING MEDICATION ANALYSIS RESULTS\n');
fprintf('========================================================\n\n');

% Save combined medication summary
all_med_vars = [available_med_ddd, available_med_freq, available_med_binary];
if ~isempty(all_med_vars)
    all_med_corr_26 = [med_ddd_corr_26; med_freq_corr_26; med_binary_corr_26];
    all_med_corr_27 = [med_ddd_corr_27; med_freq_corr_27; med_binary_corr_27];
    all_med_corr_bvftd = [med_ddd_corr_bvftd; med_freq_corr_bvftd; med_binary_corr_bvftd];
    
    medication_summary = table();
    medication_summary.Variable = all_med_vars';
    medication_summary.Type = [repmat({'DDD/Continuous'}, length(available_med_ddd), 1); ...
                              repmat({'Frequency_0_1_2'}, length(available_med_freq), 1); ...
                              repmat({'Binary_0_1'}, length(available_med_binary), 1)];
    medication_summary.Transition_26_r = all_med_corr_26(:,1);
    medication_summary.Transition_26_p = all_med_corr_26(:,2);
    medication_summary.Transition_26_Uncorrected_significant = all_med_corr_26(:,2) < 0.05;
    medication_summary.Transition_26_n = all_med_corr_26(:,3);
    medication_summary.Transition_26_CI_lower = all_med_corr_26(:,4);
    medication_summary.Transition_26_CI_upper = all_med_corr_26(:,5);
    medication_summary.bvFTD_r = all_med_corr_bvftd(:,1);
    medication_summary.bvFTD_p = all_med_corr_bvftd(:,2);
    medication_summary.bvFTD_Uncorrected_significant = all_med_corr_bvftd(:,2) < 0.05;
    medication_summary.bvFTD_n = all_med_corr_bvftd(:,3);
    medication_summary.bvFTD_CI_lower = all_med_corr_bvftd(:,4);
    medication_summary.bvFTD_CI_upper = all_med_corr_bvftd(:,5);
    
    writetable(medication_summary, [data_out_path 'Summary_Medication_Correlations_PatientsOnly_CORRECTED.csv']);
    fprintf('  Saved: Summary_Medication_Correlations_PatientsOnly_CORRECTED.csv\n');
    
    % Create visualization if there are significant results
    sig_idx_med = all_med_corr_26(:,2) < 0.05;
    if sum(sig_idx_med) > 0
        fprintf('  Creating visualization of significant medication associations...\n');
        
        figure('Position', [100 100 800 600]);
        
        sig_med_vars = all_med_vars(sig_idx_med);
        sig_med_r = all_med_corr_26(sig_idx_med, 1);
        sig_med_ci_lower = all_med_corr_26(sig_idx_med, 4);
        sig_med_ci_upper = all_med_corr_26(sig_idx_med, 5);
        
        % Apply interpretable labels
        sig_med_labels = cellfun(@(x) get_label(x), sig_med_vars, 'UniformOutput', false);
        
        n_sig_med = length(sig_med_labels);
        y_pos_med = 1:n_sig_med;
        
        for i = 1:n_sig_med
            plot([sig_med_ci_lower(i), sig_med_ci_upper(i)], [y_pos_med(i), y_pos_med(i)], ...
                'k-', 'LineWidth', 1.5);
            hold on;
        end
        
        scatter(sig_med_r, y_pos_med, 100, 'filled', 'MarkerFaceColor', [0.2 0.6 0.2]);
        plot([0 0], [0.5, n_sig_med+0.5], 'r--', 'LineWidth', 2);
        
        set(gca, 'YTick', y_pos_med, 'YTickLabel', sig_med_labels, 'FontSize', 10);
        xlabel('Correlation (r) with 95% CI', 'FontWeight', 'bold', 'FontSize', 12);
        title('Significant Medication Associations with Transition-26 (Patients Only, p<0.05)', ...
            'FontWeight', 'bold', 'FontSize', 13);
        xlim([min([sig_med_ci_lower; -0.1])-0.1, max([sig_med_ci_upper; 0.1])+0.1]);
        ylim([0.5, n_sig_med+0.5]);
        grid on;
        
        saveas(gcf, [fig_path 'Fig_Medication_Significant_Associations_PatientsOnly_CORRECTED.png']);
        % saveas(gcf, [fig_path 'Fig_Medication_Significant_Associations_PatientsOnly_CORRECTED.fig']);  % Disabled: .fig fails in remote environments
        fprintf('  Saved: Fig_Medication_Significant_Associations_PatientsOnly_CORRECTED.png/.fig\n');
    else
        fprintf('  No significant medication associations found (p<0.05)\n');
    end
end

% =========================================================================
% CURATED MEDICATION VARIABLES SUMMARY
% =========================================================================
% Rationale: Section 0B defined curated medication variables to avoid
% redundancy between binary and frequency coding of the same medications.
% This summary extracts the curated subset for valid FDR correction.
% =========================================================================
if exist('curated_medication_vars', 'var') && ~isempty(curated_medication_vars)
    fprintf('\n========================================================\n');
    fprintf('CURATED MEDICATION VARIABLES (NON-REDUNDANT SET)\n');
    fprintf('========================================================\n\n');

    fprintf('  Curated set: %s\n\n', strjoin(curated_medication_vars, ', '));

    % Extract correlations for curated variables only
    curated_med_corr_26 = [];
    curated_med_corr_bvftd = [];
    curated_med_available = {};

    for i = 1:length(curated_medication_vars)
        cur_var = curated_medication_vars{i};

        % Check if in available medication vars
        idx_in_all = find(strcmp(all_med_vars, cur_var));
        if ~isempty(idx_in_all)
            curated_med_corr_26 = [curated_med_corr_26; all_med_corr_26(idx_in_all, :)];
            curated_med_corr_bvftd = [curated_med_corr_bvftd; all_med_corr_bvftd(idx_in_all, :)];
            curated_med_available{end+1} = cur_var;
        elseif strcmp(cur_var, 'n_psychotropic_classes') && ...
               ismember('n_psychotropic_classes', analysis_data.Properties.VariableNames)
            % Handle polypharmacy count separately if not in all_med_vars
            poly_data = analysis_data.n_psychotropic_classes;
            valid_26 = ~isnan(poly_data) & ~isnan(analysis_data.Transition_26);
            valid_bv = ~isnan(poly_data) & ~isnan(analysis_data.bvFTD);

            if sum(valid_26) >= MIN_SAMPLE_SIZE
                [r, p] = corr(poly_data(valid_26), analysis_data.Transition_26(valid_26));
                n = sum(valid_26);
                ci = ci_r(r, n);
                curated_med_corr_26 = [curated_med_corr_26; r, p, n, ci(1), ci(2)];
            else
                curated_med_corr_26 = [curated_med_corr_26; NaN, NaN, sum(valid_26), NaN, NaN];
            end

            if sum(valid_bv) >= MIN_SAMPLE_SIZE
                [r, p] = corr(poly_data(valid_bv), analysis_data.bvFTD(valid_bv));
                n = sum(valid_bv);
                ci = ci_r(r, n);
                curated_med_corr_bvftd = [curated_med_corr_bvftd; r, p, n, ci(1), ci(2)];
            else
                curated_med_corr_bvftd = [curated_med_corr_bvftd; NaN, NaN, sum(valid_bv), NaN, NaN];
            end
            curated_med_available{end+1} = cur_var;
        end
    end

    if ~isempty(curated_med_available)
        fprintf('  Available curated medication variables: %d/%d\n', ...
            length(curated_med_available), length(curated_medication_vars));
        fprintf('  Variables: %s\n\n', strjoin(curated_med_available, ', '));

        % Save curated medication summary
        curated_med_summary = table();
        curated_med_summary.Variable = curated_med_available';
        curated_med_summary.Transition_26_r = curated_med_corr_26(:,1);
        curated_med_summary.Transition_26_p = curated_med_corr_26(:,2);
        curated_med_summary.Transition_26_n = curated_med_corr_26(:,3);
        curated_med_summary.bvFTD_r = curated_med_corr_bvftd(:,1);
        curated_med_summary.bvFTD_p = curated_med_corr_bvftd(:,2);
        curated_med_summary.bvFTD_n = curated_med_corr_bvftd(:,3);

        writetable(curated_med_summary, [data_out_path 'Summary_Medication_CURATED.csv']);
        fprintf('  Saved: Summary_Medication_CURATED.csv\n');
    end
end

fprintf('\nMEDICATION ANALYSIS COMPLETE (CORRECTED CODING)\n\n');

%% ==========================================================================
%  SECTION 9C: RECENCY STRATIFIED ANALYSIS (OPTION 6 - NEW!)
%  ==========================================================================
fprintf('===========================================================\n');
fprintf('|  SECTION 9C: RECENCY STRATIFIED ANALYSIS (OPTION 6)     |\n');
fprintf('===========================================================\n\n');

fprintf('OBJECTIVE: Test if symptom-brain associations differ by illness recency\n');
fprintf('         (Recent/Active vs. Remitted depression)\n\n');

results_recency = struct();

if ~isempty(available_recency_vars)
    fprintf('========================================================\n');
    fprintf('STEP 1: CREATE RECENCY STRATIFICATION VARIABLE\n');
    fprintf('========================================================\n\n');
    
    % Combine all 3 recency variables to get "most recent diagnosis"
    recency_data = NaN(height(analysis_data), 3);
    for i = 1:length(available_recency_vars)
        if ismember(available_recency_vars{i}, analysis_data.Properties.VariableNames)
            recency_data(:,i) = analysis_data.(available_recency_vars{i});
        end
    end
    
    % Get most recent diagnosis across all 3 variables (minimum recency code)
    % Lower codes = more recent (e.g., 1=past month, higher values = lifetime)
    most_recent_diagnosis = min(recency_data, [], 2, 'omitnan');
    
    % Create RECENT vs REMITTED groups
    % RECENT: codes 1-2 (past month or past 6 months)
    % REMITTED: codes 3+ (past year or lifetime but not recent)
    % Adjust these thresholds based on NESDA codebook if needed
    
    recency_group = cell(height(analysis_data), 1);
    recency_group(:) = {''};
    
    recent_idx = most_recent_diagnosis <= 2;  % Adjust threshold as needed
    remitted_idx = most_recent_diagnosis > 2;
    
    recency_group(recent_idx) = {'Recent'};
    recency_group(remitted_idx) = {'Remitted'};
    
    analysis_data.Recency_Group = recency_group;
    
    n_recent = sum(strcmp(recency_group, 'Recent'));
    n_remitted = sum(strcmp(recency_group, 'Remitted'));
    n_missing = sum(strcmp(recency_group, ''));
    
    fprintf('  Recency stratification created:\n');
    fprintf('    RECENT (active illness): %d subjects\n', n_recent);
    fprintf('    REMITTED (past illness): %d subjects\n', n_remitted);
    fprintf('    Missing/No data: %d subjects\n', n_missing);
    fprintf('    Threshold: Recency code <= 2 = Recent\n\n');
    
    if n_recent >= 30 && n_remitted >= 30
        fprintf('========================================================\n');
        fprintf('STEP 2: STRATIFIED CORRELATIONS WITH TRANSITION-26\n');
        fprintf('========================================================\n\n');
        
        % Key symptom variables to test
        symptom_vars_to_test = {};
        if ismember('aids', analysis_data.Properties.VariableNames)
            symptom_vars_to_test{end+1} = 'aids';  % IDS total
        end
        if ismember('aids_mood_cognition', analysis_data.Properties.VariableNames)
            symptom_vars_to_test{end+1} = 'aids_mood_cognition';  % Mood/Cognition subscale
        end
        if ismember('aidssev', analysis_data.Properties.VariableNames)
            symptom_vars_to_test{end+1} = 'aidssev';  % Depression severity
        end
        
        if ~isempty(symptom_vars_to_test)
            fprintf('Testing %d key symptom variables:\n', length(symptom_vars_to_test));
            fprintf('  %s\n\n', strjoin(symptom_vars_to_test, ', '));
            
            stratified_results = table();
            symptom_var_list = {};
            n_recent_list = [];
            r_recent_list = [];
            p_recent_list = [];
            ci_recent_lower_list = [];
            ci_recent_upper_list = [];
            
            n_remitted_list = [];
            r_remitted_list = [];
            p_remitted_list = [];
            ci_remitted_lower_list = [];
            ci_remitted_upper_list = [];
            
            fisher_z_list = [];
            fisher_p_list = [];
            
            for i = 1:length(symptom_vars_to_test)
                symptom_var = symptom_vars_to_test{i};
                symptom_data_var = analysis_data.(symptom_var);
                
                fprintf('--------------------------------------------------\n');
                fprintf('VARIABLE: %s\n', symptom_var);
                fprintf('--------------------------------------------------\n');
                
                % RECENT group
                recent_mask = strcmp(analysis_data.Recency_Group, 'Recent') & ...
                             ~isnan(symptom_data_var) & ~isnan(analysis_data.Transition_26);
                n_recent_valid = sum(recent_mask);
                
                if n_recent_valid >= 30
                    [r_recent, p_recent] = corr(symptom_data_var(recent_mask), ...
                        analysis_data.Transition_26(recent_mask));
                    ci_recent = ci_r(r_recent, n_recent_valid);
                    
                    fprintf('  RECENT group:\n');
                    fprintf('    n = %d\n', n_recent_valid);
                    fprintf('    r = %.3f [%.3f, %.3f]\n', r_recent, ci_recent(1), ci_recent(2));
                    fprintf('    p = %.4f', p_recent);
                    if p_recent < 0.05
                        fprintf(' ***\n');
                    else
                        fprintf('\n');
                    end
                else
                    r_recent = NaN;
                    p_recent = NaN;
                    ci_recent = [NaN, NaN];
                    fprintf('  RECENT group: n=%d (insufficient data)\n', n_recent_valid);
                end
                
                % REMITTED group
                remitted_mask = strcmp(analysis_data.Recency_Group, 'Remitted') & ...
                               ~isnan(symptom_data_var) & ~isnan(analysis_data.Transition_26);
                n_remitted_valid = sum(remitted_mask);
                
                if n_remitted_valid >= 30
                    [r_remitted, p_remitted] = corr(symptom_data_var(remitted_mask), ...
                        analysis_data.Transition_26(remitted_mask));
                    ci_remitted = ci_r(r_remitted, n_remitted_valid);
                    
                    fprintf('  REMITTED group:\n');
                    fprintf('    n = %d\n', n_remitted_valid);
                    fprintf('    r = %.3f [%.3f, %.3f]\n', r_remitted, ci_remitted(1), ci_remitted(2));
                    fprintf('    p = %.4f', p_remitted);
                    if p_remitted < 0.05
                        fprintf(' ***\n');
                    else
                        fprintf('\n');
                    end
                else
                    r_remitted = NaN;
                    p_remitted = NaN;
                    ci_remitted = [NaN, NaN];
                    fprintf('  REMITTED group: n=%d (insufficient data)\n', n_remitted_valid);
                end
                
                % Fisher's z-test to compare correlations
                if n_recent_valid >= 30 && n_remitted_valid >= 30 && ...
                   ~isnan(r_recent) && ~isnan(r_remitted)
                    
                    % Fisher's z transformation
                    z_recent = 0.5 * log((1 + r_recent) / (1 - r_recent));
                    z_remitted = 0.5 * log((1 + r_remitted) / (1 - r_remitted));
                    
                    % Standard error of difference
                    se_diff = sqrt((1/(n_recent_valid-3)) + (1/(n_remitted_valid-3)));
                    
                    % Z-test statistic
                    z_stat = (z_recent - z_remitted) / se_diff;
                    
                    % Two-tailed p-value
                    p_fisher = 2 * (1 - normcdf(abs(z_stat)));
                    
                    fprintf('\n  FISHER''S Z-TEST (comparing groups):\n');
                    fprintf('    Z = %.3f\n', z_stat);
                    fprintf('    p = %.4f', p_fisher);
                    if p_fisher < 0.05
                        fprintf(' *** GROUPS DIFFER SIGNIFICANTLY\n');
                    else
                        fprintf('\n');
                    end
                else
                    z_stat = NaN;
                    p_fisher = NaN;
                    fprintf('\n  FISHER''S Z-TEST: Cannot compute (insufficient data)\n');
                end
                
                fprintf('\n');
                
                % Store results
                symptom_var_list{end+1} = symptom_var;
                n_recent_list(end+1) = n_recent_valid;
                r_recent_list(end+1) = r_recent;
                p_recent_list(end+1) = p_recent;
                ci_recent_lower_list(end+1) = ci_recent(1);
                ci_recent_upper_list(end+1) = ci_recent(2);
                
                n_remitted_list(end+1) = n_remitted_valid;
                r_remitted_list(end+1) = r_remitted;
                p_remitted_list(end+1) = p_remitted;
                ci_remitted_lower_list(end+1) = ci_remitted(1);
                ci_remitted_upper_list(end+1) = ci_remitted(2);
                
                fisher_z_list(end+1) = z_stat;
                fisher_p_list(end+1) = p_fisher;
            end
            
            % Create summary table
            stratified_results.Variable = symptom_var_list';
            stratified_results.Recent_n = n_recent_list';
            stratified_results.Recent_r = r_recent_list';
            stratified_results.Recent_p = p_recent_list';
            stratified_results.Recent_CI_lower = ci_recent_lower_list';
            stratified_results.Recent_CI_upper = ci_recent_upper_list';
            stratified_results.Remitted_n = n_remitted_list';
            stratified_results.Remitted_r = r_remitted_list';
            stratified_results.Remitted_p = p_remitted_list';
            stratified_results.Remitted_CI_lower = ci_remitted_lower_list';
            stratified_results.Remitted_CI_upper = ci_remitted_upper_list';
            stratified_results.Fisher_Z = fisher_z_list';
            stratified_results.Fisher_p = fisher_p_list';
            
            writetable(stratified_results, [data_out_path 'Summary_Recency_Stratified_Analysis.csv']);
            fprintf('  Saved: Summary_Recency_Stratified_Analysis.csv\n\n');
            
            % Store in results structure
            results_recency.stratified_analysis = stratified_results;
            
            fprintf('========================================================\n');
            fprintf('STEP 3: CREATE VISUALIZATION\n');
            fprintf('========================================================\n\n');
            
            % Create forest plot showing Recent vs Remitted correlations
            figure('Position', [100 100 1000 500]);
            
            n_vars = length(symptom_var_list);
            y_pos = (1:n_vars) * 2;  % Space out the groups
            
            % Plot Recent group (blue)
            for i = 1:n_vars
                if ~isnan(r_recent_list(i))
                    plot([ci_recent_lower_list(i), ci_recent_upper_list(i)], ...
                        [y_pos(i)+0.2, y_pos(i)+0.2], 'b-', 'LineWidth', 2);
                    hold on;
                    scatter(r_recent_list(i), y_pos(i)+0.2, 100, 'filled', ...
                        'MarkerFaceColor', [0.2 0.4 0.8]);
                end
            end
            
            % Plot Remitted group (red)
            for i = 1:n_vars
                if ~isnan(r_remitted_list(i))
                    plot([ci_remitted_lower_list(i), ci_remitted_upper_list(i)], ...
                        [y_pos(i)-0.2, y_pos(i)-0.2], 'r-', 'LineWidth', 2);
                    hold on;
                    scatter(r_remitted_list(i), y_pos(i)-0.2, 100, 'filled', ...
                        'MarkerFaceColor', [0.8 0.2 0.2]);
                end
            end
            
            % Reference line at r=0
            plot([0 0], [0.5, max(y_pos)+0.5], 'k--', 'LineWidth', 1.5);
            
            % Apply interpretable labels
            symptom_labels_stratified = cellfun(@(x) get_label(x), symptom_var_list, 'UniformOutput', false);
            
            set(gca, 'YTick', y_pos, 'YTickLabel', symptom_labels_stratified, 'FontSize', 11);
            xlabel('Correlation (r) with Transition-26', 'FontWeight', 'bold', 'FontSize', 12);
            title('Symptom-Brain Associations by Illness Recency', 'FontWeight', 'bold', 'FontSize', 14);
            
            % Add legend
            legend({'Recent (Active)', '', 'Remitted (Past)'}, 'Location', 'best', 'FontSize', 11);
            
            xlim([min([ci_recent_lower_list, ci_remitted_lower_list, -0.2])-0.1, ...
                  max([ci_recent_upper_list, ci_remitted_upper_list, 0.2])+0.1]);
            ylim([0.5, max(y_pos)+0.5]);
            grid on;
            
            saveas(gcf, [fig_path 'Fig_9C_Recency_Stratified_Symptom_Correlations.png']);
            % saveas(gcf, [fig_path 'Fig_9C_Recency_Stratified_Symptom_Correlations.fig']);  % Disabled: .fig fails in remote environments
            fprintf('  Saved: Fig_9C_Recency_Stratified_Symptom_Correlations.png/.fig\n\n');
            
            fprintf('========================================================\n');
            fprintf('INTERPRETATION SUMMARY\n');
            fprintf('========================================================\n\n');
            
            sig_diff = sum(fisher_p_list < 0.05);
            if sig_diff > 0
                fprintf('  *** %d variable(s) show SIGNIFICANTLY DIFFERENT associations\n', sig_diff);
                fprintf('      between Recent and Remitted groups (Fisher p<0.05)\n\n');
                fprintf('  This suggests brain-symptom relationships differ by illness phase.\n');
                fprintf('  Potential interpretation:\n');
                fprintf('    - Recent/Active: Direct symptom-brain relationship\n');
                fprintf('    - Remitted: May reflect treatment history, compensation, or scar effects\n\n');
            else
                fprintf('  No significant differences found between Recent and Remitted groups.\n');
                fprintf('  Brain-symptom associations appear similar across illness phases.\n\n');
            end
            
        else
            fprintf('  WARNING: No key symptom variables (aids, aids_mood_cognition, aidssev) found\n');
            fprintf('           Cannot perform stratified analysis\n\n');
        end
        
    else
        fprintf('  WARNING: Insufficient subjects in recency groups\n');
        fprintf('           Need n>=30 in both Recent and Remitted groups\n');
        fprintf('           Skipping stratified analysis\n\n');
    end
    
else
    fprintf('  WARNING: No recency variables found in dataset\n');
    fprintf('           Cannot perform recency stratified analysis\n\n');
end

fprintf('SECTION 9C COMPLETE\n\n');

%% ==========================================================================
%  SECTION 10: PRIORITY 4.5 - COMPREHENSIVE SUMMARY (CURATED VARIABLES)
%  ==========================================================================
%  NOTE: This section now compiles results from CURATED (non-redundant)
%  variable sets, as defined in Section 0B and finalized in Section 5D.
%  This ensures valid FDR correction across all analyses.
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  PRIORITY 4.5: COMPREHENSIVE SUMMARY (CURATED)   |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('COMPILING ALL UNIVARIATE ASSOCIATIONS\n');
fprintf('  NOTE: Using curated (non-redundant) variable sets\n\n');

all_vars = {};
all_categories = {};
all_corr_26 = [];
all_corr_bvftd = [];

fprintf('  Compiling results from all analyses...\n');

if exist('r_bmi_26', 'var')
    all_vars{end+1} = 'BMI';
    all_categories{end+1} = 'Metabolic';
    n_bmi = sum(~isnan(analysis_data.abmi) & ~isnan(analysis_data.Transition_26));
    ci_26 = ci_r(r_bmi_26, n_bmi);
    ci_bvftd = ci_r(r_bmi_bvftd, n_bmi);
    all_corr_26 = [all_corr_26; r_bmi_26, p_bmi_26, n_bmi, ci_26(1), ci_26(2)];
    all_corr_bvftd = [all_corr_bvftd; r_bmi_bvftd, p_bmi_bvftd, n_bmi, ci_bvftd(1), ci_bvftd(2)];
end

if exist('symptom_corr_26', 'var')
    for i = 1:length(symptom_names_clean)
        if ~isnan(symptom_corr_26(i,1))
            all_vars{end+1} = symptom_names_clean{i};
            all_categories{end+1} = 'Symptom_Severity';
            all_corr_26 = [all_corr_26; symptom_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; symptom_corr_bvftd(i,:)];
        end
    end
end

if exist('age_onset_corr_26', 'var')
    for i = 1:length(available_age_onset_vars)
        if ~isnan(age_onset_corr_26(i,1))
            all_vars{end+1} = available_age_onset_vars{i};
            all_categories{end+1} = 'Age_of_Onset';
            all_corr_26 = [all_corr_26; age_onset_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; age_onset_corr_bvftd(i,:)];
        end
    end
end

if exist('duration_corr_26', 'var')
    for i = 1:length(illness_duration_vars)
        if i <= size(duration_corr_26, 1) && ~isnan(duration_corr_26(i,1))
            all_vars{end+1} = illness_duration_vars{i};
            all_categories{end+1} = 'Illness_Duration';
            all_corr_26 = [all_corr_26; duration_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; duration_corr_bvftd(i,:)];
        end
    end
end

% NEW: Add recency variables to comprehensive summary (OPTION 6)
if exist('recency_corr_26', 'var')
    for i = 1:length(available_recency_vars)
        if ~isnan(recency_corr_26(i,1))
            all_vars{end+1} = available_recency_vars{i};
            all_categories{end+1} = 'Recency';
            all_corr_26 = [all_corr_26; recency_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; recency_corr_bvftd(i,:)];
        end
    end
end

if exist('clinical_corr_26', 'var')
    for i = 1:length(available_clinical_history_vars)
        if ~isnan(clinical_corr_26(i,1))
            all_vars{end+1} = available_clinical_history_vars{i};
            all_categories{end+1} = 'Clinical_History';
            all_corr_26 = [all_corr_26; clinical_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; clinical_corr_bvftd(i,:)];
        end
    end
end

if exist('childhood_corr_26', 'var')
    for i = 1:length(available_childhood_vars)
        if ~isnan(childhood_corr_26(i,1))
            all_vars{end+1} = available_childhood_vars{i};
            all_categories{end+1} = 'Childhood_Adversity';
            all_corr_26 = [all_corr_26; childhood_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; childhood_corr_bvftd(i,:)];
        end
    end
end

if exist('cognitive_corr_26', 'var') && ~isempty(cognitive_corr_26)
    for i = 1:length(cognitive_vars_found)
        if ~isnan(cognitive_corr_26(i,1))
            all_vars{end+1} = cognitive_vars_found{i};
            all_categories{end+1} = 'Cognition_Functioning';
            all_corr_26 = [all_corr_26; cognitive_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; cognitive_corr_bvftd(i,:)];
        end
    end
end

if exist('demo_corr_26', 'var') && ~isempty(demo_corr_26)
    for i = 1:length(demo_vars_analyzed)
        if ~isnan(demo_corr_26(i,1))
            all_vars{end+1} = demo_vars_analyzed{i};
            all_categories{end+1} = 'Demographics';
            all_corr_26 = [all_corr_26; demo_corr_26(i,:)];
            all_corr_bvftd = [all_corr_bvftd; demo_corr_bvftd(i,:)];
        end
    end
end

% NOTE: Medication variables are NOT added to comprehensive summary
% because they were analyzed on PATIENTS ONLY subsample
% They have their own separate summary table

fprintf('  Total variables analyzed: %d\n\n', length(all_vars));

comprehensive_summary = table();
comprehensive_summary.Variable = all_vars';
comprehensive_summary.Category = all_categories';

comprehensive_summary.Trans26_r = all_corr_26(:,1);
comprehensive_summary.Trans26_p = all_corr_26(:,2);
comprehensive_summary.Trans26_n = all_corr_26(:,3);
comprehensive_summary.Trans26_CI_lower = all_corr_26(:,4);
comprehensive_summary.Trans26_CI_upper = all_corr_26(:,5);

comprehensive_summary.bvFTD_r = all_corr_bvftd(:,1);
comprehensive_summary.bvFTD_p = all_corr_bvftd(:,2);
comprehensive_summary.bvFTD_n = all_corr_bvftd(:,3);
comprehensive_summary.bvFTD_CI_lower = all_corr_bvftd(:,4);
comprehensive_summary.bvFTD_CI_upper = all_corr_bvftd(:,5);

writetable(comprehensive_summary, [data_out_path 'COMPREHENSIVE_SUMMARY_All_Associations_CURATED.csv']);
fprintf('  Saved: COMPREHENSIVE_SUMMARY_All_Associations_CURATED.csv\n\n');

fprintf('SUMMARY BY CATEGORY:\n\n');
unique_categories = unique(all_categories);

for i = 1:length(unique_categories)
    cat_idx = strcmp(all_categories, unique_categories{i});

    fprintf('  %s (%d variables):\n', unique_categories{i}, sum(cat_idx));

    sig_26 = sum(all_corr_26(cat_idx, 2) < 0.05);
    sig_bvftd = sum(all_corr_bvftd(cat_idx, 2) < 0.05);

    fprintf('    Significant associations (p<0.05):\n');
    fprintf('      Transition-26: %d/%d (%.1f%%)\n', sig_26, sum(cat_idx), 100*sig_26/sum(cat_idx));
    fprintf('      bvFTD: %d/%d (%.1f%%)\n', sig_bvftd, sum(cat_idx), 100*sig_bvftd/sum(cat_idx));

    fprintf('    Mean |r| (Transition-26): %.3f\n', mean(abs(all_corr_26(cat_idx, 1))));
    fprintf('    Mean |r| (bvFTD): %.3f\n\n', mean(abs(all_corr_bvftd(cat_idx, 1))));
end

%% ==========================================================================
%  FOREST PLOT: TRANSITION-26 SIGNIFICANT ASSOCIATIONS
%  ==========================================================================
sig_idx = all_corr_26(:,2) < 0.05;
if sum(sig_idx) > 0
    fprintf('CREATING FOREST PLOT OF SIGNIFICANT ASSOCIATIONS (TRANSITION-26)\n');
    
    figure('Position', [100 100 1000 600]);
    
    sig_vars = all_vars(sig_idx);
    sig_r = all_corr_26(sig_idx, 1);
    sig_ci_lower = all_corr_26(sig_idx, 4);
    sig_ci_upper = all_corr_26(sig_idx, 5);
    
    % Apply interpretable labels for forest plot
    sig_labels = cellfun(@(x) get_label(x), sig_vars, 'UniformOutput', false);
    
    n_sig = length(sig_labels);
    y_pos = 1:n_sig;
    
    for i = 1:n_sig
        plot([sig_ci_lower(i), sig_ci_upper(i)], [y_pos(i), y_pos(i)], 'k-', 'LineWidth', 1.5);
        hold on;
    end
    
    scatter(sig_r, y_pos, 100, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8]);
    
    plot([0 0], [0.5, n_sig+0.5], 'r--', 'LineWidth', 2);
    
    set(gca, 'YTick', y_pos, 'YTickLabel', sig_labels, 'FontSize', 9);
    xlabel('Correlation (r) with 95% CI', 'FontWeight', 'bold', 'FontSize', 12);
    title('Significant Associations with Transition-26 (p<0.05)', ...
        'FontWeight', 'bold', 'FontSize', 14);
    xlim([min([sig_ci_lower; -0.1])-0.1, max([sig_ci_upper; 0.1])+0.1]);
    ylim([0.5, n_sig+0.5]);
    grid on;
    
    saveas(gcf, [fig_path 'Fig_4_5_Forest_Plot_Significant_Associations_Transition26.png']);
    % saveas(gcf, [fig_path 'Fig_4_5_Forest_Plot_Significant_Associations_Transition26.fig']);  % Disabled: .fig fails in remote environments
    fprintf('  Saved: Fig_4_5_Forest_Plot_Significant_Associations_Transition26.png/.fig\n\n');
end

%% ==========================================================================
%  FOREST PLOT: bvFTD SIGNIFICANT ASSOCIATIONS (NEW!)
%  ==========================================================================
sig_idx_bvftd = all_corr_bvftd(:,2) < 0.05;
if sum(sig_idx_bvftd) > 0
    fprintf('CREATING FOREST PLOT OF SIGNIFICANT ASSOCIATIONS (bvFTD)\n');
    
    figure('Position', [100 100 1000 600]);
    
    sig_vars_bvftd = all_vars(sig_idx_bvftd);
    sig_r_bvftd = all_corr_bvftd(sig_idx_bvftd, 1);
    sig_ci_lower_bvftd = all_corr_bvftd(sig_idx_bvftd, 4);
    sig_ci_upper_bvftd = all_corr_bvftd(sig_idx_bvftd, 5);
    
    % Apply interpretable labels for forest plot
    sig_labels_bvftd = cellfun(@(x) get_label(x), sig_vars_bvftd, 'UniformOutput', false);
    
    n_sig_bvftd = length(sig_labels_bvftd);
    y_pos_bvftd = 1:n_sig_bvftd;
    
    for i = 1:n_sig_bvftd
        plot([sig_ci_lower_bvftd(i), sig_ci_upper_bvftd(i)], [y_pos_bvftd(i), y_pos_bvftd(i)], 'k-', 'LineWidth', 1.5);
        hold on;
    end
    
    scatter(sig_r_bvftd, y_pos_bvftd, 100, 'filled', 'MarkerFaceColor', [0.8 0.2 0.2]);
    
    plot([0 0], [0.5, n_sig_bvftd+0.5], 'r--', 'LineWidth', 2);
    
    set(gca, 'YTick', y_pos_bvftd, 'YTickLabel', sig_labels_bvftd, 'FontSize', 9);
    xlabel('Correlation (r) with 95% CI', 'FontWeight', 'bold', 'FontSize', 12);
    title('Significant Associations with bvFTD (p<0.05)', ...
        'FontWeight', 'bold', 'FontSize', 14);
    xlim([min([sig_ci_lower_bvftd; -0.1])-0.1, max([sig_ci_upper_bvftd; 0.1])+0.1]);
    ylim([0.5, n_sig_bvftd+0.5]);
    grid on;
    
    saveas(gcf, [fig_path 'Fig_4_5_Forest_Plot_Significant_Associations_bvFTD.png']);
    % saveas(gcf, [fig_path 'Fig_4_5_Forest_Plot_Significant_Associations_bvFTD.fig']);  % Disabled: .fig fails in remote environments
    fprintf('  Saved: Fig_4_5_Forest_Plot_Significant_Associations_bvFTD.png/.fig\n\n');
else
    fprintf('NO SIGNIFICANT ASSOCIATIONS FOUND FOR bvFTD (p<0.05)\n');
    fprintf('  bvFTD forest plot was not created.\n\n');
end

fprintf('PRIORITY 4.5 COMPLETE\n\n');

%% ==========================================================================
%  SECTION 10B: SESSION 2 - FEATURE 2.1: UNIVARIATE CORRELATIONS EXPORT
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  FEATURE 2.1: UNIVARIATE CORRELATIONS EXPORT    |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('Collecting all univariate correlations from all sections...\n\n');

% Helper function to create correlation row
create_corr_row = @(varname, label, r, p, n, ci_low, ci_high, p_fdr, fdr_sig) ...
    struct('Variable', varname, 'Label', label, 'r', r, 'p_uncorrected', p, ...
           'n_subjects', n, 'CI_lower', ci_low, 'CI_upper', ci_high, ...
           'p_FDR', p_fdr, 'FDR_significant', fdr_sig);

% Initialize collectors for each decision score
corr_data_26 = {};
corr_data_bvftd = {};

% ========== SYMPTOM SEVERITY ==========
if exist('symptom_names_clean', 'var') && exist('symptom_corr_26', 'var')
    for i = 1:length(symptom_names_clean)
        if ~isnan(symptom_corr_26(i, 1))
            varname = symptom_names_clean{i};
            label = get_label(varname);
            corr_data_26{end+1} = create_corr_row(varname, label, ...
                symptom_corr_26(i,1), symptom_corr_26(i,2), symptom_corr_26(i,3), ...
                symptom_corr_26(i,4), symptom_corr_26(i,5), adj_p_26(i), h_fdr_26(i));
            corr_data_bvftd{end+1} = create_corr_row(varname, label, ...
                symptom_corr_bvftd(i,1), symptom_corr_bvftd(i,2), symptom_corr_bvftd(i,3), ...
                symptom_corr_bvftd(i,4), symptom_corr_bvftd(i,5), adj_p_bvftd(i), h_fdr_bvftd(i));
        end
    end
end

% ========== CLINICAL HISTORY ==========
if exist('available_clinical_history_vars', 'var') && exist('clinical_corr_26', 'var')
    for i = 1:length(available_clinical_history_vars)
        if ~isnan(clinical_corr_26(i, 1))
            varname = available_clinical_history_vars{i};
            label = get_label(varname);
            corr_data_26{end+1} = create_corr_row(varname, label, ...
                clinical_corr_26(i,1), clinical_corr_26(i,2), clinical_corr_26(i,3), ...
                clinical_corr_26(i,4), clinical_corr_26(i,5), adj_p_clin_26(i), h_fdr_clin_26(i));
            corr_data_bvftd{end+1} = create_corr_row(varname, label, ...
                clinical_corr_bvftd(i,1), clinical_corr_bvftd(i,2), clinical_corr_bvftd(i,3), ...
                clinical_corr_bvftd(i,4), clinical_corr_bvftd(i,5), adj_p_clin_bvftd(i), h_fdr_clin_bvftd(i));
        end
    end
end

% ========== CHILDHOOD ADVERSITY ==========
if exist('available_childhood_vars', 'var') && exist('childhood_corr_26', 'var')
    for i = 1:length(available_childhood_vars)
        if ~isnan(childhood_corr_26(i, 1))
            varname = available_childhood_vars{i};
            label = get_label(varname);
            corr_data_26{end+1} = create_corr_row(varname, label, ...
                childhood_corr_26(i,1), childhood_corr_26(i,2), childhood_corr_26(i,3), ...
                childhood_corr_26(i,4), childhood_corr_26(i,5), adj_p_child_26(i), h_fdr_child_26(i));
            corr_data_bvftd{end+1} = create_corr_row(varname, label, ...
                childhood_corr_bvftd(i,1), childhood_corr_bvftd(i,2), childhood_corr_bvftd(i,3), ...
                childhood_corr_bvftd(i,4), childhood_corr_bvftd(i,5), adj_p_child_bvftd(i), h_fdr_child_bvftd(i));
        end
    end
end

% ========== DEMOGRAPHICS ==========
if exist('demo_vars_analyzed', 'var') && exist('demo_corr_26', 'var')
    for i = 1:length(demo_vars_analyzed)
        if ~isnan(demo_corr_26(i, 1))
            varname = demo_vars_analyzed{i};
            label = get_label(varname);
            corr_data_26{end+1} = create_corr_row(varname, label, ...
                demo_corr_26(i,1), demo_corr_26(i,2), demo_corr_26(i,3), ...
                demo_corr_26(i,4), demo_corr_26(i,5), adj_p_demo_26(i), h_fdr_demo_26(i));
            corr_data_bvftd{end+1} = create_corr_row(varname, label, ...
                demo_corr_bvftd(i,1), demo_corr_bvftd(i,2), demo_corr_bvftd(i,3), ...
                demo_corr_bvftd(i,4), demo_corr_bvftd(i,5), adj_p_demo_bvftd(i), h_fdr_demo_bvftd(i));
        end
    end
end

% Convert to tables and sort by p_uncorrected
if ~isempty(corr_data_26)
    univar_tbl_26 = struct2table(vertcat(corr_data_26{:}));
    univar_tbl_26 = sortrows(univar_tbl_26, 'p_uncorrected', 'ascend');
    univar_tbl_26.Decision_Score = repmat({'Transition-26'}, height(univar_tbl_26), 1);
    writetable(univar_tbl_26, [data_out_path 'Univariate_Correlations_Transition26_FDR_Sorted.csv']);
    fprintf('  ? Saved: Univariate_Correlations_Transition26_FDR_Sorted.csv (%d correlations)\n', height(univar_tbl_26));
end

if ~isempty(corr_data_bvftd)
    univar_tbl_bvftd = struct2table(vertcat(corr_data_bvftd{:}));
    univar_tbl_bvftd = sortrows(univar_tbl_bvftd, 'p_uncorrected', 'ascend');
    univar_tbl_bvftd.Decision_Score = repmat({'bvFTD'}, height(univar_tbl_bvftd), 1);
    writetable(univar_tbl_bvftd, [data_out_path 'Univariate_Correlations_bvFTD_FDR_Sorted.csv']);
    fprintf('  ? Saved: Univariate_Correlations_bvFTD_FDR_Sorted.csv (%d correlations)\n', height(univar_tbl_bvftd));
end

fprintf('\nFEATURE 2.1 COMPLETE: Univariate correlations exported\n\n');

%% ==========================================================================
%  SECTION 10C: SESSION 2 - FEATURE 2.3: COHORT-STRATIFIED BOXPLOTS
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  FEATURE 2.3: COHORT-STRATIFIED BOXPLOTS        |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('Creating decision score comparisons across diagnosis groups...\n\n');

% Check if diagnosis_group exists in full dataset (before patient filtering)
if exist('analysis_data_full', 'var') && ismember('diagnosis_group', analysis_data_full.Properties.VariableNames)
    cohort_data = analysis_data_full;
else
    fprintf('  WARNING: Using patient-only dataset (no HC for comparison)\n');
    cohort_data = analysis_data;
end

% Get unique diagnosis groups
unique_groups = unique(cohort_data.diagnosis_group);
unique_groups = unique_groups(~cellfun(@isempty, unique_groups));

fprintf('  Diagnosis groups found: %s\n', strjoin(unique_groups, ', '));

% Prepare data for each decision score
ds_names = {'Transition_26', 'bvFTD'};
ds_labels = {'Transition-26', 'bvFTD'};
ds_colors = {[0.2 0.4 0.8], [0.8 0.2 0.2]};

% Create figure
figure('Position', [100 100 1200 500]);

for ds_idx = 1:2
    ds_name = ds_names{ds_idx};
    ds_label = ds_labels{ds_idx};

    subplot(1, 2, ds_idx);
    hold on;

    % Collect data by group
    group_data = {};
    group_labels = {};
    group_n = [];

    for g = 1:length(unique_groups)
        group_name = unique_groups{g};
        group_mask = strcmp(cohort_data.diagnosis_group, group_name);
        group_scores = cohort_data.(ds_name)(group_mask);
        group_scores = group_scores(~isnan(group_scores));

        if ~isempty(group_scores)
            group_data{end+1} = group_scores;
            group_labels{end+1} = sprintf('%s (n=%d)', group_name, length(group_scores));
            group_n(end+1) = length(group_scores);
        end
    end

    % Create boxplot
    if length(group_data) >= 2
        % Combine all data for boxplot
        all_scores = [];
        all_groups = [];
        for g = 1:length(group_data)
            all_scores = [all_scores; group_data{g}];
            all_groups = [all_groups; repmat(g, length(group_data{g}), 1)];
        end

        % Boxplot
        boxplot(all_scores, all_groups, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);

        % Add individual data points with jitter
        for g = 1:length(group_data)
            x_jitter = g + 0.15 * (rand(length(group_data{g}), 1) - 0.5);
            scatter(x_jitter, group_data{g}, 30, ds_colors{ds_idx}, 'filled', ...
                'MarkerFaceAlpha', 0.3);
        end

        % Mark median and mean
        for g = 1:length(group_data)
            med_val = median(group_data{g});
            mean_val = mean(group_data{g});
            plot(g, med_val, 'r_', 'LineWidth', 3, 'MarkerSize', 15);
            plot(g, mean_val, 'gd', 'LineWidth', 2, 'MarkerSize', 8);
        end

        % One-way ANOVA if =3 groups with n=5
        valid_groups = group_n >= 5;
        if sum(valid_groups) >= 3
            % ANOVA
            [p_anova, tbl_anova, stats_anova] = anova1(all_scores, all_groups, 'off');

            % Title with ANOVA result
            if p_anova < 0.05
                title(sprintf('%s\nANOVA: F=%.2f, p=%.4f ***', ds_label, ...
                    tbl_anova{2,5}, p_anova), 'FontWeight', 'bold');

                % Post-hoc Tukey HSD
                [c, ~, ~, gnames] = multcompare(stats_anova, 'Display', 'off', 'CType', 'hsd');

                % Calculate Cohen's d for significant pairs
                fprintf('\n  %s - Significant pairwise comparisons (Tukey HSD):\n', ds_label);
                for i = 1:size(c, 1)
                    if c(i, 6) < 0.05  % Significant comparison
                        g1_idx = c(i, 1);
                        g2_idx = c(i, 2);

                        % Cohen's d
                        pooled_std = sqrt((var(group_data{g1_idx}) + var(group_data{g2_idx})) / 2);
                        cohens_d = (mean(group_data{g1_idx}) - mean(group_data{g2_idx})) / pooled_std;

                        fprintf('    %s vs %s: p=%.4f, d=%.3f\n', ...
                            unique_groups{g1_idx}, unique_groups{g2_idx}, c(i, 6), cohens_d);
                    end
                end
            else
                title(sprintf('%s\nANOVA: F=%.2f, p=%.4f', ds_label, ...
                    tbl_anova{2,5}, p_anova), 'FontWeight', 'bold');
            end
        else
            title(ds_label, 'FontWeight', 'bold');
        end

        % Labels and formatting
        set(gca, 'XTick', 1:length(group_labels), 'XTickLabel', group_labels, ...
            'XTickLabelRotation', 15);
        ylabel('Decision Score', 'FontWeight', 'bold');
        xlabel('Diagnosis Group', 'FontWeight', 'bold');
        grid on;

        % Legend for median/mean markers
        if ds_idx == 3
            legend({'Median', 'Mean'}, 'Location', 'northeast', 'FontSize', 8);
        end
    end

    hold off;
end

sgtitle('Decision Scores by Diagnosis Group', 'FontWeight', 'bold', 'FontSize', 14);

% Save figure
saveas(gcf, [fig_path 'Cohort_Stratified_Decision_Scores.png']);
% saveas(gcf, [fig_path 'Cohort_Stratified_Decision_Scores.fig']);  % Disabled: .fig fails in remote environments
fprintf('\n  ? Saved: Cohort_Stratified_Decision_Scores.png/.fig\n');

fprintf('\nFEATURE 2.3 COMPLETE: Cohort-stratified boxplots created\n\n');

%% ==========================================================================
%  SECTION 10D: SESSION 2 - FEATURE 2.4: AGE × DECISION SCORE INTERACTION
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  FEATURE 2.4: AGE × DECISION SCORE INTERACTION   |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('Investigating age moderation of brain-symptom relationships...\n\n');

% Use full dataset with all cohorts
if exist('analysis_data_full', 'var')
    interaction_data = analysis_data_full;
else
    interaction_data = analysis_data;
end

% Filter for valid age and diagnosis_group
valid_for_interaction = ~isnan(interaction_data.Age) & ...
                        ~cellfun(@isempty, interaction_data.diagnosis_group);
interaction_data = interaction_data(valid_for_interaction, :);

% Get unique groups
groups = unique(interaction_data.diagnosis_group);
groups = groups(~cellfun(@isempty, groups));

% Colors for each group
group_colors = containers.Map();
group_colors('HC') = [0.2 0.7 0.2];           % Green
group_colors('Depression') = [0.8 0.2 0.2];  % Red
group_colors('Anxiety') = [0.2 0.4 0.8];     % Blue
group_colors('Comorbid') = [0.8 0.4 0.8];    % Purple

% Process each decision score
ds_names = {'Transition_26', 'bvFTD'};
ds_labels = {'Transition-26', 'bvFTD'};
ds_filenames = {'Age_Interaction_Transition26', 'Age_Interaction_bvFTD'};

for ds_idx = 1:2
    ds_name = ds_names{ds_idx};
    ds_label = ds_labels{ds_idx};
    ds_filename = ds_filenames{ds_idx};

    fprintf('\nAnalyzing: %s\n', ds_label);

    % Filter for valid decision scores
    valid_ds = ~isnan(interaction_data.(ds_name));
    analysis_subset = interaction_data(valid_ds, :);

    if height(analysis_subset) < 20
        fprintf('  WARNING: Insufficient data (n=%d), skipping\n', height(analysis_subset));
        continue;
    end

    % Fit linear model: Decision_Score ~ Age * diagnosis_group
    % Convert diagnosis_group to categorical
    analysis_subset.diagnosis_group_cat = categorical(analysis_subset.diagnosis_group);

    % Fit model
    try
        mdl = fitlm(analysis_subset, sprintf('%s ~ Age * diagnosis_group_cat', ds_name));

        % Extract statistics
        r_squared = mdl.Rsquared.Ordinary;
        interaction_p = NaN;

        % Find interaction term p-value
        coeff_names = mdl.CoefficientNames;
        for i = 1:length(coeff_names)
            if contains(coeff_names{i}, 'Age:diagnosis_group_cat')
                interaction_p = mdl.Coefficients.pValue(i);
                break;
            end
        end

        fprintf('  Model R²: %.4f\n', r_squared);
        if ~isnan(interaction_p)
            fprintf('  Interaction p-value: %.4f', interaction_p);
            if interaction_p < 0.05
                fprintf(' ***\n');
            else
                fprintf('\n');
            end
        end

        % Create figure
        fig = figure('Position', [100 100 900 700]);
        hold on;

        % Age range for predictions
        age_range = linspace(min(analysis_subset.Age), max(analysis_subset.Age), 100)';

        % Plot for each group
        legend_entries = {};
        for g = 1:length(groups)
            group_name = groups{g};
            group_mask = strcmp(analysis_subset.diagnosis_group, group_name);

            if sum(group_mask) < 3
                continue;  % Skip groups with too few data points
            end

            % Get color
            if group_colors.isKey(group_name)
                color = group_colors(group_name);
            else
                color = [0.5 0.5 0.5];  % Default gray
            end

            % Scatter plot
            scatter(analysis_subset.Age(group_mask), analysis_subset.(ds_name)(group_mask), ...
                50, color, 'filled', 'MarkerFaceAlpha', 0.5);

            % Fit regression for this group
            group_ages = analysis_subset.Age(group_mask);
            group_scores = analysis_subset.(ds_name)(group_mask);

            % Simple linear regression with error structure for CI
            [p_fit, S] = polyfit(group_ages, group_scores, 1);
            [y_fit, delta] = polyval(p_fit, age_range, S);

            % Shaded 95% CI
            fill([age_range; flipud(age_range)], ...
                 [y_fit + 1.96*delta; flipud(y_fit - 1.96*delta)], ...
                 color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

            % Regression line
            plot(age_range, y_fit, '-', 'Color', color, 'LineWidth', 2);

            % Legend entry with regression equation
            [r_group, p_group] = corr(group_ages, group_scores);
            legend_entries{end+1} = sprintf('%s: r=%.3f, p=%.3f (y=%.3fx+%.2f)', ...
                group_name, r_group, p_group, p_fit(1), p_fit(2));
        end

        % Labels and formatting
        xlabel('Age (years)', 'FontWeight', 'bold', 'FontSize', 12);
        ylabel(sprintf('%s Score', ds_label), 'FontWeight', 'bold', 'FontSize', 12);
        title(sprintf('Age × %s Interaction', ds_label), 'FontWeight', 'bold', 'FontSize', 14);

        % Add statistics annotation
        if ~isnan(interaction_p)
            text(0.02, 0.98, sprintf('Interaction p=%.4f\nOverall R²=%.4f', interaction_p, r_squared), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        end

        % Legend
        legend(legend_entries, 'Location', 'best', 'FontSize', 9);
        grid on;
        hold off;

        % Validate and save figure
        if ishghandle(fig) && isvalid(fig)
            saveas(fig, [fig_path ds_filename '.png']);
            % saveas(fig, [fig_path ds_filename '.fig']);  % Disabled
            fprintf('  ? Saved: %s.png/.fig\n', ds_filename);
        else
            warning('Figure handle invalid, skipping save for %s', ds_filename);
        end

    catch ME
        fprintf('  ERROR fitting model: %s\n', ME.message);
    end
end

fprintf('\nFEATURE 2.4 COMPLETE: Age × Decision Score interaction analysis complete\n\n');

%% ==========================================================================
%  SECTION 10E: BINARY ANXIETY/DEPRESSION CODING & HC COMPARISONS (FULL)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  SECTION 10E: BINARY ANXIETY/DEPRESSION (FULL)  |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('CREATING BINARY DIAGNOSTIC INDICATORS (FULL IMPLEMENTATION)\n');
fprintf('  Rationale: Comorbid individuals coded as "yes" for BOTH anxiety AND depression\n');
fprintf('  Allows separate examination of each disorder while retaining comorbid cases\n');
fprintf('  WITH FDR CORRECTION across all comparisons\n\n');

% Use full dataset including HC
if exist('analysis_data_full', 'var') && ismember('diagnosis_group', analysis_data_full.Properties.VariableNames)
    binary_data = analysis_data_full;

    % Create binary indicators
    binary_data.anxiety_binary = zeros(height(binary_data), 1);
    binary_data.depression_binary = zeros(height(binary_data), 1);

    % Code anxiety (includes both "Anxiety" and "Comorbid")
    anxiety_idx = strcmp(binary_data.diagnosis_group, 'Anxiety') | ...
                  strcmp(binary_data.diagnosis_group, 'Comorbid');
    binary_data.anxiety_binary(anxiety_idx) = 1;

    % Code depression (includes both "Depression" and "Comorbid")
    depression_idx = strcmp(binary_data.diagnosis_group, 'Depression') | ...
                     strcmp(binary_data.diagnosis_group, 'Comorbid');
    binary_data.depression_binary(depression_idx) = 1;

    % Report coding
    fprintf('  Binary coding summary:\n');
    fprintf('    HC (reference): n=%d\n', sum(strcmp(binary_data.diagnosis_group, 'HC')));
    fprintf('    Anxiety=1 (Anxiety + Comorbid): n=%d\n', sum(binary_data.anxiety_binary == 1));
    fprintf('      Pure Anxiety: n=%d\n', sum(strcmp(binary_data.diagnosis_group, 'Anxiety')));
    fprintf('      Comorbid: n=%d\n', sum(strcmp(binary_data.diagnosis_group, 'Comorbid')));
    fprintf('    Depression=1 (Depression + Comorbid): n=%d\n', sum(binary_data.depression_binary == 1));
    fprintf('      Pure Depression: n=%d\n', sum(strcmp(binary_data.diagnosis_group, 'Depression')));
    fprintf('    Comorbid: n=%d (counted in both)\n\n', sum(strcmp(binary_data.diagnosis_group, 'Comorbid')));

    % =======================================================================
    % COMPREHENSIVE T-TESTS: ALL DECISION SCORES × BOTH DISORDERS
    % =======================================================================
    fprintf('COMPREHENSIVE GROUP COMPARISONS (with FDR correction)\n\n');

    % Collect all t-test results for FDR correction
    all_binary_pvals = [];
    binary_results = {};

    % Extract HC scores once
    hc_scores_26 = binary_data.Transition_26(strcmp(binary_data.diagnosis_group, 'HC'));
    hc_scores_26 = hc_scores_26(~isnan(hc_scores_26));

    hc_scores_bv = binary_data.bvFTD(strcmp(binary_data.diagnosis_group, 'HC'));
    hc_scores_bv = hc_scores_bv(~isnan(hc_scores_bv));

    % Test 1: Anxiety vs HC - Transition-26
    anx_scores_26 = binary_data.Transition_26(binary_data.anxiety_binary == 1);
    anx_scores_26 = anx_scores_26(~isnan(anx_scores_26));

    if length(hc_scores_26) >= 10 && length(anx_scores_26) >= 10
        [h, p, ci, stats] = ttest2(anx_scores_26, hc_scores_26);
        cohens_d = (mean(anx_scores_26) - mean(hc_scores_26)) / ...
                   sqrt(((length(anx_scores_26)-1)*var(anx_scores_26) + ...
                         (length(hc_scores_26)-1)*var(hc_scores_26)) / ...
                        (length(anx_scores_26) + length(hc_scores_26) - 2));

        all_binary_pvals(end+1) = p;
        binary_results{end+1} = struct('comparison', 'Anxiety_vs_HC', 'score', 'Transition-26', ...
                                       'n1', length(anx_scores_26), 'n2', length(hc_scores_26), ...
                                       'mean1', mean(anx_scores_26), 'sd1', std(anx_scores_26), ...
                                       'mean2', mean(hc_scores_26), 'sd2', std(hc_scores_26), ...
                                       'tstat', stats.tstat, 'df', stats.df, 'p', p, 'd', cohens_d);
    end

    % Test 2: Anxiety vs HC - bvFTD
    anx_scores_bv = binary_data.bvFTD(binary_data.anxiety_binary == 1);
    anx_scores_bv = anx_scores_bv(~isnan(anx_scores_bv));

    if length(hc_scores_bv) >= 10 && length(anx_scores_bv) >= 10
        [h, p, ci, stats] = ttest2(anx_scores_bv, hc_scores_bv);
        cohens_d = (mean(anx_scores_bv) - mean(hc_scores_bv)) / ...
                   sqrt(((length(anx_scores_bv)-1)*var(anx_scores_bv) + ...
                         (length(hc_scores_bv)-1)*var(hc_scores_bv)) / ...
                        (length(anx_scores_bv) + length(hc_scores_bv) - 2));

        all_binary_pvals(end+1) = p;
        binary_results{end+1} = struct('comparison', 'Anxiety_vs_HC', 'score', 'bvFTD', ...
                                       'n1', length(anx_scores_bv), 'n2', length(hc_scores_bv), ...
                                       'mean1', mean(anx_scores_bv), 'sd1', std(anx_scores_bv), ...
                                       'mean2', mean(hc_scores_bv), 'sd2', std(hc_scores_bv), ...
                                       'tstat', stats.tstat, 'df', stats.df, 'p', p, 'd', cohens_d);
    end

    % Test 3: Depression vs HC - Transition-26
    dep_scores_26 = binary_data.Transition_26(binary_data.depression_binary == 1);
    dep_scores_26 = dep_scores_26(~isnan(dep_scores_26));

    if length(hc_scores_26) >= 10 && length(dep_scores_26) >= 10
        [h, p, ci, stats] = ttest2(dep_scores_26, hc_scores_26);
        cohens_d = (mean(dep_scores_26) - mean(hc_scores_26)) / ...
                   sqrt(((length(dep_scores_26)-1)*var(dep_scores_26) + ...
                         (length(hc_scores_26)-1)*var(hc_scores_26)) / ...
                        (length(dep_scores_26) + length(hc_scores_26) - 2));

        all_binary_pvals(end+1) = p;
        binary_results{end+1} = struct('comparison', 'Depression_vs_HC', 'score', 'Transition-26', ...
                                       'n1', length(dep_scores_26), 'n2', length(hc_scores_26), ...
                                       'mean1', mean(dep_scores_26), 'sd1', std(dep_scores_26), ...
                                       'mean2', mean(hc_scores_26), 'sd2', std(hc_scores_26), ...
                                       'tstat', stats.tstat, 'df', stats.df, 'p', p, 'd', cohens_d);
    end

    % Test 4: Depression vs HC - bvFTD
    dep_scores_bv = binary_data.bvFTD(binary_data.depression_binary == 1);
    dep_scores_bv = dep_scores_bv(~isnan(dep_scores_bv));

    if length(hc_scores_bv) >= 10 && length(dep_scores_bv) >= 10
        [h, p, ci, stats] = ttest2(dep_scores_bv, hc_scores_bv);
        cohens_d = (mean(dep_scores_bv) - mean(hc_scores_bv)) / ...
                   sqrt(((length(dep_scores_bv)-1)*var(dep_scores_bv) + ...
                         (length(hc_scores_bv)-1)*var(hc_scores_bv)) / ...
                        (length(dep_scores_bv) + length(hc_scores_bv) - 2));

        all_binary_pvals(end+1) = p;
        binary_results{end+1} = struct('comparison', 'Depression_vs_HC', 'score', 'bvFTD', ...
                                       'n1', length(dep_scores_bv), 'n2', length(hc_scores_bv), ...
                                       'mean1', mean(dep_scores_bv), 'sd1', std(dep_scores_bv), ...
                                       'mean2', mean(hc_scores_bv), 'sd2', std(hc_scores_bv), ...
                                       'tstat', stats.tstat, 'df', stats.df, 'p', p, 'd', cohens_d);
    end

    % Apply FDR correction across all binary comparisons
    if ~isempty(all_binary_pvals)
        [h_fdr_binary, crit_p, adj_p_binary] = fdr_bh(all_binary_pvals, 0.05);

        fprintf('  FDR CORRECTION APPLIED: %d comparisons, q=0.05\n\n', length(all_binary_pvals));

        % Print results with FDR correction
        for i = 1:length(binary_results)
            res = binary_results{i};
            fprintf('  %s (%s):\n', res.comparison, res.score);
            fprintf('    Group: M=%.3f, SD=%.3f, n=%d\n', res.mean1, res.sd1, res.n1);
            fprintf('    HC: M=%.3f, SD=%.3f, n=%d\n', res.mean2, res.sd2, res.n2);
            fprintf('    t(%d)=%.3f, p_uncorr=%.4f', res.df, res.tstat, res.p);
            if h_fdr_binary(i)
                fprintf(', p_FDR=%.4f ***\n', adj_p_binary(i));
            else
                fprintf(', p_FDR=%.4f\n', adj_p_binary(i));
            end
            fprintf('    Cohen''s d=%.3f', res.d);
            if abs(res.d) > 0.8
                fprintf(' (large effect)\n\n');
            elseif abs(res.d) > 0.5
                fprintf(' (medium effect)\n\n');
            elseif abs(res.d) > 0.2
                fprintf(' (small effect)\n\n');
            else
                fprintf(' (negligible)\n\n');
            end
        end

        % Create comprehensive summary table (pre-allocate to avoid warnings)
        n_results = length(binary_results);
        binary_summary = table( ...
            cell(n_results, 1), ...      % Comparison
            cell(n_results, 1), ...      % DecisionScore
            zeros(n_results, 1), ...     % n_Group
            zeros(n_results, 1), ...     % n_HC
            zeros(n_results, 1), ...     % Mean_Group
            zeros(n_results, 1), ...     % SD_Group
            zeros(n_results, 1), ...     % Mean_HC
            zeros(n_results, 1), ...     % SD_HC
            zeros(n_results, 1), ...     % t_statistic
            zeros(n_results, 1), ...     % df
            zeros(n_results, 1), ...     % p_uncorrected
            zeros(n_results, 1), ...     % p_FDR
            false(n_results, 1), ...     % FDR_significant
            zeros(n_results, 1), ...     % Cohens_d
            'VariableNames', {'Comparison', 'DecisionScore', 'n_Group', 'n_HC', ...
                'Mean_Group', 'SD_Group', 'Mean_HC', 'SD_HC', 't_statistic', ...
                'df', 'p_uncorrected', 'p_FDR', 'FDR_significant', 'Cohens_d'});

        for i = 1:n_results
            res = binary_results{i};
            binary_summary.Comparison{i} = res.comparison;
            binary_summary.DecisionScore{i} = res.score;
            binary_summary.n_Group(i) = res.n1;
            binary_summary.n_HC(i) = res.n2;
            binary_summary.Mean_Group(i) = res.mean1;
            binary_summary.SD_Group(i) = res.sd1;
            binary_summary.Mean_HC(i) = res.mean2;
            binary_summary.SD_HC(i) = res.sd2;
            binary_summary.t_statistic(i) = res.tstat;
            binary_summary.df(i) = res.df;
            binary_summary.p_uncorrected(i) = res.p;
            binary_summary.p_FDR(i) = adj_p_binary(i);
            binary_summary.FDR_significant(i) = h_fdr_binary(i);
            binary_summary.Cohens_d(i) = res.d;
        end

        writetable(binary_summary, [data_out_path 'Binary_Diagnosis_Comparisons_COMPREHENSIVE_FDR.csv']);
        fprintf('  ? Saved: Binary_Diagnosis_Comparisons_COMPREHENSIVE_FDR.csv (%d comparisons)\n\n', height(binary_summary));
    end

    % Create 2×2 diagnostic matrix
    fprintf('2×2 DIAGNOSTIC MATRIX (Anxiety × Depression):\n');
    fprintf('  Neither (HC): n=%d\n', sum(binary_data.anxiety_binary == 0 & binary_data.depression_binary == 0));
    fprintf('  Anxiety only: n=%d\n', sum(binary_data.anxiety_binary == 1 & binary_data.depression_binary == 0));
    fprintf('  Depression only: n=%d\n', sum(binary_data.anxiety_binary == 0 & binary_data.depression_binary == 1));
    fprintf('  Both (Comorbid): n=%d\n\n', sum(binary_data.anxiety_binary == 1 & binary_data.depression_binary == 1));

else
    fprintf('  WARNING: diagnosis_group not available, skipping binary analysis\n');
end

fprintf('SECTION 10E COMPLETE: Binary diagnostic coding (FULL) complete\n\n');

%% ==========================================================================
%  SECTION 10F: PARTIAL CORRELATIONS (FULL IMPLEMENTATION)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  SECTION 10F: PARTIAL CORRELATIONS (FULL)       |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('COMPUTING PARTIAL CORRELATIONS FOR ALL CLINICAL VARIABLES\n');
fprintf('  Controlling for: Age, Sex, Site (whichever are available)\n');
fprintf('  Rationale: Removes spurious associations due to demographic differences\n\n');

% Check if control variables exist (use actual variable names from dataset)
has_age = ismember('Age', analysis_data.Properties.VariableNames);
has_sex = ismember('Sexe', analysis_data.Properties.VariableNames);  % Dataset uses 'Sexe' not 'Sex'
has_site = ismember('ascanloc', analysis_data.Properties.VariableNames) || ...
           ismember('site', analysis_data.Properties.VariableNames) || ...
           ismember('Site', analysis_data.Properties.VariableNames);

if has_site
    if ismember('ascanloc', analysis_data.Properties.VariableNames)
        site_var = 'ascanloc';  % Scan location variable
    elseif ismember('site', analysis_data.Properties.VariableNames)
        site_var = 'site';
    else
        site_var = 'Site';
    end
end

fprintf('  Available control variables:\n');
fprintf('    Age: %s\n', ternary(has_age, 'YES', 'NO'));
fprintf('    Sex (Sexe): %s\n', ternary(has_sex, 'YES', 'NO'));
fprintf('    Site (ascanloc): %s\n\n', ternary(has_site, 'YES', 'NO'));

% Build control matrix from available variables
control_var_list = {};
control_var_names = {};

if has_age
    control_var_list{end+1} = analysis_data.Age;
    control_var_names{end+1} = 'Age';
end

if has_sex
    control_var_list{end+1} = analysis_data.Sexe;  % Use 'Sexe' (French spelling)
    control_var_names{end+1} = 'Sex';
end

if has_site
    control_var_list{end+1} = analysis_data.(site_var);
    control_var_names{end+1} = 'Site';
end

% Proceed if we have at least one control variable
if ~isempty(control_var_list)
    % Create control matrix
    control_vars = [];
    for i = 1:length(control_var_list)
        control_vars = [control_vars, control_var_list{i}];
    end

    control_names_str = strjoin(control_var_names, ', ');
    fprintf('  Computing partial correlations controlling for: %s\n', control_names_str);
    fprintf('  Analyzing ALL %d clinical variables...\n\n', length(all_vars));

    % Create comprehensive partial correlation tables
    partial_corr_26 = [];
    partial_corr_bvftd = [];

    % Iterate through all clinical variables
    for v = 1:length(all_vars)
        var_name = all_vars{v};

        % Try to find the variable in analysis_data
        if ismember(var_name, analysis_data.Properties.VariableNames)
            var_data = analysis_data.(var_name);
        else
            % Variable not found - append NaN rows to maintain alignment
            partial_corr_26 = [partial_corr_26; NaN, NaN, NaN, NaN, NaN, 0];
            partial_corr_bvftd = [partial_corr_bvftd; NaN, NaN, NaN, NaN, NaN, 0];
            continue;
        end

        % =================================================================
        % TRANSITION-26 PARTIAL CORRELATION
        % =================================================================
        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26) & ...
                    all(~isnan(control_vars), 2);
        
        if sum(valid_idx) >= 40  % Need larger n for partial correlations
            try
                % Zero-order correlation
                [r_zero, p_zero] = corr(var_data(valid_idx), ...
                                        analysis_data.Transition_26(valid_idx));
                
                % Partial correlation
                [r_partial, p_partial] = partialcorr(var_data(valid_idx), ...
                                                      analysis_data.Transition_26(valid_idx), ...
                                                      control_vars(valid_idx, :));
                
                partial_corr_26 = [partial_corr_26; r_zero, p_zero, r_partial, p_partial, ...
                                   abs(r_zero - r_partial), sum(valid_idx)];
            catch
                partial_corr_26 = [partial_corr_26; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
            end
        else
            partial_corr_26 = [partial_corr_26; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
        end

        % =================================================================
        % bvFTD PARTIAL CORRELATION
        % =================================================================
        valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD) & ...
                    all(~isnan(control_vars), 2);
        
        if sum(valid_idx) >= 40
            try
                % Zero-order correlation
                [r_zero, p_zero] = corr(var_data(valid_idx), ...
                                        analysis_data.bvFTD(valid_idx));
                
                % Partial correlation
                [r_partial, p_partial] = partialcorr(var_data(valid_idx), ...
                                                      analysis_data.bvFTD(valid_idx), ...
                                                      control_vars(valid_idx, :));
                
                partial_corr_bvftd = [partial_corr_bvftd; r_zero, p_zero, r_partial, p_partial, ...
                                      abs(r_zero - r_partial), sum(valid_idx)];
            catch
                partial_corr_bvftd = [partial_corr_bvftd; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
            end
        else
            partial_corr_bvftd = [partial_corr_bvftd; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
        end
    end

    % Create comprehensive table
    partial_summary = table();
    partial_summary.Variable = all_vars';
    partial_summary.Category = all_categories';
    
    % Transition-26 results
    partial_summary.Trans26_r_zero = partial_corr_26(:,1);
    partial_summary.Trans26_p_zero = partial_corr_26(:,2);
    partial_summary.Trans26_r_partial = partial_corr_26(:,3);
    partial_summary.Trans26_p_partial = partial_corr_26(:,4);
    partial_summary.Trans26_delta_r = partial_corr_26(:,5);
    partial_summary.Trans26_n = partial_corr_26(:,6);
    
    % bvFTD results
    partial_summary.bvFTD_r_zero = partial_corr_bvftd(:,1);
    partial_summary.bvFTD_p_zero = partial_corr_bvftd(:,2);
    partial_summary.bvFTD_r_partial = partial_corr_bvftd(:,3);
    partial_summary.bvFTD_p_partial = partial_corr_bvftd(:,4);
    partial_summary.bvFTD_delta_r = partial_corr_bvftd(:,5);
    partial_summary.bvFTD_n = partial_corr_bvftd(:,6);

    writetable(partial_summary, [data_out_path 'Partial_Correlations_COMPREHENSIVE.csv']);
    fprintf('  ? Saved: Partial_Correlations_COMPREHENSIVE.csv (%d variables)\n\n', height(partial_summary));

    % Print summary statistics
    fprintf('  SUMMARY OF PARTIAL CORRELATION EFFECTS:\n');
    
    valid_26 = ~isnan(partial_corr_26(:,5));
    if sum(valid_26) > 0
        fprintf('    Transition-26 (n=%d variables with valid results):\n', sum(valid_26));
        fprintf('      Mean reduction in |r|: %.3f\n', mean(partial_corr_26(valid_26,5)));
        fprintf('      Median reduction: %.3f\n', median(partial_corr_26(valid_26,5)));
        fprintf('      Variables with >0.05 reduction: %d (%.1f%%)\n', ...
                sum(partial_corr_26(valid_26,5) > 0.05), ...
                100*sum(partial_corr_26(valid_26,5) > 0.05)/sum(valid_26));
    end
    
    valid_bv = ~isnan(partial_corr_bvftd(:,5));
    if sum(valid_bv) > 0
        fprintf('    bvFTD (n=%d variables with valid results):\n', sum(valid_bv));
        fprintf('      Mean reduction in |r|: %.3f\n', mean(partial_corr_bvftd(valid_bv,5)));
        fprintf('      Median reduction: %.3f\n', median(partial_corr_bvftd(valid_bv,5)));
        fprintf('      Variables with >0.05 reduction: %d (%.1f%%)\n\n', ...
                sum(partial_corr_bvftd(valid_bv,5) > 0.05), ...
                100*sum(partial_corr_bvftd(valid_bv,5) > 0.05)/sum(valid_bv));
    end

else
    fprintf('  WARNING: No control variables (Age, Sexe, Site) available\n');
    fprintf('           Cannot compute partial correlations without controls\n\n');
end

fprintf('SECTION 10F COMPLETE: Partial correlation analysis (FULL) complete\n\n');

%% ==========================================================================
%  SECTION 10G: SPEARMAN vs PEARSON CORRELATIONS (FULL IMPLEMENTATION)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  SECTION 10G: SPEARMAN CORRELATIONS (FULL)      |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('COMPUTING SPEARMAN vs PEARSON FOR ALL CLINICAL VARIABLES\n');
fprintf('  Rationale: Spearman is robust to outliers and non-linear relationships\n');
fprintf('  Provides systematic comparison across all %d variables\n\n', length(all_vars));

% Create comprehensive Spearman comparison tables
spearman_comp_26 = [];
spearman_comp_bvftd = [];

% Iterate through all clinical variables
for v = 1:length(all_vars)
    var_name = all_vars{v};

    % Try to find the variable in analysis_data
    if ismember(var_name, analysis_data.Properties.VariableNames)
        var_data = analysis_data.(var_name);
    else
        % Variable not found - append NaN rows to maintain alignment
        spearman_comp_26 = [spearman_comp_26; NaN, NaN, NaN, NaN, NaN, 0];
        spearman_comp_bvftd = [spearman_comp_bvftd; NaN, NaN, NaN, NaN, NaN, 0];
        continue;
    end

    % =================================================================
    % TRANSITION-26 CORRELATIONS
    % =================================================================
    valid_idx = ~isnan(var_data) & ~isnan(analysis_data.Transition_26);
    
    if sum(valid_idx) >= 30
        try
            % Pearson
            [r_pearson, p_pearson] = corr(var_data(valid_idx), ...
                                          analysis_data.Transition_26(valid_idx), ...
                                          'Type', 'Pearson');
            
            % Spearman
            [r_spearman, p_spearman] = corr(var_data(valid_idx), ...
                                            analysis_data.Transition_26(valid_idx), ...
                                            'Type', 'Spearman');
            
            spearman_comp_26 = [spearman_comp_26; r_pearson, p_pearson, r_spearman, p_spearman, ...
                               abs(r_pearson - r_spearman), sum(valid_idx)];
        catch
            spearman_comp_26 = [spearman_comp_26; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
        end
    else
        spearman_comp_26 = [spearman_comp_26; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
    end

    % =================================================================
    % bvFTD CORRELATIONS
    % =================================================================
    valid_idx = ~isnan(var_data) & ~isnan(analysis_data.bvFTD);
    
    if sum(valid_idx) >= 30
        try
            % Pearson
            [r_pearson, p_pearson] = corr(var_data(valid_idx), ...
                                          analysis_data.bvFTD(valid_idx), ...
                                          'Type', 'Pearson');
            
            % Spearman
            [r_spearman, p_spearman] = corr(var_data(valid_idx), ...
                                            analysis_data.bvFTD(valid_idx), ...
                                            'Type', 'Spearman');
            
            spearman_comp_bvftd = [spearman_comp_bvftd; r_pearson, p_pearson, r_spearman, p_spearman, ...
                                   abs(r_pearson - r_spearman), sum(valid_idx)];
        catch
            spearman_comp_bvftd = [spearman_comp_bvftd; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
        end
    else
        spearman_comp_bvftd = [spearman_comp_bvftd; NaN, NaN, NaN, NaN, NaN, sum(valid_idx)];
    end
end

% Create comprehensive table
spearman_summary = table();
spearman_summary.Variable = all_vars';
spearman_summary.Category = all_categories';

% Transition-26 results
spearman_summary.Trans26_r_Pearson = spearman_comp_26(:,1);
spearman_summary.Trans26_p_Pearson = spearman_comp_26(:,2);
spearman_summary.Trans26_rho_Spearman = spearman_comp_26(:,3);
spearman_summary.Trans26_p_Spearman = spearman_comp_26(:,4);
spearman_summary.Trans26_delta = spearman_comp_26(:,5);
spearman_summary.Trans26_n = spearman_comp_26(:,6);

% bvFTD results
spearman_summary.bvFTD_r_Pearson = spearman_comp_bvftd(:,1);
spearman_summary.bvFTD_p_Pearson = spearman_comp_bvftd(:,2);
spearman_summary.bvFTD_rho_Spearman = spearman_comp_bvftd(:,3);
spearman_summary.bvFTD_p_Spearman = spearman_comp_bvftd(:,4);
spearman_summary.bvFTD_delta = spearman_comp_bvftd(:,5);
spearman_summary.bvFTD_n = spearman_comp_bvftd(:,6);

writetable(spearman_summary, [data_out_path 'Spearman_vs_Pearson_COMPREHENSIVE.csv']);
fprintf('  ? Saved: Spearman_vs_Pearson_COMPREHENSIVE.csv (%d variables)\n\n', height(spearman_summary));

% Print summary statistics
fprintf('  SUMMARY OF SPEARMAN vs PEARSON DIFFERENCES:\n');

valid_26 = ~isnan(spearman_comp_26(:,5));
if sum(valid_26) > 0
    fprintf('    Transition-26 (n=%d variables with valid results):\n', sum(valid_26));
    fprintf('      Mean |?|: %.3f\n', mean(spearman_comp_26(valid_26,5)));
    fprintf('      Median |?|: %.3f\n', median(spearman_comp_26(valid_26,5)));
    fprintf('      Variables with substantial difference (|?|>0.05): %d (%.1f%%)\n', ...
            sum(spearman_comp_26(valid_26,5) > 0.05), ...
            100*sum(spearman_comp_26(valid_26,5) > 0.05)/sum(valid_26));
    
    if mean(spearman_comp_26(valid_26,5)) > 0.05
        fprintf('      ? MODERATE differences suggest non-linear relationships or outliers\n');
    else
        fprintf('      ? MINIMAL differences suggest approximately linear relationships\n');
    end
end

valid_bv = ~isnan(spearman_comp_bvftd(:,5));
if sum(valid_bv) > 0
    fprintf('    bvFTD (n=%d variables with valid results):\n', sum(valid_bv));
    fprintf('      Mean |?|: %.3f\n', mean(spearman_comp_bvftd(valid_bv,5)));
    fprintf('      Median |?|: %.3f\n', median(spearman_comp_bvftd(valid_bv,5)));
    fprintf('      Variables with substantial difference (|?|>0.05): %d (%.1f%%)\n\n', ...
            sum(spearman_comp_bvftd(valid_bv,5) > 0.05), ...
            100*sum(spearman_comp_bvftd(valid_bv,5) > 0.05)/sum(valid_bv));
    
    if mean(spearman_comp_bvftd(valid_bv,5)) > 0.05
        fprintf('      ? MODERATE differences suggest non-linear relationships or outliers\n\n');
    else
        fprintf('      ? MINIMAL differences suggest approximately linear relationships\n\n');
    end
end

fprintf('SECTION 10G COMPLETE: Spearman correlation comparison (FULL) complete\n\n');

%% ==========================================================================
%  SECTION 10H: FINAL VALIDATION - CURATED VARIABLE INDEPENDENCE
%  ==========================================================================
%  PURPOSE: Verify that the final curated variable set is sufficiently
%  independent for valid FDR correction. This is a quality control check.
%
%  SCIENTIFIC RATIONALE: Benjamini-Hochberg FDR assumes independence or
%  positive dependence among test statistics. If curated variables still
%  correlate strongly (|r| > 0.70), the FDR correction may be liberal
%  (more false positives than expected).
%
%  Reference: Benjamini & Hochberg (1995), Benjamini & Yekutieli (2001)
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|  SECTION 10H: FINAL VALIDATION (CURATED VARS)   |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('VERIFYING INDEPENDENCE OF FINAL CURATED VARIABLE SET\n\n');

% Compile final curated variables used in analysis
if exist('all_curated_vars', 'var') && ~isempty(all_curated_vars)
    fprintf('  Curated variables in final set: %d\n', length(all_curated_vars));
    fprintf('  Variables: %s\n\n', strjoin(all_curated_vars, ', '));

    % Identify which curated variables are available in analysis_data
    available_curated_final = {};
    for i = 1:length(all_curated_vars)
        if ismember(all_curated_vars{i}, analysis_data.Properties.VariableNames)
            available_curated_final{end+1} = all_curated_vars{i};
        end
    end

    if length(available_curated_final) >= 2
        fprintf('  Computing correlation matrix for %d available curated variables...\n\n', ...
            length(available_curated_final));

        % Build correlation matrix for curated variables
        n_curated = length(available_curated_final);
        curated_corr_matrix = NaN(n_curated, n_curated);

        for i = 1:n_curated
            for j = 1:n_curated
                var_i = analysis_data.(available_curated_final{i});
                var_j = analysis_data.(available_curated_final{j});

                valid_idx = ~isnan(var_i) & ~isnan(var_j);
                if sum(valid_idx) >= MIN_SAMPLE_SIZE
                    curated_corr_matrix(i, j) = corr(var_i(valid_idx), var_j(valid_idx));
                end
            end
        end

        % Find maximum off-diagonal correlation
        off_diag_mask = ~eye(n_curated);
        off_diag_corrs = curated_corr_matrix(off_diag_mask);
        off_diag_corrs = off_diag_corrs(~isnan(off_diag_corrs));

        if ~isempty(off_diag_corrs)
            [max_abs_r, max_idx] = max(abs(off_diag_corrs));

            % Find which pair has maximum correlation
            [row_max, col_max] = find(abs(curated_corr_matrix) == max_abs_r & off_diag_mask, 1);
            if ~isempty(row_max) && ~isempty(col_max)
                var1_name = available_curated_final{row_max};
                var2_name = available_curated_final{col_max};
                actual_r = curated_corr_matrix(row_max, col_max);

                fprintf('  MAXIMUM INTER-VARIABLE CORRELATION IN CURATED SET:\n');
                fprintf('    r = %.3f between %s and %s\n\n', actual_r, var1_name, var2_name);

                % Evaluate quality
                if max_abs_r > REDUNDANCY_THRESHOLD
                    fprintf('  *** WARNING: Max correlation (|r|=%.3f) EXCEEDS threshold (%.2f)\n', ...
                        max_abs_r, REDUNDANCY_THRESHOLD);
                    fprintf('      Consider reviewing these variables for potential redundancy.\n');
                    fprintf('      FDR correction may be slightly liberal.\n\n');
                else
                    fprintf('  ✓ PASS: Max correlation (|r|=%.3f) ≤ threshold (%.2f)\n', ...
                        max_abs_r, REDUNDANCY_THRESHOLD);
                    fprintf('      Curated variables are sufficiently independent.\n');
                    fprintf('      FDR correction is methodologically valid.\n\n');
                end
            end

            % Report summary statistics
            fprintf('  CORRELATION SUMMARY FOR CURATED VARIABLES:\n');
            fprintf('    Mean |r|: %.3f\n', mean(abs(off_diag_corrs)));
            fprintf('    Median |r|: %.3f\n', median(abs(off_diag_corrs)));
            fprintf('    Max |r|: %.3f\n', max_abs_r);
            fprintf('    Pairs with |r| > %.2f: %d/%d (%.1f%%)\n\n', ...
                REDUNDANCY_THRESHOLD, ...
                sum(abs(off_diag_corrs) > REDUNDANCY_THRESHOLD), ...
                length(off_diag_corrs), ...
                100*sum(abs(off_diag_corrs) > REDUNDANCY_THRESHOLD)/length(off_diag_corrs));

            % Save curated variable correlation matrix
            curated_corr_table = array2table(curated_corr_matrix, ...
                'VariableNames', available_curated_final, ...
                'RowNames', available_curated_final);
            writetable(curated_corr_table, [data_out_path 'FINAL_Curated_Variable_Correlations.csv'], ...
                'WriteRowNames', true);
            fprintf('  Saved: FINAL_Curated_Variable_Correlations.csv\n\n');
        else
            fprintf('  WARNING: No valid correlations computed\n\n');
        end
    else
        fprintf('  WARNING: Fewer than 2 curated variables available in dataset\n');
        fprintf('           Cannot compute inter-variable correlations\n\n');
    end
else
    fprintf('  WARNING: all_curated_vars not defined\n');
    fprintf('           Skipping final validation\n\n');
end

fprintf('SECTION 10H COMPLETE: Final validation of curated variable independence\n\n');

%% ==========================================================================
%  SECTION 11: SAVE COMPLETE ANALYSIS DATASET
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|         SAVING COMPLETE DATASET                  |\n');
fprintf('---------------------------------------------------\n\n');

writetable(analysis_data, [data_out_path 'COMPLETE_Analysis_Dataset_Priority_4_1_to_4_5.csv']);
fprintf('  Saved: COMPLETE_Analysis_Dataset_Priority_4_1_to_4_5.csv\n');
fprintf('    Dimensions: [%d subjects × %d variables]\n\n', ...
    height(analysis_data), width(analysis_data));

%% ==========================================================================
%  SECTION 12: FINAL SUMMARY
%  ==========================================================================
fprintf('---------------------------------------------------\n');
fprintf('|              FINAL SUMMARY                       |\n');
fprintf('---------------------------------------------------\n\n');

fprintf('ALL PRIORITIES COMPLETE!\n\n');

fprintf('ANALYSES COMPLETED:\n');
fprintf('  4.1: Metabolic Subtypes & BMI\n');
fprintf('  4.2: Symptom Severity (%d variables) + ENHANCED PCA\n', length(available_symptom_vars));
fprintf('    - PC1, PC2, PC3 correlations with ALL decision scores\n');
fprintf('  4.3: Clinical History\n');
fprintf('    - Age of Onset: %d variables\n', length(available_age_onset_vars));
fprintf('    - Illness Duration: %d variables\n', length(illness_duration_vars));
fprintf('    - Recency: %d variables (NEW - OPTION 6)\n', length(available_recency_vars));
fprintf('    - History Variables: %d variables\n', length(available_clinical_history_vars));
fprintf('    - Childhood Adversity: %d variables\n', length(available_childhood_vars));
fprintf('  4.4: Cognition & Functioning (%d variables)\n', length(cognitive_vars_found));
fprintf('  4.4B: Medication Analysis (PATIENTS ONLY) - CORRECTED CODING\n');
if exist('all_med_vars', 'var')
    fprintf('    - DDD/Continuous variables: %d\n', length(available_med_ddd));
    fprintf('    - Frequency variables (0/1/2): %d\n', length(available_med_freq));
    fprintf('    - Binary variables (0/1): %d\n', length(available_med_binary));
end
fprintf('  NEW 9C: Recency Stratified Analysis (OPTION 6)\n');
fprintf('    - Tests if symptom-brain associations differ by illness phase\n');
fprintf('  4.5: Comprehensive Statistical Summary\n');
fprintf('  Demographics Analysis (%d variables)\n', length(demo_vars_analyzed));
fprintf('  NEW 10E: Binary Anxiety/Depression Coding\n');
fprintf('    - Comorbid individuals coded as "yes" for BOTH disorders\n');
fprintf('    - Direct HC comparisons with t-tests and Cohen''s d\n');
fprintf('  NEW 10F: Partial Correlations\n');
fprintf('    - Controls for Age, Sex, and Site confounds\n');
fprintf('    - Identifies robust associations beyond demographics\n');
fprintf('  NEW 10G: Spearman vs Pearson Correlations\n');
fprintf('    - Non-parametric alternative robust to outliers\n');
fprintf('    - Identifies non-linear monotonic relationships\n');
fprintf('  NEW 10H: Final Validation - Curated Variable Independence\n');
fprintf('    - Verifies curated variables are sufficiently independent\n');
fprintf('    - Validates FDR correction methodology\n\n');

fprintf('VARIABLE REDUNDANCY ANALYSIS (SECTION 0B/5D):\n');
if exist('total_curated_final', 'var') && exist('total_original_final', 'var')
    fprintf('  Original variables: %d\n', total_original_final);
    fprintf('  Curated (non-redundant): %d\n', total_curated_final);
    fprintf('  Reduction: %.1f%%\n\n', 100*(1 - total_curated_final/total_original_final));
else
    fprintf('  See Section 5D for variable reduction summary\n\n');
end

fprintf('DECISION SCORES USED:\n');
fprintf('  - Transition-26 (OOCV-26)\n');
fprintf('  - bvFTD (OOCV-7)\n\n');

fprintf('OUTPUT FILES CREATED:\n');
fprintf('  Data Tables:\n');
fprintf('    * COMPLETE_Analysis_Dataset_Priority_4_1_to_4_5.csv\n');
fprintf('    * COMPREHENSIVE_SUMMARY_All_Associations_CURATED.csv\n');
fprintf('    * Variable_Exclusion_Decisions.csv (NEW - redundancy analysis)\n');
fprintf('    * Variable_Correlation_Matrix.csv (NEW - all clinical vars)\n');
fprintf('    * Curated_Variable_Summary.csv (NEW - reduction by domain)\n');
fprintf('    * FINAL_Curated_Variable_Correlations.csv (NEW - validation)\n');
fprintf('    * Summary_Symptom_Correlations.csv\n');
fprintf('    * Summary_Symptom_PCA_Loadings.csv\n');
fprintf('    * Summary_PCA_Correlations_ALL_Components.csv\n');
fprintf('    * Summary_Clinical_History_Correlations.csv\n');
fprintf('    * Summary_Recency_Correlations.csv (NEW - OPTION 6)\n');
fprintf('    * Summary_Recency_Stratified_Analysis.csv (NEW - OPTION 6)\n');
fprintf('    * Summary_Childhood_Adversity_Correlations.csv\n');
fprintf('    * Summary_Cognition_Functioning_Correlations.csv\n');
fprintf('    * Summary_Demographics_Correlations.csv\n');
fprintf('    * Summary_Medication_Correlations_PatientsOnly_CORRECTED.csv (CORRECTED)\n');
fprintf('    * Binary_Diagnosis_Comparisons_COMPREHENSIVE_FDR.csv (NEW - SECTION 10E - FULL)\n');
fprintf('    * Partial_Correlations_COMPREHENSIVE.csv (NEW - SECTION 10F - FULL)\n');
fprintf('    * Spearman_vs_Pearson_COMPREHENSIVE.csv (NEW - SECTION 10G - FULL)\n');

fprintf('\n  Figures:\n');
fprintf('    * Fig_0B_Correlation_Heatmap_AllVariables.png/.fig (NEW - redundancy)\n');
fprintf('    * Fig_0B_Variable_Dendrogram.png/.fig (NEW - redundancy)\n');
fprintf('    * Fig_4_1_Metabolic_Subtypes.png/.fig\n');
fprintf('    * Fig_4_1_BMI_Correlations.png/.fig\n');
fprintf('    * Fig_4_2_PCA_Comprehensive_Loadings.png/.fig\n');
fprintf('    * Fig_4_2_PCA_Correlation_Heatmap.png/.fig\n');
fprintf('    * Fig_4_2_PCA_Significant_Correlations.png/.fig\n');
fprintf('    * Fig_4_2_Symptom_Correlations_Heatmap.png/.fig\n');
fprintf('    * Fig_9C_Recency_Stratified_Symptom_Correlations.png/.fig (NEW - OPTION 6)\n');
fprintf('    * Fig_4_5_Forest_Plot_Significant_Associations_Transition26.png/.fig\n');
if sum(sig_idx_bvftd) > 0
    fprintf('    * Fig_4_5_Forest_Plot_Significant_Associations_bvFTD.png/.fig\n');
end
if exist('all_med_vars', 'var') && exist('sig_idx_med', 'var') && sum(sig_idx_med) > 0
    fprintf('    * Fig_Medication_Significant_Associations_PatientsOnly_CORRECTED.png/.fig\n');
end

fprintf('\n  Log:\n');
fprintf('    * Priority_4_1_to_4_5_Complete_Analysis_Log_OOCV26_bvFTD.txt\n');

fprintf('\nKEY FINDINGS SUMMARY:\n');
fprintf('  Total clinical variables analyzed: %d\n', length(all_vars));
fprintf('  Significant associations (p<0.05):\n');
fprintf('    - Transition-26: %d (%.1f%%)\n', sum(all_corr_26(:,2) < 0.05), ...
    100*sum(all_corr_26(:,2) < 0.05)/length(all_vars));
fprintf('    - bvFTD: %d (%.1f%%)\n', sum(all_corr_bvftd(:,2) < 0.05), ...
    100*sum(all_corr_bvftd(:,2) < 0.05)/length(all_vars));

if exist('all_med_vars', 'var')
    fprintf('\n  Medication variables analyzed (PATIENTS ONLY - CORRECTED): %d\n', length(all_med_vars));
    fprintf('    - DDD/Continuous: %d\n', length(available_med_ddd));
    fprintf('    - Frequency (0/1/2): %d\n', length(available_med_freq));
    fprintf('    - Binary (0/1): %d\n', length(available_med_binary));
    fprintf('    - Significant (Transition-26): %d (%.1f%%)\n', ...
        sum(all_med_corr_26(:,2) < 0.05), ...
        100*sum(all_med_corr_26(:,2) < 0.05)/length(all_med_vars));
end

fprintf('\nPCA COMPONENT SUMMARY:\n');
if exist('explained', 'var') && ~isempty(explained) && length(explained) >= 3
    fprintf('  PC1: %.1f%% variance\n', explained(1));
    fprintf('  PC2: %.1f%% variance\n', explained(2));
    fprintf('  PC3: %.1f%% variance\n', explained(3));
    fprintf('  Total (PC1-PC3): %.1f%% variance\n\n', sum(explained(1:3)));
else
    fprintf('  PCA not performed (insufficient data)\n\n');
end

if exist('results_recency', 'var') && isfield(results_recency, 'stratified_analysis')
    fprintf('RECENCY STRATIFIED ANALYSIS (OPTION 6):\n');
    fprintf('  Groups analyzed: Recent vs. Remitted illness\n');
    fprintf('  Variables tested: %d\n', height(results_recency.stratified_analysis));
    sig_fisher = sum(results_recency.stratified_analysis.Fisher_p < 0.05);
    if sig_fisher > 0
        fprintf('  *** %d variable(s) show DIFFERENT associations by recency group\n\n', sig_fisher);
    else
        fprintf('  No significant differences between recency groups\n\n');
    end
end

fprintf('End time: %s\n', datestr(now));
fprintf('---------------------------------------------------\n\n');

diary off;

fprintf('ALL DONE! Script executed successfully with OOCV-26 and bvFTD.\n');
fprintf('Check output folder: %s\n\n', results_path);

%% ==========================================================================
%  HELPER FUNCTIONS (SESSION 3: FEATURE 3.1)
%  ==========================================================================

function [r, CI, p, n_valid] = calculate_correlation_with_CI(x, y, alpha)
    % CALCULATE PEARSON CORRELATION WITH FISHER'S Z CONFIDENCE INTERVAL
    %
    % Computes Pearson correlation coefficient and confidence interval using
    % Fisher's Z transformation for improved accuracy with small samples
    %
    % INPUTS:
    %   x       - First variable (numeric vector)
    %   y       - Second variable (numeric vector)
    %   alpha   - Significance level for CI (default: 0.05 for 95% CI)
    %
    % OUTPUTS:
    %   r       - Pearson correlation coefficient
    %   CI      - Confidence interval [lower, upper]
    %   p       - Two-tailed p-value
    %   n_valid - Number of valid pairs (after removing NaN)
    %
    % METHOD:
    %   Fisher's Z transformation: z = atanh(r)
    %   SE(z) = 1/sqrt(n-3)
    %   CI(z) = z ± z_critical * SE(z)
    %   Back-transform: CI(r) = tanh(CI(z))
    %
    % REFERENCE:
    %   Fisher, R.A. (1915). Frequency distribution of the values of the
    %   correlation coefficient in samples from an indefinitely large population.
    %   Biometrika, 10(4), 507-521.
    %
    % Author: Claude AI Assistant (SESSION 3)
    % Date: November 8, 2025

    if nargin < 3
        alpha = 0.05;
    end

    % Remove NaN values (pairwise deletion)
    valid_idx = ~isnan(x) & ~isnan(y);
    x_clean = x(valid_idx);
    y_clean = y(valid_idx);
    n_valid = length(x_clean);

    % Handle edge cases
    if n_valid < 3
        r = NaN;
        CI = [NaN, NaN];
        p = NaN;
        return;
    end

    % Calculate Pearson correlation
    [r, p] = corr(x_clean, y_clean);

    % Fisher's Z transformation for CI
    if abs(r) >= 1
        % Perfect correlation - CI is undefined
        CI = [r, r];
    else
        z_crit = norminv(1 - alpha/2);  % Two-tailed critical value
        z_r = atanh(r);                 % Fisher's Z transform
        se_z = 1 / sqrt(n_valid - 3);   % Standard error of Z

        % CI in Z space
        z_lower = z_r - z_crit * se_z;
        z_upper = z_r + z_crit * se_z;

        % Back-transform to correlation space
        CI = [tanh(z_lower), tanh(z_upper)];
    end
end

function effect_label = interpret_effect_size(r, type)
    % INTERPRET EFFECT SIZE MAGNITUDE
    %
    % Provides standardized interpretation of effect sizes following Cohen's conventions.
    % This helps translate statistical significance into practical significance.
    %
    % INPUTS:
    %   r    - Effect size (correlation coefficient for 'r', Cohen's d for 'd')
    %   type - 'r' for correlations, 'd' for Cohen's d (default: 'r')
    %
    % OUTPUTS:
    %   effect_label - String describing effect magnitude:
    %                  'negligible', 'small', 'medium', 'large', 'very large'
    %
    % THRESHOLDS (Cohen, 1988):
    %   Correlations (r): small=0.10, medium=0.30, large=0.50
    %   Cohen's d:        small=0.20, medium=0.50, large=0.80
    %
    % REFERENCE:
    %   Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences (2nd ed.).
    %   Hillsdale, NJ: Lawrence Erlbaum Associates.
    %
    % Author: Claude AI Assistant
    % Date: November 8, 2025

    if nargin < 2
        type = 'r';  % Default to correlation
    end

    abs_r = abs(r);  % Use absolute value for magnitude interpretation

    if strcmpi(type, 'r')
        % Correlation thresholds
        if abs_r < 0.10
            effect_label = 'negligible';
        elseif abs_r < 0.30
            effect_label = 'small';
        elseif abs_r < 0.50
            effect_label = 'medium';
        elseif abs_r < 0.70
            effect_label = 'large';
        else
            effect_label = 'very large';
        end
    elseif strcmpi(type, 'd')
        % Cohen's d thresholds
        if abs_r < 0.20
            effect_label = 'negligible';
        elseif abs_r < 0.50
            effect_label = 'small';
        elseif abs_r < 0.80
            effect_label = 'medium';
        else
            effect_label = 'large';
        end
    else
        error('Unknown effect size type: %s. Use ''r'' or ''d''.', type);
    end
end

function fig = create_forest_plot(var_names, labels, correlations, CIs, p_vals, p_fdr, title_text, n_subjects, marker_color)
    % CREATE STANDARDIZED FOREST PLOT FOR CORRELATION RESULTS
    %
    % Generates a publication-ready forest plot showing correlations with
    % confidence intervals, with optional FDR significance markers
    %
    % INPUTS:
    %   var_names     - Cell array of variable names
    %   labels        - Cell array of interpretable labels for display
    %   correlations  - Vector of correlation coefficients (r values)
    %   CIs           - N×2 matrix of confidence intervals [lower, upper]
    %   p_vals        - Vector of uncorrected p-values
    %   p_fdr         - Vector of FDR-adjusted p-values (optional, [] if not available)
    %   title_text    - Plot title string
    %   n_subjects    - Vector of sample sizes per correlation
    %   marker_color  - RGB triplet for marker color (default: [0.2 0.4 0.8])
    %
    % OUTPUTS:
    %   fig          - Figure handle
    %
    % VISUAL ELEMENTS:
    %   - Horizontal error bars for 95% CI
    %   - Filled circles for correlation coefficients
    %   - Vertical reference line at r=0
    %   - FDR-significant results marked with **
    %   - Uncorrected significant results marked with *
    %
    % Author: Claude AI Assistant (SESSION 3)
    % Date: November 8, 2025

    if nargin < 9 || isempty(marker_color)
        marker_color = [0.2 0.4 0.8];  % Default blue
    end

    n_vars = length(var_names);

    % Create figure
    fig = figure('Position', [100 100 1000 max(400, 100 + 50*n_vars)]);
    hold on;

    y_pos = 1:n_vars;

    % Plot confidence intervals
    for i = 1:n_vars
        plot([CIs(i,1), CIs(i,2)], [y_pos(i), y_pos(i)], 'k-', 'LineWidth', 1.5);
    end

    % Plot correlation coefficients
    scatter(correlations, y_pos, 100, 'filled', 'MarkerFaceColor', marker_color);

    % Add significance markers to labels
    labels_with_sig = labels;
    for i = 1:n_vars
        if ~isempty(p_fdr) && ~isnan(p_fdr(i)) && p_fdr(i) < 0.05
            labels_with_sig{i} = [labels{i} ' **'];  % FDR significant
        elseif p_vals(i) < 0.05
            labels_with_sig{i} = [labels{i} ' *'];   % Uncorrected significant
        end
    end

    % Reference line at r=0
    plot([0 0], [0.5, n_vars+0.5], 'r--', 'LineWidth', 2);

    % Formatting
    set(gca, 'YTick', y_pos, 'YTickLabel', labels_with_sig, 'FontSize', 9);
    xlabel('Correlation (r) with 95% CI', 'FontWeight', 'bold', 'FontSize', 12);
    title(title_text, 'FontWeight', 'bold', 'FontSize', 14);

    % Dynamic x-axis limits
    xlim([min([CIs(:,1); -0.1])-0.1, max([CIs(:,2); 0.1])+0.1]);
    ylim([0.5, n_vars+0.5]);
    grid on;

    % Add legend for significance markers
    if ~isempty(p_fdr)
        text(0.98, 0.02, '* p<0.05   ** FDR q<0.05', ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontSize', 9, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
    else
        text(0.98, 0.02, '* p<0.05', ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontSize', 9, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
    end

    hold off;
end

function [h, crit_p, adj_p] = fdr_bh(pvals, q)
    % BENJAMINI-HOCHBERG FDR CORRECTION
    % Implements the Benjamini-Hochberg procedure for controlling
    % False Discovery Rate in multiple hypothesis testing
    %
    % INPUTS:
    %   pvals - vector of p-values to correct
    %   q     - desired FDR level (default: 0.05)
    %
    % OUTPUTS:
    %   h       - binary vector of significance flags (1=significant, 0=not)
    %   crit_p  - critical p-value threshold
    %   adj_p   - FDR-adjusted p-values (q-values)
    %
    % REFERENCE:
    %   Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery
    %   rate: A practical and powerful approach to multiple testing.
    %   Journal of the Royal Statistical Society, Series B, 57(1), 289-300.
    %
    % Author: Claude AI Assistant
    % Date: November 8, 2025

    if nargin < 2
        q = 0.05;
    end

    % Handle edge cases
    if isempty(pvals)
        h = [];
        crit_p = [];
        adj_p = [];
        return;
    end

    % Convert to column vector
    pvals = pvals(:);
    n = length(pvals);

    % Handle NaN values
    nan_mask = isnan(pvals);
    valid_pvals = pvals(~nan_mask);
    n_valid = length(valid_pvals);

    if n_valid == 0
        h = false(n, 1);
        crit_p = NaN;
        adj_p = NaN(n, 1);
        return;
    end

    % Sort p-values
    [sorted_pvals, sort_idx] = sort(valid_pvals);

    % Calculate BH threshold for each rank
    % P(i) <= (i/m) * q
    ranks = (1:n_valid)';
    bh_threshold = (ranks / n_valid) * q;

    % Find largest i where P(i) <= (i/m)*q
    significant_idx = find(sorted_pvals <= bh_threshold, 1, 'last');

    if isempty(significant_idx)
        crit_p = 0;
        h_valid = false(n_valid, 1);
    else
        crit_p = sorted_pvals(significant_idx);
        h_valid = sorted_pvals <= crit_p;
    end

    % Calculate adjusted p-values (q-values)
    % q(i) = min(1, min_{j>=i} (m/j) * P(j))
    adj_p_sorted = NaN(n_valid, 1);
    adj_p_sorted(n_valid) = sorted_pvals(n_valid);

    for i = (n_valid-1):-1:1
        adj_p_sorted(i) = min(1, min((n_valid / i) * sorted_pvals(i), adj_p_sorted(i+1)));
    end

    % Unsort to match original order
    h_unsorted = false(n_valid, 1);
    h_unsorted(sort_idx) = h_valid;

    adj_p_unsorted = NaN(n_valid, 1);
    adj_p_unsorted(sort_idx) = adj_p_sorted;

    % Insert NaN results back for invalid p-values
    h = false(n, 1);
    h(~nan_mask) = h_unsorted;

    adj_p = NaN(n, 1);
    adj_p(~nan_mask) = adj_p_unsorted;
end

function label = get_label_safe(varname, label_map)
    % SAFE LABEL GETTER WITH ERROR HANDLING
    % Retrieves human-readable label from variable_labels map with fallback
    %
    % INPUTS:
    %   varname   - Variable name (string or char)
    %   label_map - containers.Map with variable name ? label mappings
    %
    % OUTPUT:
    %   label - Readable label (or original varname if not found)

    try
        if ischar(varname) || isstring(varname)
            varname = char(varname);
            if label_map.isKey(varname)
                label = label_map(varname);
            else
                label = varname;
            end
        else
            label = varname;
        end
    catch
        % If anything fails, just return the original variable name
        label = varname;
    end
end

function cmap = redblue(m)
    % REDBLUE COLORMAP
    % Creates a blue-white-red diverging colormap for correlation matrices
    % Blue represents negative correlations, red represents positive
    %
    % INPUT:
    %   m - Number of colors (default: current figure colormap length)
    %
    % OUTPUT:
    %   cmap - m×3 RGB colormap matrix

    if nargin < 1
        m = size(get(gcf,'colormap'),1);
    end

    if mod(m,2) == 0
        m1 = m/2;
        r = [(0:m1-1)'/max(m1-1,1); ones(m1,1)];
        g = [(0:m1-1)'/max(m1-1,1); (m1-1:-1:0)'/max(m1-1,1)];
        b = [ones(m1,1); (m1-1:-1:0)'/max(m1-1,1)];
    else
        m1 = floor(m/2);
        r = [(0:m1-1)'/max(m1,1); ones(m-m1,1)];
        g = [(0:m1-1)'/max(m1,1); (m-m1-1:-1:0)'/max(m-m1-1,1)];
        b = [ones(m1,1); (m-m1-1:-1:0)'/max(m-m1-1,1)];
    end

    cmap = [r g b];
end

function result = ternary(condition, true_val, false_val)
    % TERNARY OPERATOR - inline conditional value selection
    % Simplified if-else for value assignment
    %
    % INPUTS:
    %   condition - Logical condition to evaluate
    %   true_val  - Value to return if condition is true
    %   false_val - Value to return if condition is false
    %
    % OUTPUT:
    %   result - Either true_val or false_val

    if condition
        result = true_val;
    else
        result = false_val;
    end
end
