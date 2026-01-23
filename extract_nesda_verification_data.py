#!/usr/bin/env python3
"""
NESDA Verification Data Extraction Script
==========================================
Author: Generated for VUMC Amsterdam verification
Purpose: Extract and validate NESDA subject metadata for Laura Han

This script extracts Wave 1 subjects (following NESDA_Clinical_Associations.m logic):
- SubjectID (pident)
- Diagnosis/Patient Group (diagnosis_group)
- Site
- MRI image paths
- Additional clinical variables from tabular data

Output: nesda_verification_subjects.xlsx in the standard output directory
"""

import os
import sys
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
from pathlib import Path

# =============================================================================
# CONFIGURATION - Based on existing MATLAB script paths
# =============================================================================
BASE_PATH = '/volume/projects/CV_NESDA/'
DATA_PATH = os.path.join(BASE_PATH, 'Data/tabular_data/')
WAVES_PATH = os.path.join(BASE_PATH, 'Data/NESDA_Waves/')

# Input files (from MATLAB script)
NESDA_TABULAR_FILE = os.path.join(DATA_PATH, 'NESDA_tabular_combined_data.csv')
DIAGNOSIS_HC_FILE = os.path.join(WAVES_PATH, 'Wave_1/DynStd_Preparation/NESDA_HC.csv')
DIAGNOSIS_PATIENTS_FILE = os.path.join(WAVES_PATH, 'Wave_1/DynStd_Preparation/NESDA_Patients.csv')

# Output directory (same pattern as MATLAB script)
TIMESTAMP = datetime.now().strftime('%Y-%m-%d_%H-%M')
OUTPUT_PATH = os.path.join(BASE_PATH, f'Analysis/NESDA_Verification_{TIMESTAMP}/')
OUTPUT_EXCEL = 'nesda_verification_subjects.xlsx'

# Additional clinical variables to extract from tabular data
CLINICAL_VARS_TO_ADD = [
    'aids', 'abaiscal', 'aauditsc',  # Symptom severity
    'acidep10', 'acidep11', 'acidep13',  # Depression history
    'aanxy21', 'aanxy22',  # Anxiety history
    'acontrol',  # Comorbidity status
    'abmi', 'aedu', 'amarpart',  # Demographics
    'ACTI_total', 'ACLEI',  # Childhood adversity
]


def print_header(text):
    """Print formatted section header."""
    print(f"\n{'='*60}")
    print(f"  {text}")
    print(f"{'='*60}\n")


def load_wave1_subjects():
    """Load Wave 1 subjects from HC and Patients files (primary data source)."""
    print_header("LOADING WAVE 1 SUBJECTS")

    all_subjects = []

    # Load HC data
    if os.path.exists(DIAGNOSIS_HC_FILE):
        print(f"  Loading HC: {DIAGNOSIS_HC_FILE}")
        try:
            hc_df = pd.read_csv(DIAGNOSIS_HC_FILE)
            hc_df['subject_type'] = 'HC'
            all_subjects.append(hc_df)
            print(f"    HC subjects: {len(hc_df)}")
        except Exception as e:
            print(f"    ERROR loading HC: {e}")
    else:
        print(f"  ERROR: HC file not found: {DIAGNOSIS_HC_FILE}")
        return None

    # Load Patients data
    if os.path.exists(DIAGNOSIS_PATIENTS_FILE):
        print(f"  Loading Patients: {DIAGNOSIS_PATIENTS_FILE}")
        try:
            patients_df = pd.read_csv(DIAGNOSIS_PATIENTS_FILE)
            patients_df['subject_type'] = 'Patient'
            all_subjects.append(patients_df)
            print(f"    Patient subjects: {len(patients_df)}")

            if 'diagnosis_group' in patients_df.columns:
                diag_counts = patients_df['diagnosis_group'].value_counts()
                print(f"    Diagnosis breakdown:")
                for diag, count in diag_counts.items():
                    print(f"      {diag}: {count}")
        except Exception as e:
            print(f"    ERROR loading Patients: {e}")
    else:
        print(f"  ERROR: Patients file not found: {DIAGNOSIS_PATIENTS_FILE}")
        return None

    # Combine
    if all_subjects:
        combined_df = pd.concat(all_subjects, ignore_index=True)
        print(f"\n  Total Wave 1 subjects: {len(combined_df)}")
        print(f"  Columns: {list(combined_df.columns)}")
        return combined_df

    return None


def load_tabular_data():
    """Load additional clinical variables from tabular data."""
    print_header("LOADING ADDITIONAL CLINICAL DATA")

    if not os.path.exists(NESDA_TABULAR_FILE):
        print(f"  WARNING: Tabular file not found: {NESDA_TABULAR_FILE}")
        return None

    print(f"  Loading: {NESDA_TABULAR_FILE}")

    try:
        df = pd.read_csv(NESDA_TABULAR_FILE, low_memory=False)
        print(f"  Loaded {len(df)} rows x {len(df.columns)} columns")

        # Check which clinical variables are available
        available_vars = ['pident']  # Always need ID
        for var in CLINICAL_VARS_TO_ADD:
            if var in df.columns:
                available_vars.append(var)

        print(f"  Available clinical variables: {len(available_vars)-1}/{len(CLINICAL_VARS_TO_ADD)}")

        # Return only needed columns
        return df[available_vars]
    except Exception as e:
        print(f"  ERROR loading tabular data: {e}")
        return None


def find_wave_directories():
    """Scan for available wave directories."""
    print_header("SCANNING FOR WAVE/TIMEPOINT DATA")

    waves_found = []

    if os.path.exists(WAVES_PATH):
        print(f"  Scanning: {WAVES_PATH}")
        for item in sorted(os.listdir(WAVES_PATH)):
            item_path = os.path.join(WAVES_PATH, item)
            if os.path.isdir(item_path):
                waves_found.append(item)
                print(f"    Found: {item}")
    else:
        print(f"  WARNING: Waves path not found: {WAVES_PATH}")

    return waves_found


def extract_verification_data(wave1_df, tabular_df):
    """Combine Wave 1 subjects with additional clinical data."""
    print_header("EXTRACTING VERIFICATION DATA")

    if wave1_df is None:
        print("  ERROR: No Wave 1 data available!")
        return None

    # Start with Wave 1 data as base
    result_df = wave1_df.copy()

    # Rename columns for clarity
    column_mapping = {
        'pident': 'SubjectID',
        'path': 'MRI_Path',
        'age': 'Age',
        'sex': 'Sex',
        'site': 'Site',
    }

    for old_name, new_name in column_mapping.items():
        if old_name in result_df.columns:
            result_df = result_df.rename(columns={old_name: new_name})

    # Add Wave column
    result_df['Wave'] = 'Wave_1'

    print(f"  Base data: {len(result_df)} subjects")
    print(f"  Columns: {list(result_df.columns)}")

    # Merge with tabular data for additional clinical variables
    if tabular_df is not None:
        print("\n  Merging additional clinical variables...")

        # Ensure ID types match
        result_df['SubjectID'] = result_df['SubjectID'].astype(str).str.strip()
        tabular_df['pident'] = tabular_df['pident'].astype(str).str.strip()

        # Merge
        n_before = len(result_df)
        result_df = result_df.merge(
            tabular_df,
            left_on='SubjectID',
            right_on='pident',
            how='left'
        )

        # Remove duplicate pident column
        if 'pident' in result_df.columns:
            result_df = result_df.drop(columns=['pident'])

        # Count successful merges
        clinical_cols = [c for c in CLINICAL_VARS_TO_ADD if c in result_df.columns]
        if clinical_cols:
            n_with_clinical = result_df[clinical_cols[0]].notna().sum()
            print(f"    Merged clinical data for {n_with_clinical}/{len(result_df)} subjects")
            print(f"    Added variables: {clinical_cols}")

    # Check MRI paths
    if 'MRI_Path' in result_df.columns:
        n_with_mri = result_df['MRI_Path'].notna().sum()
        print(f"\n  MRI paths available: {n_with_mri}/{len(result_df)}")

        # Check if paths exist (sample)
        sample_paths = result_df['MRI_Path'].dropna().head(3)
        print(f"  Sample MRI paths:")
        for p in sample_paths:
            exists = os.path.exists(p) if pd.notna(p) else False
            status = "EXISTS" if exists else "NOT FOUND"
            print(f"    [{status}] {p}")

    return result_df


def generate_validation_report(df):
    """Generate validation statistics and flags."""
    print_header("VALIDATION REPORT")

    report_data = []

    # 1. Sample size summary
    print("1. SAMPLE SIZE SUMMARY")
    print(f"   Total subjects: {len(df)}")
    print(f"   Unique SubjectIDs: {df['SubjectID'].nunique()}")

    # 2. Subject type (HC vs Patient)
    if 'subject_type' in df.columns:
        print("\n2. SUBJECT TYPE")
        type_counts = df['subject_type'].value_counts()
        for stype, count in type_counts.items():
            print(f"   {stype}: {count} ({100*count/len(df):.1f}%)")
            report_data.append({
                'Category': 'Subject_Type',
                'Value': stype,
                'N': count,
                'Percent': f"{100*count/len(df):.1f}%"
            })

    # 3. Diagnosis group distribution
    if 'diagnosis_group' in df.columns:
        print("\n3. DIAGNOSIS GROUP DISTRIBUTION")
        diag_counts = df['diagnosis_group'].value_counts(dropna=False)
        for diag, count in diag_counts.items():
            diag_label = diag if pd.notna(diag) else 'MISSING'
            print(f"   {diag_label}: {count} ({100*count/len(df):.1f}%)")
            report_data.append({
                'Category': 'Diagnosis',
                'Value': diag_label,
                'N': count,
                'Percent': f"{100*count/len(df):.1f}%"
            })

    # 4. Site distribution
    if 'Site' in df.columns:
        print("\n4. SITE DISTRIBUTION")
        site_counts = df['Site'].value_counts(dropna=False)
        for site, count in site_counts.items():
            site_label = site if pd.notna(site) else 'MISSING'
            print(f"   Site {site_label}: {count} ({100*count/len(df):.1f}%)")
            report_data.append({
                'Category': 'Site',
                'Value': str(site_label),
                'N': count,
                'Percent': f"{100*count/len(df):.1f}%"
            })

    # 5. Wave/Timepoint distribution
    if 'Wave' in df.columns:
        print("\n5. WAVE/TIMEPOINT DISTRIBUTION")
        wave_counts = df['Wave'].value_counts(dropna=False)
        for wave, count in wave_counts.items():
            wave_label = wave if pd.notna(wave) else 'MISSING'
            print(f"   {wave_label}: {count} ({100*count/len(df):.1f}%)")
            report_data.append({
                'Category': 'Wave',
                'Value': str(wave_label),
                'N': count,
                'Percent': f"{100*count/len(df):.1f}%"
            })

    # 6. Cross-tabulation: Diagnosis x Site
    if 'diagnosis_group' in df.columns and 'Site' in df.columns:
        print("\n6. CROSS-TABULATION: DIAGNOSIS x SITE")
        crosstab = pd.crosstab(df['diagnosis_group'], df['Site'], margins=True, margins_name='Total')
        print(crosstab.to_string())

    # 7. MRI availability
    if 'MRI_Path' in df.columns:
        print("\n7. MRI DATA AVAILABILITY")
        n_with_mri = df['MRI_Path'].notna().sum()
        print(f"   Subjects with MRI path: {n_with_mri} ({100*n_with_mri/len(df):.1f}%)")

        # Check by diagnosis group
        if 'diagnosis_group' in df.columns:
            print("   MRI by diagnosis:")
            for diag in df['diagnosis_group'].unique():
                if pd.notna(diag):
                    n_diag = len(df[df['diagnosis_group'] == diag])
                    n_mri = df[(df['diagnosis_group'] == diag) & df['MRI_Path'].notna()].shape[0]
                    print(f"     {diag}: {n_mri}/{n_diag}")

    # 8. Missing values
    print("\n8. MISSING VALUES")
    missing_summary = []
    key_cols = ['SubjectID', 'diagnosis_group', 'Site', 'Age', 'Sex', 'MRI_Path']
    for col in key_cols:
        if col in df.columns:
            n_missing = df[col].isna().sum()
            if n_missing > 0:
                pct_missing = 100 * n_missing / len(df)
                print(f"   {col}: {n_missing} missing ({pct_missing:.1f}%)")
                missing_summary.append({
                    'Variable': col,
                    'N_Missing': n_missing,
                    'Pct_Missing': f"{pct_missing:.1f}%"
                })

    if not missing_summary:
        print("   No missing values in key columns!")

    # 9. Duplicates
    print("\n9. DUPLICATE CHECK")
    n_duplicates = df['SubjectID'].duplicated().sum()
    if n_duplicates > 0:
        print(f"   WARNING: {n_duplicates} duplicate SubjectIDs found!")
        dup_ids = df[df['SubjectID'].duplicated(keep=False)]['SubjectID'].unique()
        print(f"   Duplicate IDs: {list(dup_ids[:10])}{'...' if len(dup_ids) > 10 else ''}")
    else:
        print("   No duplicate SubjectIDs found")

    # 10. Demographics summary
    print("\n10. DEMOGRAPHICS SUMMARY")
    if 'Age' in df.columns:
        age_valid = df['Age'].dropna()
        if len(age_valid) > 0:
            print(f"    Age: Mean={age_valid.mean():.1f}, SD={age_valid.std():.1f}, Range=[{age_valid.min():.0f}-{age_valid.max():.0f}]")

    if 'Sex' in df.columns:
        sex_counts = df['Sex'].value_counts(dropna=False)
        n_male = sex_counts.get(1, 0)
        n_female = sex_counts.get(2, 0)
        print(f"    Sex: Male={n_male}, Female={n_female}")
        print("    (Note: 1=Male, 2=Female per NESDA coding)")

    return pd.DataFrame(report_data), pd.DataFrame(missing_summary) if missing_summary else None


def main():
    """Main execution function."""
    print("\n" + "="*60)
    print("  NESDA VERIFICATION DATA EXTRACTION")
    print("  Wave 1 Subjects Only (matching Clinical Associations script)")
    print("  For: Laura Han (VUMC Amsterdam)")
    print(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("="*60)

    # Create output directory
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_PATH}")

    # Load Wave 1 subjects (PRIMARY DATA SOURCE)
    wave1_df = load_wave1_subjects()
    if wave1_df is None:
        print("\nFATAL: Could not load Wave 1 subject data. Exiting.")
        sys.exit(1)

    # Load additional clinical data from tabular file
    tabular_df = load_tabular_data()

    # Find wave directories (informational)
    waves = find_wave_directories()

    # Extract and combine verification data
    verification_df = extract_verification_data(wave1_df, tabular_df)

    if verification_df is None:
        print("\nFATAL: Could not extract verification data. Exiting.")
        sys.exit(1)

    # Generate validation report
    report_df, missing_df = generate_validation_report(verification_df)

    # Save to Excel with multiple sheets
    print_header("SAVING OUTPUT")

    output_file = os.path.join(OUTPUT_PATH, OUTPUT_EXCEL)

    # Count MRI availability
    n_with_mri = verification_df['MRI_Path'].notna().sum() if 'MRI_Path' in verification_df.columns else 0

    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Main data sheet
            verification_df.to_excel(writer, sheet_name='Subjects', index=False)
            print(f"  Sheet 'Subjects': {len(verification_df)} rows")

            # Summary statistics sheet
            if not report_df.empty:
                report_df.to_excel(writer, sheet_name='Summary', index=False)
                print(f"  Sheet 'Summary': {len(report_df)} rows")

            # Missing values sheet
            if missing_df is not None and not missing_df.empty:
                missing_df.to_excel(writer, sheet_name='Missing_Values', index=False)
                print(f"  Sheet 'Missing_Values': {len(missing_df)} rows")

            # Cross-tabulation sheet
            if 'diagnosis_group' in verification_df.columns and 'Site' in verification_df.columns:
                crosstab = pd.crosstab(
                    verification_df['diagnosis_group'],
                    verification_df['Site'],
                    margins=True,
                    margins_name='Total'
                )
                crosstab.to_excel(writer, sheet_name='Diagnosis_x_Site')
                print(f"  Sheet 'Diagnosis_x_Site': Cross-tabulation")

            # Metadata sheet
            metadata = pd.DataFrame({
                'Parameter': [
                    'Extraction_Date',
                    'HC_Source_File',
                    'Patients_Source_File',
                    'Total_Subjects',
                    'HC_Count',
                    'Patient_Count',
                    'Subjects_with_MRI',
                    'Waves_Available',
                    'Script_Version'
                ],
                'Value': [
                    datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    DIAGNOSIS_HC_FILE,
                    DIAGNOSIS_PATIENTS_FILE,
                    len(verification_df),
                    len(verification_df[verification_df['subject_type'] == 'HC']) if 'subject_type' in verification_df.columns else 'N/A',
                    len(verification_df[verification_df['subject_type'] == 'Patient']) if 'subject_type' in verification_df.columns else 'N/A',
                    n_with_mri,
                    ', '.join(waves) if waves else 'Wave_1',
                    '2.0'
                ]
            })
            metadata.to_excel(writer, sheet_name='Metadata', index=False)
            print(f"  Sheet 'Metadata': Extraction parameters")

        print(f"\n  OUTPUT SAVED: {output_file}")

    except Exception as e:
        print(f"\n  ERROR saving Excel: {e}")

        # Fallback: save as CSV
        csv_file = output_file.replace('.xlsx', '.csv')
        verification_df.to_csv(csv_file, index=False)
        print(f"  FALLBACK: Saved as CSV: {csv_file}")

    # Final summary
    print_header("EXTRACTION COMPLETE")
    print(f"  Total subjects extracted: {len(verification_df)}")
    print(f"  Output file: {output_file}")
    print("\n  Please review the validation report above for any issues.")
    print("  Send nesda_verification_subjects.xlsx to Laura Han (VUMC Amsterdam)")

    return verification_df


if __name__ == '__main__':
    df = main()
