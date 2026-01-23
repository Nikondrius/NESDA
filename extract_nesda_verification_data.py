#!/usr/bin/env python3
"""
NESDA Verification Data Extraction Script
==========================================
Author: Generated for VUMC Amsterdam verification
Purpose: Extract and validate NESDA subject metadata for Laura Han

This script extracts:
- SubjectID (pident)
- Diagnosis/Patient Group (diagnosis_group)
- Site (ascanloc)
- Timepoint/Wave
- Additional relevant variables
- MRI image paths (if available)

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

# Key variables to extract
ID_VARS = ['pident', 'PIDENT', 'ID', 'SubjectID', 'subject_id']
SITE_VARS = ['ascanloc', 'site', 'Site', 'SITE', 'scan_location']
DIAGNOSIS_VARS = ['diagnosis_group', 'Diagnosis', 'diagnosis', 'acontrol']
DEMOGRAPHIC_VARS = ['Age', 'Sexe', 'abmi', 'aedu', 'amarpart']
WAVE_VARS = ['wave', 'Wave', 'WAVE', 'timepoint', 'Timepoint', 'assessment']


def print_header(text):
    """Print formatted section header."""
    print(f"\n{'='*60}")
    print(f"  {text}")
    print(f"{'='*60}\n")


def find_variable(df, var_candidates):
    """Find first matching variable name from candidates."""
    for var in var_candidates:
        if var in df.columns:
            return var
    return None


def load_tabular_data():
    """Load main NESDA tabular data file."""
    print_header("LOADING NESDA TABULAR DATA")

    if not os.path.exists(NESDA_TABULAR_FILE):
        print(f"  ERROR: File not found: {NESDA_TABULAR_FILE}")
        return None

    print(f"  Loading: {NESDA_TABULAR_FILE}")

    try:
        df = pd.read_csv(NESDA_TABULAR_FILE, low_memory=False)
        print(f"  SUCCESS: Loaded {len(df)} rows x {len(df.columns)} columns")
        print(f"  First 5 columns: {list(df.columns[:5])}")
        return df
    except Exception as e:
        print(f"  ERROR loading file: {e}")
        return None


def load_diagnosis_data():
    """Load HC and Patients diagnosis files and combine."""
    print_header("LOADING DIAGNOSIS DATA")

    combined_diag = []

    # Load HC data
    if os.path.exists(DIAGNOSIS_HC_FILE):
        print(f"  Loading HC: {DIAGNOSIS_HC_FILE}")
        try:
            hc_df = pd.read_csv(DIAGNOSIS_HC_FILE, low_memory=False)
            if 'diagnosis_group' not in hc_df.columns:
                hc_df['diagnosis_group'] = 'HC'
            combined_diag.append(hc_df)
            print(f"    HC subjects: {len(hc_df)}")
        except Exception as e:
            print(f"    ERROR loading HC: {e}")
    else:
        print(f"  WARNING: HC file not found: {DIAGNOSIS_HC_FILE}")

    # Load Patients data
    if os.path.exists(DIAGNOSIS_PATIENTS_FILE):
        print(f"  Loading Patients: {DIAGNOSIS_PATIENTS_FILE}")
        try:
            patients_df = pd.read_csv(DIAGNOSIS_PATIENTS_FILE, low_memory=False)
            combined_diag.append(patients_df)
            print(f"    Patient subjects: {len(patients_df)}")

            # Show diagnosis groups
            if 'diagnosis_group' in patients_df.columns:
                print(f"    Diagnosis groups: {patients_df['diagnosis_group'].unique().tolist()}")
        except Exception as e:
            print(f"    ERROR loading Patients: {e}")
    else:
        print(f"  WARNING: Patients file not found: {DIAGNOSIS_PATIENTS_FILE}")

    if combined_diag:
        # Find common columns
        common_cols = set(combined_diag[0].columns)
        for df in combined_diag[1:]:
            common_cols &= set(df.columns)

        # Combine with common columns
        combined_df = pd.concat([df[list(common_cols)] for df in combined_diag], ignore_index=True)
        print(f"  Combined diagnosis data: {len(combined_df)} subjects")
        return combined_df

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

    # Also check for wave variable in data path
    if os.path.exists(DATA_PATH):
        print(f"\n  Scanning data path for wave info: {DATA_PATH}")
        for f in os.listdir(DATA_PATH):
            if 'wave' in f.lower() or 'timepoint' in f.lower():
                print(f"    Found file: {f}")

    return waves_found


def scan_for_mri_paths():
    """Scan for MRI image paths."""
    print_header("SCANNING FOR MRI PATHS")

    mri_paths = []
    potential_mri_dirs = [
        os.path.join(BASE_PATH, 'Data/MRI'),
        os.path.join(BASE_PATH, 'Data/Imaging'),
        os.path.join(BASE_PATH, 'Data/NIFTI'),
        os.path.join(BASE_PATH, 'MRI'),
        os.path.join(BASE_PATH, 'Imaging'),
    ]

    for mri_dir in potential_mri_dirs:
        if os.path.exists(mri_dir):
            mri_paths.append(mri_dir)
            print(f"  Found MRI directory: {mri_dir}")
            # List subdirectories
            try:
                subdirs = os.listdir(mri_dir)[:5]
                print(f"    Sample contents: {subdirs}")
            except:
                pass

    if not mri_paths:
        print("  No MRI directories found in expected locations")

    return mri_paths


def extract_verification_data(tabular_df, diagnosis_df):
    """Extract and merge verification data."""
    print_header("EXTRACTING VERIFICATION DATA")

    # Find ID variable
    id_var = find_variable(tabular_df, ID_VARS)
    if not id_var:
        print("  ERROR: No ID variable found!")
        print(f"    Searched for: {ID_VARS}")
        print(f"    Available columns: {list(tabular_df.columns[:20])}")
        return None
    print(f"  ID variable: {id_var}")

    # Find Site variable
    site_var = find_variable(tabular_df, SITE_VARS)
    if site_var:
        print(f"  Site variable: {site_var}")
    else:
        print(f"  WARNING: No site variable found. Searched: {SITE_VARS}")

    # Find Wave/Timepoint variable
    wave_var = find_variable(tabular_df, WAVE_VARS)
    if wave_var:
        print(f"  Wave variable: {wave_var}")
    else:
        print(f"  INFO: No wave variable in tabular data (may be single timepoint)")

    # Start building extraction dataframe
    extraction_cols = [id_var]

    # Add diagnosis from tabular if present
    diag_var_tabular = find_variable(tabular_df, DIAGNOSIS_VARS)
    if diag_var_tabular:
        extraction_cols.append(diag_var_tabular)
        print(f"  Diagnosis variable (tabular): {diag_var_tabular}")

    # Add site
    if site_var:
        extraction_cols.append(site_var)

    # Add wave
    if wave_var:
        extraction_cols.append(wave_var)

    # Add demographics
    for dem_var in DEMOGRAPHIC_VARS:
        if dem_var in tabular_df.columns:
            extraction_cols.append(dem_var)

    # Additional potentially useful variables
    additional_vars = ['acontrol', 'acidep10', 'acidep11', 'aanxy21', 'aanxy22',
                       'aLCAsubtype', 'ANDPBOXSX']
    for add_var in additional_vars:
        if add_var in tabular_df.columns and add_var not in extraction_cols:
            extraction_cols.append(add_var)

    # Extract subset
    extract_df = tabular_df[extraction_cols].copy()
    extract_df = extract_df.rename(columns={id_var: 'SubjectID'})

    print(f"\n  Extracted {len(extract_df)} subjects with {len(extraction_cols)} variables")

    # Merge with diagnosis data if available
    if diagnosis_df is not None and 'pident' in diagnosis_df.columns:
        print("\n  Merging with diagnosis data...")

        # Ensure ID types match
        extract_df['SubjectID'] = extract_df['SubjectID'].astype(str)
        diagnosis_df['pident'] = diagnosis_df['pident'].astype(str)

        # Merge
        if 'diagnosis_group' in diagnosis_df.columns:
            diag_subset = diagnosis_df[['pident', 'diagnosis_group']].drop_duplicates()
            extract_df = extract_df.merge(
                diag_subset,
                left_on='SubjectID',
                right_on='pident',
                how='left'
            )
            # Clean up
            if 'pident' in extract_df.columns:
                extract_df = extract_df.drop(columns=['pident'])

            print(f"    Merged diagnosis_group for {extract_df['diagnosis_group'].notna().sum()} subjects")

    # Rename site variable for clarity
    if site_var and site_var in extract_df.columns:
        extract_df = extract_df.rename(columns={site_var: 'Site'})

    # Add Wave column if not present (assume Wave 1 based on file paths)
    if 'Wave' not in extract_df.columns and wave_var is None:
        extract_df['Wave'] = 'Wave_1'
        print("  Added Wave column (defaulted to Wave_1 based on file paths)")

    return extract_df


def generate_validation_report(df):
    """Generate validation statistics and flags."""
    print_header("VALIDATION REPORT")

    report_data = []

    # 1. Sample size summary
    print("1. SAMPLE SIZE SUMMARY")
    print(f"   Total subjects: {len(df)}")
    print(f"   Unique SubjectIDs: {df['SubjectID'].nunique()}")

    # 2. Diagnosis group distribution
    if 'diagnosis_group' in df.columns:
        print("\n2. DIAGNOSIS GROUP DISTRIBUTION")
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

    # 3. Site distribution
    if 'Site' in df.columns:
        print("\n3. SITE DISTRIBUTION")
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

    # 4. Wave/Timepoint distribution
    if 'Wave' in df.columns:
        print("\n4. WAVE/TIMEPOINT DISTRIBUTION")
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

    # 5. Cross-tabulation: Diagnosis x Site
    if 'diagnosis_group' in df.columns and 'Site' in df.columns:
        print("\n5. CROSS-TABULATION: DIAGNOSIS x SITE")
        crosstab = pd.crosstab(df['diagnosis_group'], df['Site'], margins=True, margins_name='Total')
        print(crosstab.to_string())

    # 6. Missing values
    print("\n6. MISSING VALUES")
    missing_summary = []
    for col in df.columns:
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
        print("   No missing values detected!")

    # 7. Duplicates
    print("\n7. DUPLICATE CHECK")
    n_duplicates = df['SubjectID'].duplicated().sum()
    if n_duplicates > 0:
        print(f"   WARNING: {n_duplicates} duplicate SubjectIDs found!")
        dup_ids = df[df['SubjectID'].duplicated(keep=False)]['SubjectID'].unique()
        print(f"   Duplicate IDs: {list(dup_ids[:10])}{'...' if len(dup_ids) > 10 else ''}")
    else:
        print("   No duplicate SubjectIDs found")

    # 8. Demographics summary
    print("\n8. DEMOGRAPHICS SUMMARY")
    if 'Age' in df.columns:
        age_valid = df['Age'].dropna()
        print(f"   Age: Mean={age_valid.mean():.1f}, SD={age_valid.std():.1f}, Range=[{age_valid.min():.0f}-{age_valid.max():.0f}]")

    if 'Sexe' in df.columns:
        sex_counts = df['Sexe'].value_counts(dropna=False)
        print(f"   Sex distribution: {dict(sex_counts)}")
        print("   (Note: 1=Male, 2=Female per NESDA coding)")

    # 9. Consistency checks
    print("\n9. CONSISTENCY CHECKS")

    # Check acontrol vs diagnosis_group consistency
    if 'acontrol' in df.columns and 'diagnosis_group' in df.columns:
        print("   Checking acontrol vs diagnosis_group consistency...")
        inconsistent = df[(df['acontrol'].notna()) & (df['diagnosis_group'].notna())]
        # acontrol coding: typically categorical comorbidity status
        print(f"   Subjects with both variables: {len(inconsistent)}")

    return pd.DataFrame(report_data), pd.DataFrame(missing_summary) if missing_summary else None


def main():
    """Main execution function."""
    print("\n" + "="*60)
    print("  NESDA VERIFICATION DATA EXTRACTION")
    print("  For: Laura Han (VUMC Amsterdam)")
    print(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("="*60)

    # Create output directory
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_PATH}")

    # Load data
    tabular_df = load_tabular_data()
    if tabular_df is None:
        print("\nFATAL: Could not load tabular data. Exiting.")
        sys.exit(1)

    diagnosis_df = load_diagnosis_data()

    # Find wave directories
    waves = find_wave_directories()

    # Scan for MRI paths
    mri_paths = scan_for_mri_paths()

    # Extract verification data
    verification_df = extract_verification_data(tabular_df, diagnosis_df)

    if verification_df is None:
        print("\nFATAL: Could not extract verification data. Exiting.")
        sys.exit(1)

    # Generate validation report
    report_df, missing_df = generate_validation_report(verification_df)

    # Save to Excel with multiple sheets
    print_header("SAVING OUTPUT")

    output_file = os.path.join(OUTPUT_PATH, OUTPUT_EXCEL)

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
                    'Source_File',
                    'Total_Subjects',
                    'Unique_SubjectIDs',
                    'Waves_Found',
                    'MRI_Paths_Found',
                    'Script_Version'
                ],
                'Value': [
                    datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    NESDA_TABULAR_FILE,
                    len(verification_df),
                    verification_df['SubjectID'].nunique(),
                    ', '.join(waves) if waves else 'None',
                    ', '.join(mri_paths) if mri_paths else 'None',
                    '1.0'
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
