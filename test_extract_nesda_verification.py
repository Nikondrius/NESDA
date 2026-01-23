#!/usr/bin/env python3
"""
Comprehensive Test Suite for NESDA Verification Data Extraction Script
======================================================================
Tests edge cases, error handling, data validation, and performance.
"""

import os
import sys
import shutil
import tempfile
import unittest
from datetime import datetime
from io import StringIO
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np

# Import the module to test
import extract_nesda_verification_data as nesda_extract


class TestHelperFunctions(unittest.TestCase):
    """Test helper functions."""

    def test_find_variable_first_match(self):
        """Test that find_variable returns first matching column."""
        df = pd.DataFrame({'pident': [1, 2], 'Age': [30, 40]})
        result = nesda_extract.find_variable(df, ['pident', 'PIDENT', 'ID'])
        self.assertEqual(result, 'pident')

    def test_find_variable_second_match(self):
        """Test fallback to second candidate."""
        df = pd.DataFrame({'PIDENT': [1, 2], 'Age': [30, 40]})
        result = nesda_extract.find_variable(df, ['pident', 'PIDENT', 'ID'])
        self.assertEqual(result, 'PIDENT')

    def test_find_variable_no_match(self):
        """Test returns None when no match found."""
        df = pd.DataFrame({'subject': [1, 2], 'Age': [30, 40]})
        result = nesda_extract.find_variable(df, ['pident', 'PIDENT', 'ID'])
        self.assertIsNone(result)

    def test_find_variable_empty_df(self):
        """Test with empty DataFrame."""
        df = pd.DataFrame()
        result = nesda_extract.find_variable(df, ['pident', 'PIDENT'])
        self.assertIsNone(result)


class TestDataLoading(unittest.TestCase):
    """Test data loading functions."""

    def setUp(self):
        """Create temporary test directory."""
        self.test_dir = tempfile.mkdtemp()
        self.original_base_path = nesda_extract.BASE_PATH
        self.original_data_path = nesda_extract.DATA_PATH
        self.original_waves_path = nesda_extract.WAVES_PATH

        # Override paths for testing
        nesda_extract.BASE_PATH = self.test_dir
        nesda_extract.DATA_PATH = os.path.join(self.test_dir, 'Data/tabular_data/')
        nesda_extract.WAVES_PATH = os.path.join(self.test_dir, 'Data/NESDA_Waves/')
        nesda_extract.NESDA_TABULAR_FILE = os.path.join(nesda_extract.DATA_PATH, 'NESDA_tabular_combined_data.csv')
        nesda_extract.DIAGNOSIS_HC_FILE = os.path.join(nesda_extract.WAVES_PATH, 'Wave_1/DynStd_Preparation/NESDA_HC.csv')
        nesda_extract.DIAGNOSIS_PATIENTS_FILE = os.path.join(nesda_extract.WAVES_PATH, 'Wave_1/DynStd_Preparation/NESDA_Patients.csv')

        # Create directories
        os.makedirs(nesda_extract.DATA_PATH, exist_ok=True)
        os.makedirs(os.path.dirname(nesda_extract.DIAGNOSIS_HC_FILE), exist_ok=True)

    def tearDown(self):
        """Clean up test directory."""
        shutil.rmtree(self.test_dir)
        nesda_extract.BASE_PATH = self.original_base_path
        nesda_extract.DATA_PATH = self.original_data_path
        nesda_extract.WAVES_PATH = self.original_waves_path

    def test_load_tabular_data_success(self):
        """Test successful loading of tabular data."""
        df = pd.DataFrame({'pident': [1, 2, 3], 'Age': [30, 40, 50]})
        df.to_csv(nesda_extract.NESDA_TABULAR_FILE, index=False)

        result = nesda_extract.load_tabular_data()
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 3)

    def test_load_tabular_data_file_not_found(self):
        """Test handling of missing file."""
        # Don't create the file
        result = nesda_extract.load_tabular_data()
        self.assertIsNone(result)

    def test_load_tabular_data_empty_file(self):
        """Test handling of empty CSV."""
        # Create empty file with just headers
        df = pd.DataFrame(columns=['pident', 'Age'])
        df.to_csv(nesda_extract.NESDA_TABULAR_FILE, index=False)

        result = nesda_extract.load_tabular_data()
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 0)

    def test_load_diagnosis_data_both_files(self):
        """Test loading both HC and Patients files."""
        hc_df = pd.DataFrame({'pident': [1, 2], 'diagnosis_group': ['HC', 'HC']})
        hc_df.to_csv(nesda_extract.DIAGNOSIS_HC_FILE, index=False)

        patients_df = pd.DataFrame({'pident': [3, 4], 'diagnosis_group': ['Depression', 'Anxiety']})
        patients_df.to_csv(nesda_extract.DIAGNOSIS_PATIENTS_FILE, index=False)

        result = nesda_extract.load_diagnosis_data()
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 4)

    def test_load_diagnosis_data_only_hc(self):
        """Test loading only HC file when Patients missing."""
        hc_df = pd.DataFrame({'pident': [1, 2], 'diagnosis_group': ['HC', 'HC']})
        hc_df.to_csv(nesda_extract.DIAGNOSIS_HC_FILE, index=False)

        result = nesda_extract.load_diagnosis_data()
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 2)

    def test_load_diagnosis_data_adds_hc_group(self):
        """Test that HC group is added if missing."""
        hc_df = pd.DataFrame({'pident': [1, 2], 'other_col': ['a', 'b']})
        hc_df.to_csv(nesda_extract.DIAGNOSIS_HC_FILE, index=False)

        result = nesda_extract.load_diagnosis_data()
        self.assertIn('diagnosis_group', result.columns)
        self.assertTrue(all(result['diagnosis_group'] == 'HC'))


class TestDataExtraction(unittest.TestCase):
    """Test data extraction logic."""

    def test_extract_basic(self):
        """Test basic extraction with all required columns."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50],
            'Sexe': [1, 2, 1],
            'ascanloc': [1, 2, 3]
        })

        diagnosis_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'diagnosis_group': ['HC', 'Depression', 'Anxiety']
        })

        result = nesda_extract.extract_verification_data(tabular_df, diagnosis_df)

        self.assertIsNotNone(result)
        self.assertEqual(len(result), 3)
        self.assertIn('SubjectID', result.columns)
        self.assertIn('diagnosis_group', result.columns)

    def test_extract_missing_site(self):
        """Test extraction when site variable is missing."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)

        self.assertIsNotNone(result)
        self.assertNotIn('Site', result.columns)

    def test_extract_no_id_variable(self):
        """Test extraction fails gracefully without ID."""
        tabular_df = pd.DataFrame({
            'Age': [30, 40, 50],
            'Sexe': [1, 2, 1]
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)
        self.assertIsNone(result)

    def test_extract_id_type_conversion_numeric_to_string(self):
        """Test ID type conversion from numeric to string."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        diagnosis_df = pd.DataFrame({
            'pident': ['1', '2', '3'],  # String IDs
            'diagnosis_group': ['HC', 'Depression', 'Anxiety']
        })

        result = nesda_extract.extract_verification_data(tabular_df, diagnosis_df)
        self.assertEqual(result['diagnosis_group'].notna().sum(), 3)

    def test_extract_adds_wave_column(self):
        """Test that Wave column is added when not present."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)

        self.assertIn('Wave', result.columns)
        self.assertTrue(all(result['Wave'] == 'Wave_1'))


class TestValidationReport(unittest.TestCase):
    """Test validation report generation."""

    def test_report_counts_diagnoses(self):
        """Test diagnosis group counting."""
        df = pd.DataFrame({
            'SubjectID': [1, 2, 3, 4, 5],
            'diagnosis_group': ['HC', 'HC', 'Depression', 'Anxiety', None]
        })

        report_df, missing_df = nesda_extract.generate_validation_report(df)

        # Check diagnosis counts in report
        diag_rows = report_df[report_df['Category'] == 'Diagnosis']
        self.assertEqual(len(diag_rows), 4)  # HC, Depression, Anxiety, MISSING

    def test_report_detects_duplicates(self):
        """Test duplicate detection."""
        df = pd.DataFrame({
            'SubjectID': [1, 1, 2, 3],  # Duplicate ID
            'diagnosis_group': ['HC', 'Depression', 'Anxiety', 'HC']
        })

        # Capture stdout
        old_stdout = sys.stdout
        sys.stdout = StringIO()

        nesda_extract.generate_validation_report(df)

        output = sys.stdout.getvalue()
        sys.stdout = old_stdout

        self.assertIn('duplicate', output.lower())

    def test_report_missing_values(self):
        """Test missing value detection."""
        df = pd.DataFrame({
            'SubjectID': [1, 2, 3, 4],
            'Site': [1, None, 3, None],  # 2 missing
            'Age': [30, 40, None, 50]    # 1 missing
        })

        report_df, missing_df = nesda_extract.generate_validation_report(df)

        self.assertIsNotNone(missing_df)
        self.assertEqual(len(missing_df), 2)  # Site and Age

    def test_report_no_missing_values(self):
        """Test when no missing values."""
        df = pd.DataFrame({
            'SubjectID': [1, 2, 3],
            'Site': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        report_df, missing_df = nesda_extract.generate_validation_report(df)

        self.assertTrue(missing_df is None or len(missing_df) == 0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and unusual data."""

    def test_single_subject(self):
        """Test with single subject."""
        tabular_df = pd.DataFrame({
            'pident': [1],
            'Age': [30],
            'ascanloc': [1]
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)
        self.assertEqual(len(result), 1)

    def test_all_missing_diagnosis(self):
        """Test when all diagnoses are missing."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        diagnosis_df = pd.DataFrame({
            'pident': [4, 5, 6],  # No matching IDs
            'diagnosis_group': ['HC', 'Depression', 'Anxiety']
        })

        result = nesda_extract.extract_verification_data(tabular_df, diagnosis_df)
        self.assertEqual(result['diagnosis_group'].isna().sum(), 3)

    def test_special_characters_in_data(self):
        """Test handling of special characters."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50],
            'notes': ['Test äöü', 'Test <>&', 'Normal']
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)
        self.assertEqual(len(result), 3)

    def test_very_large_ids(self):
        """Test with very large ID numbers."""
        tabular_df = pd.DataFrame({
            'pident': [9999999999, 10000000000, 10000000001],
            'Age': [30, 40, 50]
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)
        self.assertEqual(len(result), 3)

    def test_negative_values(self):
        """Test handling of negative values (e.g., coding for missing)."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, -99, 50],  # -99 often means missing
            'ascanloc': [1, -1, 3]
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)
        self.assertEqual(len(result), 3)

    def test_mixed_id_formats(self):
        """Test with mixed ID formats in diagnosis data."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        # Diagnosis with string IDs including leading zeros
        diagnosis_df = pd.DataFrame({
            'pident': ['001', '002', '003'],
            'diagnosis_group': ['HC', 'Depression', 'Anxiety']
        })

        result = nesda_extract.extract_verification_data(tabular_df, diagnosis_df)
        # This may not match due to format differences - that's expected behavior
        self.assertIsNotNone(result)


class TestDataIntegrity(unittest.TestCase):
    """Test data integrity and consistency."""

    def test_no_data_loss(self):
        """Verify no subjects are lost during extraction."""
        n_subjects = 100
        tabular_df = pd.DataFrame({
            'pident': range(n_subjects),
            'Age': np.random.randint(18, 80, n_subjects),
            'ascanloc': np.random.choice([1, 2, 3], n_subjects)
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)
        self.assertEqual(len(result), n_subjects)

    def test_correct_merge(self):
        """Verify diagnosis merge is correct."""
        tabular_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        diagnosis_df = pd.DataFrame({
            'pident': [1, 2, 3],
            'diagnosis_group': ['HC', 'Depression', 'Anxiety']
        })

        result = nesda_extract.extract_verification_data(tabular_df, diagnosis_df)

        # Verify correct diagnosis for each subject
        self.assertEqual(result[result['SubjectID'] == '1']['diagnosis_group'].values[0], 'HC')
        self.assertEqual(result[result['SubjectID'] == '2']['diagnosis_group'].values[0], 'Depression')
        self.assertEqual(result[result['SubjectID'] == '3']['diagnosis_group'].values[0], 'Anxiety')

    def test_crosstab_totals(self):
        """Verify crosstab totals are consistent."""
        df = pd.DataFrame({
            'SubjectID': range(100),
            'diagnosis_group': np.random.choice(['HC', 'Depression', 'Anxiety'], 100),
            'Site': np.random.choice([1, 2, 3], 100)
        })

        crosstab = pd.crosstab(df['diagnosis_group'], df['Site'], margins=True)

        # Total should equal number of subjects
        self.assertEqual(crosstab.loc['All', 'All'], 100)


class TestPerformance(unittest.TestCase):
    """Test performance with larger datasets."""

    def test_large_dataset(self):
        """Test with 10,000 subjects."""
        n_subjects = 10000

        tabular_df = pd.DataFrame({
            'pident': range(n_subjects),
            'Age': np.random.randint(18, 80, n_subjects),
            'Sexe': np.random.choice([1, 2], n_subjects),
            'ascanloc': np.random.choice([1, 2, 3], n_subjects),
            'abmi': np.random.normal(25, 4, n_subjects),
            'aedu': np.random.choice([8, 10, 12, 14, 16, 18], n_subjects)
        })

        diagnosis_df = pd.DataFrame({
            'pident': range(n_subjects),
            'diagnosis_group': np.random.choice(['HC', 'Depression', 'Anxiety', 'Comorbid'], n_subjects)
        })

        import time
        start = time.time()
        result = nesda_extract.extract_verification_data(tabular_df, diagnosis_df)
        elapsed = time.time() - start

        self.assertEqual(len(result), n_subjects)
        self.assertLess(elapsed, 5.0, f"Extraction took too long: {elapsed:.2f}s")
        print(f"\n  Performance: 10,000 subjects extracted in {elapsed:.3f}s")


class TestOutputFormat(unittest.TestCase):
    """Test output file format and content."""

    def setUp(self):
        """Create temporary directory for output."""
        self.test_dir = tempfile.mkdtemp()
        self.original_output_path = nesda_extract.OUTPUT_PATH
        nesda_extract.OUTPUT_PATH = self.test_dir + '/'

    def tearDown(self):
        """Clean up."""
        shutil.rmtree(self.test_dir)
        nesda_extract.OUTPUT_PATH = self.original_output_path

    def test_excel_creation(self):
        """Test Excel file is created with correct sheets."""
        df = pd.DataFrame({
            'SubjectID': [1, 2, 3],
            'diagnosis_group': ['HC', 'Depression', 'Anxiety'],
            'Site': [1, 2, 3],
            'Age': [30, 40, 50]
        })

        output_file = os.path.join(self.test_dir, 'test_output.xlsx')

        report_df, missing_df = nesda_extract.generate_validation_report(df)

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Subjects', index=False)
            if not report_df.empty:
                report_df.to_excel(writer, sheet_name='Summary', index=False)

        # Verify file exists and has correct sheets
        self.assertTrue(os.path.exists(output_file))

        excel_file = pd.ExcelFile(output_file)
        self.assertIn('Subjects', excel_file.sheet_names)
        self.assertIn('Summary', excel_file.sheet_names)


class TestRealWorldScenarios(unittest.TestCase):
    """Test realistic NESDA-like data scenarios."""

    def test_nesda_typical_distribution(self):
        """Test with typical NESDA-like distribution."""
        np.random.seed(42)
        n_subjects = 500

        # NESDA typical: ~30% HC, ~25% Depression, ~20% Anxiety, ~15% Comorbid, ~10% missing
        diagnoses = np.random.choice(
            ['HC', 'Depression', 'Anxiety', 'Comorbid', np.nan],
            size=n_subjects,
            p=[0.30, 0.25, 0.20, 0.15, 0.10]
        )

        tabular_df = pd.DataFrame({
            'pident': range(10000, 10000 + n_subjects),
            'Age': np.random.normal(42, 12, n_subjects).clip(18, 75).astype(int),
            'Sexe': np.random.choice([1, 2], n_subjects, p=[0.35, 0.65]),
            'ascanloc': np.random.choice([1, 2, 3], n_subjects, p=[0.4, 0.35, 0.25]),
            'abmi': np.random.normal(25, 4, n_subjects).round(1)
        })

        diagnosis_df = pd.DataFrame({
            'pident': range(10000, 10000 + n_subjects),
            'diagnosis_group': diagnoses
        })

        result = nesda_extract.extract_verification_data(tabular_df, diagnosis_df)

        self.assertEqual(len(result), n_subjects)

        # Verify reasonable distributions
        hc_pct = (result['diagnosis_group'] == 'HC').sum() / len(result) * 100
        self.assertGreater(hc_pct, 20)  # At least 20% HC
        self.assertLess(hc_pct, 40)     # At most 40% HC

    def test_site_imbalance(self):
        """Test with imbalanced site distribution."""
        n_subjects = 300

        # Highly imbalanced: Site 1 has 80% of subjects
        tabular_df = pd.DataFrame({
            'pident': range(n_subjects),
            'Age': np.random.randint(18, 80, n_subjects),
            'ascanloc': np.random.choice([1, 2, 3], n_subjects, p=[0.8, 0.15, 0.05])
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)

        site_counts = result['Site'].value_counts()
        self.assertGreater(site_counts.get(1.0, 0), site_counts.get(2.0, 0))
        self.assertGreater(site_counts.get(1.0, 0), site_counts.get(3.0, 0))

    def test_missing_data_patterns(self):
        """Test various missing data patterns."""
        n_subjects = 100

        tabular_df = pd.DataFrame({
            'pident': range(n_subjects),
            'Age': np.where(np.random.random(n_subjects) < 0.1, np.nan,
                          np.random.randint(18, 80, n_subjects)),
            'Sexe': np.random.choice([1, 2], n_subjects),
            'ascanloc': np.where(np.random.random(n_subjects) < 0.05, np.nan,
                               np.random.choice([1, 2, 3], n_subjects)),
            'abmi': np.where(np.random.random(n_subjects) < 0.15, np.nan,
                           np.random.normal(25, 4, n_subjects))
        })

        result = nesda_extract.extract_verification_data(tabular_df, None)
        report_df, missing_df = nesda_extract.generate_validation_report(result)

        # Should detect missing values
        self.assertIsNotNone(missing_df)
        self.assertGreater(len(missing_df), 0)


def run_all_tests():
    """Run all tests and provide summary."""
    print("\n" + "="*70)
    print("  NESDA VERIFICATION SCRIPT - COMPREHENSIVE TEST SUITE")
    print("="*70 + "\n")

    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    test_classes = [
        TestHelperFunctions,
        TestDataLoading,
        TestDataExtraction,
        TestValidationReport,
        TestEdgeCases,
        TestDataIntegrity,
        TestPerformance,
        TestOutputFormat,
        TestRealWorldScenarios
    ]

    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)

    # Run tests with verbosity
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Print summary
    print("\n" + "="*70)
    print("  TEST SUMMARY")
    print("="*70)
    print(f"  Tests run: {result.testsRun}")
    print(f"  Failures: {len(result.failures)}")
    print(f"  Errors: {len(result.errors)}")
    print(f"  Skipped: {len(result.skipped)}")

    if result.failures:
        print("\n  FAILURES:")
        for test, traceback in result.failures:
            print(f"    - {test}")

    if result.errors:
        print("\n  ERRORS:")
        for test, traceback in result.errors:
            print(f"    - {test}")

    success = len(result.failures) == 0 and len(result.errors) == 0
    print(f"\n  Overall: {'PASSED' if success else 'FAILED'}")
    print("="*70 + "\n")

    return success


if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
