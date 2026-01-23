#!/bin/bash
# ==========================================================================
#  WAHN Phase 0: Data Exploration Script
# ==========================================================================
#  Author: Nikos Diederichs
#  Date: January 2025
#
#  Purpose: Explore WAHN data structure before running main preprocessing
#
#  WICHTIG: Dieses Skript ZUERST ausführen, bevor das MATLAB-Skript läuft!
#           Die Ausgabe an Claude senden, damit er das Hauptskript anpassen kann.
#
#  Usage: bash phase0_exploration.sh > exploration_output.txt 2>&1
# ==========================================================================

echo "=========================================="
echo "WAHN Phase 0: Data Exploration"
echo "Date: $(date)"
echo "=========================================="

# Configuration
WAHN_MRI_BASE="/volume/data/WAHN/MRI/19-Nov-2025/Data/Preprocessed"
WAHN_DATADUMP="/volume/data/WAHN/DataDump/19-Nov-2025"
NESDA_NM_PATH="/volume/projects/CV_NESDA/Data/NESDA_Waves/Wave_1/NM_Structures"

# ==========================================================================
# 0.1: CHECK IF MWP1 FILES EXIST
# ==========================================================================
echo ""
echo "========== 0.1: SEARCHING FOR MWP1 FILES =========="
echo ""

echo "--- Search 1: mwp1*.nii in Preprocessed folder ---"
MWP1_COUNT=$(find "$WAHN_MRI_BASE" -name "mwp1*.nii" 2>/dev/null | wc -l)
echo "Found: $MWP1_COUNT mwp1*.nii files"

if [ "$MWP1_COUNT" -gt 0 ]; then
    echo ""
    echo "First 20 mwp1 files:"
    find "$WAHN_MRI_BASE" -name "mwp1*.nii" 2>/dev/null | head -20
fi

echo ""
echo "--- Search 2: mwp*.nii (broader pattern) ---"
MWP_COUNT=$(find "$WAHN_MRI_BASE" -name "mwp*.nii" 2>/dev/null | wc -l)
echo "Found: $MWP_COUNT mwp*.nii files"

if [ "$MWP_COUNT" -gt 0 ]; then
    echo ""
    echo "First 20 mwp files:"
    find "$WAHN_MRI_BASE" -name "mwp*.nii" 2>/dev/null | head -20
fi

echo ""
echo "--- Search 3: Look for 'mri' subdirectories (CAT12 output) ---"
MRI_DIRS=$(find "$WAHN_MRI_BASE" -type d -name "mri" 2>/dev/null | wc -l)
echo "Found: $MRI_DIRS 'mri' directories"

if [ "$MRI_DIRS" -gt 0 ]; then
    echo ""
    echo "MRI directories:"
    find "$WAHN_MRI_BASE" -type d -name "mri" 2>/dev/null | head -10
fi

echo ""
echo "--- Search 4: Count T1w files as reference ---"
T1W_COUNT=$(find "$WAHN_MRI_BASE" -name "sub-*_T1w.nii" 2>/dev/null | wc -l)
echo "Found: $T1W_COUNT T1w files"

echo ""
echo "--- Search 5: Directory structure overview ---"
echo "Top-level structure in Preprocessed:"
ls -la "$WAHN_MRI_BASE" 2>/dev/null | head -20

# Check if there's a subject subdirectory structure
FIRST_SUBJ_DIR=$(ls -d "$WAHN_MRI_BASE"/sub-* 2>/dev/null | head -1)
if [ -n "$FIRST_SUBJ_DIR" ]; then
    echo ""
    echo "Example subject directory contents ($FIRST_SUBJ_DIR):"
    ls -laR "$FIRST_SUBJ_DIR" 2>/dev/null | head -50
fi

# ==========================================================================
# 0.2: CHECK EXCEL FILES
# ==========================================================================
echo ""
echo "========== 0.2: CHECK EXCEL FILES =========="
echo ""

echo "--- Excel 1 (Demographics): ---"
ls -la "$WAHN_DATADUMP/WAHN_tabular_combined_data.xlsx" 2>/dev/null || echo "NOT FOUND!"

echo ""
echo "--- Excel 2 (CAT12 Metadata): ---"
ls -la "$WAHN_DATADUMP/WAHN_baseline_data_cat12_r1207_tabular.xlsx" 2>/dev/null || echo "NOT FOUND!"

echo ""
echo "--- All files in DataDump: ---"
ls -la "$WAHN_DATADUMP/" 2>/dev/null | head -20

# ==========================================================================
# 0.3: CHECK NESDA REFERENCE FILES
# ==========================================================================
echo ""
echo "========== 0.3: CHECK NESDA REFERENCE FILES =========="
echo ""

echo "--- NM_Pat.mat: ---"
ls -la "$NESDA_NM_PATH/NM_Pat.mat" 2>/dev/null || echo "NOT FOUND!"

echo ""
echo "--- NM_HC.mat: ---"
ls -la "$NESDA_NM_PATH/NM_HC.mat" 2>/dev/null || echo "NOT FOUND!"

echo ""
echo "--- Brain Mask: ---"
ls -la "/volume/projects/CV_NESDA/Analysis/bvFTD/Mask/brainmask_T1_2mm.nii" 2>/dev/null || echo "NOT FOUND!"

# ==========================================================================
# 0.4: CRITICAL ASSESSMENT
# ==========================================================================
echo ""
echo "========== CRITICAL ASSESSMENT =========="
echo ""

if [ "$MWP1_COUNT" -eq 0 ] && [ "$MWP_COUNT" -eq 0 ]; then
    echo "!!! CRITICAL: NO MWP1 FILES FOUND !!!"
    echo ""
    echo "Action required:"
    echo "  1. Contact Maja - mwp1 files may be on her local Mac"
    echo "  2. Excel 2 shows paths like: /Users/maja/Desktop/..."
    echo "  3. Files need to be uploaded to HPC before proceeding"
    echo ""
    echo "STOP: Cannot proceed with MATLAB script until mwp1 files are available."
else
    echo "mwp1 files found. Proceed with MATLAB exploration script."
fi

echo ""
echo "=========================================="
echo "Exploration complete: $(date)"
echo "=========================================="
