########################################
# Script: 10_export_tiff_figures.py
# Purpose: Convert final locked PNG figures to 600-dpi LZW-compressed TIFF without changing plotted content.
# Input: Six final PNG files in a supplied input directory.
# Output: Fig1.tiff through Fig6.tiff.
# Software: Python + Pillow
# Version: Python 3.12.13; Pillow 12.2.0
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
from pathlib import Path
import sys
from PIL import Image

if len(sys.argv) != 3:
    raise SystemExit("Usage: python 10_export_tiff_figures.py <png_dir> <output_dir>")
png_dir, output_dir = Path(sys.argv[1]), Path(sys.argv[2])
output_dir.mkdir(parents=True, exist_ok=True)
source_names = [
    "Figure1_study_design_evidence_hierarchy_v10.png",
    "Figure2_preprocessing_donor_overview.png",
    "Figure3_GSE48350_naive_vs_donoraware.png",
    "Figure4_AD_dataset_specific_effects.png",
    "Figure5_matched_random_negative_control.png",
    "Figure6_final_evidence_and_liver_preservation_v10.png",
]
for index, name in enumerate(source_names, 1):
    with Image.open(png_dir / name) as image:
        if image.mode not in ("RGB", "L"):
            image = image.convert("RGB")
        image.save(output_dir / f"Fig{index}.tiff", format="TIFF", compression="tiff_lzw", dpi=(600, 600))
