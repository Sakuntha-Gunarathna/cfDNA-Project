# cfDNA Cancer Classification Pipeline

**Author:** Sakuntha Devaka Gunarathna  
**Published Date:** July 30, 2025  

This repository contains a machine learning pipeline using **XGBoost** for cancer classification based on cfDNA sequencing data.  
The pipeline integrates data generated from **deepTools multiBigwigSummary** with patient metadata and reference peaks.

---

## Repository Contents

- `general_model_github.py` — Main pipeline script
- `Patient_ID.csv` — Patient identifiers (used as index)
- `Reference_Peaks_ID.csv` — Reference genomic regions (coordinates for cfDNA features)
- `Cancer_TorF.csv` — Cancer status labels (`T` for Cancer, `F` for Healthy)
- `requirements.txt

---

## Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
```

---

## Input Data Format

### Patient_ID.csv
Must have a column `Patient_ID` set as the index.

### Reference_Peaks_ID.csv
Reference genomic regions (must match BED file used in deepTools).

### Cancer_TorF.csv
Cancer labels (`T` = True, `F` = False), Patient_ID must match Patient_ID.csv.

### multiBigwigSummary `.npz` File
Generated with **deepTools multiBigwigSummary** using the **same coordinates** as Reference_Peaks_ID.csv.

Example:

```bash
nohup multiBigwigSummary BED-file   -b <path_to_bigwig_group1>/*.bw <path_to_bigwig_group2>/*.bw ...   -o Output_deepTools_multiBigwigSummary.npz   --outRawCounts Output_counts.tab   --BED Reference_Peaks_ID.bed &
```

---

## Running the Pipeline

```bash
python general_model_github.py
```

Outputs include:
- Model metrics (Balanced Accuracy, AUC-PR, ROC-AUC)
- Confusion matrix
- Feature importance plots
- SHAP analysis

---

## Notes

- Change `device="cuda"` to `device="cpu"` if GPU is not available.

---
