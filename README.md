# Protein Signal Peptide Classification (LAB2 Module 2)

This project benchmarks two methods to predict **eukaryotic secretory signal peptides (SP)**, with an emphasis on robustness against **hard negatives** (N-terminal transmembrane proteins).

## Dataset (UniProtKB)
- **Positive**: Eukaryota, reviewed, non-fragment, length ≥ 40 aa, protein-level evidence, experimentally supported SP feature.
- **Negative-A**: intracellular proteins (cytoplasm/nucleus), no SP annotation.
- **Negative-B (hard negatives)**: transmembrane proteins, no SP annotation.

Redundancy removal: **MMseqs2** clustering at **30% identity** and **40% coverage**.

## Models
- **Baseline**: von Heijne statistical model (PSSM/PSWM)
- **SVM**: RBF-kernel SVM with N-terminal features (AA composition + mean hydrophobicity)

## Evaluation
Hyperparameters (C, gamma) were optimized by **GridSearchCV** with **5-fold stratified CV** using **MCC** on the training set.
Final evaluation was performed once on an independent blind test set.
## Key results (blind test set)
SVM (RBF, GridSearchCV 5-fold MCC): **MCC = 0.7292**, **ROC-AUC = 0.9447**, **Accuracy = 0.8754**  
Best params: C = 1, gamma = 0.01 (CV MCC = 0.7396)

## Repository structure
- `report/`: final PDF report
- `results/`: figures and tables (ROC/PR/CM, metrics, FP/FN lists)
- `lab2_clean_pipeline.py`: end-to-end pipeline script

pip install -r requirements.txt
## How to run
Place `train_df.csv` and `test_df.csv` under `data/processed/` (see script defaults), then run:
```bash
python lab2_clean_pipeline.py
