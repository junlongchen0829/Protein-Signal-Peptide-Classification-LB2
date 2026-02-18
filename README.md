# Protein-Signal-Peptide-Classification-LB2
A machine learning project to predict secretory signal peptides using SVM and von Heijne models, featuring rigorous sequence redundancy removal via mmseqs2 and hard-negative analysis.
# Signal Peptide Prediction (LAB2 Module 2)

This project builds and benchmarks two methods to predict **eukaryotic secretory signal peptides (SP)**, with a focus on robustness against **hard negatives** (N-terminal transmembrane proteins).

## Dataset (UniProtKB)
- **Positive set**: Eukaryota, reviewed, non-fragment, length ≥ 40 aa, protein-level evidence, experimentally supported signal peptide feature.
- **Negative-A**: intracellular proteins (cytoplasm or nucleus), no signal peptide annotation.
- **Negative-B (hard negatives)**: transmembrane proteins, no signal peptide annotation.

Redundancy was removed using **MMseqs2** clustering at **30% sequence identity** and **40% coverage**.

## Models
- **Baseline**: von Heijne statistical model (PSSM/PSWM)
- **SVM**: RBF-kernel SVM trained on N-terminal features (amino-acid composition + mean hydrophobicity)

## Evaluation
Hyperparameters (C, gamma) were optimized with **5-fold stratified cross-validation** using **MCC** on the training set.
Final evaluation was performed once on an independent blind test set.

Key outputs:
- `results/tables/svm_metrics.csv`
- `results/figures/svm_roc.png`, `svm_pr.png`, `svm_cm.png`
- `results/tables/false_positives.csv`, `false_negatives.csv`
- Final report PDF in `report/`

## How to run
```bash
python lab2_clean_pipeline.py
