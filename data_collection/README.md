# Data collection (UniProtKB)

**Objective:** retrieve positive and negative datasets of **eukaryotic proteins** for signal peptide prediction.

## A. Selection criteria

### Positive set (secretory proteins with experimentally supported SP)
1. **Taxonomy:** Eukaryota (taxon_id: 2759)
2. **Reviewed:** Yes (reviewed:true)
3. **No fragments:** fragment:false
4. **Length filter:** sequence length ≥ 40 aa
5. **Protein existence:** evidence at protein level (existence:1)
6. **Signal peptide:** presence of SIGNAL feature with experimental evidence (**ECO:0000269**)
7. **Additional filter:** signal peptide length > 14 aa (custom post-filter)

### Negative-A (intracellular proteins)
1. **Taxonomy:** Eukaryota (taxon_id: 2759)
2. **Reviewed:** Yes
3. **No fragments**
4. **Subcellular location:** Cytoplasm OR Nucleus
5. **No signal peptide annotation**

### Negative-B (hard negatives: transmembrane proteins)
1. **Taxonomy:** Eukaryota (taxon_id: 2759)
2. **Reviewed:** Yes
3. **No fragments**
4. **Transmembrane annotation present**
5. **No signal peptide annotation**

## B. Queries
The exact UniProtKB queries are reported in `uniprot_queries.txt`.

## C. Notes
Data were retrieved programmatically using the UniProtKB REST API. A custom post-filter was applied to enforce SP length > 14 aa.

