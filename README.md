# Neoantigen ranking stability

A small public-data stress test asking a narrow question: **if the same experimentally tested neoantigen candidates are ranked using different reasonable objectives, how stable is the very top of the list?**

## Data

The script uses Supplementary Data S9 from *Immunopeptidomics-guided identification of functional neoantigens in non-small cell lung cancer*:

https://github.com/ab604/lung-neoantigen-supplement

It restricts the analysis to HLA class-I rows with quantitative values for mutant binding affinity, mutant-vs-wild-type fold change, and the published composite rank.

## Result

Mean Jaccard overlap between each patient's top-3 candidates:

- binding vs agretopicity: **0.22**
- binding vs composite: **0.32**
- agretopicity vs composite: **0.40**

The top of the list changes substantially depending on the ranking objective in this small tested set.

## Reproduce

```bash
pip install -r requirements.txt
python analyze.py
```

The script downloads the public supplementary table directly and writes `top3_overlap_by_patient.csv`.

## What I would test next

A more realistic follow-up would quantify stability under predictor choice, score perturbations, and uncertainty, then identify candidates that remain high-ranked across plausible pipelines.

## Limitations

This is human NSCLC data, not a test of any company's pipeline. It compares ranking criteria rather than distinct HLA-binding model versions, and the experimentally tested set is small. The result should not be interpreted as evidence of a production bottleneck in any specific system.
