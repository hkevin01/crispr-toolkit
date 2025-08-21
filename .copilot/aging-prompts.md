# Copilot Prompt Pack: FTL1 Brain Aging Integration

## Entity Normalization
- Implement entity normalization in literature loader to map FTL1, FTL, and ferritin light chain to the same gene entity; recompute literature_evidence_score.

## CRISPRi sgRNA Design
- Add CRISPRi sgRNA design for FTL/FTL1 promoter: fetch TSS from Ensembl, design 8 sgRNAs within −50..+300 bp, compute Doench 2016 on-target and CFD off-target scores; return top 4 per species.

## Safety Model Extension
- Extend safety model with iron_homeostasis_prior: penalize systemic perturbation; lower penalty for neuron/astrocyte-restricted contexts (use GTEx/HPA expression).

## Mechanistic Feature
- In rejuvenation.predict_outcomes(), add a mechanistic feature for iron overload relief and calibrate via conformal prediction.

## GUI/Docs Integration
- Update GUI: add a provisional-evidence badge component; wire to docs/references metadata.
