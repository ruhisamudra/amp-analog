# BMAP-27 Analog Design for Selective Antimicrobial Targeting

Computational drug design study engineering novel analogs of bovine 
cathelicidin BMAP-27 to improve selective targeting of Gram-positive 
bacterial pathogens while reducing cytotoxicity toward healthy human cells.

---

## Background

Antibiotic resistance causes nearly 4.5 million deaths annually, with 
Gram-positive ESKAPE pathogens like Staphylococcus aureus and Enterococcus 
faecium increasingly resistant to frontline drugs like vancomycin and 
methicillin. Antimicrobial peptides (AMPs) are promising next-generation 
therapeutics due to their multi-target lytic mechanisms, but many — 
including the bovine cathelicidin BMAP-27 — are cytotoxic toward healthy 
human cells, limiting therapeutic use.

This study developed a reusable in silico pre-screening pipeline to engineer 
and evaluate BMAP-27 analogs before synthesis, addressing a key gap in AMP 
drug design workflows.

## Hypothesis

AMP analogs with decreased hydrophobicity and increased cationic charge will 
bind more selectively and stably to bacterial membrane proteins compared to 
healthy human cell targets.

## Repository Contents
- amp-analog-analysis.py — reusable Python toolkit (see below)
- analog pdbs/           — AlphaFold-generated 3D structures for all analogs
- BMAP-27-a0.pdb         — original BMAP-27 (PDB ID: 2KET reference)
- BMAP-27-a1.pdb         — decreased hydrophobicity, non-active region
- BMAP-27-a2.pdb         — decreased hydrophobicity, active region
- BMAP-27-a3.pdb         — decreased hydrophobicity, both regions
- BMAP-27-a4.pdb         — increased charge, non-active region
- BMAP-27-a5.pdb         — increased charge, active region
- BMAP-27-a6.pdb         — increased charge, both regions
- BMAP-27-a7.pdb         — decreased hydrophobicity + increased charge simultaneously
- BMAP-27-a8.pdb         — a7 with aromatic substitutions for enhanced amphipathicity

## Python Toolkit

`amp-analog-analysis.py` provides three reusable modules:

**Module 1 — Analog Screener**  
Calculates hydrophobicity % and net charge at pH 7.4 for any peptide 
sequence using IMGT and Prot Pi scales. Flags candidates meeting the 
design thresholds established in this study.

**Module 2 — Docking Analyzer**  
Ranks analogs by pyDock binding energy and polar contacts against any 
target, and flags which analogs improve over a reference.

**Module 3 — Selectivity Scorer**  
Computes a selectivity ratio comparing mean bacterial vs human target 
binding energy to estimate cytotoxic risk.

### Quick Start

```python
from amp-analog-analysis import screen_analog, analyze_docking, score_selectivity

# Screen a new peptide sequence
result = screen_analog("MyPeptide", "GRFKRFRKKFKK", amidated=True)

# Run with your own docking results
results = [
    {'analog': 'MyPeptide', 'binding_energy': -68.5, 'polar_contacts': 5},
]
df = analyze_docking(results, target_name="MprF")
```

### Requirements
- pandas

## Methods

**Analog Design**  
8 analogs of BMAP-27 were designed through amino acid substitutions 
targeting hydrophobicity and positive charge independently and 
simultaneously in both the N-terminus active region and C-terminus 
antimicrobial region. Validated using Peptide 2.0 and Prot Pi.

**Structure Prediction**  
AlphaFold generated 3D PDB structures for all analogs (included in 
`analog pdbs/`). Structural alignment and RMSD calculated in PyMOL.

**Molecular Docking**  
45 docking simulations across all peptide-target pairs using GRAMM 
(binding location) and pyDock (binding energy). Targets:
- Bacterial: MprF (6LVF), autolysin (1GVM), pneumolysin (5AOD)
- Human: neutrophil elastase (3Q76), erythrocytes (1HYN)
- Control: vancomycin (1AA5)

## Key Findings

- Analogs with increased positive charge (≥ +10) and hydrophobicity 
  ≤ 38% showed stronger binding to all three bacterial targets than 
  the original peptide
- **Top candidates**: BMAP-27-a5 and BMAP-27-a7 improved over the 
  original across all bacterial targets in both binding energy and 
  polar contacts
- **Design thresholds established**: hydrophobicity ≤ 38%, charge ≥ +10
- Vancomycin showed significantly weaker bacterial binding than all 
  BMAP-27 analogs, suggesting higher potency at lower concentrations
- High erythrocyte binding affinity persisted across all analogs, 
  consistent with known hemolytic activity of BMAP-27

## Tools & Databases

| Purpose | Tool |
|---|---|
| Structure prediction | AlphaFold |
| Visualization & alignment | PyMOL |
| Docking — binding location | GRAMM Web Server |
| Docking — binding energy | pyDock |
| Ligand preparation | MGLTools |
| Sequence validation | Peptide 2.0, Prot Pi |
| Structure retrieval | RCSB Protein Data Bank |

## Acknowledgements

Dr. Rakesh Kumar Tiwari (Western University of Health Sciences) for 
mentorship throughout this project and Dr. Ryan Arvidson (Case Western 
Reserve University) for computing resources for AlphaFold structure 
generation.
