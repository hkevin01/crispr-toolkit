# MAGeCK Integration for CRISPR Toolkit

## Overview

The CRISPR Toolkit now includes comprehensive **MAGeCK integration** for analyzing CRISPR screens with specialized focus on aging and senescence research. This integration provides industry-standard analysis tools specifically enhanced for studying pathways like **ReHMGB1/RAGE signaling**, **JAK/STAT**, and **NF-κB** that are critical in aging and senescence.

## 🚀 Key Features

### ✅ **Phase 1 Implementation Complete**

- **MAGeCK Wrapper**: Python interface to MAGeCK for count analysis and gene ranking
- **Screen Quality Control**: Comprehensive QC metrics and visualization
- **Pathway Enrichment**: Multi-database pathway analysis with aging focus
- **Senescence Analysis**: Specialized analysis for senescence pathways
- **ReHMGB1 Investigation**: Targeted analysis of ReHMGB1/RAGE signaling cascade

## 📁 Module Structure

```
src/crispr_toolkit/analysis/screens/
├── __init__.py                 # Main module exports
├── mageck_wrapper.py          # MAGeCK integration wrapper
├── screen_qc.py               # Quality control analysis
├── pathway_analysis.py        # Pathway enrichment analysis
└── senescence_analysis.py     # Senescence-specific analysis
```

## 🔧 Installation & Dependencies

The integration requires the following new dependencies (already added to requirements.txt):

```bash
# CRISPR Analysis Tools
mageck>=0.5.9          # Industry-standard screen analysis
crisprdesign>=1.0.0    # Guide design tools
bagel2>=1.0.0          # Bayesian screen analysis
```

Install MAGeCK separately:
```bash
# Via conda (recommended)
conda install -c bioconda mageck

# Via pip
pip install mageck
```

## 🧬 Usage Examples

### 1. MAGeCK Screen Analysis

```python
from crispr_toolkit.analysis.screens import MAGeCKAnalyzer

# Initialize analyzer
mageck = MAGeCKAnalyzer()

# Count reads from FASTQ files
count_file = mageck.count_reads(
    fastq_files=['treatment_1.fastq', 'treatment_2.fastq'],
    library_file='sgrna_library.txt',
    output_prefix='screen_analysis',
    sample_labels=['Treatment_1', 'Treatment_2']
)

# Identify essential genes
results = mageck.test_essential_genes(
    count_file=count_file,
    treatment_samples=['Treatment_1', 'Treatment_2'],
    control_samples=['Control_1', 'Control_2'],
    output_prefix='essential_genes'
)

# Analyze senescence pathways
senescence_results = mageck.senescence_pathway_analysis(results)
print(f"Found {len(senescence_results['senescence_hits'])} senescence hits")
```

### 2. Screen Quality Control

```python
from crispr_toolkit.analysis.screens import ScreenQualityControl
import pandas as pd

# Load your count matrix
count_matrix = pd.read_csv('count_matrix.txt', sep='\t', index_col=0)

# Run QC analysis
qc = ScreenQualityControl()
qc_results = qc.analyze_count_matrix(count_matrix)

# View recommendations
for rec in qc_results['recommendations']:
    print(f"- {rec}")

# Generate QC plots
plots = qc.create_qc_plots(count_matrix, output_dir='./qc_plots')

# Export QC report
qc.export_qc_report('screen_qc_report.txt')
```

### 3. Pathway Enrichment Analysis

```python
from crispr_toolkit.analysis.screens import PathwayEnrichment

# Load screen results
screen_results = pd.DataFrame({
    'gene': ['HMGB1', 'AGER', 'JAK1', 'STAT3', ...],
    'log2fc': [-1.5, -1.2, -0.8, -1.0, ...],
    'fdr': [0.001, 0.005, 0.05, 0.02, ...]
})

# Run pathway enrichment
pathway_analyzer = PathwayEnrichment()
enrichment_results = pathway_analyzer.run_enrichment_analysis(
    screen_results=screen_results,
    fdr_threshold=0.1
)

# Create summary table
summary_df = pathway_analyzer.create_enrichment_summary(enrichment_results)
print(summary_df.head())

# Export results
pathway_analyzer.export_enrichment_results(
    enrichment_results,
    'pathway_enrichment.csv'
)
```

### 4. ReHMGB1 Pathway Analysis

```python
from crispr_toolkit.analysis.screens import SenescenceScreenAnalyzer

# Initialize senescence analyzer
senescence = SenescenceScreenAnalyzer()

# Comprehensive senescence analysis
analysis_results = senescence.analyze_senescence_screen(
    screen_results=screen_results,
    fdr_threshold=0.1,
    lfc_threshold=0.5
)

# ReHMGB1 pathway specific results
rehmgb1_results = analysis_results['rehmgb1_pathway_hits']
print(f"ReHMGB1 pathway disruption score: {rehmgb1_results['pathway_disruption_score']:.3f}")

# Key regulators
for regulator in rehmgb1_results['key_regulators']:
    print(f"{regulator['gene']}: LFC={regulator['log2fc']:.2f}, FDR={regulator['fdr']:.3f}")

# Therapeutic targets
targets = analysis_results['therapeutic_targets']
print(f"Identified {len(targets['senescence_suppressors'])} potential therapeutic targets")
```

## 🧪 Senescence Pathway Coverage

The integration includes comprehensive gene sets for:

### **ReHMGB1/RAGE Signaling**
- HMGB1, AGER, MYD88, IRAK1, IRAK4, TRAF6, TAK1
- MAPK14, MAPK8, MAPK9

### **JAK/STAT Pathway**
- JAK1, JAK2, JAK3, TYK2
- STAT1, STAT2, STAT3, STAT4, STAT5A, STAT5B, STAT6

### **NF-κB Signaling**
- NFKB1, NFKB2, RELA, RELB, REL
- IKBKA, IKBKB, IKBKG, NFKBIA, NFKBIB

### **Cell Cycle Arrest**
- CDKN1A, CDKN1B, CDKN2A, CDKN2B, CDKN2C, CDKN2D
- RB1, RBL1, RBL2, TP53, TP73

### **SASP Factors**
- IL6, IL1A, IL1B, TNF, CXCL8, CCL2, CCL20
- CXCL1, CXCL2, MMP1, MMP3, MMP9

## 📊 Quality Control Metrics

The QC module provides comprehensive metrics:

- **Read Distribution**: Sample and gene-level read count statistics
- **Library Representation**: Coverage and detection rates
- **Sample Correlation**: Between-sample correlation analysis
- **Zero Count Analysis**: Missing data patterns
- **Replicate Consistency**: Coefficient of variation analysis
- **Automated Recommendations**: QC-based analysis suggestions

## 🎯 Applications for Your Research

Based on your ReHMGB1 research, this integration enables:

### **Direct Applications**
- **Screen ReHMGB1 pathway components** for senescence modulators
- **Identify JAK/STAT and NF-κB regulators** affecting senescence
- **Discover therapeutic targets** that suppress senescence
- **Validate pathway interactions** in senescence cascade

### **Analytical Capabilities**
- **Multi-condition comparison** (ReHMGB1 vs oxidized form)
- **Tissue-specific analysis** across multiple cell types
- **Time-course analysis** of senescence progression
- **Drug screening** for senescence interventions

## 🚀 Next Steps - Phase 2 Integration

```markdown
- [ ] **Phase 2: Aging Biomarker Integration** (Upcoming)
  - [ ] Add pyaging for 100+ aging clocks
  - [ ] Integrate biolearn for aging biomarker analysis
  - [ ] Enhance clinical intelligence with aging metrics
  - [ ] Create aging-specific visualizations

- [ ] **Phase 3: Advanced Clinical AI** (Future)
  - [ ] PyHealth integration for EHR processing
  - [ ] Medical code standardization
  - [ ] Enhanced clinical reporting
  - [ ] Interactive dashboards
```

## 🔬 Example Analysis Workflow

```python
# Complete ReHMGB1 pathway analysis workflow
from crispr_toolkit.analysis.screens import *

# 1. Quality Control
qc = ScreenQualityControl()
qc_results = qc.analyze_count_matrix(count_matrix)

# 2. MAGeCK Analysis
mageck = MAGeCKAnalyzer()
screen_results = mageck.test_essential_genes(
    count_file, treatment_samples, control_samples, 'rehmgb1_screen'
)

# 3. Pathway Enrichment
pathway = PathwayEnrichment()
enrichment = pathway.run_enrichment_analysis(screen_results.gene_summary)

# 4. Senescence Analysis
senescence = SenescenceScreenAnalyzer()
senescence_analysis = senescence.analyze_senescence_screen(screen_results.gene_summary)

# 5. ReHMGB1 Specific Analysis
rehmgb1_hits = senescence_analysis['rehmgb1_pathway_hits']
therapeutic_targets = senescence_analysis['therapeutic_targets']

print(f"✓ Identified {rehmgb1_hits['significant_hits']} ReHMGB1 pathway hits")
print(f"✓ Found {len(therapeutic_targets['senescence_suppressors'])} therapeutic targets")
```

## 📈 Benefits for Aging Research

This integration provides:

- **Industry-standard analysis** with MAGeCK integration
- **Aging-specific pathways** curated for longevity research
- **ReHMGB1 pathway focus** aligned with your current research
- **Therapeutic target identification** for senescence intervention
- **Quality-controlled analysis** with comprehensive metrics
- **Reproducible workflows** for collaborative research

---

**🎉 Phase 1 MAGeCK Integration: Complete!**

The CRISPR Toolkit now provides comprehensive screen analysis capabilities specifically designed for aging and senescence research. Ready to enhance your ReHMGB1 pathway investigations and accelerate discovery of senescence modulators.
