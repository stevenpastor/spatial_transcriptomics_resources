# Visium HD Spatial Transcriptomics - Slide Outline

## Slide 1: What is Spatial Transcriptomics?

**Key Message**: Spatial transcriptomics captures gene expression while preserving tissue architecture—solving the fundamental limitation of scRNA-seq.

### Content:
- **The Problem**: Single-cell RNA-seq dissociates tissue → loses spatial context. We know WHAT cells express, but not WHERE they are.
- **The Solution**: Spatial transcriptomics measures gene expression in situ, preserving tissue architecture.
- **Key Insight**: WHERE + WHAT = biological mechanism. Cell function depends on location, neighbors, and microenvironment.
- **Historical Arc**: ISH (one gene) → FISH (few genes) → Spatial transcriptomics (hundreds to whole transcriptome)
- **2020 Nature Method of the Year** - recognized as transformative technology

### Visual:
- Diagram: tissue → dissociation (loses context) vs. spatial (preserves context)
- Example: tumor microenvironment where immune cell position relative to tumor cells determines response

---

## Slide 2: Technology Landscape

**Key Message**: Two fundamental approaches—sequencing-based (whole transcriptome, lower resolution) vs. imaging-based (targeted panels, subcellular resolution)—each with distinct trade-offs.

### Content:

| Feature | Sequencing-Based | Imaging-Based |
|---------|-----------------|---------------|
| **Examples** | Visium, Visium HD, Slide-seq, HDST, Stereo-seq | MERFISH, seqFISH+, Xenium, CosMx, CODEX |
| **Gene Coverage** | Whole transcriptome (~18,000+ genes) | Targeted panels (100-6,000 genes) |
| **Resolution** | 2-55µm (Visium HD: 2µm) | Subcellular (~100nm) |
| **Throughput** | Moderate | High (millions of transcripts) |
| **Tissue Types** | FFPE + Fresh Frozen | Primarily Fresh Frozen (expanding to FFPE) |
| **Quantification** | UMI-based, digital | Spot counting, analog |

**When to use which:**
- Sequencing-based: Discovery studies, FFPE archives, when you need unbiased whole-transcriptome
- Imaging-based: Validation, subcellular localization, when you have target gene lists, multiplexed protein + RNA

---

## Slide 3: Advantages, Limitations, and Applications

**Key Message**: Spatial transcriptomics unlocks biological questions impossible with dissociative approaches, but comes with unique computational and experimental challenges.

### Advantages:
- Preserves tissue architecture and cell-cell spatial relationships
- Enables neighborhood analysis, ligand-receptor signaling in context
- Compatible with clinical FFPE specimens (Visium HD)
- Unbiased whole-transcriptome discovery (sequencing-based)
- Integrates with histology (H&E, IF)

### Limitations:
- Resolution vs. gene coverage trade-off (improving rapidly)
- 2D sections miss 3D tissue organization
- Computational demands: large datasets (millions of bins/spots)
- Cost: $1,000-3,000 per sample (reagents + sequencing)
- Batch effects between sections, tissue handling artifacts
- No temporal information (snapshot)

### Applications:
- **Tumor Microenvironment**: Immune cell infiltration patterns, therapy resistance niches
- **Neuroscience**: Brain region mapping, neurodegenerative disease spatial patterns
- **Developmental Biology**: Organogenesis, cell fate specification in spatial context
- **Clinical Pathology**: Biomarker discovery, diagnostic spatial signatures
- **Infectious Disease**: Host-pathogen spatial interactions

---

## Slide 4: Visium HD Technology Deep Dive

**Key Message**: Visium HD achieves near-single-cell resolution (2µm) with whole-transcriptome coverage, using a continuous barcoded array and computational binning/segmentation.

### Technology:
- **CytAssist instrument**: Automated tissue placement, compatible with standard histology workflows
- **Barcoded array**: Continuous 2µm barcoded oligos covering 6.5mm × 6.5mm capture area
- **Tissue types**: FFPE (primary) and Fresh Frozen
- **Staining**: H&E (standard) or IF (immunofluorescence)
- **Sequencing**: Standard Illumina short-read sequencing

### Data Processing (Space Ranger v4.0):
- **Binning**: Aggregates 2µm barcodes into square bins: 2µm, 8µm, 16µm
  - 2µm: Maximum resolution, very sparse (most bins <10 UMIs)
  - 8µm: Near single-cell (~1-3 cells per bin), good balance of resolution and signal
  - 16µm: Multi-cell (~5-15 cells), robust signal, good for initial exploration
- **Segmentation** (new in v4.0): StarDist nuclei segmentation on H&E → single-cell resolution
  - Assigns barcodes to individual nuclei
  - Output: One observation per cell (not per bin)
  - Enables true single-cell spatial analysis

### Comparison with Original Visium:
| Feature | Visium (v1) | Visium HD |
|---------|------------|-----------|
| Resolution | 55µm spots | 2µm continuous array |
| Cells per spot | ~5-30 | ~1 per 8µm bin |
| Capture area | 6.5mm × 6.5mm | 6.5mm × 6.5mm |
| Gap between spots | 100µm center-to-center | Continuous (no gaps) |
| FFPE support | Yes (v2) | Yes (primary) |
| Segmentation | No | Yes (StarDist) |

---

## Slide 5: When to Use Visium HD

**Key Message**: Visium HD is ideal for unbiased whole-transcriptome spatial discovery at near-single-cell resolution, especially from FFPE archives—but it's not the right tool for every spatial question.

### Ideal Use Cases:
- **Whole-transcriptome discovery**: No prior gene list needed → find novel spatial markers
- **FFPE tissue archives**: Unlock spatial information from clinical biobanks
- **Near-single-cell resolution needed**: 8µm bins ≈ 1-3 cells, segmentation → single cell
- **Complex tissues**: Heterogeneous tissues where cell type composition varies spatially
- **Integration with histology**: Direct H&E overlay, pathologist-friendly

### Not Ideal For:
- **Targeted validation**: If you have a gene panel, Xenium/CosMx are faster and cheaper
- **Subcellular localization**: Need imaging-based methods (MERFISH, seqFISH+)
- **Live cell dynamics**: Fixed tissue only
- **3D reconstruction**: 2D sections (serial sectioning possible but complex)
- **Large cohorts (>50 samples)**: Cost and compute scale considerations
- **Low-input tissues**: Very small biopsies may not cover capture area efficiently

### Practical Considerations:
- **Data size**: ~5-50 GB per sample (depending on sequencing depth and bin size)
- **Compute**: 16+ GB RAM for 8µm, 64+ GB for 2µm analysis
- **Runtime**: Space Ranger: 2-8 hours; downstream analysis: hours to days depending on depth
- **Cost**: ~$1,500-2,500 per sample (reagents + CytAssist + sequencing)

---

## Slide 6: Analysis Considerations & Pipeline Overview

**Key Message**: Visium HD analysis requires careful QC (spatial artifacts are common), appropriate resolution choice, and spatial-aware methods throughout the pipeline.

### Resolution Choice:
- Start with 16µm for exploration and QC
- Use 8µm for primary analysis (best balance)
- 2µm only for specific high-resolution questions (very sparse)
- Segmented data for true single-cell analysis (when available)

### QC Considerations Unique to Visium HD:
- **Spatial QC is essential**: Tissue folds, debris, edge effects create spatial artifacts
- **Moran's I on QC metrics**: High spatial autocorrelation of MT% → localized tissue damage
- **Resolution-dependent thresholds**: QC cutoffs differ between 2µm, 8µm, 16µm bins
- **Segmentation artifacts**: Over/under-segmentation, nuclear vs. cytoplasmic capture

### Analysis Pipeline Overview:
1. **QC & Filtering** → Spatial-aware quality control
2. **Normalization & HVG Selection** → Standard scanpy workflow
3. **Dimensionality Reduction** → PCA, UMAP (for visualization)
4. **Cell Type Annotation** → Deconvolution (binned) or label transfer (segmented)
5. **Spatial Analysis** → Neighborhood enrichment, spatially variable genes
6. **Cell-Cell Communication** → Ligand-receptor analysis in spatial context
7. **Validation** → Compare with known biology, orthogonal data

### Deconvolution vs. Segmentation:
- **Binned data (8/16µm)**: Multiple cells per bin → use deconvolution (cell2location, RCTD, STdeconvolve)
- **Segmented data**: One cell per observation → direct annotation (label transfer, marker-based)
- **Best practice**: Run both, compare results for concordance
