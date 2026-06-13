# Ciliary protein lists from the Human Protein Atlas — methods and data dictionary

Fetch the data to `ciliary_proteins_hpa.csv` like this:

  python3 fetch_ciliary_proteins_hpa.py

Candidate ciliary proteins were retrieved from the Human Protein Atlas (HPA;
https://www.proteinatlas.org), version 25.1, current at the time of retrieval.
The HPA assigns subcellular localizations on the basis of antibody-based
immunofluorescence confocal microscopy; localizations to the primary cilium and
basal body derive from the antibody-based spatial-proteomics expansion of the
resource described by Hansen et al. (2025), in which 715 cilium-localizing
proteins were annotated to four primary-cilium regions ("primary cilium",
"primary cilium tip", "primary cilium transition zone" and "basal body") across
three ciliated cell lines.

Genes were collected by querying the HPA programmatic download API
(`search_download.php`) once for each of six subcellular-location categories —
*Primary cilium*, *Primary cilium tip*, *Primary cilium transition zone*,
*Basal body*, *Centrosome* and *Centriolar satellite* — using the
`subcell_location` field. The latter two microtubule-organizing-centre
categories were included to capture the wider basal-body/centrosome context from
which cilia are templated. The six result sets were combined and de-duplicated
by Ensembl gene identifier, yielding 1,237 unique genes; for each gene the
compartment(s) in which it was annotated were recorded (column
*Cilium compartment(s)*). A gene may be annotated to more than one compartment,
so per-compartment counts sum to more than the number of unique genes. For each
gene, the associated annotations held by the HPA (identifiers, chromosomal
location, protein class, functional and disease annotations, subcellular
localization with reliability, and RNA/protein expression specificity) were
exported in the same query.

Gene models and chromosomal coordinates follow Ensembl version 109, i.e. the
**GRCh38 (hg38)** human genome assembly. The HPA *Position* field reports the
chromosomal location of the gene on this assembly; exact base-pair start/end and
strand can be obtained by joining the Ensembl gene identifiers against a matching
Ensembl 109 / GRCh38 annotation (e.g. a GTF).

> Note: confirm the version string and Ensembl release against the HPA release
> history (https://www.proteinatlas.org/about/releases) for the exact version
> you query, and update the assembly statement if you use a later release.

## Per-compartment gene counts

| Cilium compartment | Genes |
| --- | --- |
| Primary cilium | 467 |
| Primary cilium tip | 130 |
| Primary cilium transition zone | 96 |
| Basal body | 451 |
| Centrosome | 513 |
| Centriolar satellite | 228 |
| **Unique genes (union)** | **1,237** |

## Column dictionary

| Column | Description |
| --- | --- |
| Cilium compartment(s) | Which of the six queried subcellular categories the gene was annotated to (added by this pipeline; semicolon-separated when more than one). |
| Gene | Approved gene symbol (HGNC). |
| Gene synonym | Alternative gene symbols/aliases (semicolon-separated). |
| Ensembl | Ensembl gene identifier (ENSG…), Ensembl v109 / GRCh38. |
| Gene description | Descriptive gene/protein name. |
| Uniprot | UniProt accession(s) mapped to the gene. |
| Chromosome | Chromosome on which the gene resides (GRCh38). |
| Position | Chromosomal location of the gene on the assembly (GRCh38). |
| Protein class | HPA protein-class annotations (e.g. enzymes, transporters, predicted membrane/secreted proteins, disease-related genes). |
| Biological process | Associated biological-process annotations (GO/UniProt keywords). |
| Molecular function | Associated molecular-function annotations (UniProt keywords). |
| Disease involvement | Disease associations annotated for the gene. |
| Evidence | Protein-existence evidence level (e.g. evidence at protein level / transcript level). |
| Subcellular location | All immunofluorescence-based subcellular locations annotated for the protein (main + additional). |
| Subcellular main location | Location(s) classified as the protein's main subcellular location. |
| Subcellular additional location | Additional, non-main subcellular location(s). |
| Reliability (IF) | Reliability of the immunofluorescence-based subcellular annotation: Enhanced > Supported > Approved > Uncertain. |
| RNA tissue specificity | Consensus-RNA tissue-specificity category: Tissue enriched, Group enriched, Tissue enhanced, Low tissue specificity, or Not detected. |
| RNA tissue distribution | Breadth of detection across tissues: Detected in all / many / some / single. |
| RNA tissue specificity score | Numeric tissue-specificity score underlying the category. |
| RNA tissue specific nTPM | Tissue(s) with elevated expression and their expression levels in nTPM (normalized transcripts per million); names the tissues in which the gene is enriched/enhanced. |
| RNA single cell type specificity | Single-cell-type specificity category (analogous categories to the tissue field). |
| RNA single cell type specific nCPM | Single-cell type(s) with elevated expression and their levels in nCPM (normalized counts per million). |
| Protein cell type specificity | Protein-level cell-type specificity from Deep Visual Proteomics (DVP); the DVP cell-type panel includes ciliated cells. |

Expression units: **nTPM** = normalized transcripts per million (bulk/cluster
RNA); **nCPM** = normalized counts per million (single-cell RNA).

## Suggested citations

- Human Protein Atlas: Uhlén M, et al. Tissue-based map of the human proteome.
  *Science* 347:1260419 (2015).
- Subcellular / cilia annotation: Hansen JN, et al. Intrinsic heterogeneity of
  primary cilia revealed through spatial proteomics. *Cell* (2025),
  S0092-8674(25)01029-3.
- Data version: Human Protein Atlas version 25.1; gene models from Ensembl
  release 109 (GRCh38).