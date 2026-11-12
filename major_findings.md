# Major findings

Running summary of results from the analyses in this repo. Add new bullets under the relevant section as work proceeds. Details of scripts are in `CLAUDE.md`; early species/climate notes are also in `README.md`.

## Variables / background

- **grpc1** and **grpc2** summarize variation in a set of growth variables (see `growth_curves/rTPC/parameter_correlations.R`):
  - **grpc1**: Topt, t50low (low variability), and maximum growth rate (growth rate at Topt).
  - **grpc2**: maximum growth temperature (incl. ctmax), skew, and breadth of the temperature reaction norm.
- **MAT** (mean annual temperature) and **HDD4 growing** (growing degree days) are highly correlated climate variables, so they are not independent tests.

## Phenotype vs. climate (`growth_curves/rTPC/local_adaptation_lms.R`)

- Topt is significantly correlated with MAT and with HDD4 (growing-season GDD).
- Other results (species differences, t50low vs. non-growing-season GDD, ctmax multimodality in Nf) are in `README.md`.

## Genotype-environment/phenotype associations (LFMM, `LFMM/`)

Significant SNPs = `FDR.sig.man == "sig"` in `data/LFMM/<cov>.lfmm.txt`. Gene assignment from `LFMM/nearest_neighbor_genes.R` (in-gene, else nearest protein-coding gene; reannotated genome, `Q2T73_` IDs). Coding effect from `LFMM/snp_effects.sh` + `snp_effects.R` (`bcftools csq`; checked against an independent Biostrings calculation for 1,500 SNPs, 100% agreement).

| Covariate | Sig. SNPs | In gene | Genes hit | Genes with >1 SNP |
|---|---|---|---|---|
| topt | 5 | 1 | 5 | 0 |
| t50low | 0 | – | – | – |
| ctmax | 6 | 1 | 6 | 0 |
| grpc1 | 19 | 6 | 18 | 1 |
| grpc2 | 34 | 4 | 22 | 6 |
| mat | 261 | 44 | 159 | 39 |
| hdd4_growing | 2563 | 175 | 472 | 177 |

- Most significant SNPs are intergenic (nearest-gene distances typically a few kb). Genome-wide: 90.2% intergenic, 4.5% missense, 3.8% synonymous (424,676 SNPs).
- No significant SNP is shared between any growth covariate (topt, ctmax, grpc1, grpc2) and MAT/HDD4.
- MAT and HDD4 share 44 of MAT's 261 SNPs (50 shared genes), so the two sets overlap only partly despite the correlated variables.

### Coding-effect classes of significant SNPs vs. genome-wide background (Fisher, uncorrected)

- **MAT**: enriched for missense (10.0% vs 4.5%, p = 2e-4) and loss-of-function-type effects (3 SNPs; 1.1% vs 0.3%), depleted for intergenic (83% vs 90%).
- **HDD4**: the opposite pattern for coding classes — enriched for intergenic (93% vs 90%, p = 1e-7), depleted for synonymous (2.3% vs 3.8%) and missense (3.4% vs 4.5%). So despite the correlated variables, MAT and HDD4 do not give the same effect-class pattern; HDD4's 2,563 SNPs are many and clustered, so its p-values are inflated by non-independence.
- **grpc1**: intron/splice-region enriched (3 of 19; p = 0.002); other classes not distinguishable from background.
- topt, ctmax, grpc2: too few SNPs to say anything about effect classes.
- Files: `data/LFMM/snp_effects/sig_vs_background_effect_classes.txt`, per-SNP effects in `data/LFMM/annotations/<cov>.SNPs.gene_found.nearest_neighbors.txt`.

### MAT: genes with missense or loss-of-function SNPs (29 genes, 1 SNP each)

26 missense + 3 stop-gained. About two thirds are annotated only as "hypothetical protein"; domains/GO give a bit more. Notable groups (no enrichment test done, so treat as leads, not results):

- **Membrane trafficking / secretion**: `BET1` (ER-to-Golgi SNARE; Q2T73_000480), `SEC5` (exocyst; _000771), a Vam6/VPS39-like CNH-domain protein (_000591), a Rab-GAP TBC protein (_005037), `GEF2` chloride channel (_006146).
- **Ion / solute transport**: `CCC2` Cu2+ P-type ATPase (_000456), ABC transporter (_003021), sugar/inositol MFS transporter (_007897), and an MFS transporter with a **stop-gained** (W329*, _005433).
- **Membrane lipid**: `ERG24` C-14 sterol reductase (_002978); sterol composition affects membrane fluidity, a plausible temperature link.
- **Secreted / cell-wall enzymes (SignalP+)**: glycoside hydrolase (_001830), phosphoesterase (_002969), subtilisin-like S8 protease (_003012), S10 serine carboxypeptidase (_004546), germin/cupin (_003167), histidine acid phosphatase (_006554), PAN/Apple domain protein (_005587).
- **DNA replication/repair / transcription**: `TOP2` DNA topoisomerase 2, which also sits in an antiSMASH cluster (_008369); Rad9/Ddc1 checkpoint clamp (_004176); a DNA-binding transcription activator (_010297).
- **Other stop-gained**: AAA+ ATPase with an early stop (Q39*, _011064; likely true loss of function) and a NACHT/WD40 protein (W1149*, _006568; a late stop, so effect depends on protein length and domain position).
- Table with all 29 genes: filter `data/LFMM/annotations/mat.SNPs.gene_found.nearest_neighbors.txt` on `effect_class` in (`missense`, `LoF/other_coding`) and join to `mat.gene_SNP_hits.txt` on `geneID`.

### HDD4: genes with missense or loss-of-function SNPs (79 genes, 96 SNPs)

96 SNPs: 87 missense, 7 stop-gained, 1 splice-donor, 1 stop-lost. 13 genes carry more than one SNP (up to 4). Again mostly "hypothetical protein". Only 3 of the 79 genes are also on the MAT list: the two stop-gained genes `_005433` (MFS transporter) and `_011064` (AAA+ ATPase), plus `_010297` (transcription activator). So the overlap between MAT and HDD4 coding hits is small. Leads, not tested for enrichment:

- **Multi-SNP genes with LoF-type effects**: `_010267` (P-loop NTPase; 4 SNPs: splice-donor, stop-gained W289*, two adjacent missense), `_006507` (DEAD/DEAH helicase; stop Q20* plus two changes at codon 3, all very close to the N-terminus, so check the gene start/model), `_008739` (L-ornithine monooxygenase-like, i.e. siderophore biosynthesis-type; W141*), `_011075` (BRCT domain; two stop-gained), `_004587` (MFS transporter; stop-lost + missense).
- **Cell wall**: `FKS1` 1,3-beta-D-glucan synthase (`_003173`; 2 missense), glycoside hydrolases GH5 (`_004612`) and GH16 (`_010447`), a GT2 glycosyltransferase, pectate lyase PL1 (`_005163`), cutinase/acetylxylan esterase CE5 (`_004991`), CE1 esterase (`_004597`), AA7 FAD-linked oxidase (`_007922`), LysM-domain secreted proteins (`_005597`), a Stig1-like secreted protein, an asparaginase (`_005430`). Several of these sit on tig00000071 near the MAT/HDD4 region around `_004597`-`_005597`; see the clustering caveat below.
- **Transport**: ABC transporters (including `ATM1`, mitochondrial Fe-S transporter, `_011893`), `MAM3` (Mg2+/CNNM, `_003477`), Ca2+/H+ exchanger (`_009841`), amino-acid permeases, six more MFS transporters, sugar transporter.
- **Lipid / membrane**: `PLB1` lysophospholipase (`_004617`), `FAA4` long-chain fatty-acid-CoA ligase (`_010463`), a polyprenol/steroid reductase-like protein (`_002818`), an AMP-binding fatty-acyl-type ligase.
- **Stress / DNA / RNA**: `RAD7` (UV-damaged DNA-binding; `_004367`), `RRP6` exosome nuclease, `DED1` DEAD-box helicase, tRNA adenosine deaminase, shugoshin, `TOS4`, `VAC7`, a heterokaryon-incompatibility (HET) domain protein (`_007957`), protein kinase (`_001102`), ankyrin-repeat proteins.
- **Caveat on clustering/paralogs**: several multi-SNP "genes" have SNPs within a few codons of each other (e.g. E85K/E87K in `_011076`, codon 23 hit twice in `_001203`, 1001/1002 in `_004249`), and some hits are in tandem gene pairs (`_005596`/`_005597`, `_011075`/`_011076`). That pattern is typical of linked variants or of reads from paralogous copies mapping together, so these may be one signal rather than independent ones; check LD and read depth/heterozygosity before leaning on them.

### Functional enrichment of missense/LoF genes (`LFMM/function_enrichment.R`)

One-sided Fisher tests of GO, PFAM, InterPro, CAZyme family and secreted/protease flags, BH-adjusted within covariate x term type, terms needing >=2 genes in the set. Gene sets = genes with an in-gene missense/LoF-type significant SNP (MAT 29 genes, HDD4 79). Background = the 9,492 genes that have at least one coding SNP in the dataset (a gene with no variable coding site cannot be a hit); results with all 11,871 genes as background are in the same files and agree. Output: `data/LFMM/function_enrichment/<cov>.enrichment.txt`.

- **MAT**: significant enrichment (FDR < 0.05) for **cell wall** (GO:0005618, 2 of 29 genes vs 15 of 9,492; also fungal-type cell wall), **membrane** (7 genes, 24% vs 7%), **proteolysis** (3 genes vs 1.5%), **ATP hydrolysis activity** (3 genes) and **secreted (SignalP+) proteins** (7 genes, 24% vs 8%, FDR 0.017). This matches the by-eye impression of secreted/cell-wall enzymes (GH, subtilisin, carboxypeptidase, acid phosphatase, germin). No PFAM, InterPro or CAZyme family reaches significance.
- **HDD4**: much weaker. Only **leucine-rich repeat** (PF13516; 2 genes vs 8) passes FDR 0.05 (0.016 in the coding-SNP background); LRR superfamily, DEAD/DEAH helicase, type 1 protein exporter (ABC) and other helicase/DNA-unwinding GO terms are at FDR 0.04-0.06 but rest on 2-3 genes each. **Secreted** is not enriched for HDD4 (10 of 79; FDR 0.25). So the cell-wall/secreted signal is specific to MAT.
- Caveats: set sizes are small (the "significant" terms rest on 2-7 genes), GO terms are nested and redundant so the individual terms are not independent tests, and the HDD4 multi-SNP genes may be one signal (see clustering caveat above). Treat as hypothesis-generating.

### Promoter-region classification (`promoters/`) — baseline, before RNA-seq TSS

`promoters/classify_snps_promoter.R` places every SNP relative to gene structure, respecting neighbour orientation. Run **without** RNA-seq transcript starts (start codon used as the TSS; no UTR classes), so promoter distances are measured from the ATG and are preliminary. Pipeline to get real TSS from the 8 RNA-seq samples (strandedness check -> HISAT2 -> StringTie -> per-gene TSS) is in `promoters/remote/` (SLURM), see `promoters/README.md`.

- Genome-wide (424,676 SNPs): 78.7% distal intergenic, 9.8% within CDS span, 6.0% proximal promoter (251-1000 bp upstream), 2.2% core promoter (<=250 bp), 3.3% 3' flank.
- **MAT** significant SNPs are enriched in promoter classes (proximal 11.5% vs 6.0%, p = 6e-4; core 4.2% vs 2.2%, p = 0.03) and depleted in distal intergenic (61% vs 79%).
- **HDD4** shows the opposite: enriched in distal intergenic (86% vs 79%) and depleted in promoter classes (proximal 3.9%, core 1.4%) and in coding sequence. Same direction as the coding-effect result (HDD4 intergenic-enriched, MAT coding-enriched), reinforcing that the two climate variables are not giving the same signal. Both inherit the SNP-clustering caveat above.
- Small sets (topt, ctmax, grpc1, grpc2) have too few SNPs to test.
- Per-SNP results: `data/promoters/all_snps.region.txt` and the `region_class`/`promoter_gene`/`promoter_dist` columns in each `data/LFMM/annotations/<cov>.SNPs.gene_found.nearest_neighbors.txt`.

## To follow up

- Check that the MAT secreted/cell-wall enrichment survives accounting for SNP density and gene length (e.g. permutation that preserves number of coding SNPs per gene) and for LD clustering.
- Extend function enrichment to the other covariates' gene lists (grpc1, grpc2, topt, ctmax are too small for per-term tests; consider pooling, or testing nearest-gene sets).
- Genes near the Topt/grpc1 and ctmax/grpc2 hits.
- Account for SNP clustering (e.g. LD pruning or per-gene tests) before leaning on enrichment p-values for MAT/HDD4.
- Decide how to handle intergenic SNPs (~83-93% of significant SNPs; nearest gene is often several kb away, so the nearest-gene assignment is a weak link to function). Options: distance cutoff for "likely regulatory" (e.g. within 1-2 kb upstream), per-gene upstream/downstream flags using strand, test whether sig SNPs are closer to genes than background, or treat as candidate regions rather than genes.
- Process SweepFinder2 results (`LFMM/sweepfinder2_formating.sh`, `sweepfinder2_scan.slurm`; per-contig CLR outputs): define a CLR significance threshold, call sweep regions, and test overlap with LFMM-significant SNPs/genes.
- Run the promoter pipeline on the cluster (`promoters/README.md`), then re-run `classify_snps_promoter.R` with real TSS and re-test promoter enrichment; add motif evidence (STREME/FIMO on upstream regions) and, if wanted, conservation.
