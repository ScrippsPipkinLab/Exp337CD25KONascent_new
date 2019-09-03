# Exp337CD25KONascent

#### Contributors: <br/>
> Huitian Diao (Yolanda): Experiment & Binformatics <br/>
> Dapeng Wang: Experiment <br/>
> Matthew Pipkin: Advisor

## Manuscript progress
1. **Chromatin remodeling atlas:** [GSuite](https://drive.google.com/drive/folders/1lQkHaRpWIaQ0S_ir95EF-ODGVYXxwww4?usp=sharing); [Dropbox](https://www.dropbox.com/sh/lrswxf2msgenqcj/AADE3R-FuQcxOk59wkrtzQ5Ja?dl=0)

## External data storage <br/>

## Experiment design <br/>
- Cell types: P14 CD25KO v.s. P14 WT
- Time points (Post in vitro activation): 0h, 6h, 24h, 48h

## Analysis steps <br/>
### 0. QC <br/> 
>__Output files__: [Multi-QC results](0_QC/multiqc_report.html)

### 1. Trim & Alignment & Convert & Filter
>__Script__: [HPC codes](codes)

### 2. Strand split & HTseq count & Count compilation & Normalization
>__Script__: [HPC codes](codes)  <br/>
__Output files__: <br/>
[HTseq-count output](2_count/1_dupr_count_new/) <br/>
[Compiled count](2_count/2_compiled_csv/) <br/>
[Normalized TPM](2_count/3_tpm/)

### 3. DEseq
>__Script__: [DEseq2](codes_locol/2_2_DEseq2.R)  <br/>
__Output files__: <br/>
[ENSEMBL ID](3_DE-seq/DEseq2_out/0_original/)  <br/>
[Gene name](3_DE-seq/DEseq2_out/0.1_original_GN/)





