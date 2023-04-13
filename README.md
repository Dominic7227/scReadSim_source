[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7824380.svg)](https://doi.org/10.5281/zenodo.7824380)

# scReadSim source code
The source code for the reproduction of results in "scReadSim: a single-cell RNA-seq and ATAC-seq read simulator".

## Software scReadSim
We provide a Python-based implementation for the convenience of the community. For tutorials and other details, check [our website](http://screadsim.readthedocs.io/).


## Usage instructions

Four folders correspond to four tasks:
- `10x_scRNA-seq`: Using 10x scRNA-seq for synthetic read generation, verification, and benchmarking of deduplication tools
- `10x_scATAC-seq`: Using 10x scATAC-seq for synthetic read generation and verification
- `sci-ATAC-seq`: Using sci-ATAC-seq for synthetic read generation and verification
- `sci-ATAC-seq_designedpeaks`: Using sci-ATAC-seq for synthetic read generation with designed ground-truth peaks, verification, and benchmarking of peak calling tools

Under each folder, `Synthetic data generation` contains codes for generating synthetic reads and `Verification` contains codes for reproducing the figures in the manuscript. A detailed list is as follows.

### Task: 10X_scATAC-seq
Folder `Synthetic data generation`
|Code|Function|
|:-|:-|
|`NOINPUT_ATAC_main_20221130_GreyArea.sh`|main script| 
|`NOINPUT_ATAC_BAM2CountMatrix.py`|convert BAM file to count matrix|
|`NOINPUT_ATAC_ComplePeakFunction.py`|prepare scATAC-seq non-peaks|
|`NOINPUT_ATAC_GenerateBAMCoord_20221130_GreyArea.py`|generate synthetic read coordinates|
|`NOINPUT_ATAC_ReadError_20220421.py`|introducing sequencing errors to synthetic reads|
|`NOINPUT_ATAC_SyntheticMat.R`|synthetic count matrix generation|
|`NOINPUT_ATAC_VerifyRead_20221130.py`|calculate pseudo-bulk read coverage for each feature|

Folder `Verification`
|Code|Function|
|:-|:-|
|`Fig2bBtmS7c_ATAC_VerifyRead_Peak_MakePlot_Nversion_20230316.R`|Figs. 2b Bottom and S7c| 
|`Fig2bUpMidS7b_ATAC_VerifyRead_MakePlot_Nversion.20230316.R`|Figs. 2b Upper, Middle and S7b|
|`FigS4_ATAC_VerifyCount_20230128.R`|Fig. S4|


### Task: 10X_scRNA-seq
Folder `Synthetic data generation`
|Code|Function|
|:-|:-|
|`NOINPUT_RNA_main_20221011_UMItranscriptlevel2.sh`|main script|
|`NOINPUT_RNA_main_VaryRDandCN_20221013.sh`|main script for varying cell number and sequencing depth| 
|`NOINPUT_RNA_BAM2CountMatrix_20220928.py`|convert BAM file to count matrix|
|`NOINPUT_RNA_ComplePeakFunction_20220603.py`|prepare scRNA-seq features|
|`NOINPUT_RNA_GenerateBAMCoord_20221011_UMITranscriptLevel2.py`|generate synthetic read coordinates|
|`NOINPUT_RNA_ReadError_20220421.py`|introducing sequencing errors to synthetic reads|
|`NOINPUT_RNA_SyntheticMat.R`|synthetic count matrix generation|
|`NOINPUT_RNA_SyntheticUMIMat_VaryRDandCN_20221013.R`|synthetic count matrix generation for varying cell number and sequencing depth|
|`NOINPUT_RNA_VaryRDandCN_DownsampleBAM_20221013.R`|downsample cells for varying cell number|

Folder `Verification`
|Code|Function|
|:-|:-|
|`Fig1c_RNA_Minnow_20230306.sh`|Fig. 1c run minnow| 
|`Fig1c_RNA_Minnow_PrepareCountMatrix.R`|Fig. 1c prepare minnow input|
|`Fig1c_RNA_Minnow_PlotCoverage_20230316.py`|Fig. 1c|
|`Fig2defg_RNA_VerifyPipeline_20230315_withoutAlevin_Nversion.R`|Fig. 2defg|
|`Fig2defg_RNA_VerifyPipeline_CellRanger_20221011_UMItranscriptlevel2.sh`|Fig. 2defg run cellranger|
|`Fig2defg_RNA_VerifyPipeline_TimeComp_CellRanger_20221015.sh`|Fig. 2defg time usage benchmark of cellranger|
|`Fig2defg_RNA_VerifyPipeline_TimeComp_UMItools_20221015.sh`|Fig. 2defg time usage benchmark of UMI-tools|
|`Fig2defg_RNA_VerifyPipeline_UMITools_20221011_UMItranscriptlevel2.sh`|Fig. 2defg run UMI-tools|
|`FigS2_RNA_VerifyCount_20230128.R`|Fig. S2|
|`FigS3_RNA_VerifyRead_MakePlot_20221017.R`|Fig. S3|
|`FigS12S13_RNA_VerifyPipeline_20230315_withAlevin_Nversion.R`|Figs. S12 and S13|
|`FigS12S13_RNA_VerifyPipeline_Alvein_20221011_UMItranscriptlevel2.sh`|Figs. S12 and S13 run alevin|
|`FigS12S13_RNA_VerifyPipeline_TimeComp_Alevin_20221015.sh`|Figs. S12 and S13 time usage benchmark of alevin|


### Task: sci-ATAC-seq
Folder `Synthetic data generation`
|Code|Function|
|:-|:-|
|`NOINPUT_ATAC_main_20221130_GreyArea.sh`|main script| 
|`NOINPUT_ATAC_BAM2CountMatrix.py`|convert BAM file to count matrix|
|`NOINPUT_ATAC_ComplePeakFunction.py`|prepare scATAC-seq non-peaks|
|`NOINPUT_ATAC_GenerateBAMCoord_20221130.py`|generate synthetic read coordinates|
|`NOINPUT_ATAC_ReadError_20220421.py`|introducing sequencing errors to synthetic reads|
|`NOINPUT_ATAC_SyntheticMat.R`|synthetic count matrix generation|
|`NOINPUT_ATAC_VerifyRead_20221130.py`|calculate pseudo-bulk read coverage for each feature|

Folder `Verification`
|Code|Function|
|:-|:-|
|`Fig1deS6_ATAC_SCAN-ATAC-Sim_20230303.sh`|Figs. 1de and S6 run SCAN-ATAC-Sim| 
|`Fig1deS6_ATAC_SCAN-ATAC-Sim_CompareCountMat_20230310.R`|Figs. 1de and S6|
|`Fig1deS6_ATAC_SCAN-ATAC-Sim_CountMat_20230310.py`|Figs. 1de and S6 generate count matrices for comparison|
|`Fig2a_ATAC_SpitCellType_20230310.sh`|Fig. 2a| 
|`Fig2b_ATAC_SpitCellType_20230310.py`|Fig. 2b|
|`FigS5_ATAC_VerifyCount_20230128.R`|Fig. S5| 
|`FigS8_ATAC_VerifyRead_MakePlot_20221130.R`|Fig. S8|
|`FigS8_ATAC_VerifyRead_Peak_MakePlot_20221130.R`|Fig. S8 peak calling comparison|


### Task: sci-ATAC-seq_designedpeaks
Folder `Synthetic data generation`
|Code|Function|
|:-|:-|
|`INPUT_ATAC_main_mm9TSS_20221130_GreyArea.sh`|main script| 
|`INPUT_ATAC_BAM2CountMatrix.py`|convert BAM file to count matrix|
|`INPUT_ATAC_ComplePeakFunction.py`|prepare scATAC-seq non-peaks|
|`INPUT_ATAC_demoInputPeak.py`|design ground-truth peaks|
|`INPUT_ATAC_GenerateBAMCoord_20221130.py`|generate synthetic read coordinates|
|`INPUT_ATAC_MatchPeakFunction_20220302.py`|match ground truth (non)peaks and trust worthy (non)peaks|
|`INPUT_ATAC_SyntheticMat_withCluster_SelectCellType.R`|synthetic count matrix generation|
|`NOINPUT_ATAC_ReadError_20220421.py`|introducing sequencing errors to synthetic reads|
|`NOINPUT_ATAC_VerifyRead_20221130.py`|calculate pseudo-bulk read coverage for each feature|

Folder `Verification`
|Code|Function|
|:-|:-|
|`Fig2hiS10S14S15_ATAC_VerifyRead_TSS_ProcessPeakBed_20221130.sh`|Figs. 2hi, S10, S14 and S15 run peak calling tools| 
|`Fig2hiS10_ATAC_VerifyRead_TSS_MakePlot_Nversion_20230316.R`|Figs. 2hi and S10|
|`FigS14S15_TSS_Intervene_upset_logTransformation_20221009.R`|Figs. S14 and S15| 


## License
This pacakge is licensed under the terms
of the **MIT License**.
