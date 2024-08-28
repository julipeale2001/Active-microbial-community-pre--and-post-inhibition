# [![DOI](https://zenodo.org/badge/692083872.svg)](https://zenodo.org/doi/10.5281/zenodo.13383696)

# Active Microbial Community Analysis Pre- and Post-Inhibition

Analysis of microbial community from the hindguts and faeces of E. pulchripes and G. connexa after antibiotics and 2-bromo-ethane-sulfate treatments. The analysis involves the weight of millipedes, faecal counts, bacterial colony counts, 16S rRNA copy number, methane production, 16S rRNA sequence, mcrA copy, and RNA-SIP.

## Weight Loss Analysis

1_Epibolus_pulchripes_weight_loss_analysis.Rmd
    - 1_Epibolus_weight.csv

2_Glomeris_connexa_weight_loss_analysis.Rmd
    - 2_Glomeris_weight.csv
    - Weight_loss_cache
    - Weight_loss_figs

## Survival Analysis

3_Epibolus_pulchripes_survival_analysis.Rmd
    - 3_Epibolus_pul_weight.csv

4_Glomeris_connexa_survival_analysis.Rmd
    - 4_Glomeris_weight_surv.csv

## Faecal Counts

5_Epibolus_pulchripes_faecal_counts.Rmd
    - 5_Epibolus_Faecal_counts.csv
    - 05_SIP_diff_abund_between_DESeq2_cache
    - 05_SIP_diff_abund_between_DESeq2_figures

6_Glomeris_connexa_faecal_counts.Rmd
    - 6_Glomeris_faecal_counts.csv

## Bacterial Colony Counts

7_Epibolus_pulchripes_faecal_bacterial_colony_counts.Rmd
    - 7_Epibolus_colony_counts.csv

8_Glomeris_connexa_faecal_bacterial_colony_counts.Rmd
    - 8_Glomeris_colony_counts.csv

## 16S rRNA Copies

9_16S_rRNA_copies.Rmd
    - 9_16S_rRNA_copies.html
    - 9_16S_rRNA_copies.md
    - 9_ddPCR_Copy_Numbers.csv

9_Epibolus_pulchripes_16S_rRNA_copies_in_hindgut.Rmd
    - 9_Epibolus_Hindgut_Copy_Numbers.csv

10_Glomeris_connexa_16S_rRNA_copies_in_hindgut.Rmd
    - 10_Glomeris_Hindgut_Copy_Numbers.csv

11_Epibolus_pulchripes_16S_rRNA_copies_in_faeces_RA.Rmd
    - 11_Epibolus_pulchripes_16S_rRNA_copies_in_faeces.Rmd
    - 11_Epibolus_Faecal_Pellet_Copy_Numbers.csv

12_Glomeris_connexa_16S_rRNA_copies_in_faeces.Rmd
    - 12_Glomeris_Faecal_Pellet_Copy_Numbers.csv
    - Copy_number_cache
    - Copy_number_figs

## Methane Production Analysis

13_Epibolus_pulchripes_methane_production_analysis.Rmd
    - 13_Methane_production_Epibolus.csv

## Sequence Processing Scripts

14_cutadapt.sh

15_DADA2_16S_V8.5.R
15_run_DADA2_16S_V8.5.sh

## Decontamination

16_Decontamination.Rmd
    - 16_DADA2_seqtab_nochim_decontam.tsv
    - 16_DADA2_taxa_silva_decontam.tsv
    - 16_Metadata.csv

## Hindgut 16S Sequence Analysis

17a_Epibolus_pulchirpes_Hindgut_16S_sequence_analysis.Rmd
17b_Epibolus_pulchripes_ALDEX2_Hindgut.Rmd
17b_Epibolus_pulchripes_ANCOM_BC2_Hindgut.Rmd

18a_Glomeris_connexa_Hindgut_16S_sequence_analysis.Rmd
18b_Glomeris_connexa_ALDEX2_Hindgut.Rmd
    - ALDEX2.R
18b_Glomeris_connexa_ANCOM_BC2_Hindgut.Rmd

## Faecal 16S Sequence Analysis

19a_Epibolus_pulchirpes_Faecal_16S_sequence_analysis.Rmd
19b_Epibolus_pulchripes_ALDEX2_Faeces.Rmd
19b_Epibolus_pulchripes_ANCOM_BC2_Faeces.Rmd

20a_Glomeris_connexa_Faecal_16S_sequence_analysis.Rmd
20b_Glomeris_connexa_ALDEX2_Faeces.Rmd
20b_Glomeris_connexa_ANCOM_BC2_Faeces.Rmd
    - ANCOM_BC_Fig
    - ANCOM_BC_Table

## Diversity Analysis

21_Epibolus_pulchripes_Hindgut_Alpha_Diversity.Rmd
22_Glomeris_connexa_Faecal_Alpha_Diversity.Rmd
23_Glomeris_connexa_Hindgut_Alpha_Diversity.Rmd
24_Epibolus_pulchripes_Faecal_Alpha_Diversity.Rmd

25_Epibolus_pulchripes_Hindgut_Beta_Diversity.Rmd
26_Epibolus_pulchripes_Faeces_Beta_Diversity.Rmd
27_Glomeris_connexa_Hindgut_Beta_Diversity.Rmd
28_Glomeris_connexa_Faeces_Beta_Diversity.Rmd

## Methanogenesis Inhibition

29_Epibolus_pulchripes_methanogenesis_inhibition.Rmd
    - 29_Epibolus_pulchripes_methanogenesis_inhibition.csv

30_Epibolus_pulchripes_weight_loss_after_methanogenesis_inhibition.Rmd

## Methanogen Counts

31_Methanogens_counts.Rmd
    - 31_Methanogen_counts.ods

## RNA-SIP

32_cutadapt_v2.0.sh
33_DADA2_16S_merge_V9.1.R
33_run_DADA2_16S_V9.1.sh
    - 34_DADA2_pseudo

34_Decontamination_SIP_V1.6.Rmd
    - 34_Millipedes_SIP_metadata_decontam.csv
    - 34_Millipedes_SIP_metadata.csv

35_Filter_taxa_V1.4.Qmd

36_calc_tree_V2.0.sh

37_SIP_diff_abund_between_DESeq2.Qmd
    - 37_DADA2_pseudo
    - 37_Millipedes_SIP_metadata.csv
    - 37_SIP_diff_abund_between_DESeq2_cache
    - 37_SIP_diff_abund_between_DESeq2_files
    - 37_SIP_diff_abund_between_DESeq2.html
    - DESeq_byTime_a-0.1_0_prev_concensus.tsv
    - DESeq_byTime_a-0.1_0_prev.tsv
    - DESeq2_byTime_a-0.051_wTax.tsv
    - items
    - RDS
