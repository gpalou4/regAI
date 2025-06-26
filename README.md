# regAI

This is the repository for the manuscript "_Copy number-independent allelic imbalance in mRNA is selected in cancer and has prognostic relevance_"

Authors: Guillermo Palou-Márquez\*1, Pere Pericot-Masdevall\*1, and Fran Supek1,2,3

1 Institute for Research in Biomedicine (IRB Barcelona), The Barcelona Institute for Science and Technology (BIST), 08028 Barcelona, Spain

2 Biotech Research and Innovation Centre, Faculty of Health and Medical Sciences, University of Copenhagen, 2200 Copenhagen, Denmark

3 Catalan Institution for Research and Advanced Studies (ICREA), Barcelona, Spain

\* Authors contributed equally

Correspondence to: fran.supek@bric.ku.dk 

**Abstract**:

Allelic imbalance (AI) in levels of mRNAs that originated from the paternal and maternal copies of a gene can contribute to phenotypic variation and influence disease traits, including cancer. It is widely appreciated that AI at the DNA level, resulting from somatic copy number alterations (CNA) in tumors, generates second-hit events by which mutated tumor suppressor genes are inactivated and mutated oncogenes are activated; such somatic CNAs are also a known cause of AI at the mRNA level. Apart from CNA, other mechanisms could lead to AI of somatic mutations in mRNA expression and also drive cancer evolution. By integrating genomic and transcriptomic pan-cancer data, we show that mRNA allelic imbalance favoring the mutant allele in driver genes is subject to positive selection, generating second-hit events often independently of somatic CNA. In some cases, the somatic coding mutations could induce allele-specific expression directly, for example, with splicing-altering exonic mutations, which can be selected in various cancer genes. However, in the majority of cases, these and related somatic mutation effects (which might in principle alter transcription output via impacting promoters or intragenic enhancers) do not explain the CNA-independent mRNA-level AI, suggesting prevalent epigenetic alterations affecting alleles differently in tumors. Importantly, the mRNA AI events associate with worse overall survival across all cancer types, outperforming various other predictive markers. Our study suggests that mRNA allelic imbalances can occur independently of CNA but similarly function as second-hit events to somatic mutations, driving tumorigenesis and represent valuable prognostic biomarkers for cancer patient stratification.

<p align="center">
  <img
    src="https://pfst.cf2.poecdn.net/base/image/db187af325eade466db23177181ed8e187cac931dc0490ca920ffc485a5cc328?w=630&h=429&pmaid=408166007"
    alt="Determinants of allelic imbalance (AI) in cancer genomes"
    width="630">
</p>

# Code

There are two files:

-**AI_beta_binomial_models** --> Builds beta–binomial regression models that quantify allelic imbalance (AI) for every single-nucleotide variant (SNV) detected across the TCGA pan-cancer cohort (-> 854 ,787 SNVs). For this, it uses a pre-processed DNA- and RNA-level allele counts generated from matched WES and RNA-Seq data by Strelka2 (not provided, see Methods in the manuscript). The input data already integrates additional annotation (survival data, cancer driver status etc.) or pre-calculated predictions (DNA methylation, _Puffin-D_, _Sei_, _SpliceAI_, _Pangolin_, etc)

-**AI_cancer_analysis** --> Re-creates all manuscript figures and tables from the previous dataframe produced above. Performs AI-dN/dES selection tests, survival modelling, enrichment analyses, etc. Generates data frames used for publication-quality plots (PDF/PNG) later on.





