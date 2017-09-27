1) ampliSeq.panel.txt
from CO25560_Ion_AmpliSeq_Comprehensive_Cancer_Panel_Gene_List_final9062012.pdf

2) gene_list_trusight_cancer.txt
from gene_list_trusight_cancer.xlsx

3) cancer.panels.txt
from evernote cancer panels
# copy text from cancer panels in evernote to tmp (manually remove "and" and stuff in parentheses)
sed 's/,//g' tmp | tr ' ' '\n' | awk 'NF > 0' > cancer.panels.txt
### this file contains duplicaed lines - dealt with below

4) cancer.wtrusight.wampliseq.panels.txt
cat cancer.panels.txt gene_list_trusight_cancer.txt ampliSeq.panel.txt | sort | uniq -c | awk '{print $2"\t"$1}' > cancer.wtrusight.wampliseq.panels.txt
# therefore cancer.wtrusight.wampliseq.panels.txt has the number of sources for each gene

NOTE: CO25560_Ion_AmpliSeq_Comprehensive_Cancer_Panel_Gene_List_final9062012.pdf ignored because has nothing that isn't already in one of the above

5) TCGAbreastCancer_SMGall.txt and TCGAbreastCancer_SMGsubtypes.txt
Two parts of supplementary table two of TCGA breast cancer paper, from TCGAbreastCancer_suppTables1-4.xls. 
SMGall is the list of genes that were significantly mutated considering all samples, whereas SMGsubtypes is the list of genes that were significantly mutate in specific subtypes (but not when considering all samples)

6) multiplatformAnalysis12cancerTypes_top40SMGs.txt and multiplatformAnalysis12cancerTypes_otherSMGs.txt
From "Multiplatform analysis of 12 cancer types reveals molecular classification within and across tissue of origin"  supplementary table 2: multiplatformAnalysis12cancerTypes_suppTable2.xlsx.
top40 is, as expected, the top 40 most significantly mutated genes.
otherSMG is the entire list of significantly mutated genes minus the top 40

7) foundationOne.txt and foundationOneHeme.txt
Extracted genes from respective pdf files (only the ones for which they sequence the whole gene sequence)

8) stampv2.txt
Genes extracted from Gene sheet of stampv2.xlsx

9) ongen.excessASE.txt
The 71 genes shown to have more ASE in cancer than normals from Ongen et al. (FDR 5%)
Data are from Supplementary table 8 found in ongen.tables.xlsx

10) gene lists-for rna seq collab.xlsx
Gene lists from Jason Merker.
Split into stanford.cancer.gold.txt and ../other_genes/stanford.cardio.gold.txt

genes.ordered.txt is generated from the above files by script select_cancer_genes.R
