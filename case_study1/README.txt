Here the scripts to run the case study 1 presented in the manuscript.
The scripts should be run in the following order:
1.case1_09102018_limma.R: to perform the pre-processing of the data and the DEA with limma pipeline
2. case1_09102018_edgeR.R: to perform the DEA with EdgeR pipeline
3. correlation_logFC_09102018.R: to calculate the correlation between the top DE genes estimated by limma and EdgeR
4. upSetR_091018.R: to make a upSetR plot of the intersects between the two methods
