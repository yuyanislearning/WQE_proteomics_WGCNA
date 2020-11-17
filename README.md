# WQE_proteomics_WGCNA
By: Yu Yan, 2020/9/31

Code and results for WQE, using WGCNA analyzing proteomics data

This project aims to create a workflow for analyzing proteomics data through: performing differential expression analysis on proteomics data; WGCNA to form cluster proteins based on co-expression pattern; clusters are then associated with phenotype of interest. It also includes modules for functional annotation and comparison the proteins identified from different analysis.

Proteomics search output should be formatted to fit the input of the workflow, detail are illustrated in preprocess part.
Multiple outputs will be generated based on user's preferrence. In general, the results of differential expression analysis, WGCNA and functional annotation will be main results.

Details are shown below.

## Code 
**/prepocess**  Preprocess code for missing data imputation, data transformation
Preprocess the data to form a complete matrix by using cubic spline imputation. Differential expression analysis is also included

**/WGCNA**   code for WGCNA analysis, including module traits association, connectivity calculation and so on
WGCNA preprocess, WGCNA construction, relate WGCNA module to traits

**/Reactome**  code for Reactome parent pathway identification

**/DisGeNet**  code for retrieve disease protein association from DisGeNet API

**/ID_convertion**  code for project human protein to mouse protein or vice versa

## Results
**/Differential_expression_analysis**   list of proteins with significant change, all proteins stats test results

**/WGCNA**     module found by WGCNA and hub proteins of yellow module

**/CaseOLAP**   top proteins identified by CaseOLAP

**/Reactome**   Reactome results

**/GeneOntology**     Gene ontology results

**/Clinical**    Results of differentially expressed proteins identified in cinical study


## Dependencies
WGCNA



