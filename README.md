# WQE_proteomics_WGCNA
Code and results for WQE, using WGCNA analyzing proteomics data
This project aims to create a workflow for analyzing proteomics data through differential expression analysis and WGCNA to prioritize proteins that of interest based on experiment condition. It also includes modules for functional annotation and comparison for results coming from different analysis.

## Code 
**/prepocess**  Preprocess code for missing data imputation, data transformation

**/WGCNA**   code for WGCNA analysis, including module traits association, connectivity calculation and so on

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
