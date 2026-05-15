# Mycobacterium tuberculosis target prioritization for virtual screening

This repository contains code for the prioritization of Mtb targets to be used as a target profiling panel. Selection is based on current relevance of the proteins and essentiality and vulnerability data from [Bosch et al, 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)00824-2). Data curation for pipeline validation from ChEMBL and BindingDB can akso be found in this repository

This work is performed in collaboration with [DaltonTx](https://www.daltontx.com/).

### Outputs

1. List of curated 417 targets selected by essentiality and vulnerability for _M.tuberculosis_ (H37Rv strain)
2. List of molecules tested in ChEMBL for growth inhibition of _M.tuberculosis_ (dose-response and single-point). Data is informed by the pipeline in the [chembl-antimicrobial-tasks](https://github.com/ersilia-os/chembl-antimicrobial-tasks) repository.
3. Binding affinity data available for the 417 targets from ChEMBL and BindingDB

## About the Ersilia Open Source Initiative

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) and our [GitHub profile](https://github.com/ersilia-os).
