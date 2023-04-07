# CellHeapRelease
This repository contains supplementary materials for the paper "CellHeap: A scRNA-seq workflow for large-scale bioinformatics data analysis."

CellHeap is a flexible, extensible, and portable platform for scRNA-seq analysis, which processes hundreds of Gigabytes of scRNA-seq raw data in supercomputer environments.

CellHeap platform is composed of five phases: (i) sample curation; (ii) data pre-processing; (iii) quality control; (iv) dimensionality reduction and clustering analysis; and (v) advanced cell level and gene level analysis. 

One CellHeap’s phase can include many computational tools and couple them such that inputs and outputs consume/generate data in a flow that meets requirements for subsequent phases. It employs quality control to ensure correct results and relies on high-performance computational support for parallelizing and distributing computational tasks.

This repository includes a case study of a CellHeap workflow, with its deployment in the Brazilian Santos Dumont supercomputer, which processed COVID-19 scRNA-seq raw data considering the four initial phases. In this case study, CellHeap accessed the Gene Expression Omnibus (GEO) repository to collect samples that meet the criteria defined by the Expression Atlas of the European Molecular Biology Laboratory (EMBL) protocol.

The workflow received as input the bronchoalveolar scRNA-seq dataset PRJNA608742 from the NIH GEO repository. This dataset was generated by Liao et al. (2020) and considered the samples of 12 patients as input data.

Liao, M., Liu, Y., Yuan, J., Wen, Y., Xu, G., Zhao, J., Cheng, L., Li, J., Wang, X., Wang, F., et al.: Single-cell landscape of bronchoalveolar immune cells in patients with covid-19. Nature medicine 26(6), 842–844 (2020)

After the curation and quality control, we considered 11 patients: three healthy controls, two with mild COVID-19 symptoms, and six with severe COVID-19 symptoms.

The case study aimed to analyze the modulations of Fc receptors in lung lymphoid and myeloid cells imposed by COVID-19 by processing bronchoalveolar lavage fluid (BALF) samples.

The source code of CellHeap for Liao et al. 2020's dataset is in the Source directory, that contains four phases.
