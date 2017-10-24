# Quantifying the ease of escape of viral escape from broad and narrow antibodies to influenza hemagglutinin

This directory contains the analysis of the mutational antigenic profiling of the A/WSN/1933 (H1N1) influenza virus hemagglutinin against a collection of borad and narrow antibodies.

The study was performed by Mike Doud, Juhye Lee, and Jesse Bloom in [Bloom lab](https://research.fhcrc.org/bloom/en.html).
The analysis was performed using [dms_tools2](https://jbloomlab.github.io/dms_tools2/).

The analysis is performed by the Jupyter notebook [analysis_notebook.ipynb].

All generated results are placed in the created directory [./results/](./results).

All required input data are in the [./data/](./data/) subdirectory. Specifically, this directory includes the following:

    * [./data/WSN_HA_reference.fa](./data/WSN_HA_reference.fa) contains the wildtype sequence of the A/WSN/1933 (H1N1) HA used in this experiment.

    * [./data/samples.csv](./data/samples.csv) is a CSV file that specifies each sample, the SRA run that contains its deep sequencing data, and some other descriptive information.

    * [.data/H1toH3_renumber.csv](.data/H1toH3_renumber.csv) is a CSV file that maps the numbering from sequential (1, 2, ...) numbering of the WSN HA protein sequence to the commonly used H3 numbering scheme. The sequential number is in the *original* column, and the H3 numbering is in the *new* column.

    * [./data/Overall-WSNHA_merged_prefs_rescaled_H3numbering.csv](./data/Overall-WSNHA_merged_prefs_rescaled_H3numbering.csv) gives the site-specific amino-acid preferences for the A/WSN/1933 HA after being re-scaled to optimally fit natural sequence, as taken from the supplementary material of [Doud and Bloom (2016)](http://www.mdpi.com/1999-4915/8/6/155).
