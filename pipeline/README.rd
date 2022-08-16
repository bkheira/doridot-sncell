## Demultiplexing test

This folder aims to compare several demultiplexing method on our samples. (and to start implementing snakemake infrastructure)

Methods to compare :
- Take the max antibody
- Housemade by Ludivine
- Seurat implemented method
- Gaublomme et al. Nature 2019 (?)

Note : some of those methods also allow for the detections double celle droplets. Because of that, this test will have to go up to filtering stage.

### Input
Files with multiplexing and name of the mutliplexing antibodies

### Output
tsv file describing the classification of each cell in each method

Ideally this pipeline will be run on at least two samples
