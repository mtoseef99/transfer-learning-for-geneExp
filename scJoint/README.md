scJoint analysis
----------------

First prepare the 10xPBMC data in .h5 from RDS file(SingleCellExperiment) from
`sce_to_h5.ipynb`. The input data (RDS) is given
[here.](https://github.com/SydneyBioX/scJoint/tree/main/data_10)

For mouse primary motor cortex experiments to integrates transcriptomics,
chromatin accessibility and methylation, first download the data from this
[link](https://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scJoint/data_MOp.zip), and
place under `data_MOp` directory.

To run scJoint with different datasets and model settings, change the model
configuration in `config.py`.

To change database, change `DB` to:

-   `10x` for 10xGenomics data

-   `MOp` for Mouse Primary Motor Cortex Data

-   `db4_control` for CITE-seq and ASAP-seq PBMC data
