# T-cell receptor (TCR)-guided assignment of pathogen specificity to T cells
------------

The workflow
------------

This repository provides a Snakemake-based workflow to predict pathogen specificity of T cells from single-cell RNA-seq and TCR-seq datasets.<br/>
The pipeline is modular, reproducible, and designed for easy reuse across projects via a wrapper Snakefile.

For details on the workflow, please see [**paper cite to be updated upon publication**].


🔧 System requirements
------------

* **Operating system**: Linux or macOS (tested on Ubuntu 22.04, macOS Ventura)
* **[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)**: v7.22.0
* **[R](https://cran.r-project.org/)**: v4.2.2
* **[Python](https://cran.r-project.org/)**: v3.11.0
* **[GIANA](https://github.com/s175573/GIANA)**
* **[TCRdist3](https://tcrdist3.readthedocs.io/en/latest/)**
* **[iSMART](https://github.com/s175573/iSMART)**
* **[clusTCR](https://github.com/svalkiers/clusTCR)**
* **[GLIPH2](http://50.255.35.37:8080/)**
* Additional R packages: `optparse`, `Seurat`, `data.table`, `tidyr`, `stringr`, `fuzzyjoin` `ggplot2`, `ggpubr`, `pheatmap`, `Hmisc`, `corrplot`
* Additional Python packages: `optparse`, `numpy`, `pandas`, `sklearn.cluster`, `matplotlib`

📥 Installation
------------
### 1. Clone the repository to a permanent location
```bash
cd ~/tools
git clone https://github.com/VicenteFR/PatSpcImpByTCRs.git
```
Make sure to replace `~/tools` with the path you wish to clone the repo into, although this is the recommended location.<br>
**Note that you only need to clone once**. Each new dataset/project should live in its own folder and reference this repository through the dataset-specific config file—see below.<br/>
This step takes ~1 minute.

### 2. Download the TCR reference sets
Please fetch the files through GEO with the accession number GSExxx (**accession number to be updated upon publication**). These are relatively large files—make sure to store them in a proper location. The files should be named as follows:

```bash
Ref-1_CloneBasedInfo.csv
Ref-1_CellBasedInfo.csv
Ref-2_CloneBasedInfo.csv
Ref-2_CellBasedInfo.csv
```

**Review note**: If you are reviewer, please kindly use the access token provided on the manuscript file to fetch the files from GEO.

### 3. Common config file set-up
Edit the following file within the repo folder: `~/tools/PatSpcImpByTCRs/workflow/configs/common.yaml`. Replace `~/toools` by your selected repo folder if necessary. Please replace all instances of `/path/to/references` with the absolute path to the folder where you actually stored the reference set files per step 2.

**You're all set!** Now, time to test the installation...


🧪 Demo
------------

We provide a small demo dataset in `./tests` for quick testing.

1. Navigate to the output folder for the demo dataset:
```bash
cd ~/tools/PatSpcImpByTCRs/tests/test_out
```

2. Check the config file (config.yaml) to see the expected input/output and options structure. Please note that this will be the sole file that's meant to be changed for each dataset you apply the workflow to (more details below). If necessary, replace all path instances of `~/toools` by your alternative repo folder. Please make sure that all instances have been properly replaced before jumping into the next step.

3. Activate your Snakemake conda environment and run a dry run to check the DAG:
```bash
conda activate snakemake
snakemake -np
```
4. Indicating the proper of cores available to you (example with 2 cores only), execute the workflow on the demo dataset:
```bash
snakemake --cores 2
```
This will generate example outputs in the results/ folder in under 15 minutes when 2 cores are provided. Naturally, the greater the number of cores, the quicker the workflow will complete.<br/>
File `AllDone.txt` in the folder `test_out` indicates that the workflow ran in its entirety.<br/>
Several folders must be generated in the process, but the most relevant outputs are the following:<br/>
`./specificity_pred/EpRef-1/VDJMatch-TRUE/RefThold-1/UCMThold-X`: Output for CD8 T cells.<br/>
`./specificity_pred/EpRef-1/VDJMatch-TRUE/RefThold-2/UCMThold-X`: Output for CD4 T cells.

🚀 Running the Pipeline on Your Data
------------

1. Create a project folder outside the repo, e.g.:
```bash
mkdir ~/projects/my-dataset
cd ~/projects/my-dataset
```

2. Copy the template wrapper Snakefile:
```bash
cp ~/tools/PatSpcImpByTCRs/templates/wrapper.Snakefile Snakefile
```

3. Prepare a dataset-specific config:
```bash
mkdir configs
cp ~/tools/my-pipeline/templates/config.yaml configs/config.yaml
```
Please also make sure to modify options `PIPELINE` and `R_functions` under `paths` as necessary.

4. Dry run to confirm the workflow:
```bash
snakemake -np
```

5. Run the workflow:
```bash
snakemake --cores 8
```

Results will appear in the folder specified in `config.yaml` under options `paths` and `output`.


📄 Citation
--------------

If you use this workflow in your research, please cite:
**To update upon publication.**

📜 License
--------------
This repository is distributed under the MIT license.

Maintainers
-----------

Current maintainer:
* Vicente Fajardo Rosas (vfajardo@lji.org) 

Vijayanand Laboratory.  
La Jolla Institute for Immunology La Jolla, CA 92037, USA


Contact
-----------
Please email any of the current maintainers or feel free to open an issue in this repo.
