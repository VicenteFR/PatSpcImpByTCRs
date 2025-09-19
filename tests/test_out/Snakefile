############    ---  Pathogen specificity imputation  ---    ############
############    --------  based on TCR sequences  -------    ############
############    --------  Integrative strategy  ---------    ############

#########################################################################
# Snakemake workflow, WRAPPER Snakefile
#
# Author:    Vicente Fajardo Rosas (La Jolla Institute for Immunology)
# Contact:   vfajardo@lji.org
# Created:   2025-09
# Updated:   2025-09
#
# Description:
#   We developed three complementary strategies to assign pathogen
#   specificities to lung T cells on the basis of their TCRs:
#   * Donor-matched approach (applicable only when donor-matched TCR seqs.
#   with known pathogen specificity (e.g., from blood) are available.
#   * Reference match approach.
#   * TCR UCM approach.
#   This workflow, deeemed the integrative strategy, tactically incorporated
#   the individual results from the three approaches.
#   Details on each separate approach are documented at our GitHub repo:
#   https://github.com/VicenteFR/PatSpcImpByTCRs
#
# Usage:
#   This file should not be modified directly for dataset-specific runs.
#   Instead, create a project folder with:
#       - THIS wrapper Snakefile.
#       - a dataset-specific config.yaml
#   For more details on usage, please refer to the repo.
#
# Citation:
#   If you use this workflow in your research, please cite:
#   <Paper under review>.
#
# License:
#   This workflow is distributed under the MIT License (see LICENSE).
#########################################################################


### ------------------------- Dependencies -------------------------- ###

import yaml, os


### ------------------------------ Main ----------------------------- ###

# ---> Preflights.
configfile: './configs/config.yaml'
# Absolute or relative path to the pipeline repo
# PIPELINE = os.path.expanduser("~/tools/PatSpcImpByTCRs")
PIPELINE = os.path.expanduser("~/shared_code/PatSpcImpByTCRs")

# ---> Full config options, combination of commong and dataset-specific options.
# Load common config
with open(os.path.join(PIPELINE, "./configs/common.yaml")) as f:
    common_config = yaml.safe_load(f)
config = {**common_config, **config}

# ---> Run pipeline.
include: PIPELINE + '/Snakefile'

rule all:
    input: f'{config["paths"]["output"]}' + '/AllDone.txt'