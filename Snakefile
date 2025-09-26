############    ---  Pathogen specificity imputation  ---    ############
############    --------  based on TCR sequences  -------    ############
############    --------  Integrative strategy  ---------    ############

#########################################################################
# Snakemake workflow, general Snakefile
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
#       - a wrapper Snakefile including this workflow
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
# import os


### --------------------------- Arguments --------------------------- ###
# Config file from wrapper.


### --------------------------- Functions --------------------------- ###

def get_mem_mb_light(wildcards, attempt):
    return attempt * 40000

def get_mem_mb_heavy(wildcards, attempt):
    return attempt * 100000

def select_env(wildcards):
    return f'../envs/{wildcards.tool}.yaml'

def get_pred_outs(wildcards):
    return expand(
        f'{config["paths"]["output"]}/' + 'specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}/IntersectedPredictionDetails.csv',
        ref=wildcards.ref,
        vdj_match=config['vdj_match_opts'][wildcards.ref],
        ref_thold=config['ref_tholds'][wildcards.ref],
        ucm_thold=config['ucm_tholds'][wildcards.ref]
    )

def get_procesed_outs(wildcards):
    return expand(
        f'{config["paths"]["output"]}/' + 'clust_results/{ref}/{tool}/RefThold-{ref_thold}/ToolResults_{file}.csv',
        ref=wildcards.ref,
        tool=config['tool_opts'],
        ref_thold=wildcards.ref_thold,
        file=['Processed', 'AdjList-1', 'ClonotypeAtts-1']
    )

def get_tool_opts(wildcards):
    return ' '.join(
        [config['tool_opts'][wildcards.tool],
        config['ref_opts'][wildcards.ref][wildcards.tool]]
    )


### -------------------------- Config file -------------------------- ###
# Process paths as necessary.
config['paths'] = {key: os.path.expanduser(path) for key, path in config['paths'].items()}
for ref_id in ['EpRef-1', 'EpRef-2']:
    config['tcr_data'][ref_id] = {key: os.path.expanduser(path) for key, path in config['tcr_data'][ref_id].items()}
for ref_id in ['EpRef-1', 'EpRef-2']:
    config['gex_data'][ref_id] = {key: os.path.expanduser(val) if key in ['aggr_table', 'seurat_obj'] else val for key, val in config['gex_data'][ref_id].items()}


### ----------------------------- Rules ----------------------------- ###

include: './workflow/rules/prepare_input.smk'
include: './workflow/rules/run_tcr_ucms.smk'
include: './workflow/rules/run_preds.smk'

rule all_mock:
    input:
        expand(
            f'{config["paths"]["output"]}/' + 'best_preds/{ref}/BestPredReport.csv',
            ref=config['tcr_data']
        )
    output:
        f'{config["paths"]["output"]}/' + 'AllDone.txt'
    shell:
        'touch {output}'

# Might move to a separate file if more specific steps for determining best options come up.
rule compare_ref_options:
    input:
        get_pred_outs
    output:
        f'{config["paths"]["output"]}/' + 'best_preds/{ref}/BestPredReport.csv'
    params:
        PIPELINE = PIPELINE,
        run_path = f'{config["paths"]["output"]}/' + '/specificity_pred/{ref}',
        opts_file_1 = f'{PIPELINE}/configs/common.yaml',
        opts_file_2 = f'{config["paths"]["output"]}/configs/config.yaml',
        reports_path = f'{config["paths"]["output"]}/' + '/best_preds/{ref}',
        opts_file = f'{config["paths"]["output"]}/configs/config.yaml'
    log:
        f'{config["paths"]["output"]}/' + 'jobs_scripts/logs/compare_ref_options/{ref}.log'
    benchmark: 
        repeat(f'{config["paths"]["output"]}/'+'benchmark/compare_ref_options/{ref}.benchmark.txt', 1)
    threads: 1
    resources:
        mem_mb = get_mem_mb_light,
        runtime = '120h'
    conda:
        "../envs/r_analyses.yaml"
    shell:
        'Rscript {params.PIPELINE}/workflow/scripts/comp_opts.R --ReportsPath {params.reports_path} --RunPath {params.run_path} --OptsFile1 {params.opts_file_1} --OptsFile2 {params.opts_file_2} --RefID {wildcards.ref} > {log} 2>&1'