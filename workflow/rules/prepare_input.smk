############    ----------  Rules related to  -----------    ############
############    -------  preparing general input  -------    ############


rule get_gen_input:
    input:
        f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/CellBasedTCRData.csv'
    output:
        f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/RefThold-{ref_thold}.csv'
    params:
        # ref_data = config['ref_specs'][{wildcards.ref}]['ref_data'],
        qry_id=config['qry_id'],
        ref_data = lambda wildcards: config['ref_specs'][wildcards.ref]['clust_data'],
        filt_criteria = lambda wildcards: config['input_opts'][wildcards.ref][wildcards.ref_thold]
    threads: 1
    resources:
        mem_mb = get_mem_mb_light,
        runtime = '1h'
    shell:
        'cat <( echo -e "dataset.lab,meta.data,expansion.thold,filt.criteria" ) <( echo -e "{params.qry_id},{input},1,NA" ) <( echo -e "{wildcards.ref},{params.ref_data},1,{params.filt_criteria}" ) > {output}'


rule get_comp_tcr_info:
    output:
        f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/CellBasedTCRData.csv',
        f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/ClonotypeBasedTCRData.csv'
    params:
        PIPELINE = PIPELINE,
        reports_path = f'{config["paths"]["output"]}/' + '/gen_tcr_tables/{ref}',
        opts_file_1 = f'{PIPELINE}/configs/common.yaml',
        opts_file_2 = f'{config["paths"]["output"]}/configs/config.yaml',
        vdj_data = lambda wildcards: config['tcr_data'][wildcards.ref]['aggr_out'],
        gex_data = lambda wildcards: config['gex_data'][wildcards.ref]['seurat_obj'],
        vdj_aggr = lambda wildcards: config['tcr_data'][wildcards.ref]['aggr_table'],
        gex_aggr = lambda wildcards: config['gex_data'][wildcards.ref]['aggr_table'],
        raw_data = lambda wildcards: config['gex_data'][wildcards.ref]['aggr_out'],
        clusts_lab = lambda wildcards: config['gex_data'][wildcards.ref]['clusters_tag'],
        funct_file = f'{PIPELINE}/workflow/scripts/handy_functions.R'
    log:
        f'{config["paths"]["output"]}/' + 'jobs_scripts/logs/comp_tcr_info/{ref}.log'
    threads: 1
    resources:
        mem_mb = get_mem_mb_heavy,
        runtime = '160h'
    conda:
        "../envs/r_analyses.yaml"
    shell:
        'Rscript {params.PIPELINE}/workflow/scripts/get_comp_tcr_info.R --ReportsPath {params.reports_path} --OptsFile1 {params.opts_file_1} --OptsFile2 {params.opts_file_2} --RefID {wildcards.ref} --AggrPath {params.vdj_data} --GExData {params.gex_data} --VDJAggrTable {params.vdj_aggr} --GExAggrTable {params.gex_aggr} {params.raw_data} --ClustsLab {params.clusts_lab} > {log} 2>&1'