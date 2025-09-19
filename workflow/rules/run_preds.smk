############    ----------  Rules related to  -----------    ############
############    -------------  predictions  -------------    ############

rule intersect_modalities:
    input:
        pred_res = f'{config["paths"]["output"]}/' + 'specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}/PredictionDetails.csv',
        ref_def = f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/RefThold-{ref_thold}.csv',
        qry_clone_info = f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/ClonotypeBasedTCRData.csv'
    output:
        f'{config["paths"]["output"]}/'+ 'specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}/IntersectedPredictionDetails.csv'
    params:
        PIPELINE = PIPELINE,
        reports_path = f'{config["paths"]["output"]}/' + '/specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}',
        r_mod = config['preds_r_mod'],
        opts_file_1 = f'{PIPELINE}/configs/common.yaml',
        opts_file_2 = f'{config["paths"]["output"]}/configs/config.yaml'
    log:
        f'{config["paths"]["output"]}/' + 'jobs_scripts/logs/intersect_modalities/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}.log'
    benchmark: 
        repeat(f'{config["paths"]["output"]}/' + 'benchmark/intersect_modalities/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}.benchmark.txt', 1)
    threads: 1
    resources:
        mem_mb = get_mem_mb_light,
        runtime = '120h'
    shell:
        'module unload R\n'
        'module load {params.r_mod}\n'
        'Rscript {params.PIPELINE}/workflow/scripts/intersect_modalities.R --RefDefFile {input.ref_def} --RunPath {params.reports_path} --OptsFile1 {params.opts_file_1} --OptsFile2 {params.opts_file_2} --RefID {wildcards.ref} --CloneInfo {input.qry_clone_info} > {log} 2>&1'


rule predict_specificity:
    input:
        get_procesed_outs,
        ref_def = f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/RefThold-{ref_thold}.csv',
        ref_info = f'{config["paths"]["output"]}/' + 'tools_input/{ref}/RefThold-{ref_thold}/GeneralBetaInfo.csv'
    output:
        f'{config["paths"]["output"]}/' + 'specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}/PredictionDetails.csv'
    params:
        PIPELINE = PIPELINE,
        run_path = f'{config["paths"]["output"]}/' + '/clust_results/{ref}',
        opts_file_1 = f'{PIPELINE}/configs/common.yaml',
        opts_file_2 = f'{config["paths"]["output"]}/configs/config.yaml',
        reports_path = f'{config["paths"]["output"]}/' + '/specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}',
        r_mod = config['preds_r_mod'],
    log:
        f'{config["paths"]["output"]}/' + 'jobs_scripts/logs/specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}.log'
    benchmark: 
        repeat(f'{config["paths"]["output"]}/' + 'benchmark/specificity_pred/{ref}/VDJMatch-{vdj_match}/RefThold-{ref_thold}/UCMThold-{ucm_thold}.benchmark.txt', 1)
    threads: 1
    resources:
        mem_mb = get_mem_mb_light,
        runtime = '120h'
    shell:
        'module unload R\n'
        'module load {params.r_mod}\n'
        'Rscript {params.PIPELINE}/workflow/scripts/predict_ag_specificity.R --ReportsPath {params.reports_path} --RunPath {params.run_path} --RefDefFile {input.ref_def} --RefInfo {input.ref_info} --OptsFile1 {params.opts_file_1} --OptsFile2 {params.opts_file_2} --RefID {wildcards.ref} --ExpansionThold {wildcards.ref_thold} --UCMThold {wildcards.ucm_thold} --VDJMatch {wildcards.vdj_match} > {log} 2>&1\n'
        'touch {output}'