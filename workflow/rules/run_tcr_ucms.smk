############    ----------  Rules related to  -----------    ############
############    ----------  running TCR UCMs  -----------    ############

rule process_output:
    input:
        tool_out = f'{config["paths"]["output"]}/' + 'clust_results/{ref}/{tool}/RefThold-{ref_thold}/{tool}_Clusters.csv',
        ref_info = f'{config["paths"]["output"]}/' + 'tools_input/{ref}/RefThold-{ref_thold}/GeneralBetaInfo.csv'
    output:
        [
            f'{config["paths"]["output"]}/' + 'clust_results/{ref}/{tool}/RefThold-{ref_thold}/ToolResults_' + file + '.csv'
            for file in ['Processed', 'AdjList-1', 'ClonotypeAtts-1']
        ]
    params:
        PIPELINE = PIPELINE,
        reports_path = f'{config["paths"]["output"]}/' + '/clust_results/{ref}/{tool}/RefThold-{ref_thold}',
        run_id = '{ref}_RefThold-{ref_thold}',
        process_opts = lambda wildcards: config['process_opts'][wildcards.ref]
    log:
        f'{config["paths"]["output"]}/' + 'jobs_scripts/logs/process_output/{ref}/{tool}_RefThold-{ref_thold}.log'
    benchmark: 
        repeat(f'{config["paths"]["output"]}/' + 'benchmark/process_output/{ref}/{tool}_RefThold-{ref_thold}.benchmark.txt', 1)
    threads: 1
    resources:
        mem_mb = get_mem_mb_light,
        runtime = '160h'
    conda:
        "../envs/r_analyses.yaml"
    shell:
        'Rscript {params.PIPELINE}/workflow/scripts/process_output.R --Tool {wildcards.tool} --RunID {params.run_id} --ReportsPath {params.reports_path} --ToolOutput {input.tool_out} --RefInfo {input.ref_info} {params.process_opts} > {log} 2>&1'


rule run_tool:
    input:
        f'{config["paths"]["output"]}/' + 'tools_input/{ref}/RefThold-{ref_thold}/ToolInput.tsv'
    output:
        clusters_file = f'{config["paths"]["output"]}/' + 'clust_results/{ref}/{tool}/RefThold-{ref_thold}/{tool}_Clusters.csv'
    params:
        PIPELINE = PIPELINE,
        reports_path = f'{config["paths"]["output"]}/' + '/clust_results/{ref}/{tool}/RefThold-{ref_thold}',
        run_id = '{ref}_RefThold-{ref_thold}',
        tool_opts = get_tool_opts
    log:
        f'{config["paths"]["output"]}/' + 'jobs_scripts/logs/clust_by_sim/{ref}/{tool}_RefThold-{ref_thold}.log'
    benchmark: 
        repeat(f'{config["paths"]["output"]}/' + 'benchmark/clust_by_sim/{ref}/{tool}_RefThold-{ref_thold}.benchmark.txt', 1)
    threads: 1
    resources:
        mem_mb = get_mem_mb_heavy,
        runtime = '160h'
    conda:
        select_env
    shell:
        r'''
            cd {params.reports_path}
            if [ '{wildcards.tool}' == 'GIANA' ]; then
                python {params.PIPELINE}/workflow/scripts/giana/GIANA4.py -f {input} -o {params.reports_path} {params.tool_opts} > {log} 2>&1;
                mv {params.reports_path}/ToolInput--RotationEncodingBL62.txt {params.reports_path}/GIANA_Clusters.csv
            elif [ '{wildcards.tool}' == 'TCRdist3' ]; then
                python {params.PIPELINE}/workflow/scripts/run_tcrdist3_from_fmtd_input.py --RunID {params.run_id} --InputFile {input} --ReportsPath {params.reports_path} {params.tool_opts} > {log} 2>&1
            elif [ '{wildcards.tool}' == 'GLIPH2' ]; then
                Rscript {params.PIPELINE}/workflow/scripts/gliph2/run_gliph2_from_fmtd_input.R --RunID {params.run_id} --InputFile {input} --ReportsPath {params.reports_path} {params.tool_opts} > {log} 2>&1
            elif [ '{wildcards.tool}' == 'clusTCR' ]; then
                python {params.PIPELINE}/workflow/scripts/run_clustcr_from_fmtd_input.py --InputFile {input} --ReportsPath {params.reports_path} {params.tool_opts} > {log} 2>&1
            elif [ '{wildcards.tool}' == 'iSMART' ]; then
                cp {input} {params.reports_path}/Input_CDR3s.tsv;
                python {params.PIPELINE}/workflow/scripts/ismart/iSMARTv3.py -f {params.reports_path}/Input_CDR3s.tsv {params.tool_opts} > {log} 2>&1;
                mv {params.reports_path}/Input_CDR3s.tsv_ClusteredCDR3s_*.txt {params.reports_path}/iSMART_Clusters.csv
            else
                echo "Unproper tool definition while running directive of rule 'run_tool'." >&2
                exit 1
            fi
        '''


rule get_clust_input:
    input:
        f'{config["paths"]["output"]}/' + 'gen_tcr_tables/{ref}/RefThold-{ref_thold}.csv'
    output:
        f'{config["paths"]["output"]}/' + 'tools_input/{ref}/RefThold-{ref_thold}/ToolInput.tsv',
        f'{config["paths"]["output"]}/' + 'tools_input/{ref}/RefThold-{ref_thold}/GeneralBetaInfo.csv'
    params:
        PIPELINE = PIPELINE,
        reports_path = f'{config["paths"]["output"]}/' + '/tools_input/{ref}/RefThold-{ref_thold}',
        run_id = 'RefThold-{ref_thold}'
    log:
        f'{config["paths"]["output"]}/' + 'jobs_scripts/logs/tools_input/{ref}_RefThold-{ref_thold}.log'
    benchmark: 
        repeat(f'{config["paths"]["output"]}/' + 'benchmark/tools_input/{ref}_RefThold-{ref_thold}.benchmark.txt', 1)
    threads: 1
    resources:
        mem_mb = get_mem_mb_heavy,
        runtime = '120h'
    conda:
        "../envs/r_analyses.yaml"
    shell:
        'Rscript {params.PIPELINE}/workflow/scripts/get_clust_input.R --ReportsPath {params.reports_path} --RunID {params.run_id} --InputTable {input} > {log} 2>&1'