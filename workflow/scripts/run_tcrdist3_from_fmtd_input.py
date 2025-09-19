############    -----------   Run TCRdist3    -----------    ############
############    ---------  from formatted input  --------    ############
# By Vicente Fajardo

print('\n\n')
print('############    -----------   Run TCRdist3    -----------    ############\n')
print('############    ---------  from formatted input  --------    ############\n')

# Long version: 0.2
# Version: 0
# Subversion: 2
# Updates.
# ---> Version updates:
#   No stable version achieved yet.
# ---> Subversion updates: 0.1
#   * First subversion.
# ---> Subversion updates: 0.2
#   * We now make sure to filter out TCRs whose TCR V or J genes are not properly defined in the database used as reference by TCRdist3. Otherwise, when keeping these we'd have problems while running the clustering process.


### -------------------------- Description -------------------------- ###
# TBD.


print('\n\n')
### -------------------------- Dependencies ------------------------- ###
print('### -------------------------- Dependencies ------------------------- ###\n')
import warnings
import os.path
import numpy as np
import pandas as pd
import string
from tcrdist.repertoire import TCRrep
from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory2
from sklearn.cluster import DBSCAN #, KMeans
from optparse import OptionParser
print('Dependencies loaded...\n')


print('\n\n')
print('### --------------------------- Arguments --------------------------- ###\n')
### --------------------------- Arguments --------------------------- ###
# ---> Hardcoded.
gene_db_file = '/home/vfajardo/scripts/TCR_data_analysis/general_clustering/TCRdist3/dbs/alphabeta_gammadelta_db.tsv'
# ---> From command line.
# Set seed.
np.random.seed(1)
# Load options.
parser = OptionParser()
parser.add_option("-I", "--RunID", dest="run_id", type="string", default="Project", help="Str, ID for the comprehensive analysis to be performed.\n")
parser.add_option("-R", "--ReportsPath", dest="reports_path", type="string", default=None, help="Str, absolute path to directory to save results.\n")
parser.add_option("-F", "--InputFile", dest="input_file", type="string", default=None, help="Optional. Str, absolute path to input file listing the node attributes. Must be a csv file with the column \"general.clone.id\" and any other set of columns listing further attributes.\n")
parser.add_option("-N", "--CalcMetaclonotypes", dest="do_meta", action="store_true", default=False, help="Optional. Bool, indicates whether the script should attempt to calculate and report on metaclonotypes through the standard TCRdist3 manner or not. Although  please note that this is only possible when there is a maximum of 15,000 unique clonotypes in the dataset. Otherwise, due to memory constrains, it is not possible.\n")
parser.add_option("-C", "--Chains", dest="chains_for_clust", type="string", default="beta", help="Optional. String, indicates which chains should be used to perform clustering. Valid values are \"beta\", \"alpha\" and \"both\".\n")
parser.add_option("-D", "--MaxDist", dest="dist_radius", type="int", default=50, help="Maximum distance allowed to be saved to the TCR dissimilarity matrix. To be supplied to parameter 'radius' while calculating distances.\n")
parser.add_option("-M", "--Method", dest="clust_method", type="string", default="DBSCAN", help="Str, method to be used for clustering. Valid values are 'DBSCAN' and 'KMEANS', with 'DBSCAN' as the default.\n")
parser.add_option("-H", "--Hyper", dest="clust_hyper", type="float", default=None, help="Float, hyperparameter for clustering, whose meaning depends on the clustering method selected. In the case of DBSCAN, it refers to sklearn.cluster.DBSCAN's 'eps' and in the case of KMEANS, it refers to sklearn.cluster.kmeans's n_clusters.\n")
parser.add_option("-T", "--CPUs", dest="cpu_no", type="int", default=1, help="Int, indicates number of available CPUs.\n")
(options, args) = parser.parse_args()
# Move arguments to single variables.
run_id = options.run_id
reports_path = options.reports_path
input_file = options.input_file
do_meta = options.do_meta
chains_for_clust = options.chains_for_clust
dist_radius = options.dist_radius
clust_method = options.clust_method
clust_hyper = options.clust_hyper
cpu_no = options.cpu_no

# Sanity checks.
#       There must exist at least one proper way to save the results.
if reports_path == None:
    print('No reports path was provided, hence no pickle file with the graph to be built will be provided.\n')
else:
    to_check = os.path.exists(reports_path)
    if not to_check:
        exit('A path to save results to was provided but it does not exist. Please make sure to provide a proper path; next is the non-existing one:\n' + reports_path + '\n')
#       The input file must be properly defined.
tmp_check = os.path.exists(input_file)
if not tmp_check:
    exit('The input file provided does not exist. Please make sure to provide a proper path; next is the non-existing one:\n' + input_file + '\n')
# Present arguments
print('Run ID: ' + run_id + '\n')
print('Path to save reports to: ' + reports_path + '\n')
print('Path to input file: ' + input_file + '\n')
print('Whether metaclonotypes should be determined: ' + str(do_meta) + '\n')
print('Chains to use for clustering: ' + chains_for_clust + '\n')
print('Maximum distance allowed in dissimilarity matrix: ' + str(dist_radius) + '\n')
print('Clustering method: ' + clust_method + '\n')
print('Clustering hyperparameter: ' + str(clust_hyper) + '\n')


### --------------------------- Functions --------------------------- ###


print('\n\n')
print('### --------------------------- Load data --------------------------- ###\n')
### --------------------------- Load data --------------------------- ###
# Load input TCR data.
tcr_data = pd.read_table(input_file)
# Gene database.
gene_db = pd.read_table(gene_db_file)


print('\n\n')
print('### ---------------------- Data preprocessing ----------------------- ###\n')
### ---------------------- Data preprocessing ----------------------- ###

# ---> Input data clean up.
# Clean up: remove columns whose values are undefined across all rows.
for col in tcr_data.columns:
    if np.all(tcr_data[col].isna()):
        del tcr_data[col]
# Clean up: remove TCRs whose v and j genes are not defined in the reference database used by TCRdist3.
tmp_cols = {'trb.v', 'trb.j', 'tra.v', 'tra.j'}.intersection(tcr_data.columns)
for tmp_col in tmp_cols:
    tcr_data = tcr_data[tcr_data[tmp_col].isin(gene_db['id'])]
# Clean up: remove TCRs whose any of its CDR3 sequence is too short (<6 aas)
tmp_cols = {'cdr3b.aa.seq', 'cdr3a.aa.seq'}.intersection(tcr_data.columns)
for tmp_col in tmp_cols:
    tcr_data = tcr_data[np.logical_and(
        tcr_data[tmp_col].str.len() > 9,
        tcr_data[tmp_col].str.len() < 30
    )]
# End of clean up
tcr_data = tcr_data.reset_index()

# ---> Tool input definition.
# Add necessary columns and relabel according to what's expected by TCRdist.
cols_to_keep = {
    'donor.id' : 'subject',
    'size' : 'count',
    'tra.v' : 'v_a_gene',
    'tra.j' : 'j_a_gene',
    'cdr3a.aa.seq' : 'cdr3_a_aa',
    'trb.v' : 'v_b_gene',
    'trb.j' : 'j_b_gene',
    'cdr3b.aa.seq' : 'cdr3_b_aa',
}
tool_input = tcr_data.copy()
tool_input.rename(columns=cols_to_keep, inplace=True)
# Keep unique clonotypes in the dataframe only by combining the counts from individual donors. TCRdist fails to take into consideration the donor identity information. This step is critical, otherwise TCRdist will end up yielding clusters w/ identical CDR3 sequences coming from different donors or, in other words, public clonotypes.
tmp_cols_1 = set(['cdr3_a_aa', 'v_a_gene', 'j_a_gene', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene'])
tmp_cols_1 = list(tmp_cols_1.intersection(tool_input.columns))
tmp_cols_2 = tmp_cols_1 + ['count']
tool_input = tool_input[tmp_cols_2]
tool_input = tool_input.groupby(tmp_cols_1).sum().reset_index()
tool_input['clone_id'] = range(0, tool_input.shape[0])


print('\n\n')
print('### ------------------------- Main program -------------------------- ###\n')
### ------------------------- Main program -------------------------- ###

### -------------------- Metaclonotype definition ------------------- ###

# Perform only if both it is doable (number of unique TCRs is smaill enough) and required by user,
if do_meta and tool_input.shape[0] < 15000:
    print('\n\n')
    print('### -------------------- Metaclonotype definition ------------------- ###\n')

    # ---> Dependencies.
    import matplotlib.pyplot as plt
    from progress.bar import IncrementalBar
    from tcrdist.repertoire import TCRrep
    from tcrdist.automate import auto_pgen
    from tcrsampler.sampler import TCRsampler
    from tcrdist.background import get_stratified_gene_usage_frequency, sample_britanova, get_gene_frequencies
    from tcrdist.ecdf import distance_ecdf, _plot_manuscript_ecdfs
    from tcrdist.neighbors import bkgd_cntl_nn2
    from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius, make_motif_logo
    from tcrdist.centers import rank_centers

    # ---> Distance sparse matrix.
    tr = TCRrep(
        cell_df=tcr_data,
        organism='human',
        chains=['beta'],
        # db_file # Let's use the default reference.
        compute_distances=True,
        cpus=cpu_no
    )

    # ---> Estimation of probability of generation.
    auto_pgen(tr)
    print(f"Computing probability of generation w/ OLGA (Sethna et al 2018) for clonotypes to be used for subsequent analyses")

    # ---> Synthesize an Inverse Probability Weighted VJ Matched Background
    # Generating an appropriate set of unenriched reference TCRs is important; for each set of antigen-associated TCRs, discovered by MIRA, we created a two part background. One part consists of 100,000 synthetic TCRs whose V-gene and J-gene frequencies match those in the antigen-enriched repertoire, using the software OLGA (Sethna et al. 2019; Marcou et al. 2018). The other part consists of 100,000 umbilical cord blood TCRs sampled uniformly from 8 subjects (Britanova et al., 2017). This mix balances dense sampling of sequences near the biochemical neighborhoods of interest with broad sampling of TCRs from an antigen-naive repertoire. Importantly, we adjust for the biased sampling by using the V- and J-gene frequencies observed in the cord-blood data (see Methods for details about inverse probability weighting adjustment). Using this approach we are able to estimate the abundance of TCRs similar to a centroid TCR in an unenriched background repertoire of ~1,000,000 TCRs, using a comparatively modest background dataset of 200,000 TCRs. While this estimate may underestimate the true specificity, since some of the neighborhood TCRs in the unenriched background repertoire may in fact recognize the antigen of interest, it is useful for prioritizing neighborhoods and selecting a radius for each neighborhood that balances sensitivity and specificity. Initialize a TCRsampler -- human, beta, umbilical cord blood from 8 people.
    print(f"Using tcrsampler to construct a custom V-J-matched background.\n")
    ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
    # Stratify sample so that each subject contributes similarly to estimate of gene usage frequency
    ts = get_stratified_gene_usage_frequency(ts=ts, replace=True)
    # Synthesize an inverse probability weighted V,J gene background that matches usage in your enriched repertoire 
    vj_bg_df = tr.synthesize_vj_matched_background(ts=ts, chain='beta')

    # Get a randomly drawn stratified sampler of beta, cord blood from Britanova et al. 2016, Dynamics of Individual T Cell Repertoires: From Cord Blood to Centenarians
    df_britanova_100K = sample_britanova(size=100000)
    # Append frequency columns using, using sampler above
    df_britanova_100K = get_gene_frequencies(ts=ts, df=df_britanova_100K)
    df_britanova_100K['weights'] = 1
    df_britanova_100K['source'] = "stratified_random"
    # Combine the two parts of the background into a single DataFrame
    bg_df = pd.concat(
        [
            vj_bg_df.copy(),
            df_britanova_100K.copy()
        ],
        axis = 0
    )
    bg_df.reset_index(drop = True)
    # Assert that the backgrounds have the expected number of rows.
    assert bg_df.shape[0] == 200000
    # Save the background for future use # SKIPPED STEP. MIGHT DELETE IN FUTURE VERSIONS.
    # background_outfile = os.path.join(project_path, f"{antigen_enriched_file}.olga100K_brit100K_bkgd.csv")
    # print(f'WRITING {background_outfile}')
    # df_bkgd.to_csv(background_outfile, index = False)
    # Load the background to a TCRrep without computing pairwise distances 
    # (i.e., compute_distances = False)
    bg_tr = TCRrep(
        cell_df = bg_df,
        organism = "human", 
        chains = ['beta'], 
        compute_distances = False
    )
    # Compute rectangular distances. Those are, distances between each clone in the antigen-enriched repertoire and each TCR in the background.
    # With a single 1 CPU and < 10GB RAM, 5E2x2E5 = 100 million pairwise distances, across CDR1, CDR2, CDR2.5, and CDR3 
    # 1min 34s ± 0 ns per loop (mean ± std. dev. of 1 run, 1 loop each) 
    # %timeit -r 1 tr.compute_rect_distances(df = tr.clone_df, df2 = tr_bkdg.clone_df, store = False)

    # ---> Calculate distances.
    print(f"Compute rectangular distances.\n")
    tr.compute_sparse_rect_distances(
        df=tr.clone_df, 
        df2=bg_tr.clone_df,
        radius=50,
        chunk_size=100
    )

    # ---> Examine empirical cummulative distribution functions (ECDFs).
    # Investigate the density of neighbors to each TCR, based on expanding distance radius.
    # Compute empirical cumulative density function (ecdf)
    # Compare Antigen Enriched TCRs (against itself).
    thresholds, ecdf = distance_ecdf(
        tr.pw_beta,
        thresholds=range(0, 50, 2)
    )
    # Compute empirical cumulative density function (ecdf)
    # Compare Antigen Enriched TCRs (against) 200K probability 
    # inverse weighted background
    thresholds, bg_ecdf = distance_ecdf(
        tr.rw_beta,
        thresholds=range(0, 50, 2),
        weights= bg_tr.clone_df['weights'], 
        absolute_weight = True
    )
    # plot_ecdf similar to tcrdist3 manuscript #
    ecdf[ecdf == ecdf.min()] = 1E-10
    f1 = _plot_manuscript_ecdfs(
        thresholds, 
        ecdf, 
        ylab= 'Proportion of Antigen Enriched TCRs', 
        cdr3_len=tr.clone_df.cdr3_b_aa.str.len(), 
        min_freq=1E-10
    )
    f1.savefig(os.path.join(reports_path, f'{run_id}.ecdf_AER_plot.png'))
    f2 = _plot_manuscript_ecdfs(
        thresholds,
        bg_ecdf,
        ylab= 'Proportion of Reference TCRs',
        cdr3_len=tr.clone_df.cdr3_b_aa.str.len(),
        min_freq=1E-10
    )
    f2.savefig(os.path.join(reports_path, f'{run_id}.ecdf_BUR_plot.png'))


    # ---> Find optimal radii
    # To ascertain which meta-clonotypes are likely to be most specific, take advantage of an existing function <bkgd_cntrl_nn2>.                                                                                                                    
    level_tag = '1E5'
    centers_df = bkgd_cntl_nn2(
        tr=tr,
        tr_background=bg_tr,
        weights=bg_tr.clone_df.weights,
        ctrl_bkgd=10**-5, 
        col='cdr3_b_aa',
        add_cols=['v_b_gene', 'j_b_gene'],
        ncpus=4,
        include_seq_info=True,
        thresholds=[x for x in range(0, 50, 2)],
        generate_regex=True,
        test_regex=True,
        forced_max_radius=36
    )

    # ---> All metaclonotypes to a tsv file.
    # Save center to project_path for future use
    centers_df.to_csv(os.path.join(reports_path, f'{run_id}.centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )
        
    # Many of meta-clonotypes contain redundant information. We can winnow down to less-redundant list. We do this by ranking clonotypes from most to least specific. 
    #   <min_nsubject> is minimum publicity of the meta-clonotype,  
    #   <min_nr> is minimum non-redundancy
    # Add neighbors, k_neighbors, and nsubject columns
    centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_beta, radius_list = centers_df['radius'])
    centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
    # We determine how many <nsubjects> are in the set of neighbors 
    centers_df['nsubject']  = centers_df['neighbors'].\
            apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
    centers_df.to_csv(os.path.join(reports_path, f'{run_id}.centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )

    # ---> Output, ready to search bulk data.
    ranked_centers_df = rank_centers(
        centers_df=centers_df, 
        rank_column='chi2joint', 
        min_nsubject=2, 
        min_nr = 1
    )
    ranked_centers_df.to_csv(os.path.join(reports_path, f'{run_id}.ranked_centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )

    # ---> Output Meta-Clonotypes HTML Summary
    # Here we can make a svg logo for each NR meta-clonotype
    if ranked_centers_df.shape[0] > 0:
        cdr3_name = 'cdr3_b_aa'
        v_gene_name = 'v_b_gene'
        svgs = list()
        svgs_raw = list()
        bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
        for i,r in ranked_centers_df.iterrows():
            bar.next()
            centroid = r[cdr3_name]
            v_gene = r[v_gene_name]
            svg, svg_raw = make_motif_logo(
                tcrsampler=ts, 
                pwmat=tr.pw_beta,
                clone_df=tr.clone_df,
                centroid=centroid,
                v_gene=v_gene,
                radius=r['radius'],
                pwmat_str='pw_beta',
                cdr3_name='cdr3_b_aa',
                v_name='v_b_gene',
                gene_names=['v_b_gene','j_b_gene']
            )
            svgs.append(svg)
            svgs_raw.append(svg_raw)
        bar.next(); bar.finish()
        ranked_centers_df['svg'] = svgs
        ranked_centers_df['svg_raw'] = svgs_raw

        def shrink(s):
            return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        labels =['cdr3_b_aa','v_b_gene', 'j_b_gene', 'pgen',
                'radius', 'regex','nsubject','K_neighbors', 
                'bkgd_hits_weighted','chi2dist','chi2re','chi2joint']
        
        output_html_name = os.path.join(reports_path, f'{run_id}.ranked_centers_bkgd_ctlr_{level_tag}.html')
        with open(output_html_name, 'w') as output_handle:
            for i,r in ranked_centers_df.iterrows():
                #import pdb; pdb.set_trace()
                svg, svg_raw = r['svg'],r['svg_raw']
                output_handle.write("<br></br>")
                output_handle.write(shrink(svg))
                output_handle.write(shrink(svg_raw))
                output_handle.write("<br></br>")
                output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
                output_handle.write("<br></br>")


### -------------------- Clustering by similarity -------------------- ###

# Borrowed from Hudson et al., Immunoinformatics, 2024
def cluster_TCRDist_matrix(
    S, seqs,
    method, hyperparam=None,
    cpus=1
):
    '''Cluster distance matrix from tcrdist
    :param S: Distance matrix from tcrdist3
    :type S: array
    :param seqs: input data
    :type seqs: Pandas DataFrame
    :param method: clustering method
    :type method: str
    :param cpus: # CPUs
    :type cpus: int
    :param hyperparam: hyperparameter for clustering method
    :type hyperparam: float
    :return mapper: cluster results
    :rtype mapper: dictionary
    '''
    # Available methods
    methods = ['DBSCAN', 'KMEANS']
    method=method.upper()
    assert method in methods, r'Please choose one of the following: /n %s' % methods
    if method == 'DBSCAN':
        if not hyperparam:
            hyperparam=0.5
        clustering = DBSCAN(eps=hyperparam, min_samples=2, n_jobs=cpus).fit(S)
        labels = clustering.labels_    
    elif method == 'KMEANS':
        os.environ["OMP_NUM_THREADS"] = "%s"%(str(cpus)) # export OMP_NUM_THREADS
        if not hyperparam:
            hyperparam=500
        kmeans = KMeans(init='random',
                        n_init=10,
                        n_clusters=int(hyperparam)).fit(S)
        labels = kmeans.labels_
    return {seq: label for seq, label in zip(seqs['clone_id'].values, labels) if label!=-1}

# Adapted from Hudson et al., Immunoinformatics, 2024
def run_tcrdist3(
    df,
    chain_selection,
    radius=50,
    method='DBSCAN', hyper=None,
    chunk=True,
    cpus=1
):
    '''Run tcrdist3 clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :param cpus: # CPUs
    :type cpus: int
    :param method: clustering method
    :type method: str
    :param radius: tcrdist metaclonotype radius
    :type radius: int
    :param hyper: clustering algorithm hyperparameter
    :type hyper: float
    :return df2: cluster results
    :rtype df2: pandas DataFrame
    :return t: runtime (s)
    :rtype t: float
    '''
    # Reformat input dataframe for tcrdist
    df2 = df.copy()
    cdr3a = 'cdr3_a_aa'
    va = 'v_a_gene'
    ja = 'j_a_gene'
    cdr3b = 'cdr3_b_aa'
    vb = 'v_b_gene'
    jb = 'j_b_gene'
    # Assign constant count (n=1) to every TCR if info not provided.
    if not 'count' in df2.columns:
        df2['count']=[1]*len(df2)
    # Final data structure according to chain selection.
    if chain_selection == 'alpha':
        df_epi = df2[[cdr3a, va, ja, 'clone_id', 'count']]
    elif chain_selection == 'beta':
        df_epi = df2[[cdr3b, vb, jb, 'clone_id', 'count']]
    else:
        df_epi = df2[[va, cdr3a, ja, vb, cdr3b, jb, 'clone_id', 'count']]
    seqs = df_epi.reset_index(drop=True)
    # Final chain definition.
    if chain_selection in ['alpha', 'beta']:
        chain = [chain_selection]
    else:
        chain = ['alpha', 'beta']
    # Run tcrdist3
    print('\n*** Clustering %s %s chains with tcrdist3' % (len(seqs), chain_selection))
    # Create tcrdist object
    tr = TCRrep(cell_df=seqs,   # Input data
        organism='human',   # Organism
        chains=chain,       # Chain selection
        infer_all_genes=True,   # Infer genes
        infer_cdrs=True,        # Infer CDRs
        infer_index_cols=True,  # Infer column names
        store_all_cdr=True,     # Store CDRs
        deduplicate=False,      # Deduplicate
        compute_distances=False # Compute distances
    )
    # Compute tcrdist distances using sparse rect distances
    tr.cpus = cpus
    if chunk:
        if chain ==['alpha']:
            name = 'alpha'
        else:
            name = 'beta'
        S, _ = compute_pw_sparse_out_of_memory2(tr = tr,
            row_size      = 1000,
            pm_processes  = 1, # Our system does not allow for proper multithreading for this particular function.
            pm_pbar       = True,
            max_distance  = radius,
            reassemble    = True,
            cleanup       = True,
            assign        = False)
        S=S[name]
    else:
        tr.compute_sparse_rect_distances(radius=radius, chunk_size=500)
        if chain_selection == 'alpha':
            S = tr.rw_alpha
        elif chain_selection == 'beta':
            S = tr.rw_beta
        else:
            S = tr.rw_beta
    # Cluster tcrdist matrix
    mapper = cluster_TCRDist_matrix(S, seqs, method = method, cpus=cpus, hyperparam=hyper)
    df.loc[:,'cluster']=df['clone_id'].map(mapper)
    return df

# ---> Perform clustering.
# Determine if chunking is necessary for TCR distance matrix calculation. It is when we're dealing with >10,000 unique clonotypes.
# Source of this assertion: https://tcrdist3.readthedocs.io/en/latest/bulkdata.html
# "For paired or single-chain datasets larger than 10K clones (i.e. >10^8 pairwise comparisons) we recommend the function compute_pw_sparse_out_of_memory2."
do_chunk = tool_input.shape[0]>15000
# Carry on
# bkup = tool_input.copy()
# tool_input = bkup.copy()
# tool_input = tool_input.iloc[5000:20000]
tmp_data = run_tcrdist3(
    df=tool_input,
    chain_selection=chains_for_clust,
    radius=dist_radius,
    method=clust_method,
    hyper=clust_hyper,
    chunk=do_chunk,
    cpus=cpu_no
)

# ---> Summarize results and save
# Disregard singletons.
tmp_data = tmp_data[tmp_data['cluster'].notna()]
tmp_data = tmp_data.sort_values(by=['cluster', 'count'], ascending=[True, False])
tmp_data['cluster'] = tmp_data['cluster'].astype(int)
# Merge info with 
cols_to_keep = {
    'v_a_gene' : 'tra.v',
    'j_a_gene' : 'tra.j',
    'cdr3_a_aa' : 'cdr3a.aa.seq',
    'v_b_gene': 'trb.v',
    'j_b_gene' : 'trb.j',
    'cdr3_b_aa': 'cdr3b.aa.seq',
}
tmp_data.rename(columns=cols_to_keep, inplace=True)
# Save results.
tmp_file_name = '/'.join([reports_path, 'TCRdist3_Clusters.csv'])
tmp_data.to_csv(tmp_file_name, index=False)

print('\n\n\n')
print('Process complete!\n\n\n')