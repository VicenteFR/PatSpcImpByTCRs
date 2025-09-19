############    -----------   Run ClusTCR    ------------    ############
############    ---------  from formatted input  --------    ############
# By Vicente Fajardo

print('\n\n')
print('############    -----------   Run ClusTCR    ------------    ############\n')
print('############    ---------  from formatted input  --------    ############\n')

# Long version: 0.1
# Version: 0
# Subversion: 1
# Updates.
# ---> Version updates:
#   No stable version achieved yet.
# ---> Subversion updates:
#   * First subversion.


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
from optparse import OptionParser
from clustcr import Clustering
print('Dependencies loaded...\n')


print('\n\n')
print('### --------------------------- Arguments --------------------------- ###\n')
### --------------------------- Arguments --------------------------- ###
# Load options.
parser = OptionParser()
parser.add_option("-R", "--ReportsPath", dest="reports_path", type="string", default=None, help="Str, absolute path to directory to save results.\n")
parser.add_option("-I", "--InputFile", dest="input_file", type="string", default=None, help="Optional. Str, absolute path to input file listing the node attributes. Must be a csv file with the column \"general.clone.id\" and any other set of columns listing further attributes.\n")
parser.add_option("-M", "--Method", dest="clust_method", type="string", default="two-step", help="Str, method to be used for clustering. Valid values taken by clusTCR are 'two-step', 'faiss' and 'mcl', with 'two-step' as the default.\n")
parser.add_option("-E", "--Expansion", dest="exp_hyper", type="float", default=2, help="Int, value for MCL expansion hyperparameter.\n")
parser.add_option("-N", "--Inflation", dest="inf_hyper", type="float", default=1.2, help="Int, value for MCL inflation hyperparameter.\n")
parser.add_option("-A", "--UseAlpha", dest="use_alpha", default=True, help="Optional. Bool, indicates whether the alpha chain CDR3 aminoacid sequence has to be considered during clustering or not.\n")
(options, args) = parser.parse_args()
# Move arguments to single variables.
reports_path = options.reports_path
input_file = options.input_file
clust_method = options.clust_method
exp_hyper = options.exp_hyper
inf_hyper = options.inf_hyper
use_alpha = options.use_alpha
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
print('Path to save reports to: ' + reports_path + '\n')
print('Path to input file: ' + input_file + '\n')
print('Clustering method: ' + clust_method + '\n')
print('Use alpha: ' + str(use_alpha) + '\n')


### --------------------------- Functions --------------------------- ###


print('\n\n')
print('### --------------------------- Load data --------------------------- ###\n')
### --------------------------- Load data --------------------------- ###
# Load input TCR data.
tcr_data = pd.read_table(input_file)


# print('\n\n')
# print('### ---------------------- Data preprocessing ----------------------- ###\n')
### ---------------------- Data preprocessing ----------------------- ###


print('\n\n')
print('### ------------------------- Main program -------------------------- ###\n')
### ------------------------- Main program -------------------------- ###

# ---> Perform clustering.
alpha_col = 'cdr3a.aa.seq' if use_alpha else None # Whether alpha chain should be considered during clusering or not.
mcl_hypers = [inf_hyper, exp_hyper]
clust_obj = Clustering(
    method=clust_method,
    n_cpus='all',
    mcl_params=mcl_hypers
)
output = clust_obj.fit(
    data=tcr_data,
    include_vgene=True,
    cdr3_col='cdr3b.aa.seq',
    alpha=alpha_col,
    v_gene_col='trb.v'
)

# ---> Summarize results and save
tmp_data_1 = output.clusters_df.copy()
tmp_data_2 = output.summary()
tmp_data_2['cluster'] = tmp_data_2.index
tmp_data = tmp_data_1.merge(
    tmp_data_2,
    on='cluster', how='outer'
)
tmp_data = tmp_data.sort_values(by=['cluster', 'junction_aa'])
tmp_file_name = '/'.join([reports_path, 'clusTCR_Clusters.csv'])
tmp_data.to_csv(tmp_file_name, index=False)

print('\n\n\n')
print('Process complete!\n\n\n')