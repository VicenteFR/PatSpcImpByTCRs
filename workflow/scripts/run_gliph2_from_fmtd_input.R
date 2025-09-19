############    ------------   Run GLIPH2    ------------    ############
############    ---------  from formatted input  --------    ############

# By Vicente Fajardo

# Long version: 0.1
# Version: 0
# Subversion: 1
# Updates.
# ---> Version updates:
#   No stable version achieved yet.
# ---> Subversion updates:
#   * First subversion.


### -------------------------- Description -------------------------- ###
# TBD


cat('\n\n')
### --------------------------- Libraries --------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
cat('Importing libraries...\n\n')
library(stringr)
library(data.table)
library(optparse)
source('/home/vfajardo/scripts/functions/R_handy_functions.0.4.R')
source('/home/vfajardo/scripts/functions/R_visualization_functions.1.5.R')
cat('Libraries imported!\n')


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
    # Server to work on.
    make_option(opt_str="--Server", type="character", default='slurm', dest="server.type", help="Char, indicates sort of server to work on. Available valid options are 'torque' and 'slurm'.\n"),
    # Input- and output-related paths and files.
    make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Char, absolute path to directory to save DE genes table.\n"),
    make_option(opt_str="--RunID", type="character", default=NULL, dest="run.id", help="ID for the comprehensive analysis to be performed.\n"),
    make_option(opt_str="--InputFile", type="character", default=NULL, dest="input.file", help="Char, absolute path formatted input file.\n"),
    make_option(opt_str="--HLAData", type="character", default=NULL, dest="hla.data.file", help="Char, absolute path to file listing the HLA data from all donors listed across the datasets provided in 'InputTable'.\n"),
    # @ Tool-specific parameters
    make_option(opt_str="--RefSpecies", type="character", default='human', dest="ref.species", help="Char, reference species, either human or mouse.'\n"),
    make_option(opt_str="--RefVersion", type="character", default='1.0', dest="ref.version", help="Char, reference version to be used during the GLIPH2 analysis. Valid values are '1.0' and '2.0'\n"),
    make_option(opt_str="--RefCellType", type="character", default='CD4', dest="ref.cell.type", help="Char, reference cell type to be used during the GLIPH2 analysis. Valid values are 'CD4' and 'CD8'.\n"),
    make_option(opt_str="--LocalMinP", type="numeric", default=0.001, dest="local.min.p", help="Numeric, value for local minimum P-value hyperparameter.\n"),
    make_option(opt_str="--KMin", type="numeric", default=3, dest="k.min.depth", help="Numeric, value for K-mer minimum depth hyperparameter.\n"),
    make_option(opt_str="--JobMem", type="character", default='50', dest="job.mem", help="Char, memory to be allotted to the tool analysis job.\n"),
    make_option(opt_str="--TimeLimit", type="integer", default=60, dest="time.limit", help="Int, maximum time (in minutes) to allow tool to run. If tool run is not complete before this limit, the program will end with an error, most likely because the tool did not run properly.\n"),
    make_option(opt_str="--TMP", type="character", default='/mnt/hpcscratch', dest="tmp.path", help="Character, indicates absolute path to a directory in your system to be used as a temporary directory. It must have a short string length overall for GLIPH2 to work properly.\n"),
    # @ Bonafide specificity group definition.
    make_option(opt_str="--PttnLength", type="integer", default=3, dest="pttn.len", help="Int, minimum number of characters that a bonafide specificity group's pattern must have.\n"),
    make_option(opt_str="--FisherPThold", type="numeric", default=0.05, dest="fisher.p.thold", help="Float, Fisher score (P-value) threshold for specificity group filtering.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
server.type <- opt$server.type
reports.path <- opt$reports.path
run.id <- opt$run.id
input.file <- opt$input.file
hla.data.file <- opt$hla.data.file
ref.species <- opt$ref.species
ref.version <- opt$ref.version
ref.cell.type <- opt$ref.cell.type
local.min.p <- opt$local.min.p
k.min.depth <- opt$k.min.depth
job.mem <- opt$job.mem
time.limit <- opt$time.limit
tmp.path <- opt$tmp.path
pttn.len <- opt$pttn.len
fisher.p.thold <- opt$fisher.p.thold
# ---> Define constant vars.
slurm.header.pttn <- '#!/bin/sh
#SBATCH --job-name={JOB_NAME}
#SBATCH --output={OUT_FILE}
#SBATCH --error={ERR_FILE}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem={MEM}g
#SBATCH --time={TIME}'
torque.header.pttn <- '#PBS -N {JOB_NAME}
#PBS -o {OUT_FILE}
#PBS -e {ERR_FILE}
#PBS -m abe
#PBS -q default
#PBS -l nodes=1:ppn=2
#PBS -l mem={MEM}gb
#PBS -l walltime={TIME}'
job.header.pttn <- if(server.type=='slurm') slurm.header.pttn else torque.header.pttn


### --------------------------- Functions --------------------------- ###


# cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
# cat('### ------------------------- Data Loading ------------------------- ###\n')


# cat('\n\n')
### ---------------------- Data Preprocessing ---------------------- ###
# cat('### ---------------------- Data Preprocessing ---------------------- ###\n')


cat('\n\n')
### --------------------------- Run tool --------------------------- ###
cat('### --------------------------- Run tool --------------------------- ###\n\n')

# ---> Hardcoded parameters.
tool.path <- '/home/vfajardo/bin/GLIPH2/irtools.centos'
ref.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/general_data/gliph2'
# ---> Prepare output directory.
# @ Jobs path.
jobs.path <- paste0(reports.path, '/jobs_scripts'); if(!dir.exists(jobs.path)) dir.create(jobs.path)
jobs.path <- paste0(jobs.path, '/', run.id)
if(!dir.exists(jobs.path)) dir.create(jobs.path)
# @ Temporary output path.
tmp.user <- system(command='echo $USER', intern=TRUE)
tmp.path <- paste0(tmp.path, '/', tmp.user, '/tmp')
if(!dir.exists(tmp.path)) dir.create(tmp.path)
if(!dir.exists(tmp.path)) stop('Could not create temporary output directory.\n')
tmp.path <- paste0(tmp.path, '/gliph2_analysis_', run.id, str_replace(string=Sys.time(), pattern=' ', replacement='_'), '_', paste0(sample(x=0:9, size=10), collapse='-'))
if(!dir.exists(tmp.path)) dir.create(tmp.path)
if(!dir.exists(tmp.path)) stop('Could not create temporary output directory.\n')
# ---> Prepare input files.
# @ Reference files.
ref.path <- paste0(ref.path, '/', ref.species, '_v', ref.version)
ref.seq.file <- paste0(ref.path, '/ref_', ref.cell.type, '_v', ref.version, '.txt')
ref.v.file <- paste0(ref.path, '/ref_V_', ref.cell.type, '_v', ref.version, '.txt')
ref.len.file <- paste0(ref.path, '/ref_L_', ref.cell.type, '_v', ref.version, '.txt')
ref.files <- c(ref.seq.file, ref.v.file, ref.len.file)
tmp.check <- all(file.exists(ref.files))
if(!tmp.check){
    ref.files <- ref.files[!file.exists(ref.files)]
    tmp.err <- paste0('Following reference files could not be found:\n', paste0(ref.files, collapse='\n'), '\nMake sure to provide a valid reference species (either human or mouse), reference version (1.0 and 2.0 for human and 1.0 for mouse) and cell type (either CD4, CD8 or a combination of both, \'CD48\').\n')
    stop(tmp.err)
}
# @ Input CDR3 seq. file.
cdr3.file <- paste0(tmp.path, '/Input_CDR3s.tsv')
file.copy(from=input.file, to=cdr3.file)
# @ Input HLA data.
if(!is.null(hla.data.file)){
    hla.file <- paste0(tmp.path, '/Input_HLA-Data.tsv')
    file.copy(from=input.hla.file, to=hla.file)
}
# @ Parameters file.
# Define content.
param.content <- paste0(
    '# Parameters file for project: ', run.id, '\n',
    'out_prefix=GLIPH2\n',
    'cdr3_file=', cdr3.file, '\n'
)
if(!is.null(hla.data.file)) param.content <- paste0(param.content, 'hla_file=', hla.file, '\n')
param.content <- paste0(
    param.content,
    'refer_file=', ref.seq.file, '\n',
    'v_usage_freq_file=', ref.v.file, '\n',
    'cdr3_length_freq_file=', ref.len.file, '\n',
    'local_min_pvalue=', local.min.p, '\n',
    'p_depth=1000\n',
    'global_convergence_cutoff=1\n',
    'simulation_depth=1000\n',
    'kmer_min_depth=', k.min.depth, '\n',
    'local_min_OVE=10\n',
    'algorithm=GLIPH2\n',
    'all_aa_interchangeable=1'
)
# Output.
param.file <- paste0(tmp.path, '/ParamsFile.txt')
write(x=param.content, file=param.file)
# ---> Define job script.
# @ Job header.
job.sh.file <- paste0(jobs.path, '/', run.id, '.job.sh')
job.out.file <- paste0(jobs.path, '/', run.id, '.out.txt')
job.err.file <- paste0(jobs.path, '/', run.id, '.err.txt')
job.header <- str_replace(string=job.header.pttn, pattern='\\{JOB_NAME\\}', replacement=run.id)
job.header <- str_replace(string=job.header, pattern='\\{OUT_FILE\\}', replacement=job.out.file)
job.header <- str_replace(string=job.header, pattern='\\{ERR_FILE\\}', replacement=job.err.file)
job.header <- str_replace(string=job.header, pattern='\\{MEM\\}', replacement=as.character(job.mem))
job.header <- str_replace(string=job.header, pattern='\\{TIME\\}', replacement=as.character(time.limit))
# @ Job body.
job.body <- paste0(
    'echo -e "Preparing to run analysis.\\n"\n',
    'cp ', tool.path, ' ', tmp.path, '\n',
    'cd ', tmp.path, '\n',
    'echo -e "Start of analysis.\\n"\n',
    './', basename(tool.path), ' -c ', param.file, '\n',
    'echo -e "Analysis complete!\\n"\n',
    'echo -e "Job finished.\\n"\n'
)
# @ Output job script.
job.content <- paste0(job.header, '\n', job.body)
write(file=job.sh.file, x=job.content)
# @ If everything went fine, proceed to run the job script.
submit.cmmd <- if(server.type=='slurm') paste0('sbatch ', job.sh.file) else paste0('qsub ', job.sh.file)
if(file.exists(job.sh.file)) system(command=submit.cmmd)

# ---> Wait until the process is complete.
tmp.file.name <- paste0(
  tmp.path,
  '/GLIPH2_cluster.csv'
)
time.count <- 0
while(!file.exists(tmp.file.name)){
  time.count <- time.count + 1
  if(time.count > time.limit){
    tmp.error <- paste0('Tool did not run properly within the allotted time. Please check whether any relevant messages from the tool were provided in the following path: ', analysis.path, '/jobs_scripts', '\n')
    stop(tmp.error)
  }
  Sys.sleep(time=60)
}
cat('GLIPH2 output ready. Putting output files in order...\n')
# Wait one more minute to make sure final output is fully written to the corresponding file.
Sys.sleep(time=60)
# Put output files in order.
tmp.files <- list.files(path=tmp.path, full.names=TRUE)
file.copy(from=tmp.files, to=reports.path)
# Remove temporary directory.
tmp.cmmd <- paste0('rm -r ', tmp.path)
system(command=tmp.cmmd)

# ---> Internal output preprocessing.
# Load results.
tmp.file.name <- paste0(
  reports.path,
  '/GLIPH2_cluster.csv'
)
tool.res <- fread(file=tmp.file.name, blank.lines.skip=TRUE)
# Preprocessing.
# - 1. Pattern/Motif length larger than or equal to 'pattern.length'.
tool.res <- tool.res[str_length(string=pattern)>=pttn.len, ]
# - 2. Reported as statistically significant by GLIPH2.
tool.res <- tool.res[Fisher_score<=fisher.p.thold]
# Save results.
out.file.name <- paste0(
  reports.path,
  '/GLIPH2_Clusters.csv'
)
fwrite(file=out.file.name, x=tool.res)


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
