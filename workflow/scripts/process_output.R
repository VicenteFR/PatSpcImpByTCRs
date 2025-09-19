############    --------  Process output from  --------    ############
############    ---  TCR clustering by similarity  ----    ############

# By Vicente Fajardo

# Version: 0.1
# Version updates:
#   No stable version yet.
# Subversion: 1
# Subversion updates:
#     0.1
#       First subversion. Template taken: TCRdist3 script v0.2.



### -------------------------- Description -------------------------- ###
# TBD


cat('\n\n')
### --------------------------- Libraries --------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
cat('Importing libraries...\n\n')
library(stringr)
library(data.table)
library(optparse)
cat('Libraries imported!\n')


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
    # Server to work on.
    make_option(opt_str="--Server", type="character", default='slurm', dest="server.type", help="Char, indicates sort of server to work on. Available valid options are 'torque' and 'slurm'.\n"),
    # @ Input- and output-related paths and files.
    make_option(opt_str="--Tool", type="character", default=NULL, dest="tool.name", help="Char, tool name among GLIPH2, TCRdist3, clusTCR, GIANA, iSMART and ALICE.\n"),
    make_option(opt_str="--RunID", type="character", default=NULL, dest="run.id", help="ID for the comprehensive analysis to be performed.\n"),
    make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Char, absolute path to directory to save DE genes table.\n"),
    make_option(opt_str="--ToolOutput", type="character", default=NULL, dest="tool.output.file", help="Char, absolute path to input file, which should correspond to an output file from any of the supported tools for TCR clustering by similarity (so far, GLIPH2, TCRdist3, clusTCR, GIANA, iSMART, ALICE).\n"),
    make_option(opt_str="--RefInfo", type="character", default=NULL, dest="ref.info.file", help="Char, absolute path to reference info file as produced while preprocessing individual TCR files.\n"),
    # @ Filtering-related parameters for the tool input.
    make_option(opt_str="--TagsOfInt", type="character", default=NULL, dest="tags.of.int", help="Char, a set of tags that should be carried during the whole analysis. The tags provided must be properly defined in the metadata of all datasets. The format required is a string with the names of the tags separated by semicolons; for example, if we are interested in both tags 'virus.tag' and 'pr.tag', we must provide the single string 'virus.tag;pr.tag'.\n"),
    # @ Bonafide specificity group definition.
    make_option(opt_str="--DonorMinCount", type="integer", default=3, dest="donor.min.count", help="Int, minimum number of donors that the TCRs from a bonafide specificity groups must come from.\n"),
    # @ NDEx cloud-specific information.
    # If user would like their clonograph to be created and saved to the cloud, valid NDEx user-specific information must be provided.
    make_option(opt_str="--NDExServer", type="character", default=NULL, dest="ndex.server", help="Char, NDEx server to connect to the cloud. If provided, in order to establish a connection with NDEx, the program also requires the complementary information.\n"),
    make_option(opt_str="--NDExUser", type="character", default=NULL, dest="ndex.user", help="Char, NDEx user to connect to the cloud. If provided, in order to establish a connection with NDEx, the program also requires the complementary information.\n"),
    make_option(opt_str="--NDExPswd", type="character", default=NULL, dest="ndex.pswd", help="Char, NDEx password to connect to the cloud. If provided, in order to establish a connection with NDEx, the program also requires the complementary information.\n"),
    make_option(opt_str="--JobMem", type="character", default='50', dest="job.mem", help="Char, memory to be allotted to the tool analysis job.\n"),
    make_option(opt_str="--JobCPUs", type="character", default='2', dest="job.cpus", help="Char, CPUs to be allotted to the tool analysis job.\n"),
    make_option(opt_str="--TimeLimit", type="integer", default=60, dest="time.limit", help="Int, maximum time (in minutes) to allow tool to run. If tool run is not complete before this limit, the program will end with an error, most likely because the tool did not run properly.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
server.type <- opt$server.type
tool.name <- opt$tool.name
reports.path <- opt$reports.path
run.id <- opt$run.id
tool.output.file <- opt$tool.output.file
ref.info.file <- opt$ref.info.file
tags.of.int <- opt$tags.of.int
if(!is.null(tags.of.int)) tags.of.int <- str_split(string=tags.of.int, pattern=';', simplify=TRUE)[1, ]
if('donor.id.tag' %in% tags.of.int) tags.of.int <- setdiff(x=tags.of.int, y='donor.id.tag') # Sanity check. The donor identity will be mandatorily conserved during the process, so no need to carry it as a tag of interest.
donor.min.count <- opt$donor.min.count
job.mem <- opt$job.mem
job.cpus <- opt$job.cpus
time.limit <- opt$time.limit
ndex.server <- opt$ndex.server
ndex.user <- opt$ndex.user
ndex.pswd <- opt$ndex.pswd
tmp.check <- (is.null(ndex.server) & is.null(ndex.server) & is.null(ndex.server)) | (!is.null(ndex.server) & !is.null(ndex.server) & !is.null(ndex.server))
if(!tmp.check) stop('When information for NDEx is to be provided, user must provide all three items to successfully establish a connection to NDEx\'s cloud, specifically: server, user and password. Otherwise, please leave all three flags out of your code line.\n')
# ---> Define constant vars.
slurm.header.pttn <- '#!/bin/sh
#SBATCH --job-name={JOB_NAME}
#SBATCH --output={OUT_FILE}
#SBATCH --error={ERR_FILE}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem={MEM}g
#SBATCH --time={TIME}'
torque.header.pttn <- '#PBS -N {JOB_NAME}
#PBS -o {OUT_FILE}
#PBS -e {ERR_FILE}
#PBS -m abe
#PBS -q default
#PBS -l nodes=1:ppn={CPUS}
#PBS -l mem={MEM}gb
#PBS -l walltime={TIME}'
job.header.pttn <- if(server.type=='slurm') slurm.header.pttn else torque.header.pttn


### --------------------------- Functions --------------------------- ###


cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')

# ---> Input and output of tool for TCR clustering by similarity.
tool.res <- if(tool.name == 'GIANA') fread(file=tool.output.file, skip=2, header=FALSE, fill=TRUE) else fread(file=tool.output.file, blank.lines.skip=TRUE)
# ---> Reference info.
ref.info <- fread(file=ref.info.file)
# Define clonotype ID.
ref.info[, general.clone.id:=1:.N]
max.clone.len <- ref.info[, max(str_length(string=general.clone.id))]
ref.info[, tmp.val:=max.clone.len - str_length(string=general.clone.id)]
ref.info$tmp.val <- sapply(X=ref.info[, tmp.val], FUN=function(x) return(paste0(rep(x=0, times=x), collapse='')))
ref.info[, general.clone.id:=paste0(tmp.val, general.clone.id)]
ref.info[, tmp.val:=NULL]
ref.info[, general.clone.id:=paste0('gen.clone.', general.clone.id)]


cat('\n\n')
### ---------------------- Data Preprocessing ---------------------- ###
cat('### ---------------------- Data Preprocessing ---------------------- ###\n')

# ---> Reference-specific modifications.
ref.info[, `:=`(trb.v=paste0(trb.v, '*01'), trb.j=paste0(trb.j, '*01'))]

# ---> Tool output-specific modifications.

# NOTE: Being cautious and removing the size column to avoid mismatches between reference and tool output.

# GLIPH2-specific modifications.
if(tool.name == 'GLIPH2'){
    tool.res[, donor.id.tag:=str_replace(string=Sample, pattern=':Exp\\d+', replacement='')]
    tool.res[, Sample:=NULL]
}

# GIANA-specific modifications.
if(tool.name == 'GIANA'){
    # NOTE: THESE COLUMN NAMES HEAVILY DEPEND ON THE INPUT CREATED FROM THE BEGINNING. IF THE COLUMNS WERE TO PROVIDED IN A DIFFERENT ORDER, THIS WOULD HAVE TO CHANGE IN CONSEQUENCE.
    tmp.cols <- c('cdr3b.aa.seq', 'cluster.id', 'trb.v', 'trb.j', 'cdr3b.aa.seq', 'sample', 'size')
    colnames(tool.res) <- tmp.cols
    tool.res[, donor.id.tag:=str_replace(string=sample, pattern=':Exp\\d+', replacement='')]
    tool.res[, `:=`(sample=NULL, size=NULL)]
}

# TCRdist3-specific modifications.
if(tool.name == 'TCRdist3'){
    tool.res[, count:=NULL] # Removed because info is redundant.
}

# clusTCR-specific modifications.
if(tool.name == 'clusTCR'){
    tool.res[, size:=NULL] # Removed because info is redundant.
}

# iSMART-specific modifications
if(tool.name == 'iSMART'){
    tool.res[, donor.id.tag:=str_replace(string=donor.id, pattern=':Exp\\d+', replacement='')]
    tool.res[, `:=`(size=NULL, donor.id=NULL)]
}

cat('\n\n')
### ---------------------- Results processing ---------------------- ###
cat('### ---------------------- Results processing ---------------------- ###\n\n')


# ---> Process results.

# -- @ Set standard names to relevant columns.
out.col.names <- list(
    `ALICE` = c(
        `cluster.id`='D',
        `cdr3b.aa.seq`='CDR3.amino.acid.sequence',
        `trb.v`='bestVGene',
        `trb.j`='bestJGene'
    ),
    `clusTCR` = c(
        `cluster.id`='cluster',
        `cdr3b.aa.seq`='junction_aa',
        `trb.v`='v_call'
    ),
    `TCRdist3` = c(
        `cluster.id`='cluster'
    ),
    `GIANA` = c(
        `cluster.id`='cluster.id' # Mock
    ),
    `GLIPH2` = c(
        `cluster.id`='index',
        `cdr3b.aa.seq`='TcRb',
        `trb.v`='V',
        `trb.j`='J'
    ),
    `iSMART` = c(
        `cluster.id`='Group'
    )
)
# Column renaming according to tool-specific output.
tmp.cols <- out.col.names[[tool.name]]
cols.to.keep <- c(tmp.cols, setdiff(x=colnames(tool.res), y=tmp.cols))
names(cols.to.keep) <- c(names(tmp.cols), setdiff(x=colnames(tool.res), y=tmp.cols))
tool.res <- tool.res[, ..cols.to.keep]; colnames(tool.res) <- names(cols.to.keep)

# -- @ Include in-tool process discarded information according to tool.
merge.cols <- list(
    ALICE = c('cdr3b.aa.seq', 'trb.v', 'trb.j'),
    clusTCR = c('cdr3b.aa.seq', 'trb.v'),
    TCRdist3 = c('cdr3b.aa.seq', 'trb.v', 'trb.j'),
    GIANA = c('cdr3b.aa.seq', 'trb.v', 'trb.j', 'donor.id.tag'),
    GLIPH2 = c('cdr3b.aa.seq', 'trb.v', 'trb.j', 'donor.id.tag'),
    iSMART = c('cdr3b.aa.seq', 'trb.v', 'trb.j', 'donor.id.tag')
)
tmp.cols <- c('general.clone.id', 'cdr3b.aa.seq', 'trb.v', 'trb.j', 'size', 'donor.id.tag')
tmp.data <- ref.info[, ..tmp.cols]
tool.res <- merge(x=tmp.data, y=tool.res, by=merge.cols[[tool.name]], all.x=FALSE, all.y=TRUE, sort=FALSE)

# -- @ Further preprocessing.
# Keep only bonafide specificity groups as defined by the following features:
# - 0. Must not be reported as single.
tmp.data <- tool.res[, .(seq.count=uniqueN(cdr3b.aa.seq)), by=cluster.id]
tmp.data <- tmp.data[seq.count>1]
tool.res <- tool.res[cluster.id %in% tmp.data[, cluster.id]]
# - 1. Must contain >=1 clonotypes from larger than or equal to 'donor.min.count' donors.
tmp.data <- tool.res[, .(donor.count=uniqueN(donor.id.tag)), by=.(cluster.id)]
tmp.data <- tmp.data[donor.count>=donor.min.count, cluster.id]
tool.res <- tool.res[cluster.id %in% tmp.data]

# -- @ Include remaining important information.
tmp.cols <- c('cdr3b.aa.seq', 'trb.v', 'trb.j', 'donor.id.tag', 'clonotype.tag', tags.of.int)
if(any(tags.of.int %in% colnames(tool.res))) tmp.cols <- setdiff(x=tmp.cols, y=tags.of.int)
tmp.data <- ref.info[, ..tmp.cols]
tool.res <- merge(x=tmp.data, y=tool.res, by=c('cdr3b.aa.seq', 'trb.v', 'trb.j', 'donor.id.tag'), all.x=FALSE, all.y=TRUE, sort=FALSE)

# -- @ Processed results table.
tmp.cols <- c('cluster.id', setdiff(x=colnames(tool.res), y='cluster.id'))
tool.res <- tool.res[, ..tmp.cols]
setorderv(x=tool.res, cols=c('cluster.id', 'size'), order=c(1, -1))
tmp.file.name <- paste0(
    reports.path,
    '/ToolResults_Processed.csv'
)
fwrite(file=tmp.file.name, x=tool.res, na=NA, quote=FALSE)


# ---> Create network adjacency list.
uniq.clones <- sort(unique(tool.res[, general.clone.id]))
adj.list <- lapply(X=1:length(uniq.clones), FUN=function(clone.id){
    # Define clone IDs that were considered in previous iterations as well as the clone ID for the current iteration.
    prev.clone.ids <- uniq.clones[1:(clone.id-1)]
    clone.id <- uniq.clones[clone.id]
    # Specificity groups where that the clone ID belongs in.
    sgs <- tool.res[general.clone.id==clone.id, unique(cluster.id)]
    # Define interactions with clone IDs that have not been considered in previous iterations.
    clone.ids.to.discard <- c(prev.clone.ids, clone.id)
    tmp.data <- tool.res[cluster.id %in% sgs & !general.clone.id %in% clone.ids.to.discard, .(weight=.N), by=.(clone.id.2=general.clone.id)]
    return(tmp.data)
})
names(adj.list) <- uniq.clones
adj.list <- rbindlist(l=adj.list, use.names=TRUE, idcol='clone.id.1')
adj.list.file <- paste0(
    reports.path,
    '/ToolResults_AdjList-1.csv'
)
fwrite(file=adj.list.file, x=adj.list, na=NA, quote=FALSE)


# ---> Define clonotype-specific (node) attributes.
if(is.null(tags.of.int)){
    tags.of.int <- 'donor.id.tag'
}else{
    tags.of.int <- c('donor.id.tag', tags.of.int)
}
node.atts <- lapply(X=tags.of.int, FUN=function(tmp.tag){
    tool.res[, tmp.tag:=NULL]
    tool.res[, tmp.tag:=get(tmp.tag)]
    node.atts <- tool.res[,
        .(
            tcr.b.aa=unique(cdr3b.aa.seq),
            v.gene=unique(trb.v),
            j.gene=unique(trb.j),
            clonotype.tags=paste0(sort(unique(clonotype.tag)), collapse=';'),
            clone.size=.SD[, .(size=unique(size)), by=clonotype.tag][, sum(size)],
            tmp.tag=.SD[tmp.tag!='', paste0(sort(unique(tmp.tag)), collapse=';')]
        ),
        by=.(
            general.clone.id
        )
    ]
    colnames(node.atts)[colnames(node.atts)=='tmp.tag'] <- tmp.tag
    return(node.atts)
})
if(length(tags.of.int) == 1){
    tags.of.int <- NULL
    node.atts <- node.atts[[1]]
}else{
    node.atts <- Reduce(x=node.atts, f=function(x, y){
        merge(x=x, y=y, by=c('general.clone.id', 'tcr.b.aa', 'v.gene', 'j.gene', 'clone.size', 'clonotype.tags'))
    })
}
# Sanity check: at this point, we should have unique entries for clonotypes.
tmp.check <- node.atts[, .N]==node.atts[, uniqueN(general.clone.id)]
if(!tmp.check) stop('Unexpected error 1.\n')
# Output attributes.
node.atts.file <- paste0(
    reports.path,
    '/ToolResults_ClonotypeAtts-1.csv'
)
fwrite(file=node.atts.file, x=node.atts, na=NA, quote=FALSE)

# ---> Create graph and load to NDEx's cloud through a separate python script.
# This process will be carried out only if all necessary information to establish a successful connection to the cloud has been provided.
tmp.check <- !is.null(ndex.server) & !is.null(ndex.user) & !is.null(ndex.pswd)
if(tmp.check){
    # @ Hardcoded parameters.
    python.script <- '/home/vfajardo/scripts/TCR_data_analysis/general_clustering/GLIPH2/load_gliph2_network_to_ndex.0.1.py'
    tmp.user <- system(command='echo $USER', intern=TRUE)
    # @ Jobs path.
    jobs.path <- paste0(tool.path, '/jobs_scripts'); if(!dir.exists(jobs.path)) dir.create(jobs.path)
    jobs.path <- paste0(jobs.path, '/', run.id)
    if(!dir.exists(jobs.path)) dir.create(jobs.path)
    # @ Acual reports path.
    output.path <- paste0(tool.path, '/', run.id)
    # @ Define job script.
    # Job header.
    job.sh.file <- paste0(jobs.path, '/Tool-NDEx_', run.id, '.job.sh')
    job.out.file <- paste0(jobs.path, '/Tool-NDEx_', run.id, '.out.txt')
    job.err.file <- paste0(jobs.path, '/Tool-NDEx_', run.id, '.err.txt')
    job.header <- str_replace(string=job.header.pttn, pattern='\\{JOB_NAME\\}', replacement=paste0(run.id, '-NDEx'))
    job.header <- str_replace(string=job.header, pattern='\\{OUT_FILE\\}', replacement=job.out.file)
    job.header <- str_replace(string=job.header, pattern='\\{ERR_FILE\\}', replacement=job.err.file)
    job.header <- str_replace(string=job.header, pattern='\\{MEM\\}', replacement=as.character(job.mem))
    job.header <- str_replace(string=job.header, pattern='\\{TIME\\}', replacement=as.character(time.limit))
    # Job body.
    job.body <- paste0(
        'source activate env_for_nx\n',
        'echo -e "\\n\\n"\n',
        'echo -e "Running script to load tool graph to NDEx...\\n"\n',
        paste('python', python.script, '--ReportsPath', output.path, '--AdjListFile', adj.list.file, '--NodeAttsFile', node.atts.file, '--NDExServer', ndex.server, '--NDExUser', ndex.user, '--NDExPswd', ndex.pswd, '--NDExGraphID', paste0('\'', run.id, '\''), sep=' '), '\n',
        'echo -e "Job done!\\n"\n'
    )
    # Output job script.
    job.content <- paste0(job.header, '\n', job.body)
    write(file=job.sh.file, x=job.content)
    # If everything went fine, proceed to run the job script.
    submit.cmmd <- if(server.type=='slurm') paste0('sbatch ', job.sh.file) else paste0('qsub ', job.sh.file)
    if(file.exists(job.sh.file)) system(command=submit.cmmd)
}


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
