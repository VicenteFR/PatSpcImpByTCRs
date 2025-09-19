############    -------   Specificity prediction    -------    ############
cat('\n\n')
cat('############    -------   Specificity prediction    -------    ############\n')

# By: Vicente Fajardo-Rosas

### -------------------------- Description -------------------------- ###


cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
library(data.table)
library(tidyr)
library(fuzzyjoin)
library(stringr)
library(ggplot2)
library(optparse)


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
    make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Char, absolute path to directory to save results to.\n"),
    make_option(opt_str="--RunPath", type="character", default=NULL, dest="run.path", help="Char, absolute path to directory where clustering results have been saved to.\n"),
    make_option(opt_str="--RefDefFile", type="character", default=NULL, dest="input.table.file", help="Char, indicates the reference ID for the results to be evaluated.\n"),
    make_option(opt_str="--RefInfo", type="character", default=NULL, dest="ref.info.file", help="Char, absolute path to reference info file as produced while preprocessing individual TCR files.\n"),
    make_option(opt_str="--OptsFile1", type="character", default=NULL, dest="yaml.file.1", help="Char, absolute path to yaml file for workflow; general options.\n"),
    make_option(opt_str="--OptsFile2", type="character", default=NULL, dest="yaml.file.2", help="Char, absolute path to yaml file for workflow; project-specific options.\n"),
    make_option(opt_str="--RefID", type="character", default=NULL, dest="ref.name", help="Char, indicates the reference ID for the results to be evaluated.\n"),
    make_option(opt_str="--ExpansionThold", type="numeric", default=NULL, dest="expansion.thold", help="Int, expansion threshold that was set for clonality.\n"),
    make_option(opt_str="--UCMThold", type="numeric", default=NULL, dest="ucm.thold", help="Int, weight threshold for specificity prediction.\n"),
    make_option(opt_str="--VDJMatch", type="logical", default=FALSE, dest="vdj.match", help="Logical, indicates whether VDJ match criterion should be taken into account for matches between inquiry dataset and the reference and experimental datasets.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
reports.path <- opt$reports.path
run.path <- opt$run.path
input.table.file <- opt$input.table.file
ref.info.file <- opt$ref.info.file
yaml.file.1 <- opt$yaml.file.1
yaml.file.2 <- opt$yaml.file.2
ref.name <- opt$ref.name
expansion.thold <- opt$expansion.thold
ucm.thold <- opt$ucm.thold
vdj.match <- opt$vdj.match
# ---> Define constant vars.
set.seed(seed=1)

# ---> Options file.
# Define files from general input table.
input.table <- fread(file=input.table.file, quote='\'')
tcr.meta.file <- input.table[dataset.lab!=ref.name, meta.data]
# tcr.meta.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/HIPC/paper_developments/HIPC-Lung/preliminary_process/Batches-1-to-8_2024-10-07/CellBasedTCRData_ExVivo.csv'
cell.clons.info.file <- if('filtered_contig_annotations' %in% colnames(input.table)) input.table[dataset.lab!=ref.name, filtered_contig_annotations] else NA
# Load yaml files. Priority is given to dataset-specific options while hanlding potentially redundant entries.
if(!file.exists(yaml.file.1)) stop(paste0('Attempted definition of yaml file (general options) failed. Next file could not be found:\n', yaml.file.1, '\n'))
if(!file.exists(yaml.file.2)) stop(paste0('Attempted definition of yaml file (project-specific options) failed. Next file could not be found:\n', yaml.file.2, '\n'))
tmp.data.1 <- yaml::read_yaml(file=yaml.file.1)
tmp.data.2 <- yaml::read_yaml(file=yaml.file.2)
gen.wflow.opts <- c(tmp.data.2, tmp.data.1[!names(tmp.data.1) %in% names(tmp.data.2)])
# Define files and general options.
tool.names <- names(gen.wflow.opts$tool_opts)
w.vector <- unlist(gen.wflow.opts$weight_vector)
qry.donor.tag <- gen.wflow.opts$qry_specs$donor_id_tag
ref.meta.file <- gen.wflow.opts$ref_specs[[ref.name]]$ref_data
t.subset <- gen.wflow.opts$ref_specs[[ref.name]]$t_subset
clone.id.tag <- gen.wflow.opts$ref_specs[[ref.name]]$clonotype_id_tag
score.tag <- gen.wflow.opts$ref_specs[[ref.name]]$score_tag
ag.spc.tag <- gen.wflow.opts$ref_specs[[ref.name]]$ag_spc_tag
cdr3b.aa.tag <- gen.wflow.opts$ref_specs[[ref.name]]$cdr3b_aa_tag
cdr3a.aa.tag <- gen.wflow.opts$ref_specs[[ref.name]]$cdr3a_aa_tag
trb.v <- gen.wflow.opts$ref_specs[[ref.name]]$trb_v
trb.j <- gen.wflow.opts$ref_specs[[ref.name]]$trb_j
tra.v <- gen.wflow.opts$ref_specs[[ref.name]]$tra_v
tra.j <- gen.wflow.opts$ref_specs[[ref.name]]$tra_j
tags.of.int <- gen.wflow.opts$ref_specs[[ref.name]]$tags_of_int
tags.of.int <- str_split(string=tags.of.int, pattern=';', simplify=TRUE)[1, ]
ags.to.keep <- gen.wflow.opts$ref_specs[[ref.name]]$ags_to_keep
ags.to.keep <- str_split(string=ags.to.keep, pattern=',', simplify=TRUE)[1, ]
exp.res.file <- gen.wflow.opts$exp_data

# ---> Hardcoded arguments.
# Info for reference match approach
rma.chain.cols <- c('RMB'='cdr3b.aa.seq', 'RMA'='cdr3a.aa.seq') # 'rm' stands for "reference match"
rma.chain.vdj <- list(
    'RMB'=c(`v`='trb.v', `j`='trb.j'),
    'RMA'=c(`v`='tra.v', `j`='tra.j')
)
# Info for donor-matched approach
dma.chain.cols <- c('DMB'='cdr3b.aa.seq', 'DMA'='cdr3a.aa.seq')


### --------------------------- Functions --------------------------- ###

# ---------------------------------------->
# Name: Separate clonotypes w/ multiple same-type chains into individual rows taking vdj info into consideration.
# Description:
# Provided a TCR metadata table, this function will split the info of a multi-chain clonotype (for either or both types of chains) into multiple rows, one per chain. In the flight, the function will take the vdj info into consideration. That is, separation occurs in a gene-aware manner.
# Arguments ------------------------------>
# trc.info - data table, lists TCR metadata. Following columns are expected (though other metadata may be included): 'cdr3b.aa.seq', 'cdr3a.aa.seq', 'trb.v', 'trb.j', 'tra.v' and 'tra.j'. Although not explicitly required (i.e., function will not check for this, but is rather expected), every row should represent an individual combination of a clonotype ID and a donor ID, and there should be a number of unique rows as there are unique clonotype-donor ID relationships.
# Value ------------------------------>
# data table, modified TCR metadata (see description).

# Function:

sep.tcr.row.info <- function(tcr.inf){
    # Basic check. When multiple sequences are provided for a chain, multiple V and J have to be provided as well, one VJ gene set for each sequence.
    tmp.check.1 <- tcr.inf[
        str_detect(string=cdr3a.aa.seq, pattern=';'),
        all(
            str_count(string=cdr3a.aa.seq, pattern=';')==str_count(string=tra.v, pattern=';') &
            str_count(string=cdr3a.aa.seq, pattern=';')==str_count(string=tra.j, pattern=';')
        )
    ]
    tmp.check.2 <- tcr.inf[
        str_detect(string=cdr3b.aa.seq, pattern=';'),
        all(
            str_count(string=cdr3b.aa.seq, pattern=';')==str_count(string=trb.v, pattern=';') &
            str_count(string=cdr3b.aa.seq, pattern=';')==str_count(string=trb.j, pattern=';')
        )
    ]
    tmp.check <- tmp.check.1 & tmp.check.2
    if(!tmp.check) stop('Faulty reference data. Different number of V-J genes provided compared w/ the number of unique sequences provided for certain clonotypes.\n')
    # Proceed to separate the pieces collectively for clonotypes w/ multiple sequences for each chain.
    #       For alpha chain.
    tmp.check <- tcr.inf[, .N] + tcr.inf[, sum(str_count(string=cdr3a.aa.seq, pattern=';'), na.rm=TRUE)]
    tcr.inf <- as.data.table(separate_rows(data=tcr.inf, cdr3a.aa.seq, tra.v, tra.j, sep=';'))
    tmp.check <- tmp.check==tcr.inf[, .N]
    if(!tmp.check) stop('Unexpected error while tidying data for reference clonotypes w/ multiple alpha chain sequences.\n')
    #       For beta chain.
    tmp.check <- tcr.inf[, .N] + tcr.inf[, sum(str_count(string=cdr3b.aa.seq, pattern=';'), na.rm=TRUE)]
    tcr.inf <- as.data.table(separate_rows(data=tcr.inf, cdr3b.aa.seq, trb.v, trb.j, sep=';'))
    tmp.check <- tmp.check==tcr.inf[, .N]
    if(!tmp.check) stop('Unexpected error while tidying data for reference clonotypes w/ multiple beta chain sequences.\n')
    return(tcr.inf)
}


cat('\n\n')
### --------------------------- Load data --------------------------- ###
cat('### --------------------------- Load data --------------------------- ###\n')

# ---> TCR reference metadata.
# Input used for clustering by similarity analysis.
ref.meta <- fread(file=ref.meta.file)
# NOTE: Next is a temporary solution to problem. Ideally, this should be solved at the reference compilation (upstream) step.
# Simplify PMID column when a full link to GitHub is provided. Otherwise, it'll mess the output up.
if('pmid' %in% colnames(ref.meta)){
    ref.meta[
        !is.na(pmid) & !pmid %like% '^\\d{7,8}$' & str_detect(string=pmid, pattern='github.com/antigenomics/vdjdb-db/issues'),
        pmid:=paste0(
            'vdj-db-issue;',
            basename(pmid)
        )
    ]
    ref.meta[
        pmid=='https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#',
        pmid:='10x-page'   
    ]
    ref.meta[
        pmid=='https://doi.org/10.1101/2020.05.04.20085779',
        pmid:='unpublished'
    ]
    # ref.meta[
    #     !is.na(pmid) & !pmid %like% '^\\d{7,8}$',
    #     unique(pmid)
    # ]
}

# ---> Specific reference info.
ref.info <- fread(file=ref.info.file)
# Define reference TCRs in processed reference info.
ref.clones <- unlist(str_split(
    string=ref.info[, clonotype.tag],
    pattern=';'
))
tmp.pttn <- paste0('\\|', ref.name, '$')
ref.clones <- ref.clones[str_detect(string=ref.clones, pattern=tmp.pttn)]
ref.clones <- str_replace(
    string=ref.clones, pattern=tmp.pttn, replacement=''
)
# Subset reference metadata to keep info only for clonotypes in processed reference.
ref.meta <- ref.meta[get(clone.id.tag) %in% ref.clones]
# If score tag is provided, filter accordingly.
# NOTE: We noticed that filtering based on clonotype IDs is not perfect because of the putative public clonotypes. Therefore, this second step of cleaning is needed. Naturally, this step is quite specific for the in-house class-II reference.
if(ref.name=='EpRef-2' & !is.null(score.tag) & score.tag=='freq.conf.lvl'){
    ref.meta <- ref.meta[get(score.tag)>=expansion.thold]
}

# ---> Query TCR metadata
if(!is.na(tcr.meta.file)){
    # When a processed file has already been provided.
    tcr.meta <- fread(file=tcr.meta.file)
    if(!qry.donor.tag %in% colnames(tcr.meta)) stop(paste0('Donor ID column must be defined for query dataset.\nValue provided: ', qry.donor.tag, '\n'))
    tmp.cols <- c(qry.donor.tag, 'barcode', 'clonotype.tag', 'cdr3b.aa.seq', 'cdr3a.aa.seq', 'trb.v', 'trb.j','tra.v', 'tra.j')
    names(tmp.cols) <- tmp.cols
    names(tmp.cols)[1] <- 'donor.id.tag'
    # if('donor.id.tag' %in% colnames(tcr.meta)) tmp.cols <- c('donor.id.tag', tmp.cols)
    cells.clons.info <- unique(tcr.meta[, ..tmp.cols])
    colnames(cells.clons.info)[colnames(cells.clons.info)=='clonotype.tag'] <- 'qry.clone.id'
}else{
    # When only the raw outputs from cellranger (or an aggr) have been provided.
    cell.clons.info <- fread(file=cell.clons.info.file)
    cell.clons.info <- cell.clons.info[, .(
        barcode, raw_clonotype_id, chain,
        cdr3, v_gene, j_gene
    )]
    cell.clons.info <- unique(cell.clons.info)
    tmp.data.1 <- unique(cell.clons.info[, .(barcode, qry.clone.id=raw_clonotype_id)])
    tmp.data.2 <- unique(cell.clons.info[, .(raw_clonotype_id, chain, cdr3, v_gene, j_gene)])
    tmp.data.2 <- tmp.data.2[,
        .(
            cdr3b.aa.seq=.SD[chain=='TRB', paste0(sort(unique(cdr3)), collapse=';')],
            cdr3a.aa.seq=.SD[chain=='TRA', paste0(sort(unique(cdr3)), collapse=';')],
            trb.v=.SD[chain=='TRB', paste0(sort(unique(v_gene)), collapse=';')],
            trb.j=.SD[chain=='TRB', paste0(sort(unique(j_gene)), collapse=';')],
            tra.v=.SD[chain=='TRA', paste0(sort(unique(v_gene)), collapse=';')],
            tra.j=.SD[chain=='TRA', paste0(sort(unique(j_gene)), collapse=';')]
        ),
        by=.(qry.clone.id=raw_clonotype_id)
    ]
    cells.clons.info <- merge(x=tmp.data.1, y=tmp.data.2, by='qry.clone.id')
}
# Unique info on a query clonotype basis.
tmp.cols <- c(
    'qry.clone.id',
    'cdr3b.aa.seq', 'cdr3a.aa.seq',
    'trb.v', 'trb.j', 'tra.v', 'tra.j'
)
if('donor.id.tag' %in% colnames(cells.clons.info)) tmp.cols <- c('donor.id.tag', tmp.cols)
tcr.meta <- unique(cells.clons.info[
    !is.na(qry.clone.id),
    ..tmp.cols
])
if(!'donor.id.tag' %in% colnames(tcr.meta)){
    tmp.check <- tcr.meta[, .N==uniqueN(qry.clone.id)]
    if(!tmp.check) stop('Something went wrong while attempting to define query TCR metadata on a clonotype basis. (1)\n')
    # Add mock donor ID.
    tcr.meta[, donor.id.tag:='Mock']
}else{
    tcr.meta[!is.na(donor.id.tag)]
}


# ---> Experimental validation data.
exp.res <- if(!is.null(exp.res.file)) fread(file=exp.res.file) else NULL

cat('All necessary input data in session.\n')


cat('\n\n')
### ------------------------- Main program -------------------------- ###
cat('### ------------------------- Main program -------------------------- ###\n')

cat('\n\n')
### -------------- Reference chain-specific clonotypes -------------- ###
cat('### -------------- Reference chain-specific clonotypes -------------- ###\n')

# ---
# Clonotype calling for the sequences in our references when performed at the chain level. That is, a unique set of CDR3 *aminoacid* sequence and v and j genes (for either of the chains) will define a unique clonotype and, accordingly, unique and unambiguous metadata entries should be defined. Ambiguous metadata entries might come up as a matter of clonotype calling being performed considering both clonotype-related items and donor ID. To overcome this, priority is given to clonotypes according to score provided by user (if any), so that reference clonotypes with the same defining set are ranked and only the top set reference clonotype is considered. If still multiple ambiguous reference clonotypes come up as a top tie for a given set, that set will be discarded.
# ---

# ---> Retrieve reference TCR data.
# Retrieve basic reference metadata.
tmp.cols <- c(
    clone.id.tag, cdr3b.aa.tag, cdr3a.aa.tag, trb.v, trb.j, tra.v, tra.j, score.tag
)
tmp.data.1 <- ref.meta[, ..tmp.cols]
colnames(tmp.data.1)[colnames(tmp.data.1)==clone.id.tag] <- 'ref.clone.id'
colnames(tmp.data.1)[colnames(tmp.data.1)==cdr3b.aa.tag] <- 'cdr3b.aa.seq'
colnames(tmp.data.1)[colnames(tmp.data.1)==cdr3a.aa.tag] <- 'cdr3a.aa.seq'
colnames(tmp.data.1)[colnames(tmp.data.1)==trb.v] <- 'trb.v'
colnames(tmp.data.1)[colnames(tmp.data.1)==trb.j] <- 'trb.j'
colnames(tmp.data.1)[colnames(tmp.data.1)==tra.v] <- 'tra.v'
colnames(tmp.data.1)[colnames(tmp.data.1)==tra.j] <- 'tra.j'
if(!is.null(score.tag)) colnames(tmp.data.1)[colnames(tmp.data.1)==score.tag] <- 'score.tag'
tmp.cols <- unique(c(
    ag.spc.tag, tags.of.int
))
tmp.cols <- setdiff(x=tmp.cols, y=colnames(tmp.data.1))
tmp.data.2 <- ref.meta[, ..tmp.cols]
match.inf.ref <- unique(cbind(tmp.data.1, tmp.data.2))
# Ag subset if required.
if(!is.null(ags.to.keep)) match.inf.ref <- match.inf.ref[get(ag.spc.tag) %in% ags.to.keep]
# Score column (if not provided, give the same score to all clonotypes).
if(!is.null(score.tag)){
    match.inf.ref[, score.tag:=as.numeric(score.tag)]
}else{
    match.inf.ref[, score.tag:=1]
}
# Separate consideration for multiple chain sequences defined within the same clonotype.
match.inf.ref <- sep.tcr.row.info(tcr.inf=match.inf.ref)

# ---> Define reference chain-specific clonotypes.
# Define chain-specific clonotypes for each chain.
match.inf.ref <- lapply(X=names(rma.chain.cols), FUN=function(chain.lab){
    chain.col <- rma.chain.cols[chain.lab]
    vdj.genes <- rma.chain.vdj[[chain.lab]]
    merge.cols <- c(chain.col, vdj.genes)
    names(merge.cols)[1] <- chain.col
    tmp.data <- match.inf.ref[!is.na(get(chain.col))]
    # Order according to clonotype score.
    setorderv(x=tmp.data, cols='score.tag', order=-1)
    # Set chain-defining columns as the key. That is, defined clonotypes on a chain-basis manner.
    tmp.cols <- setdiff(x=colnames(match.inf.ref), y=chain.col)
    tmp.data.1 <- lapply(X=tmp.cols, FUN=function(tmp.col){
        tmp.cols <- c(merge.cols, tmp.col)
        tmp.data <- tmp.data[, ..tmp.cols]
        colnames(tmp.data)[ncol(tmp.data)] <- 'tmp.col'
        tmp.data <- tmp.data[,
            .(
                tmp.col=paste0(tmp.col, collapse='|')
            ),
            by=merge.cols
        ]
        colnames(tmp.data) <- c(names(merge.cols), tmp.col)
        return(tmp.data)
    })
    tmp.data.1 <- Reduce(x=tmp.data.1, f=function(x, y){
        tmp.data <- merge(x=x, y=y, by=names(merge.cols))
        return(tmp.data)
    })
    # Set consensus antigen specificity for chain-specific clonotypes. This is done by considering the score of individual reference clonotypes. Only clonotypes w/ the greatest of the scores will be considered.
    tmp.cols <- c(merge.cols, 'ref.clone.id', ag.spc.tag, 'score.tag')
    tmp.data.2 <- tmp.data[, ..tmp.cols]
    tmp.data.2 <- tmp.data.2[,
        .SD[score.tag==max(score.tag)],
        by=merge.cols
    ]
    #       Bug-specific ammendment. To make sure column were named as expected. Note: Different data.table versions might yield different column names.
    tmp.check <- all(names(merge.cols) %in% colnames(tmp.data.2))
    if(!tmp.check){
        tmp.vals <- setdiff(x=colnames(tmp.data.2), merge.cols)
        to.fix <- c(
            names(merge.cols),
            tmp.vals
        )
        names(to.fix) <- c(
            merge.cols,
            tmp.vals
        )
        colnames(tmp.data.2) <- to.fix[colnames(tmp.data.2)]
    }
    #       Continue w/ process.
    tmp.cols <- setdiff(x=colnames(tmp.data.2), y=names(merge.cols))
    tmp.data.2 <- lapply(X=tmp.cols, FUN=function(tmp.col){
        tmp.cols <- c(names(merge.cols), tmp.col)
        tmp.data <- tmp.data.2[, ..tmp.cols]
        tmp.cols <- names(merge.cols)
        tmp.data <- tmp.data[,
            .(
                tmp.col=paste0(unique(get(tmp.col)), collapse='|')
            ),
            by=tmp.cols
        ]
        colnames(tmp.data) <- c(names(merge.cols), tmp.col)
        return(tmp.data)
    })
    tmp.data.2 <- Reduce(x=tmp.data.2, f=function(x, y){
        tmp.data <- merge(x=x, y=y, by=names(merge.cols))
        return(tmp.data)
    })
    tmp.vals.1 <- paste0('ref.consensus.', c('clone.id', 'ag', 'score'))
    names(tmp.vals.1) <- c('ref.clone.id', ag.spc.tag, 'score.tag')
    tmp.vals.2 <- names(merge.cols); names(tmp.vals.2) <- tmp.vals.2
    tmp.vals <- c(tmp.vals.1, tmp.vals.2)
    tmp.check <- all(colnames(tmp.data.2) %in% names(tmp.vals))
    if(!tmp.check) stop(paste0('Unexpected error during chain-specific complete matching for chain', chain.lab, '.\n'))
    colnames(tmp.data.2) <- tmp.vals[colnames(tmp.data.2)]
    # Chain-specific clonotypes w/ multiple specificities will be labelled "Multiple"
    tmp.data.2[
        str_detect(string=ref.consensus.ag, pattern='\\|'),
        ref.consensus.ag:='Multiple'
    ]
    # Merge reference's pieces together.
    tmp.check <- tmp.data.1[, uniqueN(get(chain.col))] == tmp.data.2[, uniqueN(get(chain.col))]
    if(!tmp.check) stop(paste0('Unexpected error during chain-specific complete matching for chain', chain.lab, '.\n'))
    tmp.data <- merge(
        x=tmp.data.1, y=tmp.data.2,
        by=names(merge.cols)
    )
    # Final format and return 
    tmp.data[['ref.consensus.score']] <- as.numeric(tmp.data[['ref.consensus.score']])
    return(tmp.data)
})
names(match.inf.ref) <- names(rma.chain.cols)


cat('\n\n')
### -------------------- Load clustering results -------------------- ###
### ---------------------- Load TCR UCM results --------------------- ###
cat('### ---------------------- Load TCR UCM results --------------------- ###\n')

# ---> Process per UCM.
cat('Loading and preprocessing results for next tools:\n')
ucm.res <- lapply(X=tool.names, FUN=function(tool.name){
    cat(tool.name, '\n')

    # ---> File definition.
    # Define output path.
    ucm.res.path <- paste0(run.path, '/', tool.name, '/RefThold-', expansion.thold, '')
    # Adjacency list.
    adj.list.file <- paste0(ucm.res.path,  '/ToolResults_AdjList-1.csv')
    # Attributes of clonotypes.
    clone.atts.file <- paste0(ucm.res.path,  '/ToolResults_ClonotypeAtts-1.csv')
    # Processed results from tool
    # sgs.info.file <- paste0(ucm.res.path,  '/ToolResults_Processed.csv')
    # Sanity checks
    to.check <- c(
        adj.list.file,
        clone.atts.file
        # sgs.info.file,
    )
    to.check <- to.check[!file.exists(to.check)]
    if(length(to.check)>0){
        to.check <- paste0('Next files do not exist:\n', paste0(to.check, collapse='\n'), '\n')
        stop(to.check)
    }

    # ---> Load results.
    # Adjacency list.
    adj.list <- fread(file=adj.list.file)
    # Attributes of clonotypes.
    clone.atts <- fread(file=clone.atts.file)
    # Processed results from GLIPH2
    # sgs.info <- fread(file=sgs.info.file)

    # ---> Tidy adjacency list 
    # Replace clonotype clustering clonotype ID by the source clonotype IDs.
    clone.atts <- separate_rows(data=clone.atts, clonotype.tags, sep=';')
    clone.atts <- as.data.table(separate(data=clone.atts, col=clonotype.tags, into=c('set.clone.id', 'data.set'), sep='\\|'))
    clone.atts <- clone.atts[!is.na(set.clone.id)] # RFI: This was necessary because there are some clonotypes that aren't properly defined for the reference dataset we used to develop this script. This must be fixed upstream and then this line should be removed.
    clone.atts[, clone.type:=ifelse(
        test=data.set==ref.name,
        yes='ref', no='qry'
    )]
    tmp.check <- setdiff(x=clone.atts[, unique(clone.type)], y=c('ref', 'qry'))
    if(length(tmp.check)>0) stop('Unexpected error. More types of dataset defined than anticipated. Should be only either query or reference.\n')
    # Merge w/ adjacency list for proper tidying.
    tmp.data <- clone.atts[, .(general.clone.id, set.clone.id, clone.type)]
    adj.list <- merge(x=adj.list, y=tmp.data, by.x='clone.id.1', by.y='general.clone.id', all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE)
    adj.list <- merge(x=adj.list, y=tmp.data, by.x='clone.id.2', by.y='general.clone.id', all.x=TRUE, all.y=FALSE, suffixes=c('.1', '.2'), allow.cartesian=TRUE)
    adj.list[, `:=`(clone.id.1=set.clone.id.1, clone.id.2=set.clone.id.2)]
    adj.list <- adj.list[, .(clone.id.1, clone.id.2, clone.type.1, clone.type.2, weight)]
    # Remove duplicates generated by "cartesian" option.
    adj.list <- unique(adj.list) 
    # Define interaction type.
    adj.list[clone.type.1=='qry' & clone.type.2=='qry', int.type:='UQ']
    adj.list[clone.type.1=='ref' & clone.type.2=='ref', int.type:='UR']
    adj.list[(clone.type.1=='ref' & clone.type.2=='qry') | (clone.type.1=='qry' & clone.type.2=='ref'), int.type:='NU']
    adj.list <- adj.list[!is.na(int.type)] # RFI: This was necessary because there are some clonotypes that aren't properly defined for the reference dataset we used to develop this script. This must be fixed upstream and then this line should be removed.

    # ---> A^{X} matrix creation
    # Identify and keep heterogeneous interactions.
    tmp.data <- adj.list[int.type=='NU']
    # Define query and reference-related clonotypes in that order.
    # For query clonotype.
    tmp.data[, qry.clone.id:='']
    tmp.data[clone.type.1=='qry', qry.clone.id:=clone.id.1]
    tmp.data[clone.type.2=='qry', qry.clone.id:=clone.id.2]
    tmp.check <- tmp.data[, all(qry.clone.id!='')] # Sanity check.
    if(!tmp.check) stop('Unexpected error 1.\n')
    # For reference clonotype.
    tmp.data[, ref.clone.id:='']
    tmp.data[clone.type.1=='ref', ref.clone.id:=clone.id.1]
    tmp.data[clone.type.2=='ref', ref.clone.id:=clone.id.2]
    tmp.check <- tmp.data[, all(ref.clone.id!='')] # Sanity check.
    if(!tmp.check) stop('Unexpected error 1.\n')
    # New table with relevant info (specifically, reactivity info from the reference dataset)
    tmp.data.1 <- tmp.data[, .(qry.clone.id, ref.clone.id)]
    tmp.data.2 <- match.inf.ref[['RMB']]
    # colnames(tmp.data.2)[colnames(tmp.data.2)=='v'] <- trb.v # NOTE: CAREFUL
    # colnames(tmp.data.2)[colnames(tmp.data.2)=='j'] <- trb.j
    colnames(tmp.data.2)[colnames(tmp.data.2)=='ref.clone.id'] <- 'ref.clone.ids'
    tmp.data.2 <- as.data.table(separate_rows(data=tmp.data.2, `ref.consensus.clone.id`, sep='\\|'))
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by.x='ref.clone.id', by.y='ref.consensus.clone.id', all.x=TRUE, all.y=FALSE, sort=FALSE)
    # Sanity check. After this process, we should recover an antigen specificity for at lest 50% of the reference clonotypes. Otherwise, ask the user to revisit their reference.
    tmp.check <- tmp.data[, .SD[!is.na(ref.consensus.ag), .N]/.N]>0.5
    if(is.na(tmp.check) | !tmp.check) warning(paste0('Something went wrong while defining specificities of clustered reference clonotypes for tool ', tool.name, '.\n'))
    tmp.data <- tmp.data[!is.na(ref.consensus.ag)]
    # ---> Output.
    return(tmp.data)
})
names(ucm.res) <- tool.names


cat('\n\n')
### ----------------- Anitgen specifciity prediction ---------------- ###
### ----------------------- TCR UCM approach ------------------------ ###
cat('### ----------------- Anitgen specifciity prediction ---------------- ###\n')
cat('### ----------------------- TCR UCM approach ------------------------ ###\n')

# ---> Function to define A^{X} matrix per antigen X.
get.a.x.mat <- function(ucm.res, ag.x){
    # Retrieve matrix. See definition somewhere else.
    tmp.data <- lapply(X=names(ucm.res), FUN=function(tool.name){
        tmp.data <- ucm.res[[tool.name]]
        tmp.data[['ag']] <- tmp.data[['ref.consensus.ag']]; tmp.data[['ref.consensus.ag']] <- NULL
        tmp.data[['score']] <- tmp.data[['ref.consensus.score']]; tmp.data[['ref.consensus.score']] <- NULL
        tmp.data <- unique(tmp.data[ag==ag.x])
        # Weight calculation.
        tmp.data.1 <- tmp.data[,
            .(
                weight=sum(score)
            ), # NOTE: We went back to unique reference entries after calling chain-specific reference clonotypes.
            by=qry.clone.id
        ]
        colnames(tmp.data.1)[colnames(tmp.data.1)=='weight'] <- tool.name
        # Carry tags of interest
        tmp.data.2 <- lapply(X=tags.of.int, FUN=function(tmp.tag){
            tmp.data <- tmp.data[,
                .(
                    tmp.col=paste0(
                        get(tmp.tag), collapse='/'
                    )
                ),
                by=qry.clone.id
            ]
            colnames(tmp.data)[2] <- tmp.tag
            return(tmp.data)
        })
        tmp.data.2 <- Reduce(x=tmp.data.2, f=function(x, y){
            tmp.data <- merge(x=x, y=y, by='qry.clone.id')
            return(tmp.data)
        })
        colnames(tmp.data.2)[colnames(tmp.data.2)!='qry.clone.id'] <- paste(
            tool.name,
            colnames(tmp.data.2)[colnames(tmp.data.2)!='qry.clone.id'],
            sep='|'
        )
        # Merge and return
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='qry.clone.id')
        return(tmp.data)
    })
    tmp.data.1 <- Reduce(x=tmp.data, f=function(x, y){
        tmp.data <- merge(x=x, y=y, by='qry.clone.id', all=TRUE)
        return(tmp.data)
    })
    # Prepare A^{X} matrix as a separate object.
    tmp.data.2 <- as.matrix(tmp.data.1[, ..tool.names])
    row.names(tmp.data.2) <- tmp.data.1[, qry.clone.id]
    tmp.data.2[is.na(tmp.data.2)] <- 0
    # Return
    to.return <- list(
        `a.x.mat`=tmp.data.2,
        `score.meta`=tmp.data.1
    )
    return(to.return)
}

# ---> Function to define B^{X} matrix per antigen X.
get.b.x.mat <- function(ucm.res, ag.x, w=w.vector){
    # Retrieve A^{X} and score metadata. Separate accordingly.
    score.meta <- get.a.x.mat(ucm.res=ucm.res, ag.x=ag.x)
    a.x.mat <- score.meta[['a.x.mat']]
    score.meta <- score.meta[['score.meta']]
    # Weigth observations according to tool.
    for(tool.name in names(ucm.res)){
        a.x.mat[, tool.name] <- a.x.mat[, tool.name] * w[tool.name]
    }
    # Return.
    to.return <- list(
        `b.x.mat`=a.x.mat,
        `score.meta`=score.meta
    )
    return(to.return)
}

# ---> Function to define C^{X} vector per antigen X.
get.c.x.vec <- function(ucm.res, ag.x){
    # Retrieve B^{X} and score metadata. Separate accordingly.
    score.meta <- get.b.x.mat(ucm.res=ucm.res, ag.x=ag.x)
    b.x.mat <- score.meta[['b.x.mat']]
    score.meta <- score.meta[['score.meta']]
    # Condense general support for antigen X per clonotype.
    c.x.vec <- rowSums(b.x.mat)
    c.x.vec <- data.table(
        qry.clone.id=names(c.x.vec),
        ag.x=c.x.vec
    )
    # Return.
    if(nrow(c.x.vec)>0){
        to.return <- list(
            `c.x.vec`=c.x.vec,
            `score.meta`=score.meta
        )
    }else{
        to.return <- NA
    }
    return(to.return)
}

# ---> Define D matrix
# Specify antigens to give predictions for.
tmp.ags <- sort(unique(ref.meta[, get(ag.spc.tag)]))
if(!is.null(ags.to.keep)) tmp.ags <- tmp.ags[tmp.ags %in% ags.to.keep] # Custom exclusion accorsing to reference.

# Retrieve C^{X} and score metadata. Separate accordingly.
tmp.data <- lapply(X=tmp.ags, FUN=function(ag.x){
    tmp.data <- get.c.x.vec(ucm.res=ucm.res, ag.x=ag.x)
    return(tmp.data)
})
names(tmp.data) <- tmp.ags
tmp.data <- tmp.data[!is.na(tmp.data)] # Remove entries for ags with empty info.
# Update ag list accordingly. Checkpoint.
tmp.ags <- names(tmp.data)
tmp.check <- length(tmp.ags)>0
if(!tmp.check) stop('Something went wrong while calculating D matrix. No scores were retrieved.\n')
# Separate results accordingly.
c.x.vecs <- lapply(X=tmp.data, FUN=function(x) return(x$c.x.vec))
score.meta <- lapply(X=tmp.data, FUN=function(x) return(x$score.meta))
# Merge C^{X} vectors into single D matrix.
tmp.data <- rbindlist(l=c.x.vecs, use.names=TRUE, idcol='ag.name')
tmp.data <- as.data.table(spread(data=tmp.data, key=ag.name, value=ag.x, fill=0))
d.matrix <- as.matrix(tmp.data[, ..tmp.ags])
row.names(d.matrix) <- tmp.data[, qry.clone.id]
# Combined score metadata.
score.meta <- rbindlist(l=score.meta, use.names=TRUE, idcol='ag.name')
colnames(score.meta)[3:ncol(score.meta)] <- paste0(
    'UCM|M|',
    colnames(score.meta)[3:ncol(score.meta)]
)
# Set an upper bound in character count for metadata entries.
meta.entry.thold <- 3000
tmp.cols <- setdiff(x=colnames(score.meta), y=c('ag.name', 'qry.clone.id'))
for(tmp.col in tmp.cols){
    tmp.flag <- score.meta[, !is.na(get(tmp.col)) & str_count(string=get(tmp.col))>meta.entry.thold]
    tmp.vals <- score.meta[, get(tmp.col)]
    tmp.vals[tmp.flag] <- paste0(
        str_extract(
            string=tmp.vals[tmp.flag],
            pattern=paste0('.{', meta.entry.thold, '}')
        ),
        '[BOUNDED]'
    )
    set(x=score.meta, j=tmp.col, value=tmp.vals)
}

# ---> Ratio distribution exploration.
# Define ratio
ratio.num <- apply(X=d.matrix, MARGIN=1, FUN=max)
ratio.den <- apply(X=d.matrix, MARGIN=1, FUN=function(x){
    return(max(setdiff(x=x, y=max(x))))
})
ratio.vals <- ratio.num / ratio.den
ratio.vals[is.infinite(ratio.vals)] <- ratio.num[is.infinite(ratio.vals)]
# Define antigens w/ max value.
max.ags <- apply(X=d.matrix, MARGIN=1, FUN=function(x){
    tmp.vals <- colnames(d.matrix)[x==max(x)]
    tmp.vals <- paste0(sort(tmp.vals), collapse=';')
    return(tmp.vals)
})
# Summarize.
pred.summ <- data.table(
    qry.clone.id=row.names(d.matrix),
    UCM.ratio=ratio.vals,
    UCM.max.ags=max.ags
)
# Define general prediction outcome.
pred.summ[, 
    UCM.ag.class:=ifelse(
        test=str_detect(string=UCM.max.ags, pattern=';'),
        yes='Multiple',
        no='Single'
    )
]
# Explore ratio.
tmp.ggplot.1 <- ggplot(data=pred.summ, aes(x=UCM.ratio)) +
    geom_density(alpha=0, color='black', linewidth=2) +
    geom_vline(xintercept=ucm.thold, linewidth=2, color='red', linetype='dashed') +
    scale_y_continuous(expand=c(0, 0)) +
    labs(title='UCM ratio', x='Ratio', y='Density') +
    theme_bw()
tmp.ggplot.2 <- ggplot(data=pred.summ, aes(x=UCM.ratio)) +
    geom_density(alpha=0, color='black', linewidth=2) +
    geom_vline(xintercept=ucm.thold, linewidth=2, color='red', linetype='dashed') +
    scale_y_continuous(expand=c(0, 0)) +
    scale_x_log10(expand=c(0, 0)) +
    labs(title='UCM ratio', x='Ratio (log10)', y='Density') +
    theme_bw()
tmp.ggplot.3 <- ggplot(data=pred.summ, aes(x=UCM.ratio, color=UCM.ag.class)) +
    geom_density(alpha=0, linewidth=2) +
    geom_vline(xintercept=ucm.thold, linewidth=2, color='red', linetype='dashed') +
    scale_y_continuous(expand=c(0, 0)) +
    scale_color_viridis_d(option='viridis') +
    labs(title='UCM ratio according to general potential prediction antigen class', x='Ratio', y='Density', color='General\nantigen\nclass') +
    theme_bw()
tmp.ggplot.4 <- ggplot(data=pred.summ, aes(x=UCM.ratio, color=UCM.ag.class)) +
    geom_density(alpha=0, linewidth=2) +
    geom_vline(xintercept=ucm.thold, linewidth=2, color='red', linetype='dashed') +
    scale_y_continuous(expand=c(0, 0)) +
    scale_x_log10(expand=c(0, 0)) +
    scale_color_viridis_d(option='viridis') +
    labs(title='UCM ratio according to general potential prediction antigen class', x='Ratio', y='Density', color='General\nantigen\nclass') +
    theme_bw()
tmp.file.name <- paste0(reports.path, '/RatioEmpDist.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
dev.off()

# ---> Final clustering-related prediction and gather prediction-associated details.
# Final prediction accordingly
pred.summ[
    UCM.ratio>=ucm.thold & UCM.ag.class!='Multiple',
    UCM.pred:=UCM.max.ags
]
# General prediction class.
pred.summ[,
    UCM.clust.class:=ifelse(
        test=is.na(UCM.pred),
        yes='Clustered',
        no=UCM.ag.class
    )
]
# Include donor ID
tmp.data <- unique(tcr.meta[, .(qry.clone.id, donor.id.tag)])
pred.summ <- merge(
    x=pred.summ, y=tmp.data,
    by='qry.clone.id',
    all.x=TRUE, all.y=FALSE
)
# Merge w/ D matrix.
tmp.data <- as.data.table(d.matrix)
colnames(tmp.data) <- paste0('score.', colnames(tmp.data))
tmp.data[, qry.clone.id:=row.names(d.matrix)]
pred.summ <- merge(x=pred.summ, y=tmp.data, by='qry.clone.id', all.x=TRUE, all.y=FALSE)
# Bring relevant columns to the front.
tmp.cols <- c('UCM.pred', 'UCM.max.ags', 'UCM.ag.class', 'UCM.clust.class', 'UCM.ratio')
tmp.cols <- c(tmp.cols, setdiff(x=colnames(pred.summ), y=tmp.cols))
pred.summ <- pred.summ[, ..tmp.cols]
# Merge w/ score metadata for sucessful final "single" predictions.
pred.summ <- merge(
    x=pred.summ, y=score.meta,
    by.x=c('qry.clone.id', 'UCM.pred'),
    by.y=c('qry.clone.id', 'ag.name'),
    all.x=TRUE, all.y=FALSE
)
# bkup <- copy(pred.summ)


cat('\n\n')
### ----------------- Anitgen specifciity prediction ---------------- ###
### -------------------- Reference match approach ------------------- ###
cat('### ----------------- Anitgen specifciity prediction ---------------- ###\n')
cat('### -------------------- Reference match approach ------------------- ###\n')

# ---> Retrieve inquiry TCR data.
# Retrieve basic inquiry metadata.
tmp.cols <- c('donor.id.tag', 'qry.clone.id', 'cdr3b.aa.seq', 'cdr3a.aa.seq', 'trb.v', 'trb.j', 'tra.v', 'tra.j')
match.inf.qry <- unique(tcr.meta[, ..tmp.cols])
# Separate consideration for multiple chain sequences defined within the same clonotype.
match.inf.qry <- sep.tcr.row.info(tcr.inf=match.inf.qry)

# ---> Function to set inquiry clonotype-specific info after chain-specific matches
# To be used both for complete and incomplete matches.
set.qry.key <- function(
    match.info,
    chain.lab, chain.col
){
    # Set combination of query clonotype ID and donor ID as key.
    tmp.cols <- setdiff(x=colnames(match.info), y=c('donor.id.tag', 'qry.clone.id'))
    tmp.data.1 <- lapply(X=tmp.cols, FUN=function(tmp.col){
        tmp.cols <- c('donor.id.tag', 'qry.clone.id', tmp.col)
        tmp.data <- match.info[, ..tmp.cols]
        tmp.data <- tmp.data[,
            .(tmp.col=paste0(get(tmp.col), collapse='/')),
            by=.(donor.id.tag, qry.clone.id)
        ]
        colnames(tmp.data) <- c('donor.id.tag', 'qry.clone.id', tmp.col)
        return(tmp.data)
    })
    tmp.data.1 <- Reduce(x=tmp.data.1, f=function(x, y){
        tmp.data <- merge(x=x, y=y, by=c('donor.id.tag', 'qry.clone.id'))
        return(tmp.data)
    })
    # Set consensus antigen specificity for inquiry clonotypes. This is done by considering the priority of individual reference chain-specific clonotypes. Only clonotypes w/ the greatest of the priorities will be considered.
    tmp.cols <- c(
        'qry.clone.id', 'donor.id.tag',
        paste0(chain.lab, '|ref.consensus.ag'),
        paste0(chain.lab, '|ref.consensus.clone.id'),
        paste0(chain.lab, '|ref.consensus.score')
    )
    names(tmp.cols) <- c(
        'qry.clone.id', 'donor.id.tag',
        'ref.consensus.ag',
        'ref.consensus.clone.id',
        'ref.consensus.score'
    )
    tmp.data.2 <- match.info[, ..tmp.cols]
    colnames(tmp.data.2) <- names(tmp.cols)
    tmp.data.2[['ref.consensus.score']] <- tmp.data.2[['ref.consensus.score']]
    # Reference clonotypes w/ an ambiguous specificity will be discarded.
    tmp.data.2 <- tmp.data.2[ref.consensus.ag!='Multiple']
    tmp.data.2 <- tmp.data.2[,
        .(
            consensus.ag=.SD[
                ref.consensus.score==max(ref.consensus.score),
                ifelse(
                    test=uniqueN(ref.consensus.ag)>1,
                    yes='Multiple',
                    no=unique(ref.consensus.ag)
                )
            ],
            consensus.clone.id=.SD[
                ref.consensus.score==max(ref.consensus.score),
                ifelse(
                    test=uniqueN(ref.consensus.ag)>1,
                    yes='Multiple',
                    no=paste0(unique(ref.consensus.clone.id), collapse='/')
                )
            ]
        ),
        by=.(qry.clone.id, donor.id.tag)
    ]
    tmp.data.2[
        consensus.ag=='Multiple',
        `:=`(
            consensus.ag=NA,
            consensus.clone.id=NA
        )
    ]
    colnames(tmp.data.2)[colnames(tmp.data.2)=='consensus.ag'] <- paste0(chain.lab, '|consensus.ag')
    colnames(tmp.data.2)[colnames(tmp.data.2)=='consensus.clone.id'] <- paste0(chain.lab, '|consensus.clone.id')
    # Sanity checking as we merge.
    tmp.data <- merge(
        x=tmp.data.1, y=tmp.data.2,
        by=c('qry.clone.id', 'donor.id.tag'),
        all=TRUE
    )
    tmp.check.1 <- match.info[, .N, by=.(qry.clone.id, donor.id.tag)][, .N]
    tmp.check.2 <- tmp.data[, .N, by=.(qry.clone.id, donor.id.tag)][, .N]
    tmp.check <- tmp.check.1==tmp.check.2
    if(!tmp.check) stop(paste0('Unexpected error during inquiry clonotype-specific complete matching for chain', chain.lab, '.\n'))
    # Final sanity check. There should be the same amount of chain-specific unique reference clonotypes as there are entries for the info regarding their VDJ gene match and scoring.
    tmp.check <- tmp.data[,
        .(
            to.check.1=str_count(
                string=get(paste(chain.lab, chain.col, sep='|')),
                pattern='/'
            ),
            to.check.2=str_count(
                string=get(paste(chain.lab, 'vj.match', sep='|')),
                pattern='/'
            )
        )
    ]
    tmp.check <- tmp.check[,
        all(to.check.1==to.check.2)
    ]
    if(!tmp.check) stop(paste0('VDJ gene match assessment failed. Different number of info entries compared w/ the number of chain-specific clonotypes in dataset for chain ', chain.lab, '.\n'))
    # Final return 
    return(tmp.data)
}

# ---> Chain-specific complete matches.
# Merge on a chain basis.
comp.match.info <- lapply(X=names(rma.chain.cols), FUN=function(chain.lab){
    # ---> Merge reference and inquiry dataset.
    # Retrieve basic info and reference chain-specific clonotypes.
    chain.col <- rma.chain.cols[chain.lab]
    vdj.genes <- rma.chain.vdj[[chain.lab]]
    merge.cols <- c(chain.col, vdj.genes)
    names(merge.cols)[1] <- chain.col
    tmp.data.1 <- match.inf.ref[[chain.lab]]
    # Unique chain-specific labels.
    colnames(tmp.data.1) <- paste(chain.lab, colnames(tmp.data.1), sep='|')
    tmp.data.1[['v']] <- tmp.data.1[[paste0(chain.lab, '|v')]]; tmp.data.1[[paste0(chain.lab, '|v')]] <- NULL
    tmp.data.1[['j']] <- tmp.data.1[[paste0(chain.lab, '|j')]]; tmp.data.1[[paste0(chain.lab, '|j')]] <- NULL
    # Merge w/ query info.
    tmp.cols <- c('qry.clone.id', chain.col, vdj.genes, 'donor.id.tag')
    names(tmp.cols) <- c('qry.clone.id', chain.col, names(vdj.genes), 'donor.id.tag')
    tmp.data.2 <- match.inf.qry[, ..tmp.cols]
    colnames(tmp.data.2) <- names(tmp.cols)
    tmp.data <- merge(
        x=tmp.data.1, y=tmp.data.2,
        by.x=paste(chain.lab, chain.col, sep='|'),
        by.y=chain.col,
        all=FALSE
    )
    # Add vdj gene match info and discard individual genes.
    tmp.data.1 <- tmp.data[,
        .(
            v.match=v.x==v.y,
            j.match=j.x==j.y,
            vj.match=(v.x==v.y & j.x==j.y)
        )
    ]
    colnames(tmp.data.1) <- paste(chain.lab, colnames(tmp.data.1), sep='|')
    tmp.data <- cbind(tmp.data, tmp.data.1)
    cols.to.rm <- paste(c('v', 'j'), rep(x=c('x', 'y'), each=2), sep='.')
    if(!all(cols.to.rm %in% colnames(tmp.data))) stop(paste0('Unexpected error during chain-specific complete matching for chain', chain.lab, '.\n'))
    tmp.cols <- setdiff(x=colnames(tmp.data), y=cols.to.rm)
    tmp.data <- tmp.data[, ..tmp.cols]
    # Maintain only matches w/ a proper VDJ gene match (if required).
    # if(vdj.match) tmp.data <- tmp.data[get(paste(chain.lab, 'vj.match', sep='|'))==TRUE]
    if(vdj.match) tmp.data <- tmp.data[get(paste(chain.lab, 'v.match', sep='|'))==TRUE]
    # ---> Inquiry clonotype-specific info. 
    tmp.data <- set.qry.key(
        match.info=tmp.data,
        chain.lab=chain.lab, chain.col=chain.col
    )
    # Final return 
    return(tmp.data)
})
# Define general type of CDR3 sequence match.
comp.match.info <- merge(x=comp.match.info[[1]], y=comp.match.info[[2]], by=c('qry.clone.id', 'donor.id.tag'), all=TRUE)
comp.match.info[,
    RM.chain.class:=ifelse(
        test=!is.na(`RMB|consensus.ag`) & !is.na(`RMA|consensus.ag`),
        yes='Two-chain perfect match',
        no=ifelse(
            test=!is.na(`RMB|consensus.ag`),
            yes='Beta-only perfect match',
            no=ifelse(
                test=!is.na(`RMA|consensus.ag`),
                yes='Alpha-only perfect match',
                no=NA
            )
        )
    )
]

# ---> Incomplete matches for complete single-chain matches.
# Merge on a chain basis.
icpt.match.info <- lapply(X=names(rma.chain.cols), FUN=function(chain.lab){
    # ---> Merge reference and inquiry dataset.
    # Retrieve basic info and reference chain-specific clonotypes.
    chain.col <- rma.chain.cols[chain.lab]
    vdj.genes <- rma.chain.vdj[[chain.lab]]
    merge.cols <- c(chain.col, vdj.genes)
    names(merge.cols)[1] <- chain.col
    tmp.data.1 <- match.inf.ref[[chain.lab]]
    colnames(tmp.data.1) <- paste(chain.lab, colnames(tmp.data.1), sep='|')
    tmp.data.1[['v']] <- tmp.data.1[[paste0(chain.lab, '|v')]]; tmp.data.1[[paste0(chain.lab, '|v')]] <- NULL
    tmp.data.1[['j']] <- tmp.data.1[[paste0(chain.lab, '|j')]]; tmp.data.1[[paste0(chain.lab, '|j')]] <- NULL

    # Identify query clonotypes w/ single-chain complete matches for the alternative chain under consideration.
    match.type <- if(chain.lab=='RMA') 'Beta-only perfect match' else 'Alpha-only perfect match'
    tmp.clones <- comp.match.info[
        RM.chain.class==match.type,
        qry.clone.id, 
    ]
    tmp.cols <- c(
        'qry.clone.id', chain.col,
        vdj.genes, 'donor.id.tag'
    )
    names(tmp.cols) <- c(
        'qry.clone.id', paste(chain.lab, chain.col, sep='|'),
        names(vdj.genes), 'donor.id.tag'
    )
    tmp.data.2 <- match.inf.qry[
        qry.clone.id %in% tmp.clones & !is.na(get(chain.col)), 
        ..tmp.cols
    ]
    colnames(tmp.data.2) <- names(tmp.cols)
    # Merge w/ query info only for clonotypes w/ incomplete single-chain matches when referring to the opposite chain under evaluation.
    tmp.data <- as.data.table(stringdist_inner_join(
        x=tmp.data.1, y=tmp.data.2,
        by=paste(chain.lab, chain.col, sep='|'),
        max_dist=1, method='hamming'
    ))
    # Keep chains equence columns only for the reference CDR3 sequence, since the query will still be defined in final output.
    tmp.cols <- setdiff(x=colnames(tmp.data), y=paste0(chain.lab, '|', chain.col, '.y'))
    tmp.data <- tmp.data[, ..tmp.cols]
    colnames(tmp.data)[colnames(tmp.data)==paste0(chain.lab, '|', chain.col, '.x')] <- paste0(chain.lab, '|', chain.col)
    # Add vdj gene match info and discard individual genes.
    tmp.data.1 <- tmp.data[,
        .(
            v.match=v.x==v.y,
            j.match=j.x==j.y,
            vj.match=(v.x==v.y & j.x==j.y)
        )
    ]
    colnames(tmp.data.1) <- paste(chain.lab, colnames(tmp.data.1), sep='|')
    tmp.data <- cbind(tmp.data, tmp.data.1)
    cols.to.rm <- paste(c('v', 'j'), rep(x=c('x', 'y'), each=2), sep='.')
    if(!all(cols.to.rm %in% colnames(tmp.data))) stop(paste0('Unexpected error during chain-specific complete matching for chain', chain.lab, '.\n'))
    tmp.cols <- setdiff(x=colnames(tmp.data), y=cols.to.rm)
    tmp.data <- tmp.data[, ..tmp.cols]
    # Maintain only matches w/ a proper VDJ gene match (if required).
    if(vdj.match) tmp.data <- tmp.data[get(paste(chain.lab, 'vj.match', sep='|'))==TRUE]    
    # ---> Inquiry clonotype-specific info. 
    # Merge between inquiry and reference datasets.
    tmp.data <- set.qry.key(
        match.info=tmp.data,
        chain.lab=chain.lab, chain.col=chain.col
    )
    # Merge info w/ that corresponding to the chain that matched completely.
    tmp.clones <- tmp.data[, qry.clone.id]
    tmp.data.1 <- comp.match.info[qry.clone.id %in% tmp.clones]
    tmp.cols <- c('qry.clone.id', 'donor.id.tag', setdiff(x=colnames(tmp.data.1), y=colnames(tmp.data)))
    tmp.data.1 <- tmp.data.1[, ..tmp.cols]
    tmp.data <- merge(
        x=tmp.data.1,
        y=tmp.data,
        by=c('qry.clone.id', 'donor.id.tag')
    )
    tmp.data[, RM.chain.class:=NULL]
    tmp.data[, RM.chain.class:=if(chain.lab=='RMA') 'Perfect beta match with alpha quasi match' else 'Perfect alpha match with beta quasi match']
    # Final return 
    return(tmp.data)
})
icpt.match.info <- rbindlist(l=icpt.match.info, use.names=TRUE)

# ---> Full match info.
# Join info from complete and incomplete matches while avoiding redundancy in the final info.
ref.match.info <- list(
    comp.match.info[!qry.clone.id %in% icpt.match.info[, qry.clone.id]],
    icpt.match.info
)
ref.match.info <- rbindlist(l=ref.match.info, use.names=TRUE)
tmp.check <- ref.match.info[, uniqueN(qry.clone.id)]==comp.match.info[, uniqueN(qry.clone.id)]
if(!tmp.check) stop('Something went wrong while merging complete and incomplete match info (1).\n')

# ---> General match class.
# Required info.
tmp.data.1 <- str_split(
    string=ref.match.info[, `RMB|consensus.clone.id`],
    pattern='\\||/'
)
tmp.data.2 <- str_split(
    string=ref.match.info[, `RMA|consensus.clone.id`],
    pattern='\\||/'
)
# tmp.check <- length(tmp.data.1)==length(tmp.data.2)
tmp.data <- unlist(lapply(X=1:length(tmp.data.1), FUN=function(i){
    tmp.data <- intersect(x=tmp.data.1[[i]], y=tmp.data.2[[i]])
    tmp.data <- tmp.data[!is.na(tmp.data)]
    tmp.data <- if(length(tmp.data)>0) paste0(tmp.data, collapse='|') else NA
    return(tmp.data)
}))
ref.match.info[['RMG|consensus.clone.id']] <- tmp.data
#       Perfect matches.
tmp.classes <- c('Two-chain perfect match', 'Perfect beta match with alpha quasi match', 'Perfect alpha match with beta quasi match')
tmp.vals <- ifelse(
    test=ref.match.info[, RM.chain.class%in%tmp.classes & !is.na(`RMG|consensus.clone.id`)],
    yes='Perfect', no=NA
)
ref.match.info[, RM.match.class:=tmp.vals]
# set(x=ref.match.info, j='RM.match.class', value=tmp.vals)
#       Imperfect matches.
#           Unambiguous in terms of specificity.
ref.match.info[
    RM.chain.class%in%tmp.classes & is.na(`RMG|consensus.clone.id`) & `RMA|consensus.ag`==`RMB|consensus.ag`,
    RM.match.class:='Unambiguous'
]
#           Ambiguous in terms of specificity.
ref.match.info[
    RM.chain.class%in%tmp.classes & is.na(RM.match.class),
    RM.match.class:='Ambiguous'
]
#       Single-chain matches.
ref.match.info[
    (RM.chain.class=='Alpha-only perfect match' | RM.chain.class=='Beta-only perfect match'),
    RM.match.class:='Unambiguous'
]

# ---> Reference-matched consensus antigen
# Set consensus antigen. Priority is given to the beta chain when available (i.e., for single-chain beta and double-chain matches).
ref.match.info[,
    RM.pred:=`RMB|consensus.ag`
]
ref.match.info[
    is.na(RM.pred),
    RM.pred:=`RMA|consensus.ag`
]

# ---> Final format.
# Bring relevant columns to the front.
rmb.rel.cols <- c('RMB|consensus.ag')
rma.rel.cols <- c('RMA|consensus.ag')
tmp.cols <- c(
    'qry.clone.id', 'donor.id.tag',
    'RM.pred', 'RM.chain.class', 'RM.match.class', 'RMG|consensus.clone.id', rmb.rel.cols, rma.rel.cols,
    setdiff(x=colnames(ref.match.info)[str_detect(string=colnames(ref.match.info), pattern='^RMB\\|')], y=rmb.rel.cols),
    setdiff(x=colnames(ref.match.info)[str_detect(string=colnames(ref.match.info), pattern='^RMA\\|')], y=rma.rel.cols)
)
tmp.check <- length(setdiff(x=colnames(ref.match.info), y=tmp.cols))==0 & length(setdiff(x=tmp.cols, y=colnames(ref.match.info)))==0
if(!tmp.check) stop('Unexpected error. Faulty column set for reference match approach.\n')
ref.match.info <- ref.match.info[, ..tmp.cols]


cat('\n\n')
### ----------------- Anitgen specifciity prediction ---------------- ###
### --------------------- Donor-matched approach -------------------- ###
cat('### ----------------- Anitgen specifciity prediction ---------------- ###\n')
cat('### --------------------- Donor-matched approach -------------------- ###\n')

# ---> To perform only when this type of data is available and a donor ID has been provided.
if(!is.null(exp.res) & !is.null(qry.donor.tag)){
    # ---> Clean experimental validation data.
    # Retrieve results.
    match.inf.exp <- exp.res[
      cell.subset==t.subset & !is.na(specificity.class.tag) & specificity.class.tag %in% ags.to.keep,
      .(
        total.abs.freq=sum(total.abs.freq),
        freq.conf.lvl=max(freq.conf.lvl)
      ),
      by=.(
        tcr.chain, donor.id.tag=donor.id,
        cdr3.aa.seq, trv, trj,
        DM.pred=specificity.class.tag
      )
    ]
    # @ Unique entries per clonotype.
    # That is, define unique clonotypes based on their aminoacid sequence.
    # ---
    # NOTE: There might not necessarily be unique entries here because we defined the clonotypes based on the nucleotide sequence and there might exist cases of multiple different nucleotide sequences within the same donor giving rise to the same aminoacid sequence. For these cases, we will keep the nucleotide sequence with the greates UMI support within the donor.
    # ---
    # If multiple specificities were assigned to the same aminoacid sequence, we keep the one(s) w/ the greatest UMI support.
    clone.def.cols <- c('tcr.chain', 'donor.id.tag', 'cdr3.aa.seq', 'trv', 'trj')
    match.inf.exp <- match.inf.exp[,
        .(
            DM.pred=.SD[total.abs.freq==max(total.abs.freq), DM.pred],
            freq.conf.lvl=.SD[total.abs.freq==max(total.abs.freq), freq.conf.lvl],
            total.abs.freq=max(total.abs.freq)
        ),
      by=clone.def.cols
    ]
    # Identify and remove clonotypes that still remain repeated. These repeats won't be considered a problem unless they represent over 0.5% of the dataset.
    to.discard <- match.inf.exp[, .N, by=clone.def.cols][N>1]
    tmp.check.1 <- to.discard[, .N]/match.inf.exp[, .N]
    tmp.check.2 <- tmp.check.1<0.005
    tmp.check.1 <- round(x=tmp.check.1*100, digits=3)
    if(!tmp.check.2) stop(paste0('Problem found with the exp. validation data for reference ', ref.name, '.\n', tmp.check.1, '% of the clonotypes represent aminoacid sequences that were given rise by multiple unequal nucleotide sequences and that were ultimately assigned to different antigen specificities. When this fraction is less than 0.5%, they are just disregarded, but this case was flagged problematic.\nPlease inspect your dataset and, if found necessary, request/perform changes in/to the workflow.\n'))
    to.discard <- merge(x=match.inf.exp, y=to.discard, by=clone.def.cols, all=FALSE)
    to.discard[, discard:=TRUE]
    tmp.cols <- c(clone.def.cols, 'discard'); to.discard <- to.discard[, ..tmp.cols]
    match.inf.exp <- merge(x=match.inf.exp, y=to.discard, by=clone.def.cols, all=TRUE, allow.cartesian=TRUE)
    match.inf.exp <- match.inf.exp[is.na(discard)]; match.inf.exp[, discard:=NULL]
    # Final sanity check.
    tmp.check <- match.inf.exp[, .N, by=clone.def.cols][N>1, .N==0]
    if(!tmp.check) stop('Failed to define unique clonotypes based on the aminoacid sequence for data from the prolif. assay experiments.\n')
    # For clonotypes w/ multiple V and J gene instances, separate them into multiple rows.
    tmp.cols <- c('trv', 'trj')
    for(tmp.col in tmp.cols){
        match.inf.exp <- as.data.table(separate_rows(data=match.inf.exp, tmp.col, sep=';'))
    }
    # ---> Retrieve query TCR metadata
    tmp.cols <- c('donor.id.tag', 'qry.clone.id', 'cdr3b.aa.seq', 'cdr3a.aa.seq')
    match.inf.qry <- unique(tcr.meta[, ..tmp.cols])
    match.inf.qry <- as.data.table(separate_rows(data=match.inf.qry, cdr3b.aa.seq, sep=';'))
    match.inf.qry <- as.data.table(separate_rows(data=match.inf.qry, cdr3a.aa.seq, sep=';'))
    # ---> Merge on a chain basis.
    exp.match.info <- lapply(X=names(dma.chain.cols), FUN=function(chain.lab){
        # Retrieve necessary data for the process.
        chain.col <- dma.chain.cols[chain.lab]
        chain.val <- if(chain.col=='cdr3b.aa.seq') 'TCRB' else 'TCRA'
        tmp.data.1 <- match.inf.exp[tcr.chain==chain.val]; tmp.data.1[, tcr.chain:=NULL]
        colnames(tmp.data.1)[colnames(tmp.data.1)!='donor.id.tag'] <- paste(
            chain.lab,
            colnames(tmp.data.1)[colnames(tmp.data.1)!='donor.id.tag'],
            sep='|'
        )
        tmp.data.2 <- match.inf.qry[!is.na(get(chain.col)), .(qry.clone.id, get(chain.col), donor.id.tag)]
        colnames(tmp.data.2) <- c('qry.clone.id', chain.col, 'donor.id.tag')
        # Merge data.
        tmp.data <- merge(
            x=tmp.data.1, y=tmp.data.2,
            by.x=c(paste(chain.lab, 'cdr3.aa.seq', sep='|'), 'donor.id.tag'),
            by.y=c(chain.col, 'donor.id.tag'),
            all=FALSE
        )
        # Set combination of query clonotype ID and donor ID as the key.
        tmp.cols <- setdiff(x=colnames(tmp.data), y=c('donor.id.tag', 'qry.clone.id'))
        tmp.data <- lapply(X=tmp.cols, FUN=function(tmp.col){
            tmp.cols <- c('donor.id.tag', 'qry.clone.id', tmp.col)
            tmp.data <- tmp.data[, ..tmp.cols]
            tmp.data <- tmp.data[,
                .(tmp.col=paste0(get(tmp.col), collapse='|')),
                by=.(donor.id.tag, qry.clone.id)
            ]
            colnames(tmp.data) <- c('donor.id.tag', 'qry.clone.id', tmp.col)
            return(tmp.data)
        })
        tmp.data <- Reduce(x=tmp.data, f=function(x, y){
            tmp.data <- merge(x=x, y=y, by=c('donor.id.tag', 'qry.clone.id'))
            return(tmp.data)
        })
        # Set consensus specificity according to chain.
        tmp.col <- paste(chain.lab, 'DM.pred', sep='|')
        tmp.col <- str_split(
            string=tmp.data[, get(tmp.col)], pattern='\\|'
        )
        tmp.col <- unlist(lapply(X=tmp.col, FUN=function(x) paste0(sort(unique(x)), collapse='|')))
        tmp.data[[paste(chain.lab, 'consensus.ag', sep='|')]] <- tmp.col
        return(tmp.data)
    })
    # ---> Consensus experimental annotations.
    # Define general type of CDR3 sequence match.
    exp.match.info <- merge(x=exp.match.info[[1]], y=exp.match.info[[2]], by=c('qry.clone.id', 'donor.id.tag'), all=TRUE)
    exp.match.info[,
        DM.chain.class:=ifelse(
            test=!is.na(`DMB|consensus.ag`) & !is.na(`DMA|consensus.ag`),
            yes='Two-chain perfect match',
            no=ifelse(
                test=!is.na(`DMB|consensus.ag`),
                yes='Beta-only perfect match',
                no='Alpha-only perfect match'
            )
        )
    ]
    # Set consensus antigen. Priority is given to the beta chain when available (i.e., for single-chain beta and double-chain matches).
    exp.match.info[,
        DM.pred:=`DMB|consensus.ag`
    ]
    exp.match.info[
        is.na(DM.pred),
        DM.pred:=`DMA|consensus.ag`
    ]
    exp.match.info[
        str_detect(string=DM.pred, pattern='\\|'),
        DM.pred:='Multiple'
    ]
    # General match class.
    #       Required info.
    tmp.data.1 <- str_split(
        string=exp.match.info[, `DMB|consensus.ag`],
        pattern='\\|'
    )
    tmp.data.2 <- str_split(
        string=exp.match.info[, `DMA|consensus.ag`],
        pattern='\\|'
    )
    # tmp.check <- length(tmp.data.1)==length(tmp.data.2)
    tmp.data <- unlist(lapply(X=1:length(tmp.data.1), FUN=function(i){
        tmp.data <- intersect(x=tmp.data.1[[i]], y=tmp.data.2[[i]])
        tmp.data <- tmp.data[!is.na(tmp.data)]
        tmp.data <- if(length(tmp.data)>0) paste0(tmp.data, collapse='|') else NA
        return(tmp.data)
    }))
    exp.match.info[['DMG|consensus.ag']] <- tmp.data
    #       Perfect matches.
    tmp.classes <- c('Two-chain perfect match')
    exp.match.info[
        DM.chain.class%in%tmp.classes & DM.pred!='Multiple' & !is.na(`DMG|consensus.ag`),
        DM.match.class:='Unambiguous'
    ]
    #       Umambiguous matches when matched to both chains.
    exp.match.info[
        DM.chain.class%in%tmp.classes & DM.pred!='Multiple' & is.na(`DMG|consensus.ag`),
        DM.match.class:='Unambiguous'
    ]
    #       Ambiguous matches when matched to both chains.
    exp.match.info[
        is.na(DM.match.class) & DM.chain.class%in%tmp.classes & DM.pred=='Multiple',
        DM.match.class:='Ambiguous'
    ]
    #       Umambiguous matches when matched to either chain.
    exp.match.info[
        is.na(DM.match.class) & (DM.chain.class=='Alpha' | DM.chain.class=='Beta') & DM.pred!='Multiple',
        DM.match.class:='Unambiguous'
    ]
    #       Ambiguous matches when matched to either chain.
    exp.match.info[
        is.na(DM.match.class) & (DM.chain.class=='Alpha' | DM.chain.class=='Beta') & DM.pred=='Multiple',
        DM.match.class:='Ambiguous'
    ]
    # ---> Final format.
    # Bring relevant columns to the front.
    tmp.cols <- c('DM.pred', 'DM.chain.class', 'DM.match.class', 'DMG|consensus.ag', 'DMB|consensus.ag', 'DMA|consensus.ag')
    tmp.cols <- c(tmp.cols, setdiff(x=colnames(exp.match.info), y=tmp.cols))
    exp.match.info <- exp.match.info[, ..tmp.cols]
}


cat('\n\n')
### ------------------- Consensus of predictions -------------------- ###
cat('### ------------------- Consensus of predictions -------------------- ###\n')

# ---> Merge prediction data from multiple approaches.
# By clustering by similarity.
# pred.summ <- copy(bkup)
# By match w/ reference.
pred.summ <- merge(
    x=pred.summ,
    y=ref.match.info,
    by=c('qry.clone.id', 'donor.id.tag'),
    all=TRUE
)
# By match w/ experimental data (if available).
if(!is.null(exp.res) & !is.null(qry.donor.tag)) pred.summ <- merge(
    x=pred.summ,
    y=exp.match.info,
    by=c('qry.clone.id', 'donor.id.tag'),
    all=TRUE
)

# ---> Consensus annotations.
# @ Consensus prediction & Consensus approach class.
# Priority is given in the following order: 1. Donor-matched approach (if donor-matched ag-specific TCR data are available); 2. Reference match approach; 3. Unbiased clustering model (UCM) approach.
pred.summ[, consensus.pred:=as.character(x=NA)]
# 1. Donor-matched approach (if experimental data are available)
if(!is.null(exp.res) & !is.null(qry.donor.tag)){
    tmp.vals <- c('Two-chain perfect match', 'Alpha-only perfect match', 'Beta-only perfect match')
    pred.summ[
        is.na(consensus.pred) & !is.na(DM.pred) & DM.pred!='Multiple' & DM.chain.class%in%tmp.vals,
        `:=`(
            consensus.pred=DM.pred,
            consensus.approach='Donor-matched'
        )
    ]
}
# 2. Reference match approach
tmp.vals <- c('Two-chain perfect match', 'Beta-only perfect match', 'Perfect beta match with alpha quasi match', 'Perfect alpha match with beta quasi match')
pred.summ[
    is.na(consensus.pred) & !is.na(RM.pred) & RM.pred!='Multiple' & RM.chain.class%in%tmp.vals,
    `:=`(
        consensus.pred=RM.pred,
        consensus.approach='Reference match'
    )
]
# 3. Unbiased clustering model (UCM) approach.
pred.summ[
    is.na(consensus.pred) & !is.na(UCM.pred) & UCM.pred!='Multiple',
    `:=`(
        consensus.pred=UCM.pred,
        consensus.approach='UCM'
    )
]

# @ General consensus class. DEPRECATED.
# if(is.null(exp.res) | is.null(qry.donor.tag)) pred.summ[, DM.pred:=NA]
# tmp.data <- pred.summ[,
#     .(DM.pred, RM.pred, UCM.pred)
# ]
# tmp.data[!is.na(DM.pred), DM.pred:='DM']
# tmp.data[!is.na(RM.pred), RM.pred:='RM']
# tmp.data[!is.na(UCM.pred), UCM.pred:='UCM']
# tmp.data <- apply(X=tmp.data, MARGIN=1, FUN=function(x){
#     x <- x[!is.na(x)]
#     if(length(x)==0) return(NA)
#     x <- paste0(x, collapse='|')
#     return(x)
# })
# pred.summ[['consensus.gen.class']] <- tmp.data
# pred.summ[
#     !is.na(consensus.pred),
#     consensus.gen.class:='Single',
# ]


# ---> Final format.
# @ Final column name fix.
colnames(pred.summ)[colnames(pred.summ) %like% '^score\\.'] <- paste0(
    'UCM|',
    colnames(pred.summ)[colnames(pred.summ) %like% '^score\\.']
)
# @ Bring relevant columns to the front.
# General columns.
tmp.cols.1 <- c(
    'qry.clone.id', 'donor.id.tag',
    'consensus.pred', 'consensus.approach'
)
if(!is.null(exp.res)) tmp.cols.1 <- c(
    tmp.cols.1,
    'DM.pred', 'DM.chain.class', 'DM.match.class'
)
tmp.cols.1 <- c(
    tmp.cols.1,
    'RM.pred', 'RM.chain.class', 'RM.match.class',
    'UCM.pred', 'UCM.ag.class', 'UCM.max.ags', 'UCM.clust.class', 'UCM.ratio'
)
# Donor-matched approach columns.
tmp.cols.2 <- setdiff(x=colnames(pred.summ), y=tmp.cols.1)
tmp.cols.2 <- c(
    tmp.cols.2[tmp.cols.2 %like% '^DMG\\|'],
    tmp.cols.2[tmp.cols.2 %like% '^DMB\\|'],
    tmp.cols.2[tmp.cols.2 %like% '^DMA\\|']
)
tmp.cols.1 <- c(tmp.cols.1, tmp.cols.2)
# Reference match approach columns.
tmp.cols.2 <- setdiff(x=colnames(pred.summ), y=tmp.cols.1)
tmp.cols.2 <- c(
    tmp.cols.2[tmp.cols.2 %like% '^RMG\\|'],
    tmp.cols.2[tmp.cols.2 %like% '^RMB\\|'],
    tmp.cols.2[tmp.cols.2 %like% '^RMA\\|']
)
tmp.cols.1 <- c(tmp.cols.1, tmp.cols.2)
# UCM approach columns.
tmp.cols.2 <- setdiff(x=colnames(pred.summ), y=tmp.cols.1)
tmp.cols.2 <- c(
    tmp.cols.2[tmp.cols.2 %like% '^UCM\\|score'],
    tmp.cols.2[tmp.cols.2 %like% '^UCM\\|M\\|']
)
tmp.cols.1 <- c(tmp.cols.1, tmp.cols.2)
# Done!
tmp.check <- length(setdiff(x=colnames(pred.summ), y=tmp.cols.1))==0
if(!tmp.check) stop('Unexpected error during final column ordering.\n')
pred.summ <- pred.summ[, ..tmp.cols.1]
# Final row order.
setorderv(x=pred.summ, cols=c('qry.clone.id', 'donor.id.tag', 'consensus.pred', 'consensus.approach'))
# Small things to consider.
if(is.null(qry.donor.tag)) pred.summ[, donor.id.tag:=NULL]
# Save all prediction details.
tmp.file.name <- paste0(reports.path, '/PredictionDetails.csv')
fwrite(file=tmp.file.name, x=pred.summ, na=NA, quote=TRUE)


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
