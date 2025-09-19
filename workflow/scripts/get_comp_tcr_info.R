cat('\n\n')
############    ------   Retrieve comprehensive   -------    ############
############    ----------   TCR information   ----------    ############
cat('############    ------   Retrieve comprehensive   -------    ############\n')
cat('############    ----------   TCR information   ----------    ############\n')

# By: Vicente Fajardo
# By: Donaldo Sosa


### -------------------------- Description -------------------------- ###


cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### -------------------------- Dependencies ------------------------- ###\n')
library(optparse)
library(Seurat)
library(data.table)
library(tidyr)
library(stringr)


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')

# ---> Arguments from command line.
option.list <- list(
    make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Character, absolute path to reports directory.\n"),
    make_option(opt_str="--OptsFile1", type="character", default=NULL, dest="yaml.file.1", help="Char, absolute path to yaml file for workflow; general options.\n"),
    make_option(opt_str="--OptsFile2", type="character", default=NULL, dest="yaml.file.2", help="Char, absolute path to yaml file for workflow; project-specific options.\n"),
     make_option(opt_str="--RefID", type="character", default=NULL, dest="ref.name", help="Char, indicates the reference ID for the results to be evaluated.\n"),
    make_option(opt_str="--AggrPath", type="character", default=NULL, dest="aggr.res.path", help="Character, absolute path to vdj aggregation directory.\n"),
    make_option(opt_str="--GExData", type="character", default=NULL, dest="seurat.obj.file", help="Character, absolute path to seurat object file.\n"),
    make_option(opt_str="--VDJAggrTable", type="character", default=NULL, dest="aggr.table.1.file", help="Character, absolute path to vdj aggr table file.\n"),
    make_option(opt_str="--GExAggrTable", type="character", default=NULL, dest="aggr.table.2.file", help="Character, absolute path to GEx aggr table file.\n"),
    make_option(opt_str="--RawData", type="character", default=NULL, dest="raw.data.file", help="Absolute path to original 10x output directory (same as the one used to create seurat object), where features.tsv file must exist. Necessary if an output describing both ENSEMBL ID and gene name (regardless of the ID type used for the seurat object) must be provided.\n"),
    make_option(opt_str="--FeatureID", type="character", default='name', dest="feature.id", help="Character defining the feature ID type that was taken from 10X cellranger output (either from features.tsv file or from the h5 matrix) to originally create the seurat object provided. Default is 'name' (features.tsv column 2), though possible values are 'name' and 'ensembl', the last one corresponding to ENSEMBL ID. This option applies only if a cellranger output path is provided.\n"),
    make_option(opt_str="--ClustsLab", type="character", default=NULL, dest="clusts.lab", help="Character, label for column defining clusters in metadata.\n"),
    # ---> Filtering.
    # Pre-DEA filtering.
    make_option(opt_str="--FiltCriteria", type="character", default=NULL, dest="filt.criteria.file", help="Character, absolute path to csv file listing the criteria to consider in oder to filter the expression data (seurat object) previous to apply the differential analysis. File must include the following specifications: 1) first column (unnamed) must define the row names, and 2) two rows must be included, 'keep' and 'discard' (described below; any other row name may be read but will not be considered at all); 3) there's not any limit for the number of columns, however, 4) each column name should match that of a tag already defined in the input seurat object's metadata; 5) the values for each column must list the set of valid, existing tag values that are to be kept (as defined in the row 'keep') or discarded (as defined in the row 'discard'); 6) tag values defined, other than being valid and defined for a given tag, must be separated by a semicolon if more than one provided (e.g., 'val.1;val.2') and 7) if NA values must be considered (whether to be kept or discarded), you must add 'NA' to the list (e.g., 'val.1;'NA';val.2'); finally, 8) you must define either only the values of a tag to be kept or discarded (they're mutually exclusive criteria), then 9) leaving as NA (beware, NA value not defined in between quotation marks and therefore read as bonafide NA values in R) the non-defined cell.\n"),
    # ---> Rules to create new tags and transform the current ones.
    make_option(opt_str="--NewTags", type="character", default=NULL, dest="new.tags.file", help="Character, absolute path to file describing the rules to add new tags based on the combination of the ones already defined. See more details in the description above.\n"),
    make_option(opt_str="--NewMeta", type="character", default=NULL, dest="new.meta.file", help="Absolute path to file listing new metadata that are relevant for analysis. Details about the format can be found in module 'add_new_tags'.\n"),
    make_option(opt_str="--TransformRules", type="character", default=NULL, dest="transform.rules.file", help="Character, absolute path to file describing the rules to transform any tags to be transferred. See any further documentation about the format of this file in this script.\n"),
    # ---> Specific metadata
    make_option(opt_str="--MetaCols", type="character", default=NULL, dest="donor.meta.cols", help="Character, set of columns from full metadata that should be included in the final output. Taken into account only if an extra file of donors metadata is provided as well, otherwise ignored.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser <- OptionParser(option_list=option.list);
opt <- parse_args(opt.parser);
# Moving options to their own variables
reports.path <- opt$reports.path
yaml.file.1 <- opt$yaml.file.1
yaml.file.2 <- opt$yaml.file.2
ref.name <- opt$ref.name
aggr.res.path <- opt$aggr.res.path
seurat.obj.file <- opt$seurat.obj.file
aggr.table.1.file <- opt$aggr.table.1.file
aggr.table.2.file <- opt$aggr.table.2.file
raw.data.file <- opt$raw.data.file
feature.id <- opt$feature.id
clusts.lab <- opt$clusts.lab
filt.criteria.file <- opt$filt.criteria.file
new.tags.file <- opt$new.tags.file
new.meta.file <- opt$new.meta.file
transform.rules.file <- opt$transform.rules.file
donor.meta.cols <- opt$donor.meta.cols
if(!is.null(donor.meta.cols)) donor.meta.cols <- unique(unlist(str_split(string=donor.meta.cols, pattern=';')))

# ---> Options file.
# Load yaml files. Priority is given to dataset-specific options while hanlding potentially redundant entries.
if(!file.exists(yaml.file.1)) stop(paste0('Attempted definition of yaml file (general options) failed. Next file could not be found:\n', yaml.file.1, '\n'))
if(!file.exists(yaml.file.2)) stop(paste0('Attempted definition of yaml file (project-specific options) failed. Next file could not be found:\n', yaml.file.2, '\n'))
tmp.data.1 <- yaml::read_yaml(file=yaml.file.1)
tmp.data.2 <- yaml::read_yaml(file=yaml.file.2)
gen.wflow.opts <- c(tmp.data.2, tmp.data.1[!names(tmp.data.1) %in% names(tmp.data.2)])
# Define files and general options.
srt.donor.id <- gen.wflow.opts$gex_data[[ref.name]]$donor_id_tag
if(is.null(srt.donor.id)) srt.donor.id <- 'donor.id.tag'
srt.library.id <- gen.wflow.opts$gex_data[[ref.name]]$lib_id_tag
if(is.null(srt.library.id)) srt.library.id <- 'library.id.tag'
mod.sig.path <- gen.wflow.opts$mod_sig_path
exp.set.file <- gen.wflow.opts$exp_set_file
stat.cell.thold <- gen.wflow.opts$stat_cell_thold
if(is.null(stat.cell.thold)) stat.cell.thold <- 10
functions.file <- gen.wflow.opts$paths$R_functions
if(!file.exists(functions.file)) stop('Faulty R functions file.\n')

# ---> Hardcoded arguments.
# ---> Color defintions.
# @ Other terms' colors.
# signatures.col.scale <- c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# ---> File definitions.
# TCR data.
tcr.clones.file <- paste0(aggr.res.path, '/clonotypes_aggr.csv')
tcr.anns.file <- paste0(aggr.res.path, '/filtered_contig_annotations_aggr.csv')
# ---> Check directories and files.
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(
  tcr.clones.file, tcr.anns.file,
  seurat.obj.file,
  aggr.table.1.file, aggr.table.2.file,
  new.meta.file, new.tags.file, transform.rules.file
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### --------------------------- Functions --------------------------- ###
source(functions.file)


cat('\n\n')
### --------------------------- Load data --------------------------- ###
cat('### --------------------------- Load data --------------------------- ###\n')

# ---> Feature info.
if(!is.null(raw.data.file)){
  feature.info.file <- unique(list.files(path=raw.data.file, pattern='features.tsv', full.names=TRUE))
  if(length(feature.info.file)!=1) stop('File with gene names -features.tsv- not properly defined.\n')
  is.zipped <- grepl(x=feature.info.file, pattern='\\.gz')
  if(is.zipped){
    # Then, read relatiosnships between feature names and ENSEMBL IDs.
    tmp.dir <- paste0(raw.data.file, '/tmp_', paste0(sample(x=LETTERS, size=10, replace=TRUE), collapse=''))
    tmp.cmmnd <- paste0('mkdir ',  tmp.dir, ' && cp ', feature.info.file, ' ', tmp.dir)
    system(command=tmp.cmmnd)
    tmp.file.name <- paste0(tmp.dir, '/features.tsv')
    tmp.cmmnd <- paste0('gunzip ', tmp.file.name, '.gz')
    system(tmp.cmmnd)
    feature.info <- read.delim(file=tmp.file.name, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
    tmp.cmmnd <- paste0('rm -r ', tmp.dir)
    system(tmp.cmmnd)
  }else{
    feature.info <- read.delim(file=feature.info.file, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
    tmp.cmmnd <- paste0('rm -r ', tmp.dir)
  }
}else{
  feature.info <- NULL
}

# ---> TCR data.
tcr.clones <- fread(file=tcr.clones.file)
tcr.anns <- fread(file=tcr.anns.file)

# ---> Aggr table.
id.cols <- c('sample_id', 'library_id')
aggr.table.1 <- fread(file=aggr.table.1.file)
tmp.check <- sum(id.cols %in% colnames(aggr.table.1))==1
if(!tmp.check) stop('Sample ID not defined in TCR aggr table.\n')
tmp.col <- id.cols[id.cols %in% colnames(aggr.table.1)]
aggr.table.1[, library_id:=get(tmp.col)]
aggr.table.1[, library_id:=str_replace(string=library_id, pattern='_TCR$', replacement='')]
aggr.table.2 <- fread(file=aggr.table.2.file)
tmp.check <- sum(id.cols %in% colnames(aggr.table.2))==1
if(!tmp.check) stop('Sample ID not defined in GEx aggr table.\n')
tmp.col <- id.cols[id.cols %in% colnames(aggr.table.2)]
aggr.table.2[, library_id:=get(tmp.col)]
aggr.table.2[, library_id:=str_replace(string=library_id, pattern='_Gex$|_GEX$|_GEx$', replacement='')]
aggr.table <- merge(x=aggr.table.1, y=aggr.table.2, by='library_id')
# Sanity check.
tmp.check <- aggr.table.1[, .N] == aggr.table.2[, .N] & aggr.table.1[, .N] == aggr.table[, .N]
if(!tmp.check) stop('Faulty aggr tables provided. Base library names (disregarding data modality suffix) have to match between individual GEx and TCR tables, and the library universe should captured for both.\n')

# ---> Set of genes for expression
# Load only if necessary.
if(!is.null(exp.set.file)){
    exp.set <- fread(file=exp.set.file)
    exp.set <- exp.set[[1]]
    names(exp.set) <- if(feature.id=='ensembl') names(translate.ids(ids=exp.set, ensembl=FALSE)) else exp.set
    # Remove genes without values defined in reference.
    tmp.check <- any(is.na(names(exp.set)))
    if(tmp.check){
        tmp.check <- paste0(exp.set[is.na(names(exp.set))], collapse='\n')
        exp.set <- exp.set[is.na(names(exp.set))]
        if(length(exp.set)==0) stop('Gene set for expression summary was provided, but none in the list could be found in reference.\n')
        tmp.check <- paste0('From gene set for expression summary, next symbols could not be found in reference:\n', tmp.check, '\n')
        warning(tmp.check)
    }
}else{
    exp.set <- NULL
}

# ---> Seurat object.
# Load data.
seurat.obj <- readRDS(file=seurat.obj.file)

# ---> Tags of interest.
transform.rules <- if(!is.null(transform.rules.file)) read.csv(file=transform.rules.file, stringsAsFactors=FALSE) else NULL
new.tags.rules <- if(!is.null(new.tags.file)) read.csv(file=new.tags.file, stringsAsFactors=FALSE) else NULL
new.meta.data <- if(!is.null(new.meta.file)) read.csv(file=new.meta.file, stringsAsFactors=FALSE) else NULL
cat('All rules for manipulation of tags loaded (if any)...\n')

# ---> Subset criteria
if(!is.null(filt.criteria.file)){
  if(!file.exists(filt.criteria.file)) stop(paste0('A file for pre-DEA filtering purposes was provided, but it does not exist. If you wish no file to be considered for this step, simply drop the flag.\nHere is the file:\n', filt.criteria.file, '\n'))
  filt.criteria <- read.csv(file=filt.criteria.file, stringsAsFactors=FALSE, row.names=1, na.strings=c("NA", ""), colClasses='character')
  cat('Subset criteria loaded...\n')
}else{
  filt.criteria <- NULL
}

# ---> Module scoring
# Retrieve signatures in provided folder. NOTE: These are needed to discern between signatures that may have already been added to the seurat before and the ones requested to be added.
if(!is.null(mod.sig.path)){
    signs.of.int <- list.files(path=mod.sig.path, pattern='_signature.csv$')
    signs.of.int <- str_replace(
        string=signs.of.int,
        pattern='_signature.csv$',
        replacement='.score'
    )
    signs.of.int <- str_replace_all(
        string=signs.of.int,
        pattern='_',
        replacement='.'
    )
}


cat('\n\n')
### ------------------------- Preprocessing ------------------------- ###
cat('### ------------------------- Preprocessing ------------------------- ###\n')

# ---> New tags.
# Proceed only if required.
if(!is.null(new.tags.rules)){
    # Get new tags.
    seurat.obj <- add.new.tags(seurat.obj=seurat.obj, new.tags.rules=new.tags.rules)
    # Add newly added tags to the ones of interest.
    tags.of.int <- unique(c(tags.of.int, new.tags.rules$tag.name))
}
if(!is.null(new.meta.data)){
    # Get new tags.
    seurat.obj <- add.new.tags(seurat.obj=seurat.obj, new.meta.data=new.meta.data)
    # Add newly added tags to the ones of interest. NOT NECESSARY HERE. REMOVE CHUNK IN FOLLOWING VERSION UNLESS FOUND NECESSARY IN PARTICULAR CASES AT SOME POINT (IN WHICH, THE CHUNK WOULD NEED AMMENDING).
    # tags.of.int <- unique(c(
    #     tags.of.int,
    #     setdiff(x=colnames(new.meta.data), y='barcode')
    # ))
}
# ---> Tranform currently existing tags.
if(!is.null(transform.rules)){
    seurat.obj <- transform.tags(seurat.obj=seurat.obj, transform.rules=transform.rules)
}

# ---> Filtering according to tags criteria.
# Apply filtering if criteria were provided.
if(!is.null(filt.criteria)){
  cat('\n')
  cat('# ---> Filtering according to tags criteria.\n')
  cat('Subsetting seurat object according to the criteria supplied...\n')
  # Check tags and features definition, as well as their possible values.
  if(!all(colnames(filt.criteria) %in% colnames(seurat.obj@meta.data))) stop('Some criteria term may not be properly defined in seurat object.')
  # Check rownames for the criteria terms.
  if(!all(rownames(filt.criteria)=='discard'|rownames(filt.criteria)=='keep')) stop('In tags criteria, terms not properly defined as being criteria to discard or keep cells.')
  # Subset.
  seurat.obj <- subset.seurat.obj(seurat.obj=seurat.obj, tags.criteria=filt.criteria, feats.criteria=NULL)
  cat('Seurat object subsetted!\n')
}

# ---> Module scoring
# Annotate seurat object w/ requested signatures.
if(!is.null(mod.sig.path)){
    seurat.obj <- get.module.scores(
        seurat.obj=seurat.obj,
        module.feats.dir=mod.sig.path,
        is.ensembl=feature.id=='ensembl',
        reports.path=NULL
    )
    # Confirm that annotation went well.
    tmp.check <- all(signs.of.int %in% colnames(seurat.obj@meta.data))
    if(!tmp.check) stop('Something went wrong during module scoring...\n')
}


# ---> Seurat object.
if(!srt.donor.id %in% colnames(seurat.obj@meta.data)) stop('Faulty donor ID metadata column provided. Not defined in seurat object.\n')
seurat.obj@meta.data[, 'donor.id.tag'] <- seurat.obj@meta.data[, srt.donor.id]
if(!srt.library.id %in% colnames(seurat.obj@meta.data)) stop('Faulty library ID metadata column provided. Not defined in seurat object.\n')
seurat.obj@meta.data[, 'library.id.tag'] <- seurat.obj@meta.data[, srt.library.id]

# ---> Metadata
# @ Basic metadata.
meta.data <- as.data.table(seurat.obj@meta.data)
meta.data[, barcode:=Cells(seurat.obj)]
meta.data[, clusters.tag:=as.character(get(clusts.lab))]
# rm(seurat.obj)
# Add sample ID info.
tmp.data <- aggr.table[, .(sample.sffx=as.character(sample.sffx), sample.id=library_id)]
meta.data[, sample.sffx:=str_extract(string=barcode, pattern='\\d+$')]
meta.data <- merge(x=meta.data, y=tmp.data, by='sample.sffx', all.x=TRUE, all.y=FALSE, sort=FALSE)
meta.data[, sample.sffx:=NULL]
# Tissue label. If none provided, we will assumme, data from healthy tissue was provided.
if(!'tissue.tag' %in% colnames(meta.data)){
    meta.data[, tissue.tag:='normal.lung']
    single.tissue <- TRUE
}else{
    tmp.check <- meta.data[!is.na(tissue.tag), uniqueN(tissue.tag)]
    if(tmp.check==1){
        meta.data[, tissue.tag:='normal.lung']
        single.tissue <- TRUE
    }else{
        single.tissue <- FALSE
    }
}
# Change tissue labels.
tmp.vals <- c(`normal.lung`='healthy', `tumor.lung`='tumor')
meta.data[, tissue.tag:=tmp.vals[tissue.tag]]



### ------------------------- Main program -------------------------- ###
cat('\n')
### ------------------------- Preprocessing ------------------------- ###
cat('### ------------------------- Preprocessing ------------------------- ###\n')


# ---> TCR data preflights.
# Find all suffixes
all.sffxs <- tcr.anns[, str_extract(string=barcode, pattern='-.+$')]
all.sffxs <- unique(all.sffxs)
all.sffxs <- gtools::mixedsort(str_replace(string=all.sffxs, pattern='^-', replacement=''))

# ---> Barcode-clone associations.
cell.clone.info <- tcr.anns[is_cell & high_confidence & full_length & productive]
cols.to.keep <- c(
    `barcode`='barcode',
    `chain`='chain', `v.gene`='v_gene', `j.gene`='j_gene',
    `cdr3.nt.seq`='cdr3_nt',
    `cdr3.aa.seq`='cdr3',
    `umi.count`='umis',
    `clonotype.tag`='raw_clonotype_id'
)
cell.clone.info <- cell.clone.info[, ..cols.to.keep]
colnames(cell.clone.info) <- names(cols.to.keep)
cell.clone.info <- unique(cell.clone.info)

# ---> vdj gene usage info.
# General vdj gene usage info.
# We obtain unique v and j genes per clonotype. When cellranger has defined multiple v and j genes for a given clonotype, we take the one that is most supported. The most supported is the gene that has the largest amount of entries (contigs) in the table below.
#       We first split the dataset into clonotypes w/ a single V gene and a single J gene definition from the rest, for which we need to pick one per the criterion established right above.
tmp.data <- cell.clone.info[,
    .(
        v.gene.count=uniqueN(v.gene),
        j.gene.count=uniqueN(j.gene)
    ),
    by=.(
        clonotype.tag,
        chain, cdr3.nt.seq
    )
]
tmp.clonotypes <- tmp.data[v.gene.count>1 | j.gene.count>1, unique(clonotype.tag)]
#          For well-defined clonotypes, collect genes.
vdj.gene.info.1 <- cell.clone.info[
    !clonotype.tag %in% tmp.clonotypes,
    .(
        v.gene=unique(v.gene),
        j.gene=unique(j.gene)
    ),
    by=.(
        clonotype.tag,
        chain, cdr3.nt.seq
    )
]
#          For the rest of the clonotypes, select best gene according to the established criterion.
vdj.gene.info.2 <- cell.clone.info[
    clonotype.tag %in% tmp.clonotypes,
    .(
        v.gene=.SD[,
            .(support=sum(umi.count)),
            by=v.gene
        ][support==max(support), paste0(unique(v.gene), collapse='|')],
        j.gene=.SD[,
            .(support=sum(umi.count)),
            by=j.gene
        ][support==max(support), paste0(unique(j.gene), collapse='|')]
    ),
    by=.(
        clonotype.tag,
        chain, cdr3.nt.seq
    )
]
#           Out of this second group, for those clonotypes with multiple v and j genes that were equally supported according to our criteorion, make a guess and keep only one gene per set of clonotype, chain and CDR3 nt seq.
vdj.gene.info.2[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
vdj.gene.info.2[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
#           Concatenate info from individual approaches.
vdj.gene.info <- rbindlist(l=list(vdj.gene.info.1, vdj.gene.info.2), use.names=TRUE)
# Sanity check. Confirm we have unique entries for each combination of clonotype ID, chain and CDR3 sequence.
tmp.check <- vdj.gene.info[, .N, by=.(clonotype.tag, chain, cdr3.nt.seq)][N>1, .N==0]
if(!tmp.check) stop('Unexpected error post vdj data retrieval.\n')

# ---> Unique sequences per clonotype.
tmp.data <- unique(cell.clone.info[, .(clonotype.tag, chain, cdr3.nt.seq, cdr3.aa.seq)])
# Add vdj gene info.
tmp.data <- merge(
    x=tmp.data, y=vdj.gene.info,
    by=c('clonotype.tag', 'chain', 'cdr3.nt.seq'),
    all=TRUE
)
tmp.check <- tmp.data[is.na(v.gene) | is.na(j.gene), .N==0]
if(!tmp.check) stop('Error while defining unique sequences per clonotype (1).\n')
# Collapse chain-specific info per clonotype. That is, if several, collapse info from multiple CDR3 and V/J gene sequences into single entries while keeping track of the correspondance between CDR3 and gene sequences.
tmp.data <- tmp.data[,
    .(
        cdr3.nt.seq=paste0(cdr3.nt.seq, collapse=';'),
        cdr3.aa.seq=paste0(cdr3.aa.seq, collapse=';'),
        v.gene=paste0(v.gene, collapse=';'),
        j.gene=paste0(j.gene, collapse=';')
    ),
    by=.(clonotype.tag, chain)
]
# For clonotypes w/ multiple sequences of either chain, confirm that we have complete info for all four relevant pieces: CDR3 nt and aa sequences and V and J genes.
tmp.check <- tmp.data[,
    .(
        cdr3.nt.seq=str_count(string=cdr3.nt.seq, pattern=';'),
        cdr3.aa.seq=str_count(string=cdr3.aa.seq, pattern=';'),
        v.gene=str_count(string=v.gene, pattern=';'),
        j.gene=str_count(string=j.gene, pattern=';')
    )
]
tmp.check <- apply(X=tmp.check, MARGIN=1, FUN=function(x) length(unique(x)))
tmp.check <- all(tmp.check==1)
if(!tmp.check) stop('Error while defining unique sequences per clonotype (2).\n')
# Spread chain-specific info into multiple columns.
tmp.cols <- c('cdr3.nt.seq', 'cdr3.aa.seq', 'v.gene', 'j.gene')
these.col.names <- list(
    `cdr3.nt.seq`=c('TRA'='cdr3a.nt.seq', 'TRB'='cdr3b.nt.seq'),
    `cdr3.aa.seq`=c('TRA'='cdr3a.aa.seq', 'TRB'='cdr3b.aa.seq'),
    `v.gene`=c('TRA'='tra.v', 'TRB'='trb.v'),
    `j.gene`=c('TRA'='tra.j', 'TRB'='trb.j')
)
these.col.names <- lapply(X=these.col.names, FUN=function(x) c(c(`clonotype.tag`='clonotype.tag'), x))
tmp.data <- lapply(X=tmp.cols, FUN=function(tmp.col){
    tmp.cols <- c('clonotype.tag', 'chain', tmp.col)
    tmp.data <- tmp.data[, ..tmp.cols]
    tmp.data <- as.data.table(spread(data=tmp.data, key=chain, value=tmp.col, fill=NA))
    colnames(tmp.data) <- these.col.names[[tmp.col]][colnames(tmp.data)]
    return(tmp.data)
})
clone.info <- Reduce(x=tmp.data, f=function(x, y) merge(x=x, y=y, by='clonotype.tag', all=TRUE)); rm(tmp.data)
# Final sanity check.
tmp.check <- clone.info[, uniqueN(clonotype.tag)]==cell.clone.info[, uniqueN(clonotype.tag)]
if(!tmp.check) stop('Error while defining unique sequences per clonotype (3).\n')


cat('\n')
### ------------------------ Cell-based info ------------------------ ###
cat('### ------------------------ Cell-based info ------------------------ ###\n')

# ---> Piece parts together.
# Find associations between barcodes and clonotypes.
tmp.data <- unique(cell.clone.info[, .(barcode, clonotype.tag)])
# Identify barcodes in TCR data that weren't found in the gene expression metadata.
# tmp.data[!barcode %in% meta.data[, barcode], uniqueN(barcode)]
meta.data <- merge(x=meta.data, y=tmp.data, by='barcode', all=FALSE)
meta.data <- merge(x=meta.data, y=clone.info, by='clonotype.tag')
# ---> Output results tables.
# @ Full info set.
tmp.file.name <- paste0(reports.path, '/CellBasedTCRData.csv')
fwrite(file=tmp.file.name, x=meta.data, quote=TRUE, na=NA)
# @ Basic info set.
cols.to.keep <- c(
    'library.id.tag', 'barcode',
    'tissue.tag',
    'donor.id.tag',
    'clusters.tag',
    colnames(clone.info)
)
tmp.check <- all(cols.to.keep %in% colnames(meta.data))
if(!tmp.check) stop('Unexpected error.\n')
tmp.file.name <- paste0(reports.path, '/CellBasedTCRData_Clean.csv')
tmp.data <- meta.data[!is.na(donor.id.tag), ..cols.to.keep]
fwrite(file=tmp.file.name, x=tmp.data, quote=TRUE, na=NA)


cat('\n')
### --------------------- Clonotype-based info ---------------------- ###
cat('### --------------------- Clonotype-based info ---------------------- ###\n')

# ---> General information.
tmp.data.1 <- unique(meta.data[
    !is.na(donor.id.tag),
    .(
        cdr3b.nt.seq, cdr3a.nt.seq,
        cdr3b.aa.seq, cdr3a.aa.seq,
        trb.v, trb.j,
        tra.v, tra.j
    ),
    by=clonotype.tag
])

# ---> Identify potentially public TCRs.
tmp.data.2 <- meta.data[
    !is.na(donor.id.tag),
    .(
        public.status=ifelse(
            test=uniqueN(donor.id.tag)>1,
            yes='likely public',
            no='private'
        )
    ),
    by=.(clonotype.tag)
]

# ---> Counts per donor.
tmp.data.3 <- meta.data[
    !is.na(donor.id.tag),
    .(abs.freq=.N),
    by=.(clonotype.tag, donor.id.tag, tissue.tag)
]
tmp.tmp <- meta.data[
    !is.na(donor.id.tag),
    .(total.freq=.N),
    by=.(donor.id.tag, tissue.tag)
]
tmp.data.3 <- merge(x=tmp.data.3, y=tmp.tmp, by=c('donor.id.tag', 'tissue.tag'))
tmp.data.3[, rel.freq:=abs.freq/total.freq]
tmp.check <- tmp.data.3[, .(total.freq=sum(rel.freq)), by=.(donor.id.tag, tissue.tag)][, all(round(x=total.freq, digits=1)==1)] # 10 should be enough rounding, but it might need to be further adjusted for some instances.
if(!tmp.check) stop('Unexpected error while calculating relative frequencies.\n')
tmp.data.3[, total.freq:=NULL]
# Spread according to tissue.
tmp.vals <- c('abs.freq', 'rel.freq')
tmp.data.3 <- lapply(X=tmp.vals, FUN=function(tmp.col){
    tmp.cols <- c('donor.id.tag', 'tissue.tag', 'clonotype.tag', tmp.col)
    tmp.data <- tmp.data.3[, ..tmp.cols]
    tmp.data <- as.data.table(spread(data=tmp.data, key=tissue.tag, value=tmp.col, fill=0))
    return(tmp.data)
})
names(tmp.data.3) <- tmp.vals
tmp.data.3 <- merge(
    x=tmp.data.3[[tmp.vals[1]]],
    y=tmp.data.3[[tmp.vals[2]]],
    by=c('clonotype.tag', 'donor.id.tag'),
    all=TRUE,
    suffixes=paste0('.', tmp.vals)
)
# Relative frequency ratio (if applies)
tmp.check <- all(c('tumor.rel.freq', 'healthy.rel.freq') %in% colnames(tmp.data.3))
if(tmp.check) tmp.data.3[, rel.freq.ratio:=tumor.rel.freq/healthy.rel.freq]


# ---> Counts per cluster per donor.
tmp.data.4 <- meta.data[
    !is.na(donor.id.tag),
    .(cluster.abs.freq=.N),
    by=.(clonotype.tag, donor.id.tag, cluster=clusters.tag, tissue.tag)
]
tmp.tmp <- meta.data[
    !is.na(donor.id.tag),
    .(total.freq=.N),
    by=.(donor.id.tag, tissue.tag)
]
tmp.data.4 <- merge(x=tmp.data.4, y=tmp.tmp, by=c('donor.id.tag', 'tissue.tag'))
tmp.data.4[, cluster.rel.freq:=cluster.abs.freq/total.freq]
tmp.check <- tmp.data.4[, .(total.freq=sum(cluster.rel.freq)), by=.(donor.id.tag, tissue.tag)][, all(round(x=total.freq, digits=1)==1)] # 10 should be enough rounding, but it might need to be further adjusted for some instances.
if(!tmp.check) stop('Unexpected error while calculating relative frequencies.\n')
tmp.data.4[, total.freq:=NULL]
# Spread according to tissue.
tmp.vals <- c('cluster.abs.freq', 'cluster.rel.freq')
tmp.data.4 <- lapply(X=tmp.vals, FUN=function(tmp.col){
    tmp.cols <- c('donor.id.tag', 'tissue.tag', 'clonotype.tag', 'cluster', tmp.col)
    tmp.data <- tmp.data.4[, ..tmp.cols]
    tmp.data <- as.data.table(spread(data=tmp.data, key=tissue.tag, value=tmp.col, fill=0))
    return(tmp.data)
})
names(tmp.data.4) <- tmp.vals
tmp.data.4 <- merge(
    x=tmp.data.4[[tmp.vals[1]]],
    y=tmp.data.4[[tmp.vals[2]]],
    by=c('clonotype.tag', 'donor.id.tag', 'cluster'),
    all=TRUE,
    suffixes=paste0('.', tmp.vals)
)
# Spread according to cluster.
tmp.vals <- c(
    `healthy.abs`='healthy.cluster.abs.freq',
    `tumor.abs`='tumor.cluster.abs.freq',
    `healthy.rel`='healthy.cluster.rel.freq',
    `tumor.rel`='tumor.cluster.rel.freq'
)
tmp.vals <- tmp.vals[tmp.vals %in% colnames(tmp.data.4)]
tmp.data.4 <- lapply(X=names(tmp.vals), FUN=function(tmp.col.name){
    tmp.col <- tmp.vals[tmp.col.name]
    tmp.cols <- c('clonotype.tag', 'donor.id.tag', 'cluster', tmp.col)
    tmp.data <- tmp.data.4[, ..tmp.cols]
    tmp.data <- as.data.table(spread(data=tmp.data, key=cluster, value=tmp.col, fill=0))
    to.exclude <- c('clonotype.tag', 'donor.id.tag')
    colnames(tmp.data)[!colnames(tmp.data) %in% to.exclude] <- paste(
        tmp.col.name, 'C',
        colnames(tmp.data)[!colnames(tmp.data) %in% to.exclude],
        sep='.'
    )
    return(tmp.data)
})
tmp.data.4 <- Reduce(x=tmp.data.4, f=function(x, y){
    merge(
        x=x, y=y,
        by=c('clonotype.tag', 'donor.id.tag')
    )
})

# ---> Put pieces together.
tmp.data.1 <- merge(x=tmp.data.1, y=tmp.data.2, by='clonotype.tag', all=TRUE)
tmp.data.2 <- merge(x=tmp.data.3, y=tmp.data.4, by=c('clonotype.tag', 'donor.id.tag'), all=TRUE)
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('clonotype.tag'), all=TRUE)
# Sanity check
to.check <- tmp.data[, .(clonotype.tag, donor.id.tag)]
to.check.1 <- to.check[, .N]
to.check.2 <- unique(to.check); to.check.2 <- to.check.2[, .N]
tmp.check <- to.check.1==to.check.2
if(!tmp.check) stop('Unexpected error. There should be as many lines as there are unique pairs of clonotypes and donor IDs.\n')

# Final format.
clone.meta <- copy(tmp.data)
to.order <- c('healthy.rel.freq', 'tumor.rel.freq')
to.order <- to.order[to.order %in% colnames(clone.meta)]
setorderv(x=clone.meta, col=to.order, order=c(-1))

# Simplify if there's a single tissue in the dataset.
if(single.tissue){
    colnames(clone.meta) <- str_replace(string=colnames(clone.meta), pattern='healthy\\.', replacement='')
}

# @ Whole dataset.
tmp.file.name <- paste0(reports.path, '/ClonotypeBasedTCRData.csv')
fwrite(file=tmp.file.name, x=clone.meta, quote=TRUE, na=NA)
# @ For expanded clonotypes only
to.filter <- c('healthy.abs.freq', 'tumor.abs.freq', 'abs.freq')
to.filter <- to.filter[to.filter %in% colnames(clone.meta)]
to.filter <- lapply(X=to.filter, FUN=function(tmp.col){
    tmp.data <- data.table(opt=clone.meta[[tmp.col]]>1)
    colnames(tmp.data) <- tmp.col
    return(tmp.data)
})
to.filter <- Reduce(x=to.filter, f=cbind)
to.filter <- as.matrix(to.filter)
to.filter <- rowSums(to.filter)>0
tmp.data <- clone.meta[to.filter]
tmp.file.name <- paste0(reports.path, '/ClonotypeBasedTCRData_ExpandedOnly.csv')
fwrite(file=tmp.file.name, x=tmp.data, quote=TRUE, na=NA)


cat('\n')
### ----------------------- Module signatures ----------------------- ###
cat('### ----------------------- Module signatures ----------------------- ###\n')

tmp.check <- any(colnames(meta.data) %like% 'score$')
if(tmp.check){
    sig.cols <- colnames(meta.data)[colnames(meta.data) %like% 'score$']
    merge.vals <- c('clonotype.tag', 'donor.id.tag')
    # Retrieve summary statistics for signature scores.
    tmp.data <- lapply(X=sig.cols, FUN=function(this.sig){
        tmp.cols <- c('clonotype.tag', 'donor.id.tag', 'tissue.tag', this.sig)
        names(tmp.cols) <- tmp.cols
        names(tmp.cols)[names(tmp.cols)==this.sig] <- 'sig.score'
        tmp.data <- meta.data[, ..tmp.cols]
        colnames(tmp.data) <- names(tmp.cols)
        # Signature summary per donor per tissue.
        tmp.data.1 <- tmp.data[,
            .(
                abs.freq=.N,
                both=median(sig.score)
            ),
            by=.(clonotype.tag, donor.id.tag)
        ]
        tmp.data.1[
            abs.freq<stat.cell.thold,
            both:=NA
        ]
        tmp.data.1[, abs.freq:=NULL]
        # Signature summary per donor per tissue.
        tmp.data.2 <- tmp.data[,
            .(
                abs.freq=.N,
                median.score=median(sig.score)
            ),
            by=.(
                clonotype.tag, donor.id.tag,
                tissue.tag
            )
        ]
        tmp.data.2[
            abs.freq<stat.cell.thold,
            median.score:=NA
        ]
        tmp.data.2[, abs.freq:=NULL]
        tmp.data.2 <- as.data.table(spread(data=tmp.data.2, key=tissue.tag, value=median.score))
        # Merge and final format.
        tmp.data <- merge(
            x=tmp.data.1,
            y=tmp.data.2,
            by=merge.vals,
            all=TRUE
        )
        colnames(tmp.data)[!colnames(tmp.data)%in%merge.vals] <- paste0(
            this.sig, '.',
            colnames(tmp.data)[!colnames(tmp.data)%in%merge.vals]
        )
        return(tmp.data)
    })
    tmp.data <- Reduce(x=tmp.data, f=function(x, y){
        tmp.data <- merge(
            x=x, y=y,
            by=merge.vals,
            all=TRUE
        )
        return(tmp.data)
    })
    # Simplify if there's a single tissue in the dataset.
    if(single.tissue){
        tmp.cols <- colnames(tmp.data)[str_detect(string=colnames(tmp.data), patter='score.both')]
        tmp.cols <- setdiff(x=colnames(tmp.data), y=tmp.cols)
        tmp.data <- tmp.data[, ..tmp.cols]
        colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='.healthy$', replacement='')
    }
    # Add info to metadata.
    clone.meta <- merge(
        x=clone.meta,
        y=tmp.data,
        by=merge.vals,
        all=TRUE
    )
}


cat('\n')
### ------------------------ Gene expression ------------------------ ###
cat('### ------------------------ Gene expression ------------------------ ###\n')

tmp.check <- !is.null(exp.set)
if(tmp.check){
    merge.vals <- c('clonotype.tag', 'donor.id.tag')
    # Sanity check. Genes in seurat object.
    tmp.check <- all(names(exp.set) %in% row.names(seurat.obj))
    if(!tmp.check){
        tmp.check <- names(exp.set)[!names(exp.set) %in% row.names(seurat.obj)]
        tmp.err <- paste0('Following genes were not found in Seurat object:\n', paste0(tmp.check, collapse='\n'), '\n')
        warning(tmp.err)
        tmp.check <- sum(names(exp.set) %in% row.names(seurat.obj))/length(exp.set) > 0.8
        if(!tmp.check){
            stop('>20% of expression gene set not in seurat object.\n')
        }else{
            exp.set <- exp.set[names(exp.set) %in% row.names(seurat.obj)]
        }
    }
    # Retrieve gene expression values for selected genes.
    exp.data <- GetAssayData(object=seurat.obj, slot='data', assay='RNA')
    exp.data <- exp.data[names(exp.set), ]
    row.names(exp.data) <- exp.set
    exp.data <- as.data.frame(t(as.matrix(exp.data)))
    exp.data$barcode <- row.names(exp.data)
    # Basic data to compute desired values.
    tmp.cols <- c('barcode', 'clonotype.tag', 'donor.id.tag', 'tissue.tag')
    tmp.data <- meta.data[, ..tmp.cols]
    exp.data <- merge(
        x=tmp.data,
        y=exp.data,
        by='barcode',
        all.x=TRUE, all.y=FALSE
    )
    # Retrieve expression summary statistics.
    tmp.data <- lapply(X=exp.set, FUN=function(this.gene){
        tmp.cols <- c('clonotype.tag', 'donor.id.tag', 'tissue.tag', this.gene)
        names(tmp.cols) <- tmp.cols
        names(tmp.cols)[names(tmp.cols)==this.gene] <- 'gene.exp'
        tmp.data <- exp.data[, ..tmp.cols]
        colnames(tmp.data) <- names(tmp.cols)
        # Signature summary per donor per tissue.
        tmp.data.1 <- tmp.data[,
            .(
                abs.freq=.N,
                mean.both=mean(gene.exp),
                fract.both=sum(gene.exp>1)/.N
            ),
            by=.(clonotype.tag, donor.id.tag)
        ]
        tmp.data.1[
            abs.freq<stat.cell.thold,
            `:=`(mean.both=NA, fract.both=NA)
        ]
        tmp.data.1[, abs.freq:=NULL]
        # Expression summary per donor per tissue.
        tmp.data.2 <- tmp.data[,
            .(
                abs.freq=.N,
                 mean=mean(gene.exp),
                 fract=sum(gene.exp>1)/.N
            ),
            by=.(
                clonotype.tag, donor.id.tag,
                tissue.tag
            )
        ]
        tmp.data.2[
            abs.freq<stat.cell.thold,
            `:=`(mean=NA, fract=NA)
        ]
        tmp.data.2[, abs.freq:=NULL]
        tmp.vals <- c('mean', 'fract')
        tmp.data.2 <- lapply(X=tmp.vals, FUN=function(this.stat){
            tmp.cols <- c(merge.vals, 'tissue.tag', this.stat)
            tmp.data <- tmp.data.2[, ..tmp.cols]
            tmp.data <- as.data.table(spread(data=tmp.data, key=tissue.tag, value=this.stat))
            colnames(tmp.data)[!colnames(tmp.data)%in%merge.vals] <- paste0(
                this.stat, '.',
                colnames(tmp.data)[!colnames(tmp.data)%in%merge.vals]
            )
            return(tmp.data)
        })
        tmp.data.2 <- Reduce(x=tmp.data.2, f=function(x, y){
            tmp.data <- merge(
                x=x, y=y,
                by=merge.vals,
                all=TRUE
            )
            return(tmp.data)
        })
        # Merge and final format.
        tmp.data <- merge(
            x=tmp.data.1,
            y=tmp.data.2,
            by=merge.vals,
            all=TRUE
        )
        colnames(tmp.data)[!colnames(tmp.data)%in%merge.vals] <- paste0(
            'exp.', this.gene, '.',
            colnames(tmp.data)[!colnames(tmp.data)%in%merge.vals]
        )
        return(tmp.data)
    })
    tmp.data <- Reduce(x=tmp.data, f=function(x, y){
        tmp.data <- merge(
            x=x, y=y,
            by=merge.vals,
            all=TRUE
        )
        return(tmp.data)
    })
    # Simplify if there's a single tissue in the dataset.
    if(single.tissue){
        tmp.cols <- colnames(tmp.data)[str_detect(string=colnames(tmp.data), patter='mean.both$|fract.both')]
        tmp.cols <- setdiff(x=colnames(tmp.data), y=tmp.cols)
        tmp.data <- tmp.data[, ..tmp.cols]
        colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='.healthy$', replacement='')
    }
    # Add info to metadata.
    clone.meta <- merge(
        x=clone.meta,
        y=tmp.data,
        by=merge.vals,
        all=TRUE
    )
}


### ------------------------- Final output -------------------------- ###

#---> Required only if extra info was added.
tmp.check <- !is.null(exp.set) | is.null(mod.sig.path)
if(tmp.check){
    # Donor metadata columns to maintain in final output.
    if(!is.null(donor.meta.cols)){
        tmp.check <- all(donor.meta.cols %in% colnames(meta.data))
        if(!tmp.check) stop('Not all requested donor meta columns were found in actual metadata.\n')
        donor.meta.cols <- unique(c('donor.id.tag', donor.meta.cols))
        tmp.data <- unique(meta.data[!is.na(donor.id.tag), ..donor.meta.cols])
        tmp.check <- tmp.data[, .N, by=donor.id.tag][N>1, .N==0]
        if(!tmp.check) stop('Some requested donor meta columns might not represent donor-specific columns. Multiple values found for them across donors:\n')
        # Proceed if cool
        tmp.cols <- c(colnames(clone.meta), setdiff(x=donor.meta.cols, y='donor.id.tag'))
        clone.meta <- merge(
            x=clone.meta, y=tmp.data,
            by='donor.id.tag',
            all.x=TRUE, all.y=FALSE,
            sort=FALSE
        )
        clone.meta <- clone.meta[, ..tmp.cols]
    }
    # Final format.
    to.order <- c('rel.freq', 'healthy.rel.freq', 'tumor.rel.freq')
    to.order <- to.order[to.order %in% colnames(clone.meta)]
    setorderv(x=clone.meta, col=to.order, order=c(-1))
    # @ For expanded clonotypes only
    to.filter <- c('abs.freq', 'healthy.abs.freq', 'tumor.abs.freq')
    to.filter <- to.filter[to.filter %in% colnames(clone.meta)]
    to.filter <- lapply(X=to.filter, FUN=function(tmp.col){
        tmp.data <- data.table(opt=clone.meta[[tmp.col]]>1)
        colnames(tmp.data) <- tmp.col
        return(tmp.data)
    })
    to.filter <- Reduce(x=to.filter, f=cbind)
    to.filter <- as.matrix(to.filter)
    to.filter <- rowSums(to.filter)>0
    tmp.data <- clone.meta[to.filter]
    tmp.file.name <- paste0(reports.path, '/ClonotypeBasedTCRData_ExtendedInfo.csv')
    fwrite(file=tmp.file.name, x=tmp.data, quote=TRUE, na=NA)
}


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('Program finished!\n\n')