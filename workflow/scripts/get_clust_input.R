############    -----  Prepare general input for  -----    ############
############    ---  TCR clustering by similarity  ----    ############

# By Vicente Fajardo-Rosas


### -------------------------- Description -------------------------- ###
# TBD


cat('\n\n')
### --------------------------- Libraries --------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
cat('Importing libraries...\n\n')
library(Seurat)
library(stringr)
library(tidyr)
library(data.table)
library(ggplot2)
library(optparse)
cat('Libraries imported!\n')


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
    # Input- and output-related paths and files.
    make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Char, absolute path to directory to save DE genes table.\n"),
    make_option(opt_str="--RunID", type="character", default=NULL, dest="run.id", help="ID for the comprehensive analysis to be performed.\n"),
    make_option(opt_str="--InputTable", type="character", default=NULL, dest="input.table.file", help="Char, absolute path to input table listing, for a set of datasets, the following items (with the specified names within brackets as columns): 1) A dataset label ('dataset.lab'); 2) absolute path to a seurat object file to retrieve clonotype-related annotations from ('seurat.obj'); 3) absolute path to a metadata file to retrieve clonotype-related annotations from ('meta.data'); 4) absolute path to 'clonotypes.csv' file as output by cellranger vdj or a related tool ('clonotypes'); 5) absolute path to 'filtered_contig_annotations.csv' file as output by cellranger vdj or a related tool ('filtered_contig_annotations'). You might either at least one or both items between a seurat object and a metadata file to retrieve clonotype annotations from; however, if both are provided, priority will be given to the seurat object.\n"),
    # Filtering-related parameters for the tool input.
    make_option(opt_str="--TagsOfInt", type="character", default=NULL, dest="tags.of.int", help="Char, a set of tags that should be carried during the whole analysis. The tags provided must be properly defined in the metadata of all datasets. The format required is a string with the names of the tags separated by semicolons; for example, if we are interested in both tags 'virus.tag' and 'pr.tag', we must provide the single string 'virus.tag;pr.tag'.\n"),
    make_option(opt_str="--FiltCriteria", type="character", default=NULL, dest="filt.criteria", help="Char, take the following as an example: 'pr.tag:pR,non-pR;virus.tag:SARS-CoV-2,FLU'.\n"),
    make_option(opt_str="--ExpansionThold", type="integer", default=1, dest="expansion.thold", help="Int, expansion threshold to set for clonality. Clonotypes with a expansion lesser than this threshold will be excluded from the analysis.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
reports.path <- opt$reports.path
run.id <- opt$run.id
input.table.file <- opt$input.table.file
tags.of.int <- opt$tags.of.int
if(!is.null(tags.of.int)) tags.of.int <- str_split(string=tags.of.int, pattern=';', simplify=TRUE)[1, ]
if('donor.id.tag' %in% tags.of.int) tags.of.int <- setdiff(x=tags.of.int, y='donor.id.tag') # Sanity check. The donor identity will be mandatorily conserved during the process, so no need to carry it as a tag of interest.
filt.criteria <- opt$filt.criteria
# filt.criteria <- 'pr.tag:pR,non-pR;virus.tag:SARS-CoV-2,FLU' # TAKE THIS AS AN EXAMPLE.
expansion.thold <- opt$expansion.thold

# ---> Define input files.
# Define files from general input table.
input.table <- fread(file=input.table.file, quote='\'')
dataset.labs <- input.table$dataset.lab
objs.files.list <- if('seurat.obj' %in% colnames(input.table)) input.table$seurat.obj else rep(x=NA, times=nrow(input.table))
names(objs.files.list) <- dataset.labs
meta.files.list <- if('meta.data' %in% colnames(input.table)) input.table$meta.data else rep(x=NA, times=nrow(input.table))
names(meta.files.list) <- dataset.labs
clons.info.list <- if('clonotypes' %in% colnames(input.table)) input.table$clonotypes else rep(x=NA, times=nrow(input.table))
names(clons.info.list) <- dataset.labs
cells.clons.info.list <- if('filtered_contig_annotations' %in% colnames(input.table)) input.table$filtered_contig_annotations else rep(x=NA, times=nrow(input.table))
names(cells.clons.info.list) <- dataset.labs
if('expansion.thold' %in% colnames(input.table)){
    expansion.tholds <- input.table$expansion.thold
}else{
    expansion.tholds <- unlist(lapply(X=input.table$dataset.lab, FUN=function(x) return(expansion.thold)))
}
names(expansion.tholds) <- dataset.labs
process.filt.criteria <- function(this.criteria){
    if(is.na(this.criteria)) return(NA)
    this.criteria <- str_split(string=this.criteria, pattern=';', simplify=TRUE)[1, ]
    this.criteria <- lapply(X=this.criteria, FUN=function(x){
        tmp.data <- str_split(string=x, pattern=':', simplify=TRUE)[1, ]
        if(length(tmp.data)!=2) stop('Unexpected error while parsing filtering criteria.\n')
        to.return <- str_split(string=tmp.data[2], pattern=',', simplify=FALSE)
        names(to.return) <- tmp.data[1]
        return(to.return)
    })
    this.criteria <- Reduce(x=this.criteria, f=c)
    return(this.criteria)
}
if('filt.criteria' %in% colnames(input.table)){
    filt.criteria <- input.table$filt.criteria
    filt.criteria <- lapply(X=filt.criteria, FUN=process.filt.criteria)
}else{
    if(!is.null(filt.criteria)){
        filt.criteria <- process.filt.criteria(this.criteria=filt.criteria)
        filt.criteria <- lapply(X=input.table$dataset.lab, FUN=function(x) return(filt.criteria))
    }else{
        filt.criteria <- lapply(X=input.table$dataset.lab, FUN=function(x) return(NA))
    }
}
names(filt.criteria) <- dataset.labs
# Check files.
essential.files <- c(
  unlist(objs.files.list),
  unlist(meta.files.list),
  unlist(cells.clons.info.list),
  unlist(clons.info.list)
)
essential.files <- essential.files[!is.na(essential.files)]
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### --------------------------- Functions --------------------------- ###


cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# ---> Gene expression info.
srt.objs.list <- lapply(X=objs.files.list, FUN=function(x) if(!is.na(x)) readRDS(file=x) else NA)
# ---> Metadata files.
meta.objs.list <- lapply(X=meta.files.list, FUN=function(x) if(!is.na(x)) fread(file=x) else NA)
# ---> TCR data files
# Clonotypes info
clons.info.list <- lapply(X=clons.info.list, FUN=function(x) if(!is.na(x)) fread(file=x) else NA)
names(clons.info.list) <- dataset.labs
# Cells-clonotypes relationships info.
cells.clons.info.list <- lapply(X=cells.clons.info.list, FUN=function(x) if(!is.na(x)) fread(file=x) else NA)
names(cells.clons.info.list) <- c(dataset.labs)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')


cat('\n\n')
### ---------------------- Data Preprocessing ---------------------- ###
cat('### ---------------------- Data Preprocessing ---------------------- ###\n')

# ---> Cells-clonotypes relationships info (filtering).
# Filter out non-productive or not-of-interest contigs from the cells-clonotypes info.
# Conditions:
# * For barcodes called as cells.
# * For clonotype contigs marked with high confidence.
# * For chains ultimately defined as TRA or TRB.
# * For productive clonotypes (see 10X support website for a broader definition of 'productive').
# Also, we'll take out this info since it has been taken into consideration already.
cells.clons.info.list <- lapply(X=cells.clons.info.list, FUN=function(x){
    if(is.na(x)) return(NA)
    # Filter for good quality cells and alpha/beta TCRs.
    these.chains <- c('TRA', 'TRB')
    x <- x[is_cell & high_confidence & productive & chain%in%these.chains]
    # Keep relevant column information for downstream processes.
    cols.to.keep <- c(
        'barcode', 'chain', 'v_gene', 'j_gene', 'cdr3_nt', 'cdr3', 'umis', 'raw_clonotype_id'
    )
    x <- x[, ..cols.to.keep]
    x <- unique(x)
    return(x)
})


cat('\n\n')
### ----------------- Create individual tool input ----------------- ###
cat('### ----------------- Create individual tool input ----------------- ###\n')

# ---> Assess clonal expansion distribution across unique clonotypes for the datasets of interest in order to pick cutoffs.
tmp.file.name <- paste0(reports.path, '/CloneSizeDensitiesPerSet.pdf')
pdf(file=tmp.file.name)
# Get plot per dataset.
for(data.set in dataset.labs){
  if(all(is.na(srt.objs.list[[data.set]]))){
    meta.data <- meta.objs.list[[data.set]]
  }else{
    meta.data <- as.data.table(srt.objs.list[[data.set]]@meta.data)
  }
  tmp.data <- meta.data[, .(clone.size=.N), by=clonotype.tag]
  tmp.thold <- tmp.data[, quantile(x=clone.size, probs=0.99)]
  tmp.data[clone.size>=tmp.thold, clone.size:=tmp.thold]
  tmp.caption <- paste0('Cutoff was set at the 99% quantile which is ', tmp.thold, '.\nThat is, all values greater than this were set to this cutiff for better viz.\n')
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=clone.size)) +
    geom_histogram(bins=tmp.thold, fill='indianred', colour='black', alpha=0.6) +
    geom_vline(xintercept=5, linetype='dashed') +
    scale_x_continuous(expand=c(0, 0)) + 
    scale_y_continuous(expand=c(0, 0)) +
    labs(title=data.set, x='Clone size', y='Frequency', caption=tmp.caption) +
    theme_bw()
  print(tmp.ggplot)
}
dev.off()

# ---> Program as a function.
get.input.1 <- function(
  # Gene expression-related information.
  seurat.obj=NULL, meta.data=NULL,
  # Clonotype-related information.
  cells.clons.info, clons.info,
  # Tags of interest.
  tags.of.int=NULL,
  # Filtering criteria
  filt.criteria=NULL, keep=TRUE, expansion.thold=1
){
    # ---> Determine if inference of CDR3 and v and j gene sequences is needed.
    mandatory.cols <- unique(c("clonotype.tag", "donor.id.tag", "cdr3a.aa.seq", "cdr3b.aa.seq", "cdr3a.nt.seq", "cdr3b.nt.seq", "tra.v", "tra.j", "trb.v", "trb.j", tags.of.int, names(filt.criteria)))
    tmp.check <- all(mandatory.cols %in% colnames(meta.data))
    if(!tmp.check){
        ### ----------- Inference of relevant sequences is needed ---------- ###
        # ---> Fill cell-clonotypes relationships with tags of interest info.
        # Take the overall metadata if necessary
        if(!is.null(seurat.obj)){
            meta.data <- seurat.obj@meta.data
            meta.data$barcode <- rownames(meta.data)
            meta.data <- as.data.table(meta.data)
        }
        # Then, combine clones data with other meta data and update seurat object.
        tmp.data <- cells.clons.info[, .(barcode, raw_clonotype_id)]
        tmp.data <- merge(x=tmp.data, y=meta.data, by='barcode', all.x=FALSE, all.y=TRUE)
        # Check previous clonotype annotations match the ones here. Else, break loop.
        tmp.check <- all(tmp.data$clonotype.tag==tmp.data$raw_clonotype_id, na.rm=TRUE)
        if(!tmp.check) stop(paste0('Unexpected base error. Previous clonotype tags do not match the ones obtained based on the vdj data provided.\n'))
        gex.clons.info <- merge(x=tmp.data, y=clons.info[, c('clonotype_id', 'cdr3s_aa', 'cdr3s_nt')], by.x='clonotype.tag', by.y='clonotype_id', all=FALSE)
        gex.clons.info <- as.data.table(gex.clons.info)

        # ---> vdj gene usage info.
        # General vdj gene usage info.
        # We obtain unique v and j genes per clonotype. When cellranger has defined multiple v and j genes for a given clonotype, we take the one that is most supported. The most supported is the gene that has the largest amount of entries (contigs) in the table below.
        cells.clons.info <- cells.clons.info[barcode %in% meta.data[, barcode]]
        vdj.gene.info <- cells.clons.info[,
            .(
                v.gene=.SD[,
                .(count=.N),
                by=v_gene
                ][count==max(count), paste0(unique(v_gene), collapse='|')],
                j.gene=.SD[,
                .(count=.N),
                by=j_gene
                ][count==max(count), paste0(unique(j_gene), collapse='|')]
            ),
            by=.(
                raw_clonotype_id,
                chain
            )
        ]
        # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
        vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
        vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
        # Spread values according to clonotype ID (i.e., to disregard chain information)
        vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, sep=',')]
        vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
        # Separate values according to gene type for each chain.
        vdj.gene.info <- separate(data=vdj.gene.info, col=TRA, into=c('tra.v', 'tra.j'), sep=',', convert=FALSE)
        vdj.gene.info <- separate(data=vdj.gene.info, col=TRB, into=c('trb.v', 'trb.j'), sep=',', convert=FALSE)
        vdj.gene.info <- as.data.table(vdj.gene.info)
        # Sanity check. Confirm we have unique entries for the amount of unique clonotypes.
        tmp.check <- vdj.gene.info[, uniqueN(raw_clonotype_id)==.N]
        if(!tmp.check) stop('Unexpected error post vdj data retrieval.\n')

        # ---> Info per clonotype (i.e., disregarding GEx and extended information).
        # We must keep only the clonotypes that meet certain criteria (as described below) at the time that we keep track of the proportions clonotypes that are kept after each filtering step.
        # Sanity check. All tags from filtering criteria and tags of interest should be properly defined in the clone information at this point.
        if(!is.null(filt.criteria)){
            tags.of.int <- if(!is.null(tags.of.int)) unique(c(tags.of.int, names(filt.criteria))) else names(filt.criteria)
        }
        tmp.check <- all(tags.of.int %in% colnames(gex.clons.info))
        if(!tmp.check) stop('Some of the tags of interest or tags to filter based on are not defined in the gene expression metadata.\n')
        # Preprocess data.
        if(!is.null(tags.of.int)){
            clons.info <- lapply(X=tags.of.int, FUN=function(tmp.tag){
                clons.info <- gex.clons.info[
                    !is.na(donor.id.tag),
                    .(
                        cdr3a.aa.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_aa), pattern='TRA:\\w+[^;\\n]', simplify=TRUE), pattern='TRA:|;', replacement=''), collapse=';'),
                        cdr3b.aa.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_aa), pattern='TRB:\\w+[^;\\n]', simplify=TRUE), pattern='TRB:|;', replacement=''), collapse=';'),
                        cdr3a.nt.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_nt), pattern='TRA:\\w+[^;\\n]', simplify=TRUE), pattern='TRA:|;', replacement=''), collapse=';'),
                        cdr3b.nt.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_nt), pattern='TRB:\\w+[^;\\n]', simplify=TRUE), pattern='TRB:|;', replacement=''), collapse=';'),
                        size=.N,
                        tmp.tag=paste0(sort(unique(get(tmp.tag))), collapse=';')
                    ),
                    by=.(
                        clonotype.tag,
                        donor.id.tag
                    )
                ]
                colnames(clons.info)[colnames(clons.info)=='tmp.tag'] <- tmp.tag
                return(clons.info)
            })
            clons.info <- Reduce(x=clons.info, f=function(x, y){
                merge(x=x, y=y, by=c('clonotype.tag', 'cdr3a.aa.seq', 'cdr3b.aa.seq', 'cdr3a.nt.seq', 'cdr3b.nt.seq', 'size', 'donor.id.tag'))
            })
        }else{
            clons.info <- gex.clons.info[
                !is.na(donor.id.tag),
                .(
                    cdr3a.aa.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_aa), pattern='TRA:\\w+[^;\\n]', simplify=TRUE), pattern='TRA:|;', replacement=''), collapse=';'),
                    cdr3b.aa.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_aa), pattern='TRB:\\w+[^;\\n]', simplify=TRUE), pattern='TRB:|;', replacement=''), collapse=';'),
                    cdr3a.nt.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_nt), pattern='TRA:\\w+[^;\\n]', simplify=TRUE), pattern='TRA:|;', replacement=''), collapse=';'),
                    cdr3b.nt.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_nt), pattern='TRB:\\w+[^;\\n]', simplify=TRUE), pattern='TRB:|;', replacement=''), collapse=';'),
                    size=.N
                ),
                by=.(
                    clonotype.tag,
                    donor.id.tag
                )
            ]
        }
        # NOTE: Here, we won't have unique entries per clonotype since it might be expressed by multiple donors and, hence, have different values for extra tags across donors.
        # Final formatting.
        clons.info[cdr3a.aa.seq=='', cdr3a.aa.seq:=NA]
        clons.info[cdr3b.aa.seq=='', cdr3b.aa.seq:=NA]
        clons.info[cdr3a.nt.seq=='', cdr3a.nt.seq:=NA]
        clons.info[cdr3b.nt.seq=='', cdr3b.nt.seq:=NA]
        clons.info <- merge(x=clons.info, y=vdj.gene.info, by.x='clonotype.tag', by.y='raw_clonotype_id', all.x=TRUE, all.y=FALSE)
    }else{
        ### --------- Inference of relevant sequences is NOT needed -------- ###
        clons.info <- unique(meta.data[!is.na(donor.id.tag), ..mandatory.cols])
        # Determine clone size on a cellranger's clonotype ID and donor ID basis as it would be calculated when inferring sequences. Done for consistency.
        tmp.data <- meta.data[,
            .(size=.N),
            by=.(clonotype.tag, donor.id.tag)
        ]
        clons.info <- merge(x=clons.info, y=tmp.data, by=c('clonotype.tag', 'donor.id.tag'))
        # Finally, rename metadata for consistency.
        gex.clons.info <- meta.data
    }

    # ---> Further cleaning of TCR data.
    # To keep track of proportions.
    filt.track <- clons.info[, uniqueN(clonotype.tag)]
    #   Filtering rule #1: Clone size must be equal to or greater than the provided threshold.
    clons.info <- clons.info[size>=expansion.thold]
    filt.track[2] <- clons.info[, uniqueN(clonotype.tag)]
    #   Filtering rule #2: Remove ambiguously annotated clonotype annotations for donors where they're not expanded.
    tmp.data.1 <- clons.info[,
        .(donor.count=uniqueN(donor.id.tag)),
        by=clonotype.tag
    ]
    tmp.data.1 <- tmp.data.1[donor.count>1, clonotype.tag] # Potentially public TCRs.
    tmp.data.2 <- gex.clons.info[!is.na(donor.id.tag) & clonotype.tag %in% tmp.data.1, .N, by=.(donor.id.tag, clonotype.tag)] # Donors w/ expanded clonotype from ambiguously assigned clonotypes.
    tmp.data.2 <- tmp.data.2[N>=1, .(clonotype.tag, donor.id.tag)] # For ambiguous clonotypes, keep track of their presence only in donors where it is expanded.
    tmp.data.3 <- gex.clons.info[!is.na(donor.id.tag) & !clonotype.tag %in% tmp.data.1, .(clonotype.tag, donor.id.tag)] # Keep track of all non-ambiguous clonotypes.
    tmp.data.1 <- unique(rbind(tmp.data.2, tmp.data.3))
    clons.info <- merge(x=clons.info, y=tmp.data.1, by=c('clonotype.tag', 'donor.id.tag'))
    filt.track[3] <- clons.info[, uniqueN(clonotype.tag)] # Keep track of proportions.
    #   Filtering rule #3: Clonotypes must have at least a single properly defined beta sequence.
    clons.info <- clons.info[!is.na(cdr3b.aa.seq)]
    filt.track[4] <- clons.info[, uniqueN(clonotype.tag)] # Keep track of proportions.
    #   Filtering rule #4: Clonotypes must have defined beta chain V and J genes.
    clons.info <- clons.info[!is.na(trb.v) & !is.na(trb.j)]
    filt.track[5] <- clons.info[, uniqueN(clonotype.tag)] # Keep track of proportions.
    #   Filtering rule #5: Clonotypes must have **properly defined** V and J genes.
    clons.info <- clons.info[
        str_detect(string=trb.v, pattern='^TRBV') &
        str_detect(string=trb.j, pattern='^TRBJ') &
        (is.na(tra.v) | str_detect(string=tra.v, pattern='^TRAV')) &
        (is.na(tra.j) | str_detect(string=tra.j, pattern='^TRAJ'))
    ]
    filt.track[6] <- clons.info[, uniqueN(clonotype.tag)] # Keep track of proportions.
    #   Filtering rule #6: CDR3 beta sequences must be properly defined.
    # Make sure that all aa sequences are written w/ upper case letters.
    clons.info[, cdr3b.aa.seq:=str_to_upper(string=cdr3b.aa.seq)]
    # Were only valid aminoacids one-letter abbreviations used to define CDR3 sequences?
    valid.aas <- c(
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
    )
    tmp.pttn <- paste0('[^', paste0(valid.aas, collapse=''), ']')
    clons.info <- clons.info[!str_detect(string=cdr3b.aa.seq, pattern=tmp.pttn)]
    filt.track[7] <- clons.info[, uniqueN(clonotype.tag)] # Keep track of proportions.
    #   Filtering rule #7: Custom filtering criteria.
    if(!is.null(filt.criteria)){
        for(criteria.term in names(filt.criteria)){
            clons.info[, filt.criteria.term:=get(criteria.term)]
            filt.vals <- filt.criteria[[criteria.term]]
            clons.info <- if(keep) clons.info[filt.criteria.term %in% filt.vals] else clons.info[!filt.criteria.term %in% filt.vals]
            clons.info[, filt.criteria.term:=NULL]
        }
    }
    filt.track[8] <- clons.info[, uniqueN(clonotype.tag)]

    # ---> Pack and go.
    names(filt.track) <- c('original', 'expansion.thold', 'ambiguous', 'beta.seq', 'beta.vj.genes', 'valid.vj.genes', 'valid.cdr3.aa', 'filt.criteria')
    to.return <- list(
        clons.info,
        filt.track
    )
    return(to.return)
}

# ---> Run for every dataset.
final.tool.list <- lapply(X=dataset.labs, FUN=function(data.set){
    cat(data.set, '\n')
    # Define input.
    seurat.obj <- srt.objs.list[[data.set]]
    seurat.obj <- if(!is.na(seurat.obj)) seurat.obj else NULL
    meta.data <- if(is.null(seurat.obj)) meta.objs.list[[data.set]] else NULL
    expansion.thold <- expansion.tholds[data.set]
    this.criteria <- if(all(is.na(filt.criteria[[data.set]]))) NULL else filt.criteria[[data.set]]
    # Run.
    to.return <- get.input.1(
        seurat.obj=seurat.obj, meta.data=meta.data,
        cells.clons.info=cells.clons.info.list[[data.set]], clons.info=clons.info.list[[data.set]],
        tags.of.int=tags.of.int, filt.criteria=this.criteria, keep=TRUE, expansion.thold=expansion.thold
    )
    return(to.return)
})
# Prepare report.
tmp.report <- lapply(X=final.tool.list, FUN=function(x) return(x[[2]]))
if(length(tmp.report)==1){
    tmp.report <- as.data.frame(as.matrix(tmp.report[[1]]), stringsAsFactors=FALSE)
}else{
    tmp.report <- t(Reduce(x=tmp.report, f=rbind))
}
colnames(tmp.report) <- dataset.labs
tmp.file.name <- paste0(reports.path, '/ClonotypeSummarySetpwisely.csv')
write.csv(file=tmp.file.name, x=tmp.report, quote=FALSE)
# Keep tool input only.
final.tool.list <- lapply(X=final.tool.list, FUN=function(x) return(x[[1]]))
names(final.tool.list) <- dataset.labs


cat('\n\n')
### -------------------- Get actual tool inputs -------------------- ###
cat('### -------------------- Get actual tool inputs -------------------- ###\n\n')

# ---> TCR input, program as a function.
# @ Program as a function.
get.input.2 <- function(tool.list, reports.path, run.id=NA){
    # Get general data.
    tmp.data <- rbindlist(l=tool.list, use.names=TRUE, idcol='data.set', fill=TRUE)
    tmp.data <- as.data.table(separate_rows(data=tmp.data, cdr3b.aa.seq, sep=';'))
    # ---> Split info between chain types.
    # @ Beta TCR chains.
    beta.info <- tmp.data[!is.na(cdr3b.aa.seq)]
    beta.info[, cdr3a.aa.seq:=NULL]
    get.beta.info <- function(tmp.tag){
        beta.info[, tmp.tag:=NULL]
        beta.info[, tmp.tag:=get(tmp.tag)]
        tmp.data <- beta.info[,
            .(
            clonotype.tag=paste0(
                unique(
                sort(paste(clonotype.tag, data.set, sep='|'))
                ), collapse=';'
            ),
            size=sum(size),
            tmp.tag=.SD[,
                .(size=sum(size)), by=tmp.tag][
                size==max(size),
                paste0(sort(unique(tmp.tag)), collapse=';')
            ]
            ),
            by=.(cdr3b.aa.seq, donor.id.tag, trb.v, trb.j)
        ]
        colnames(tmp.data)[colnames(tmp.data)=='tmp.tag'] <- tmp.tag
        return(tmp.data)
    }
    if(is.null(tags.of.int)){
        tags.of.int <- 'mock'
        beta.info[, mock:='mock']
    }
    beta.info <- lapply(X=tags.of.int, FUN=get.beta.info)
    if(tags.of.int == 'mock'){
        tags.of.int <- NULL
        beta.info <- beta.info[[1]]
        beta.info[, mock:=NULL]
    }else{
        beta.info <- Reduce(x=beta.info, f=function(x, y){
            merge(x=x, y=y, by=c('clonotype.tag', 'cdr3b.aa.seq', 'trb.v', 'trb.j', 'size', 'donor.id.tag'))
        })
    }
    # ---> Clone size assessment.
    # Get plot.
    tmp.thold <- beta.info[, quantile(x=size, prob=.99)]
    tmp.data <- beta.info
    tmp.data[size>=tmp.thold, size:=tmp.thold]
    tmp.caption <- paste0('Cutoff was set at the 99% quantile which is ', tmp.thold, '.\nThat is, all values greater than this were set to this cutiff for better viz.\n')
    tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=size)) +
        geom_histogram(bins=tmp.thold, fill='indianred', colour='black', alpha=0.6) +
        geom_vline(xintercept=5, linetype='dashed') +
        scale_x_continuous(expand=c(0, 0)) + 
        scale_y_continuous(expand=c(0, 0)) +
        labs(title='Combined datasets', x='Clone size', y='Frequency', caption=tmp.caption) +
        theme_bw()
    # Output plot.
    tmp.file.name <- paste0(reports.path, '/CloneSizeDensitiesAcrossSets.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.1)
    dev.off()
    # Pick expansion threshold and indicate it below.
    # ---> Actual inputs.
    # @ Keep only clonotypes with clone size over certain value.
    tmp.data <- beta.info[, .(cdr3b.aa.seq, trb.v, trb.j, cdr3a.aa.seq=NA, donor.id=paste(donor.id.tag, 'Exp1', sep=':'), size)]
    # Add mock gene allele.
    tmp.data[, trb.v:=paste0(trb.v, '*01')]
    tmp.data[, trb.j:=paste0(trb.j, '*01')]
    # Output.
    if(!is.na(run.id)){
        tmp.file.name <- paste0(reports.path, '/ToolInput.tsv')
        fwrite(file=tmp.file.name, x=tmp.data, sep='\t', na=NA, quote=FALSE)
    }
    # Return results
    return(beta.info)
}

# @ Obtain input.
final.beta.info <- get.input.2(
    tool.list=final.tool.list,
    reports.path=reports.path, run.id=run.id
)
# Save results for further use in general clustering pipeline.
tmp.file.name <- paste0(reports.path, '/GeneralBetaInfo.csv')
fwrite(file=tmp.file.name, x=final.beta.info, na=NA)

cat('Data has been preprocessed!\n\n')


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')