############    -----   Data modality intersection    -----    ############
cat('\n\n')
cat('############    -----   Data modality intersection    -----    ############\n')

# By: Vicente Fajardo-Rosas

### -------------------------- Description -------------------------- ###

cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
library(data.table)
library(tidyr)
library(stringr)
library(gtools)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(optparse)
library(Hmisc)
library(corrplot)


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
    make_option(opt_str="--RefDefFile", type="character", default=NULL, dest="input.table.file", help="Char, indicates the reference ID for the results to be evaluated.\n"),
    make_option(opt_str="--RunPath", type="character", default=NULL, dest="run.path", help="Char, absolute path to directory where clustering results have been saved to.\n"),
    make_option(opt_str="--OptsFile1", type="character", default=NULL, dest="yaml.file.1", help="Char, absolute path to yaml file for workflow; general options.\n"),
    make_option(opt_str="--OptsFile2", type="character", default=NULL, dest="yaml.file.2", help="Char, absolute path to yaml file for workflow; project-specific options.\n"),
    make_option(opt_str="--RefID", type="character", default=NULL, dest="ref.name", help="Char, indicates the reference ID for the results to be evaluated.\n"),
    make_option(opt_str="--CloneInfo", type="character", default=NULL, dest="clone.info.file", help="Char, path to clonotype-based TCR info.\n"),
    make_option(opt_str="--HLAInfo", type="character", default=NULL, dest="qry.hla.info.file", help="Character, absolute path to HLA typing data for query dataset. If not available, drop the flag.\n"),
    make_option(opt_str="--BgSize", type="integer", default=10000, dest="cell.bg.size", help="Integer, number of cells left without an antigen specificity prediction to take as the background set.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
input.table.file <- opt$input.table.file
run.path <- opt$run.path
yaml.file.1 <- opt$yaml.file.1
yaml.file.2 <- opt$yaml.file.2
ref.name <- opt$ref.name
clone.info.file <- opt$clone.info.file
qry.hla.info.file <- opt$qry.hla.info.file
cell.bg.size <- opt$cell.bg.size
# ---> Define constant vars.
set.seed(seed=1)

# ---> Options file.
# Define files from general input table.
input.table <- fread(file=input.table.file, quote='\'')
tcr.meta.file <- input.table[dataset.lab!=ref.name, meta.data]
cell.clons.info.file <- if('filtered_contig_annotations' %in% colnames(input.table)) input.table[dataset.lab!=ref.name, filtered_contig_annotations] else NA
# Load yaml files. Priority is given to dataset-specific options while hanlding potentially redundant entries.
if(!file.exists(yaml.file.1)) stop(paste0('Attempted definition of yaml file (general options) failed. Next file could not be found:\n', yaml.file.1, '\n'))
if(!file.exists(yaml.file.2)) stop(paste0('Attempted definition of yaml file (project-specific options) failed. Next file could not be found:\n', yaml.file.2, '\n'))
tmp.data.1 <- yaml::read_yaml(file=yaml.file.1)
tmp.data.2 <- yaml::read_yaml(file=yaml.file.2)
gen.wflow.opts <- c(tmp.data.2, tmp.data.1[!names(tmp.data.1) %in% names(tmp.data.2)])
# Define files and general options.
preds.file <- paste0(run.path, '/PredictionDetails.csv')
tool.names <- names(gen.wflow.opts$tool_opts)
gex.meta <- gen.wflow.opts$gex_data[[ref.name]]
ref.hla.info.file <- gen.wflow.opts$ref_specs[[ref.name]]$hla_info
ref.donor.id.tag <- gen.wflow.opts$ref_specs[[ref.name]]$donor_id_tag
ag.groups <- gen.wflow.opts$ref_specs[[ref.name]]$ag_groups
seurat.obj.file <- gex.meta$seurat_obj
donor.id.tag <- gex.meta$donor_id_tag
clusters.tag <- gex.meta$clusters_tag
cluster.defs <- unlist(gex.meta$cluster_defs)
limit.of.detect <- gen.wflow.opts$limit_of_detect
if(is.null(limit.of.detect)) limit.of.detect <- 10
donor.meta.file <- gen.wflow.opts$donor_meta
# Define antigen groups to evaluate (if any).
if(!is.null(ag.groups)){
    ag.groups <- lapply(X=ag.groups, FUN=function(x) str_split(string=x, pattern=',')[[1]])
}
functions.file <- gen.wflow.opts$paths$R_functions
if(!file.exists(functions.file)) stop('Faulty R functions file.\n')

# ---> General aesthetic definitions.
# @ General distinction between original and bounded values.
tmp.shapes <- c('ori'=1, `bounded`=16)
# @ General dot color.
gen.dot.col <- '#929292'
# @ General dot link color and width.
gen.link.col <- '#BFBFBF'
gen.link.width <- 0.7

# ---> Color definitions.
# @ For clusters.
cluster.cols <- unlist(gen.wflow.opts$gex_data[[ref.name]][['cluster_cols']])
# @ For populations.
pop.cols <- unlist(gen.wflow.opts$gex_data[[ref.name]][['pop_cols']])
# @ For antigen specificities.
# Default attempted definition.
def.ag.spc.cols <- c(
    `CMV`='#42a2cc',
    `EBV`='#005746',
    `IAV`='#00947F',
    `FLU`='#00947F',
    `SARS-CoV-2`='#FE0002',
    `MPV`='#a76400',
    `RSPV`='#EEA101',
    `RSV`='#a5492a',    
    `HAV`='#e69500',
    `HBV`='#cc8400',
    `HCV`='#b37400',
    `YFV`='#C1779B',
    `SARS-CoV`='#A24B4A',
    `HIV`='#720000',
    `B. pertussis Vax`='#800080',
    `B. pertussis Rest`='#cd00cd',
    `Aspergillus`='#ffc0cb',
    `Alternaria`='#e0c5ca',
    `Homo sapiens`='#779BC1',
    `Cross-reactive`='#cccc00',
    `Multiple`='#808080'
)
if(!is.null(gen.wflow.opts$ag_spc_cols)){
    ag.spc.cols <- unlist(gen.wflow.opts$ag_spc_cols)
}else{
    ag.spc.cols <- def.ag.spc.cols
}


### --------------------------- Functions --------------------------- ###
source(functions.file)


### --------------------------- Load data --------------------------- ###
### ----------------------- and preprocessing ----------------------- ###

# ---> Prediction details.
pred.summ <- fread(file=preds.file)
pred.summ[, donor.id.tag:=as.character(donor.id.tag)]

# ---> Query GEx meta.
seurat.obj <- readRDS(file=seurat.obj.file)
cell.emb.cols <- paste0('UMAP_', 1:2)
if(any(cell.emb.cols %in% colnames(seurat.obj@meta.data))){
    tmp.cols <- setdiff(x=colnames(seurat.obj@meta.data), cell.emb.cols)
    seurat.obj@meta.data <- seurat.obj@meta.data[, tmp.cols]
}
gex.meta <- as.data.table(cbind(
    seurat.obj@meta.data,
    seurat.obj@reductions$umap@cell.embeddings
))
gex.meta[, barcode:=Seurat::Cells(seurat.obj)]
rm(seurat.obj)
if(!'donor.id.tag' %in% colnames(gex.meta)){
    gex.meta[, donor.id.tag:=get(donor.id.tag)]
}else{
    if(donor.id.tag!='donor.id.tag'){
        gex.meta[, donor.id.tag:=NULL]
        gex.meta[, donor.id.tag:=get(donor.id.tag)]
    }
}
gex.meta[, donor.id.tag:=as.character(donor.id.tag)]
# if(!'cluster.tag' %in% colnames(gex.meta)) gex.meta[, cluster.tag:=as.character(get(clusters.tag))]
gex.meta[, cluster.tag:=as.character(get(clusters.tag))]
if(!is.null(cluster.defs)){
    # At least 90% of the cells should have a defined identity for the labels to be useful.
    tmp.check <- gex.meta[, .SD[cluster.tag %in% names(cluster.defs), .N]/.N] > 0.9
    if(tmp.check){
        gex.meta[, pop.tag:=cluster.defs[cluster.tag]]
        gex.meta[is.na(pop.tag), pop.tag:='Undefined']
    }else{
        gex.meta[, pop.tag:=cluster.tag]
    }
}else{
    gex.meta[, pop.tag:=cluster.tag]
}
# Check that provided colors (if any) are meaningful.
tmp.check <- all(gex.meta[, cluster.tag] %in% names(cluster.cols))
if(!tmp.check){
    tmp.vals <- gex.meta[, sort(unique(cluster.tag))]
    cluster.cols <- scales::hue_pal()(length(tmp.vals))
    names(cluster.cols) <- tmp.vals
}
if(!'Undefined' %in% names(pop.cols)) pop.cols <- c(pop.cols, `Undefined`='#a6a6a6')
tmp.check <- all(gex.meta[, pop.tag] %in% names(pop.cols))
if(!tmp.check){
    tmp.vals <- gex.meta[, sort(unique(pop.tag))]
    pop.cols <- scales::hue_pal()(length(tmp.vals))
    names(pop.cols) <- tmp.vals
}

# ---> Query TCR metadata
# @ Cell-based info.
if(!is.na(tcr.meta.file)){
    # When a processed file has already been provided.
    tcr.meta <- fread(file=tcr.meta.file)
    tmp.cols <- c('barcode', 'clonotype.tag', 'cdr3b.aa.seq', 'cdr3a.aa.seq', 'trb.v', 'trb.j','tra.v', 'tra.j')
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
# @ Clonotype-based info.
clone.info <- fread(file=clone.info.file)
clone.info[, donor.id.tag:=as.character(donor.id.tag)]

# ---> HLA info.
# For query data.
qry.hla.info <- if(!is.null(qry.hla.info.file)) fread(file=qry.hla.info.file) else NULL
if(!is.null(qry.hla.info)){
    tmp.cols <- c('donor.id.tag', 'hla.gene.allele')
    tmp.check <- all(tmp.cols %in% colnames(qry.hla.info))
    if(!tmp.check) stop('HLA type info file provided for query dataset does not follow the required format. It must contain columns "donor.id.tag" and "hla.gene.allele".')
    qry.hla.info <- qry.hla.info[, ..tmp.cols]
    tmp.donors <- qry.hla.info[, sort(unique(donor.id.tag))]
    qry.hla.info <- lapply(X=tmp.donors, FUN=function(x) return(qry.hla.info[donor.id.tag==x, unique(hla.gene.allele)]))
    names(qry.hla.info) <- tmp.donors
}
# For reference data.
if(!is.null(ref.hla.info.file)){
    if(file.exists(ref.hla.info.file)){
        ref.hla.info <- fread(file=ref.hla.info.file)
        tmp.cols <- c('donor.id.tag', 'hla.gene.allele')
        tmp.check <- all(tmp.cols %in% colnames(ref.hla.info))
        if(!tmp.check) stop('HLA type info file provided for query dataset does not follow the required format. It must contain columns "donor.id.tag" and "hla.gene.allele".')
        ref.hla.info <- ref.hla.info[, ..tmp.cols]
        tmp.donors <- ref.hla.info[, sort(unique(donor.id.tag))]
        ref.hla.info <- lapply(X=tmp.donors, FUN=function(x) return(ref.hla.info[donor.id.tag==x, unique(hla.gene.allele)]))
        names(ref.hla.info) <- tmp.donors
    }else{
        ref.hla.info <- ref.hla.info.file
    }
}else{
    ref.hla.info <- NULL
}

# ---> Donor metadata
donor.meta <- if(!is.null(donor.meta.file)) fread(file=donor.meta.file) else NULL
if(!is.null(donor.meta)){
    tmp.check.1 <- 'donor.id.tag' %in% colnames(donor.meta)
    donor.meta[, donor.id.tag:=as.character(donor.id.tag)]
    tmp.check.2 <- gex.meta[!is.na(donor.id.tag), all(donor.id.tag %in% donor.meta[, donor.id.tag])]
    tmp.check <- tmp.check.1 & tmp.check.2
    if(!tmp.check) stop('Faulty donor metadata file. Must contain column \'donor.id.tag\' and info for all donors described in GEx input must be provided.\n')
}


### ------------------------- Main program -------------------------- ###


### --------------------------- Preflights --------------------------- ###

# @ Merge individual items into a single table with all details.
plots.data <- gex.meta[,
    .(
        barcode, donor.id.tag,
        cluster.tag, pop.tag, UMAP_1, UMAP_2
    )
]
plots.data <- merge(
    x=plots.data, y=cells.clons.info,
    by='barcode',
    all=FALSE 
)
merge.cols <- if('donor.id.tag' %in% colnames(pred.summ)) c('qry.clone.id', 'donor.id.tag') else 'qry.clone.id'
plots.data <- merge(
    x=plots.data, y=pred.summ,
    by=merge.cols,
    all.x=TRUE, all.y=FALSE # Throughout, the metadata only for cells w/ available TCR info are kept.
)
# @ Remove unnecessary material.
tmp.cols <- colnames(plots.data)[!str_detect(string=colnames(plots.data), pattern='\\|')]
plots.data <- plots.data[, ..tmp.cols]
# @ Set general factors.
plots.data$donor.id.tag <- factor(x=plots.data$donor.id.tag, levels=sort(unique(plots.data$donor.id.tag)))
tmp.lvls <- names(ag.spc.cols)[names(ag.spc.cols) %in% plots.data[, consensus.pred]]
tmp.check <- all(plots.data[!is.na(consensus.pred), consensus.pred] %in% tmp.lvls)
if(!tmp.check) stop('Unexpected error while setting levels for consensus prediction.\n')
plots.data$consensus.pred <- factor(x=as.character(plots.data$consensus.pred), levels=tmp.lvls)
tmp.lvls <- c('Donor-matched', 'Reference match', 'UCM')
tmp.check <- all(plots.data[!is.na(consensus.approach), consensus.approach] %in% tmp.lvls)
if(!tmp.check) stop('Unexpected error while setting levels for consensus approach.\n')
plots.data$consensus.approach <- factor(x=plots.data$consensus.approach, levels=tmp.lvls)
tmp.lvls <- gtools::mixedsort(unique(as.character(plots.data$cluster.tag)))
plots.data$cluster.tag <- factor(x=as.character(plots.data$cluster.tag), levels=tmp.lvls)
if(!is.null(cluster.defs) & !is.null(pop.cols)){
    plots.data$pop.tag <- factor(x=plots.data$pop.tag, levels=names(pop.cols))
}else{
    tmp.lvls <- gtools::mixedsort(unique(as.character(plots.data$pop.tag)))
    plots.data$pop.tag <- factor(x=plots.data$pop.tag, levels=tmp.lvls)
}



### ------------------------ General summaries ----------------------- ###

summs.path <- paste0(run.path, '/gen_summaries')
if(!dir.exists(summs.path)) dir.create(summs.path)

# @ General summary, across the board.
# Cell count per prediction type.
tmp.ggplot.1 <- ggplot(data=plots.data[!is.na(consensus.pred)], aes(x=consensus.approach)) +
    geom_bar(width=0.6, alpha=0.6, color='black', fill='lightblue', linewidth=0.6) +
    scale_y_continuous(expand=c(0, 0)) +
    labs(title='Cell frequency accounted for by each approach.', x='Consensus approach', y='Cell count') +
    theme_bw() + theme(axis.text.x=element_text(angle=-90))
# Number of cells per single specificity.
tmp.ggplot.2 <- ggplot(data=plots.data[!is.na(consensus.pred)], aes(x=consensus.pred)) +
    geom_bar(width=0.6, alpha=0.6, color='black', fill='lightblue', linewidth=0.6) +
    scale_y_continuous(expand=c(0, 0)) +
    labs(title='Cell frequency accounted for by each organism.', x='Consensus organism', y='Cell count') +
    theme_bw() + theme(axis.text.x=element_text(angle=-90))
# Clonotype count per prediction type.
tmp.cols <- c('qry.clone.id', 'consensus.approach', 'consensus.pred')
tmp.report <- unique(plots.data[!is.na(consensus.pred), ..tmp.cols])
tmp.ggplot.3 <- ggplot(data=tmp.report, aes(x=consensus.approach)) +
    geom_bar(width=0.6, alpha=0.6, color='black', fill='lightblue', linewidth=0.6) +
    scale_y_continuous(expand=c(0, 0)) +
    labs(title='Clonotype count per prediction type.', x='Prediction type', y='Clonotype count') +
    theme_bw() + theme(axis.text.x=element_text(angle=-90))
# Clonotype count per single reactivity.
tmp.ggplot.4 <- ggplot(data=tmp.report, aes(x=consensus.pred)) +
    geom_bar(width=0.6, alpha=0.6, color='black', fill='lightblue', linewidth=0.6) +
    scale_y_continuous(expand=c(0, 0)) +
    labs(title='Clonotype count per prediction type.', x='Prediction type', y='Clonotype count') +
    theme_bw() + theme(axis.text.x=element_text(angle=-90))
tmp.file.name <- paste0(summs.path, '/PredictionSummary_General.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
dev.off()

# @ General summary per item.
tags.of.int <- c(`Donor ID`='donor.id.tag', `Cell cluster`='cluster.tag', `Cell population`='pop.tag')
for(tag.of.int in names(tags.of.int)){
    tmp.val <- tags.of.int[tag.of.int]
    # Predictions per item.
    tmp.report <- plots.data[,
        .(
            cell.fraction=.SD[!is.na(consensus.pred), .N]/.N,
            clone.fraction=.SD[!is.na(consensus.pred), uniqueN(qry.clone.id)]/uniqueN(qry.clone.id)
        ),
        by=.(tmp.val=get(tmp.val))
    ]
    tmp.note <- paste0('Thick, horizontal line indicates median across individual values.')
    tmp.ggplot.1 <- ggplot(data=tmp.report, aes(x=tmp.val, y=clone.fraction)) +
        geom_bar(stat='identity', width=0.6, alpha=0.6, color='black', fill='lightblue', linewidth=0.6) +
        geom_hline(yintercept=tmp.report[, median(clone.fraction)], color='#0000b3', linewidth=1.2) +
        scale_y_continuous(expand=c(0, 0)) +
        labs(title='Fraction of imputed clonotypes', x=tag.of.int, y='Clonotype fraction', caption=tmp.note) +
        theme_bw() + theme(axis.text.x=element_text(angle=-90))
    tmp.ggplot.2 <- ggplot(data=tmp.report, aes(x=tmp.val, y=cell.fraction)) +
        geom_bar(stat='identity', width=0.6, alpha=0.6, color='black', fill='lightblue', linewidth=0.6) +
        geom_hline(yintercept=tmp.report[, median(cell.fraction)], color='#0000b3', linewidth=1.2) +
        scale_y_continuous(expand=c(0, 0)) +
        labs(title='Fraction of imputed cells', x=tag.of.int, y='Cell fraction', caption=tmp.note) +
        theme_bw() + theme(axis.text.x=element_text(angle=-90))
    # Association between cell and clonotype fractions for single-reactivity predictions.
    tmp.report <- as.data.table(gather(data=tmp.report, key='fraction.type', value='value', -`tmp.val`))
    tmp.vals <- c(`clone.fraction`='Clonotypes', `cell.fraction`='Cells')
    tmp.report[, fraction.type:=tmp.vals[fraction.type]]
    tmp.report$fraction.type <- factor(x=as.character(tmp.report$fraction.type), levels=tmp.vals)
    tmp.ggplot.3 <- ggplot(data=tmp.report, aes(x=fraction.type, y=value, group=tmp.val)) +
        geom_point(size=1.2, color='black') +
        geom_line(linewidth=0.3, color='black') +
        scale_y_continuous(limits=c(0, 1)) +
        labs(title='Comparison between cells\' and clonotypes\' fractions', x='', y='Fraction') +
        theme_bw() + theme(axis.text.x=element_text(angle=-90))
    tmp.file.name <- paste0(summs.path, '/PredictionSummary_', str_replace(string=tag.of.int, pattern=' ', replacement=''), '.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    print(tmp.ggplot.3)
    dev.off()
}


### ---------------- Antigen-specific TCR repertoires --------------- ###

reps.path <- paste0(run.path, '/ag_spc_tcr_reps')
if(!dir.exists(reps.path)) dir.create(reps.path)

# ---> Prevalence of pathogen-specific T cells in the lung. Donor-specific info.
# Fetch data to plot.
tmp.data.1 <- plots.data[
    !is.na(donor.id.tag) & !is.na(qry.clone.id) & !is.na(consensus.pred),
    .(cell.count=.N),
    by=.(donor.id.tag, consensus.pred)
]
tmp.data.2 <- plots.data[
    !is.na(donor.id.tag) & !is.na(qry.clone.id),
    .(total.count=.N),
    by=.(donor.id.tag)
]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag', all.x=TRUE, all.y=FALSE)
tmp.data[, cell.fract:=cell.count/total.count]
tmp.data <- tmp.data[, .(donor.id.tag, consensus.pred, cell.fract)]
tmp.thold <- tmp.data[, quantile(x=cell.fract, probs=0.95)]
if(tmp.thold<0.05) tmp.thold <- 0.05
tmp.data[cell.fract>tmp.thold, cell.fract:=tmp.thold]
tmp.data <- spread(data=tmp.data, key=consensus.pred, value=cell.fract, drop=FALSE)
row.names(tmp.data) <- tmp.data$donor.id.tag; tmp.data$donor.id.tag <- NULL
tmp.data <- as.matrix(tmp.data)
# Set metadata.
col.metadata <- data.frame(
    row.names=colnames(tmp.data),
    `Specificity`=colnames(tmp.data)
)
# Define colors for metadata tracks.
tmp.cols <- ag.spc.cols[names(ag.spc.cols) %in% colnames(tmp.data)]
ann.colors <- list(
    `Specificity`=tmp.cols
)
# Set color scale and breaks for heatmap.
break.no <- 200
col.breaks <- seq(from=min(range(tmp.data, na.rm=TRUE)), to=max(range(tmp.data, na.rm=TRUE)), length.out=break.no)
mid.point <- .005
mid.point <- which.min(abs(col.breaks-mid.point))
hmap.col.scale.1 <- colorRampPalette(c('#414487', '#2A788E'))(mid.point)
hmap.col.scale.2 <- colorRampPalette(c('#2A788E', '#22A884', '#FDE725', '#FFEA00'))(break.no-(mid.point+1))
hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
# File specs.
tmp.height <- (nrow(tmp.data) * 6.68/40) + 0.32
tmp.width.c <- (ncol(tmp.data) * 7)/10
tmp.width.b <- (ncol(tmp.data) * 1.8)/10
# @ Version with hierarchical clustering.
# To complete only if possible.
# tmp.check <- any(c(
#     any(rowSums(is.na(tmp.data))==ncol(tmp.data)),
#     any(colSums(is.na(tmp.data))==nrow(tmp.data))
# ))
# if(!tmp.check){
#     # Complete version.
#     tmp.file.name <- paste0(reps.path, '/Donor_Reactivity_CellFract_Clust.C.pdf')
#     pheatmap(
#         mat=tmp.data, scale='none',
#         color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
#         cluster_rows=TRUE, cluster_cols=TRUE,
#         annotation_col=col.metadata, annotation_colors=ann.colors,
#         show_colnames=FALSE, show_rownames=TRUE,
#         legend=TRUE, annotation_legend=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE,
#         filename=tmp.file.name, heigh=tmp.height, width=tmp.width.c
#     )
#     # Blank version.
#     tmp.file.name <- paste0(reps.path, '/Donor_Reactivity_CellFract_Clust.B.pdf')
#     pheatmap(
#         mat=tmp.data, scale='none',
#         color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
#         cluster_rows=TRUE, cluster_cols=TRUE,
#         annotation_col=col.metadata, annotation_colors=ann.colors,
#         show_colnames=FALSE, show_rownames=FALSE,
#         legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
#         filename=tmp.file.name, heigh=tmp.height, width=tmp.width.b+1
#     )
# }
# @ Version w/out hierarchical clustering.
# Complete version.
tmp.file.name <- paste0(reps.path, '/Donor_Reactivity_CellFract_NoClust.C.pdf')
pheatmap(
    mat=tmp.data, scale='none',
    color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
    cluster_rows=FALSE, cluster_cols=FALSE,
    annotation_col=col.metadata, annotation_colors=ann.colors,
    show_colnames=FALSE, show_rownames=TRUE,
    legend=TRUE, annotation_legend=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE,
    filename=tmp.file.name, heigh=tmp.height, width=tmp.width.c
)
# Blank version.
tmp.file.name <- paste0(reps.path, '/Donor_Reactivity_CellFract_NoClust.B.pdf')
pheatmap(
    mat=tmp.data, scale='none',
    color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
    cluster_rows=FALSE, cluster_cols=FALSE,
    annotation_col=col.metadata, annotation_colors=ann.colors,
    show_colnames=FALSE, show_rownames=FALSE,
    legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
    filename=tmp.file.name, heigh=tmp.height, width=tmp.width.b
)
# @ Color legend only.
tmp.data <- as.data.frame(tmp.data)
tmp.data$donor.id.tag <- row.names(tmp.data)
tmp.data <- gather(data=tmp.data, key='spc', value='freq.rel', -`donor.id.tag`)
tmp.data <- tmp.data[!is.na(tmp.data$freq.rel), ]
# W/ labels
tmp.ggplot <- ggplot(data=tmp.data, aes(x=freq.rel, y=freq.rel, col=freq.rel)) + 
    geom_point() +

    scale_color_gradientn(colors=hmap.col.scale) +
    theme(
        legend.ticks=element_line(color='black', linewidth=0.6),
        legend.ticks.length=unit(0.22, "cm"),
        legend.frame=element_rect(color='black', linewidth=0.6)
    )
tmp.ggplot <- get_legend(p=tmp.ggplot)
tmp.file.name <- paste0(reps.path, '/Donor_Reactivity_CellFract.L1.pdf')
pdf(file=tmp.file.name, height=2, width=2)
print(as_ggplot(tmp.ggplot))
dev.off()
# W/out labels
tmp.ggplot <- ggplot(data=as.data.frame(tmp.data), aes(x=freq.rel, y=freq.rel, col=freq.rel)) + 
    geom_point() +
    scale_color_gradientn(colors=hmap.col.scale, name=NULL, labels=NULL) +
    theme(
        legend.ticks=element_line(color='black', linewidth=0.6),
        legend.ticks.length=unit(0.22, "cm"),
        legend.frame=element_rect(color='black', linewidth=0.6)
    )
tmp.ggplot <- get_legend(p=tmp.ggplot)
tmp.file.name <- paste0(reps.path, '/Donor_Reactivity_CellFract.L2.pdf')
pdf(file=tmp.file.name, height=1.25, width=0.3)
print(as_ggplot(tmp.ggplot))
dev.off()

# ---> Prevalence of pathogen-specific T cells in the lung. Summary across donors
# @ Data fetch
# Fetch data to plot.
tmp.data.1 <- plots.data[
    !is.na(consensus.pred),
    .(cell.count=.N),
    by=.(donor.id.tag, consensus.pred)
]
tmp.data.2 <- plots.data[
    !is.na(qry.clone.id),
    .(total.count=.N),
    by=.(donor.id.tag)
]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag', all.x=TRUE, all.y=FALSE)
tmp.data[, cell.fract:=cell.count/total.count]
# Spread to capture a donor-wide distribution.
tmp.data <- tmp.data[, .(donor.id.tag, consensus.pred, cell.fract)]
tmp.data <- spread(data=tmp.data, key=consensus.pred, value=cell.fract, fill=0)
tmp.data <- as.data.table(gather(data=tmp.data, key=consensus.pred, value=cell.fract, -`donor.id.tag`))
# Determine basic stats.
tmp.caption <- tmp.data[, paste0(
    'Median: ', round(median(cell.fract), 5), '. ',
    'Min: ', round(min(cell.fract), 5), '. ',
    'Max: ', round(max(cell.fract), 5)
)]
# Set upper bound.
tmp.data[, bounded.cell.fract:=cell.fract]
tmp.data[, point.status:='ori']
tmp.bound <- tmp.data[, quantile(x=cell.fract, probs=0.95)]
tmp.data[cell.fract>tmp.bound, `:=`(bounded.cell.fract=tmp.bound, point.status='bounded')]
# @ Set general factors.
tmp.lvls <- names(ag.spc.cols)[names(ag.spc.cols) %in% tmp.data[, consensus.pred]]
tmp.data$consensus.pred <- factor(x=as.character(tmp.data$consensus.pred), levels=tmp.lvls)
# ---> Across pathogens.
tmp.width <- (0.8*14/0.53) * tmp.data[, uniqueN(consensus.pred)]/10
# @ Plot.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=consensus.pred)) +
    geom_jitter(aes(y=bounded.cell.fract, shape=point.status), stroke=4, size=8, color=gen.dot.col, width=0.2, height=0) +
    geom_boxplot(aes(color=consensus.pred, fill=consensus.pred, y=cell.fract), alpha=0.4, width=0.7, linewidth=3, fatten=4, outlier.shape=NA) +
    # geom_hline(yintercept=0.01, linewidth=2, linetype='dashed', color='black') +
    coord_cartesian(ylim=c(0, tmp.bound)) +
    # scale_x_discrete(drop=FALSE) +
    scale_x_discrete(drop=TRUE) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_shape_manual(values=tmp.shapes) +
    scale_color_manual(values=ag.spc.cols) +
    scale_fill_manual(values=ag.spc.cols) +
    labs(x='Reactivity', y='% of T cells', color='', caption=tmp.caption)
tmp.lab <- paste0('/CellFract_Reactivity_Donor_Summ')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=reps.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=14
)

# ---> Number of cells per donor per clonotype, bins filled according to reactivity.
# Stacked clone info.
tmp.data.1 <- plots.data[
    !is.na(donor.id.tag) & !is.na(consensus.pred),
    .(clone.size=.N),
    by=.(
        donor.id.tag, qry.clone.id, consensus.pred
    )
]
to.order <- 1:length(ag.spc.cols); names(to.order) <- names(ag.spc.cols)
tmp.data.1[, order.tag:=to.order[consensus.pred]]
setorderv(x=tmp.data.1, cols=c('order.tag', 'clone.size'), order=c(1, -1))
tmp.data.1[, order.tag:=paste(donor.id.tag, qry.clone.id, consensus.pred, sep=';')]
tmp.data.1$order.tag <- factor(x=tmp.data.1$order.tag, levels=tmp.data.1$order.tag)
# Number of unique clonotypes per donor and cell type.
tmp.data.2 <- plots.data[
    !is.na(donor.id.tag) & !is.na(consensus.pred),
    .(clone.count=uniqueN(qry.clone.id)),
    by=.(donor.id.tag)
]
tmp.data.2[clone.count<100, clone.count:=NA]
# To force all donors to appear on the plot.
tmp.data.3 <- data.table(
    donor.id.tag=levels(tmp.data.2$donor.id.tag)[!levels(tmp.data.2$donor.id.tag) %in% tmp.data.2$donor.id.tag]
)
tmp.data.1 <- rbindlist(
    l=list(tmp.data.1, tmp.data.3),
    use.names=TRUE, fill=TRUE
)
# Define file width according to donor amount.
# tmp.width <- tmp.data.1[, uniqueN(donor.id.tag)]*14.62/40 + 0.38
# tmp.width <- tmp.data[, uniqueN(donor.id.tag)] * 15/40
tmp.width <- (tmp.data.1[, uniqueN(donor.id.tag)] * 6.68/40) + 0.32
tmp.width <- (tmp.width * 5)/1.8
# Plot.
tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=donor.id.tag)) +
    geom_bar(aes(y=clone.size, group=order.tag, color=consensus.pred), stat='identity', position='stack', width=0.8, linewidth=1.4, alpha=0) +
    # geom_point(data=tmp.data.2, aes(y=clone.count), size=7, shape=21, fill='white', color='black') +
    scale_color_manual(values=ag.spc.cols) +
    scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=3)) +
    labs(x='Donor ID', y='Cell count', fill='Reactivity')
tmp.lab <- paste0('/CellCount_Donor_ConsReactivity-StackedClones')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=reps.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.7.2, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
)

# ---> Percentage of clones accounted for by each pathogen among all cells w/ specificity assignments. Summary across donors.
tmp.data <- plots.data[!is.na(consensus.pred)]
tmp.data[, clone.id:=paste(qry.clone.id, donor.id.tag, ';')]
tmp.vals <- c('Cell', 'Clone')
for(tmp.item in tmp.vals){
    if(tmp.item=='Cell'){
        to.plot <- tmp.data[,
            .(freq.abs=.N),
            by=.(cluster=cluster.tag)
        ]
    }else{
        to.plot <- tmp.data[,
            .(freq.abs=uniqueN(clone.id)),
            by=.(cluster=cluster.tag)
        ]
    }
    tmp.ggplot <- ggplot(data=to.plot, aes(x='', y=freq.abs, fill=cluster)) +
        geom_bar(stat='identity', linewidth=0.5, color='black') +
        coord_polar(theta='y') +
        scale_fill_manual(values=cluster.cols) +
        labs(x='', y='', fill='Specificity') +
        theme_void()
    tmp.lab <- paste0('/', tmp.item, 'Fract_Cluster-Summ')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=reps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.4, do.legend=FALSE, do.rotate=FALSE, width=14, height=14
    )
}


### --------------- Antigen-specific T cell phenotypes -------------- ###

phens.path <- paste0(run.path, '/ag_spc_tcr_phens')
if(!dir.exists(phens.path)) dir.create(phens.path)

### -------------------- Quantitative assessment -------------------- ###
### ------------------------- First attempt ------------------------- ###

### -------------------- Quantitative assessment -------------------- ###

# ---> Quantitative enrichment for ag groups within T-cell subsets.

# ---> Preflights
# @ Retrieve data.
quant.data <- plots.data[, .(
    barcode, donor.id.tag,
    clusters.tag=cluster.tag,
    clonotype.tag=qry.clone.id, consensus.pred
)]
tmp.groups <- plots.data[!is.na(consensus.pred), as.character(unique(consensus.pred))]
names(tmp.groups) <- tmp.groups
tmp.groups <- c(
    as.list(tmp.groups),
    ag.groups
)

#   @ Fractions of ag set-specific cells/clonotypes in the group (cluster) out of all ag set-specific cells/clonotypes.
# Preflights
cat.var <- 'clusters.tag'
quant.data <- lapply(X=names(tmp.groups), FUN=function(ag.group){
    tmp.vals <- tmp.groups[[ag.group]]
    tmp.data.1 <- quant.data[
        !is.na(donor.id.tag),
        .(
            cell.freq.abs=.SD[consensus.pred %in% tmp.vals, .N],
            clon.freq.abs=.SD[consensus.pred %in% tmp.vals, uniqueN(clonotype.tag)]
        ),
        by=.(donor.id.tag, cat.tag=get(cat.var))
    ]
    tmp.data.2 <- quant.data[
        !is.na(donor.id.tag) &
        consensus.pred %in% tmp.vals,
        .(
            cell.freq.total=.N,
            clon.freq.total=uniqueN(clonotype.tag)
        ),
        by=.(donor.id.tag)
    ]
    # Apply limit of detection filtering step.
    tmp.data.2 <- tmp.data.2[cell.freq.total>=limit.of.detect]
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor.id.tag'), all=FALSE)
    tmp.data[, `:=`(
        cell.freq.rel=cell.freq.abs/cell.freq.total,
        clon.freq.rel=clon.freq.abs/clon.freq.total
    )]
    tmp.cols <- setdiff(x=colnames(tmp.data), y=c('cell.freq.total', 'clon.freq.total'))
    tmp.data <- tmp.data[, ..tmp.cols]
    tmp.check <- tmp.data[, .(check=round(x=sum(cell.freq.rel))), by=.(donor.id.tag)][, all(check==1)]
    if(!tmp.check) stop('Unexpected error!\n')
    return(tmp.data)
})
names(quant.data) <- names(tmp.groups)
quant.data <- rbindlist(l=quant.data, use.names=TRUE, idcol='ag.group')

# ---> Donor-specific information.
tmp.reports.path <- paste0(phens.path, '/donor-spc_info')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

uniq.pops <- quant.data[, sort(unique(as.character(cat.tag)))]
uniq.pps <- names(ag.spc.cols)[names(ag.spc.cols) %in% quant.data$ag.group]
# Fetch data to plot (bar plots, absolute cell frequency).
brpt.data <- quant.data[ag.group %in% uniq.pps]
brpt.data <- brpt.data[, .(cell.count=sum(cell.freq.abs)), by=.(ag.group, donor.id.tag)]
tmp.vals <- factor(x=brpt.data$ag.group, levels=uniq.pps)
set(x=brpt.data, j='ag.group', value=tmp.vals)
# Fetch data to plot (heatmap).
hmap.data <- quant.data[, .(donor.id.tag=as.character(donor.id.tag), group=paste(cat.tag, ag.group, sep='.'), cell.freq.rel)]
hmap.data <- spread(data=hmap.data, key=group, value=cell.freq.rel, drop=FALSE)
row.names(hmap.data) <- hmap.data$donor.id.tag; hmap.data$donor.id.tag <- NULL
hmap.data <- as.matrix(hmap.data)
tmp.cols <- paste(rep(x=uniq.pops, each=length(uniq.pps)), uniq.pps, sep='.')
hmap.data <- hmap.data[, tmp.cols]
# Hierarchical clustering (if possible)
clust.obj <- dist(x=hmap.data, method='euclidean')
clust.obj <- tryCatch(
    expr={hclust(d=clust.obj, method='ward.D')},
    error=function(e) NULL
)
# Set metadata.
col.meta <- data.frame(
    row.names=colnames(hmap.data),
    `Cluster`=str_extract(
        string=colnames(hmap.data),
        pattern='^\\d+'
    ),
    `Specificity`=str_replace(
        string=colnames(hmap.data),
        pattern='^\\d+\\.', replacement=''
    )
)
# Define colors for metadata tracks.
ann.colors <- list(
    `Cluster`=cluster.cols[names(cluster.cols) %in% col.meta$Cluster],
    `Specificity`=ag.spc.cols[names(ag.spc.cols) %in% col.meta$Specificity]
)
# Set color scale and breaks for heatmap.
break.no <- 200
col.breaks <- seq(from=min(range(hmap.data, na.rm=TRUE)), to=max(range(hmap.data, na.rm=TRUE)), length.out=break.no)
mid.point <- (max(col.breaks) - min(col.breaks)) / 2
mid.point <- which.min(abs(col.breaks-mid.point))
hmap.col.scale.1 <- colorRampPalette(c('#414487', '#355F8A', '#2A788E', '#22A884'))(mid.point)
hmap.col.scale.2 <- colorRampPalette(c('#22A884', '#7AD151', '#FDE725', '#FFEA00'))(break.no-(mid.point+1))
hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
na.col <- '#000000'
# File specs.
tmp.height <- (nrow(hmap.data) * 0.18) + 0.3
tmp.width <- ncol(hmap.data) * 0.18
# @ Version with hierarchical clustering.
if(!is.null(clust.obj)){
    # Complete version.
    tmp.file.name <- paste0(
        tmp.reports.path, '/',
        'Donor_Reactivity_CellFract_WClust.C.pdf'
    )
    pheatmap(
        mat=hmap.data, scale='none',
        color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
        cluster_rows=clust.obj, cluster_cols=FALSE,
        annotation_row=NULL, annotation_col=col.meta, annotation_colors=ann.colors,
        show_colnames=FALSE, show_rownames=TRUE, legend=TRUE, annotation_legend=TRUE,
        filename=tmp.file.name, height=tmp.height+3, width=tmp.width+3
    )
    # Blank version.
    tmp.file.name <- paste0(
        tmp.reports.path, '/',
        'Donor_Reactivity_CellFract_WClust.B.pdf'
    )
    pheatmap(
        mat=hmap.data, scale='none',
        color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
        cluster_rows=clust.obj, cluster_cols=FALSE,
        annotation_col=col.meta, annotation_colors=ann.colors,
        show_colnames=FALSE, show_rownames=FALSE,
        legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
        filename=tmp.file.name, height=tmp.height+1, width=tmp.width+1
    )
    # @ Number of cells per specificity per donor.
    tmp.lvls <- rev(x=clust.obj$labels[clust.obj$order])
    tmp.vals <- factor(x=as.character(brpt.data[, donor.id.tag]), levels=tmp.lvls)
    set(x=brpt.data, j='donor.id.tag', value=tmp.vals)
    tmp.ggplot <- ggplot(data=brpt.data, aes(x=donor.id.tag, y=cell.count, fill=ag.group)) +
        geom_bar(stat='identity', position='stack', width=0.6, color='black', linewidth=1.3) +
        scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
        scale_fill_manual(values=ag.spc.cols) +
        labs(x='Donor ID', y='Cell count', fill='Specificity')
    tmp.lab <- paste0('SpcCounts_DonorID_Clusters_HClustOrder')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE,
        width=tmp.height*5, height=7
    )
}
# @ Version w/out hierarchical clustering.
# Complete version.
tmp.file.name <- paste0(
    tmp.reports.path, '/',
    'Donor_Reactivity_CellFract_NoClust.C.pdf'
)
pheatmap(
    mat=hmap.data, scale='none',
    color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
    cluster_rows=FALSE, cluster_cols=FALSE,
    annotation_row=NULL, annotation_col=col.meta, annotation_colors=ann.colors,
    show_colnames=FALSE, show_rownames=TRUE, legend=TRUE, annotation_legend=TRUE,
    filename=tmp.file.name, height=tmp.height+2, width=tmp.width+2
)
# Blank version.
tmp.file.name <- paste0(
    tmp.reports.path, '/',
    'Donor_Reactivity_CellFract_NoClust.B.pdf'
)
pheatmap(
    mat=hmap.data, scale='none',
    color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
    cluster_rows=FALSE, cluster_cols=FALSE,
    annotation_col=col.meta, annotation_colors=ann.colors,
    show_colnames=FALSE, show_rownames=FALSE,
    legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
    filename=tmp.file.name, height=tmp.height, width=tmp.width
)
# @ Color legend only.
# W/ labels
# tmp.breaks <- if(tmp.group %in% names(ag.group.breaks)) NULL else ag.group.breaks[[tmp.group]]
tmp.breaks <- NULL
tmp.data <- as.data.table(gather(data=as.data.frame(hmap.data), key=key, value=value))
tmp.data <- tmp.data[!is.na(value)]
tmp.ggplot <- ggplot(data=tmp.data, aes(x=value, y=value, col=value)) + 
    geom_point() +
    scale_color_gradientn(colors=hmap.col.scale) +
    theme(
        legend.ticks=element_line(color='black', linewidth=0.6),
        legend.ticks.length=unit(0.22, "cm"),
        legend.frame=element_rect(color='black', linewidth=0.6)
    )
tmp.ggplot <- get_legend(p=tmp.ggplot)
tmp.file.name <- paste0(
    tmp.reports.path, '/',
    'Donor_Reactivity_CellFract_Legend.L1.pdf'
)
pdf(file=tmp.file.name, height=2, width=2)
print(as_ggplot(tmp.ggplot))
dev.off()
# W/out labels
tmp.ggplot <- ggplot(data=tmp.data, aes(x=value, y=value, col=value)) + 
    geom_point() +
    scale_color_gradientn(colors=hmap.col.scale, name=NULL, labels=NULL) +
    theme(
        legend.ticks=element_line(color='black', linewidth=0.6),
        legend.ticks.length=unit(0.22, "cm"),
        legend.frame=element_rect(color='black', linewidth=0.6)
    )
tmp.ggplot <- get_legend(p=tmp.ggplot)
tmp.file.name <- paste0(
    tmp.reports.path, '/',
    'Donor_Reactivity_CellFract_Legend.L2.pdf'
)
pdf(file=tmp.file.name, height=1.25, width=0.3)
print(as_ggplot(tmp.ggplot))
dev.off()

# ---> Summarization across donors.
tmp.reports.path <- paste0(phens.path, '/summ_info')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# General cluster-wise fraction of all T cells.
tmp.data.1 <- gex.meta[
    !is.na(donor.id.tag),
    .(freq.abs=.N),
    by=.(donor.id.tag, cat.tag=cluster.tag)
]
tmp.data.2 <- gex.meta[
    !is.na(donor.id.tag),
    .(freq.total=.N),
    by=.(donor.id.tag)
]
gen.plot.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
gen.plot.data[, cell.freq.rel:=freq.abs/freq.total]
# Apply consistent absolute frequency threshold.
gen.plot.data <- gen.plot.data[freq.total>limit.of.detect]
gen.plot.data[, `:=`(freq.abs=NULL, freq.total=NULL)]
gen.plot.data[, ag.group:='ALL T']
# Process per cluster.
tmp.pops <- as.character(quant.data[, unique(cat.tag)])
for(tmp.pop in tmp.pops){
    # Fetch data to plot.
    tmp.data <- quant.data[cat.tag==tmp.pop, ]
    tmp.data <- tmp.data[ag.group %in% names(ag.spc.cols)]
    # ---> Stats.
    # Kruskal Wallis test.
    tmp.test.1 <- kruskal.test(
        x=tmp.data,
        formula=cell.freq.rel~ag.group
    )
    tmp.test.1 <- tmp.test.1$p.value
    # Multiple pairwise Wilcoxon signed-rank test.
    tmp.test.2 <- pairwise.wilcox.test(
        x=tmp.data$cell.freq.rel, g=tmp.data$ag.group,
        p.adjust.method='BH'
    )
    tmp.test.2 <- tmp.test.2$p.value
    tmp.lab <- paste0('Pop-', tmp.pop, '_CellFract_Pat_Donor_MultPairwiseComps.csv')
    tmp.file.name <- paste0(tmp.reports.path, '/', tmp.lab)
    write.csv(file=tmp.file.name, x=tmp.test.2, quote=TRUE, row.names=TRUE)
    # Continue w/ preparation.
    tmp.thold <- 1
    # Add cell fractions for all T cells.
    tmp.data <- list(
        tmp.data,
        gen.plot.data
    )
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, fill=TRUE)
    # Set factors
    tmp.lvls <- names(ag.spc.cols)[names(ag.spc.cols) %in% tmp.data[, ag.group]]
    tmp.lvls <- c('ALL T', tmp.lvls)
    tmp.vals <- factor(x=tmp.data[, ag.group], levels=tmp.lvls)
    set(x=tmp.data, j='ag.group', value=tmp.vals)
    # Set upper bound.
    tmp.data[, bounded.cell.freq.rel:=cell.freq.rel]
    tmp.data[, point.status:='ori']
    tmp.data[bounded.cell.freq.rel>tmp.thold, `:=`(bounded.cell.freq.rel=tmp.thold, point.status='bounded')]
    # Set bound.
    if(tmp.thold==1) tmp.thold <- NA
    # Dimensions and theme
    tmp.height <- 14
    tmp.width <- ((tmp.data[, uniqueN(ag.group)]-1)*2)
    tmp.blank.comp <- blank.complement.3.2
    # ---> Version w/ ag-specific cell fractions only.
    tmp.ggplot <- ggplot(data=tmp.data[ag.group!='ALL T'], aes(x=ag.group)) +
        geom_jitter(aes(y=bounded.cell.freq.rel, shape=point.status), stroke=5, size=8, color=gen.dot.col, width=0, height=0) +
        geom_boxplot(aes(y=cell.freq.rel, color=ag.group, fill=ag.group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        # coord_cartesian(ylim=c(0, tmp.thold)) +
        scale_y_continuous(
            limits=c(0, 1),
            breaks=scales::pretty_breaks(n=3)
        ) +
        scale_shape_manual(values=tmp.shapes) +
        scale_color_manual(values=ag.spc.cols) +
        scale_fill_manual(values=ag.spc.cols) +
        labs(
            x='Group',
            y=paste0('Percentage ag-specific cells in indicated groups'),
            color='',
            caption=paste0(
                'Each dot represents an individual donor.\n',
                'Kruskal-Wallis Rank Sum Test P-value: ',
                tmp.test.1
            )
        )
    tmp.lab <- paste0('Pop-', tmp.pop, '_CellFract_Pat_Donor')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
        height=tmp.height, width=tmp.width,
        blank.comp=tmp.blank.comp, do.legend=FALSE, do.rotate=TRUE
    )
    # ---> Version including ALL T cell group.
    tmp.cols <- c(
        `ALL T`='#000000',
        ag.spc.cols
    )
    tmp.width <- tmp.width + 2
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=ag.group)) +
        geom_jitter(aes(y=bounded.cell.freq.rel, shape=point.status), stroke=5, size=8, color=gen.dot.col, width=0, height=0) +
        geom_boxplot(aes(y=cell.freq.rel, color=ag.group, fill=ag.group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        # coord_cartesian(ylim=c(0, tmp.thold)) +
        scale_y_continuous(
            limits=c(0, 1),
            breaks=scales::pretty_breaks(n=3)
        ) +
        scale_shape_manual(values=tmp.shapes) +
        scale_color_manual(values=tmp.cols) +
        scale_fill_manual(values=tmp.cols) +
        labs(
            x='Group',
            y=paste0('Percentage ag-specific cells in indicated groups'),
            color='',
            caption=paste0(
                'Each dot represents an individual donor.\n',
                'Kruskal-Wallis Rank Sum Test P-value: ',
                tmp.test.1
            )
        )
    tmp.lab <- paste0('Pop-', tmp.pop, '_CellFract_Pat-and-ALL_Donor')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
        height=tmp.height, width=tmp.width,
        blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE
    )
}

if(!is.null(donor.meta)){
    # ----> Associations, preflights.
    # @ Population-specific frequencies
    tmp.data.1 <- gex.meta[,
        .(freq.abs=.N),
        by=.(
            donor.id.tag,
            pop.tag=paste0('C', cluster.tag)
        )
    ]
    tmp.data.2 <- gex.meta[,
        .(freq.total=.N),
        by=.(donor.id.tag)
    ]
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
    tmp.data[, freq.rel:=freq.abs/freq.total]
    tmp.data <- tmp.data[, .(donor.id.tag, pop.tag, freq.rel)]
    tmp.data <- as.data.table(spread(data=tmp.data, key=pop.tag, value=freq.rel, fill=0))
    # Merge w/ antigen-specific fractions.
    plot.data <- merge(
        x=quant.data, y=tmp.data,
        by='donor.id.tag', all.x=TRUE, all.y=FALSE
    )
    # Merge w/ donor metadata
    plot.data <- merge(x=plot.data, y=donor.meta, by='donor.id.tag', all.x=TRUE)

    # ---> Association between predictions and continuous variables.
    tmp.reports.path <- paste0(phens.path, '/assoc_vars-cont')
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    # Define continuous variables.
    cont.vars <- sapply(X=colnames(donor.meta), FUN=function(tmp.col){
        is.numeric(donor.meta[[tmp.col]]) &
        tmp.col!='donor.id.tag'
    })
    cont.vars <- c(
        names(donor.meta)[cont.vars],
        paste0('C', gex.meta[, unique(cluster.tag)])
    )
    names(cont.vars) <- str_to_upper(str_replace_all(
        string=str_replace_all(string=cont.vars, pattern='\\.tag$|^donor\\.', replacement=''),
        pattern='\\.|_', replacement=' '
    ))
    # Process per continuous variable.
    uniq.spcs <- plot.data[, unique(ag.group)]
    for(cont.var in names(cont.vars)){
        for(tmp.spc in uniq.spcs){
            tmp.data <- plot.data[
                ag.group==tmp.spc,
                .(
                    cont.var=get(cont.vars[cont.var]),
                    freq.var=cell.freq.rel,
                    tmp.var=cat.tag
                )
            ]
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=cont.var, y=freq.var, color=tmp.var)) +
                geom_point(shape=1, stroke=5, size=12, color='black') +
                geom_smooth(method='lm', formula=y~x, se=FALSE, fullrange=TRUE, linewidth=5) +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
                scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                scale_color_manual(values=cluster.cols) +
                labs(x='Continuous variable', y='% of ag-specific T cells')
            tmp.lab <- paste0(
                '/Association_', tmp.spc, '_',  cont.var
            )
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                stat.cor=TRUE, cor.group='tmp.var',
                blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
            )
        }
    }

    # ---> Association between predictions and discrete variables.
    tmp.reports.path <- paste0(phens.path, '/assoc_vars-disc')
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    # Define discrete variables.
    disc.vars <- sapply(X=colnames(donor.meta), FUN=function(tmp.col){
        (is.character(donor.meta[[tmp.col]]) | is.factor(donor.meta[[tmp.col]])) &
        tmp.col!='donor.id.tag'
    })
    disc.vars <- colnames(donor.meta)[disc.vars]
    names(disc.vars) <- str_to_upper(str_replace_all(
        string=str_replace_all(string=disc.vars, pattern='\\.tag$|^donor\\.', replacement=''),
        pattern='\\.|_', replacement=' '
    ))
    # Process per cell type.
    uniq.spcs <- plot.data[, unique(ag.group)]
    for(tmp.spc in uniq.spcs){
        for(disc.var in names(disc.vars)){
            tmp.vals <- sort(as.character(plot.data[, unique(cat.tag)]))
            tmp.file.name <- paste0(
                tmp.reports.path, '/Comparison_', tmp.spc, '_',  disc.var, '.pdf'
            )
            pdf(file=tmp.file.name)
            for(tmp.val in tmp.vals){
                tmp.data <- plot.data[
                    ag.group==tmp.spc &
                    cat.tag==tmp.val,
                    .(
                        disc.var=get(disc.vars[disc.var]),
                        freq.var=cell.freq.rel
                    )
                ]
                # Stats.
                if(tmp.data[, uniqueN(disc.var)==2]){
                    tmp.groups <- tmp.data[, unique(disc.var)]
                    tmp.test <- wilcox.test(
                        x=tmp.data[disc.var==tmp.groups[1], freq.var],
                        y=tmp.data[disc.var==tmp.groups[2], freq.var],
                        paired=FALSE
                    )
                    tmp.caption <- paste0('P-value is ', round(x=tmp.test$p.value, digits=6))
                }else{
                    tmp.caption <- ''
                }
                # Set upper bound.
                tmp.bound <- tmp.data[, max(freq.var)]
                tmp.data[, bounded.var:=freq.var]
                tmp.data[, point.status:='ori']
                tmp.data[freq.var>tmp.bound, `:=`(bounded.var=tmp.bound, point.status='bounded')]
                # Plot.
                tmp.ggplot <- ggplot(data=tmp.data, aes(x=disc.var)) +
                    geom_boxplot(aes(y=freq.var, color=disc.var), width=0.7, linewidth=4, outlier.shape=NA) +
                    geom_jitter(aes(y=bounded.var, shape=point.status), stroke=5, size=10, color='black', width=0.22) +
                    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                    scale_shape_manual(values=tmp.shapes) +
                    # scale_color_manual(values=tmp.cols) +
                    coord_cartesian(ylim=c(0, tmp.bound)) +
                    labs(
                        title=paste0('Cluster ', tmp.val),
                        x='', y='% of ag-specific T cells',
                        color='', caption=tmp.caption
                    ) +
                    theme(legend.position='none') +
                    theme_bw()
                print(tmp.ggplot)
            }
            dev.off()
        }
    }

    # ----> Associations, summary.
    tmp.reports.path <- paste0(phens.path, '/assoc_vars-cont-summ')
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    # Data retrieval
    #       Response variables.
    tmp.data.1 <- quant.data[
        ag.group %in% names(ag.spc.cols),
        .(ag.group, donor.id.tag, cluster=paste0('C', cat.tag), cell.freq.rel)
    ]
    tmp.data.1[, key:=paste(ag.group, cluster, sep='.')]
    tmp.data.1 <- tmp.data.1[, .(donor.id.tag, key, cell.freq.rel)]
    tmp.data.1 <- as.data.table(spread(data=tmp.data.1, key=key, value=cell.freq.rel, fill=NA))
    var.set.1 <- setdiff(x=colnames(tmp.data.1), y='donor.id.tag')
    #       Explanatory variables.
    tmp.data.2.1 <- gex.meta[
        !is.na(donor.id.tag),
        .(freq.abs=.N),
        by=.(
            donor.id.tag,
            cluster=paste0('C', cluster.tag)
        )
    ]
    tmp.data.2.2 <- gex.meta[
        !is.na(donor.id.tag),
        .(freq.total=.N),
        by=.(donor.id.tag)
    ]
    tmp.data.2 <- merge(x=tmp.data.2.1, y=tmp.data.2.2, by='donor.id.tag')
    tmp.data.2[, freq.rel:=freq.abs/freq.total]
    tmp.data.2 <- tmp.data.2[, .(
        donor.id.tag, cluster, freq.rel
    )]
    tmp.data.2 <- as.data.table(spread(data=tmp.data.2, key=cluster, value=freq.rel, fill=0))
    #       Other clinical and demographical variables.
    tmp.cols <- c('donor.id.tag', cont.vars[cont.vars %in% colnames(donor.meta)])
    tmp.data.2 <- merge(
        x=tmp.data.2, y=donor.meta[, ..tmp.cols],
        by='donor.id.tag', all=TRUE
    )
    var.set.2 <- setdiff(x=colnames(tmp.data.2), y='donor.id.tag')
    # @ Full set of vars.
    to.plot <- merge(x=tmp.data.2, y=tmp.data.1, by='donor.id.tag')
    to.plot <- as.data.frame(to.plot)
    row.names(to.plot) <- to.plot$donor.id.tag; to.plot$donor.id.tag <- NULL
    to.plot <- as.matrix(to.plot)
    # @ Pairwise distance and significance calculation
    corr.data <- rcorr(to.plot, type='spearman')
    r.mat <- corr.data$r
    r.mat <- r.mat[
        var.set.2,
        var.set.1
    ]
    p.mat <- corr.data$P
    p.mat <- p.mat[
        var.set.2,
        var.set.1
    ]
    # @ Correlation plot, across the board
    tmp.width <- ncol(r.mat)/2
    tmp.height <- nrow(r.mat)/2
    tmp.file.name <- paste0(
        tmp.reports.path, '/',
        'CorrPlot_Dist-Spearman_ALL',
        '.C.pdf'
    )
    pdf(file=tmp.file.name, width=tmp.width+2, height=tmp.height+2)
    corrplot(
        corr=r.mat,
        method = 'circle', type='full',
        col=rev(x=COL2('RdBu', 200)),
        diag=TRUE,
        p.mat=p.mat, sig.level=0.01, insig='blank',
        outline=TRUE,
        addgrid.col='black',
        tl.pos='lt',
        tl.cex=1, tl.col='black',
        cl.pos='r'
    )
    dev.off()
    # @ Correlation plot, one plot per pathogen specificity.
    uniq.spcs <- plots.data[!is.na(consensus.pred), as.character(unique(consensus.pred))]
    for(tmp.spc in uniq.spcs){
        tmp.cols <- colnames(r.mat)[colnames(r.mat) %like% paste0('^', tmp.spc)]
        r.plot.mat <- r.mat[, tmp.cols]
        if(all(is.na(r.plot.mat))) next
        p.plot.mat <- p.mat[, tmp.cols]
        tmp.width <- ncol(r.plot.mat)/2
        tmp.height <- nrow(r.plot.mat)/2
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            'CorrPlot_Dist-Spearman_Spc-', tmp.spc,
            '.C.pdf'
        )
        pdf(file=tmp.file.name, width=tmp.width+2, height=tmp.height+2)
        corrplot(
            corr=r.plot.mat,
            method = 'circle', type='full',
            col=rev(x=COL2('RdBu', 200)),
            diag=TRUE,
            p.mat=p.plot.mat, sig.level=0.01, insig='blank',
            outline=TRUE,
            addgrid.col='black',
            tl.pos='lt',
            tl.cex=1, tl.col='black',
            cl.pos='r'
        )
        dev.off()
    }
}


### ----------------- Dimensionality reduction plots ---------------- ###

dim.red.path <- paste0(run.path, '/dim_reduction')
if(!dir.exists(dim.red.path)) dir.create(dim.red.path)

# ---> Predicted TCRs on UMAP
# @ Fetch data.
tmp.data <- copy(plots.data)
# @ Plot.
# Define scale values.
alpha.scale <- c(
    `no-spc`=0,
    `spc`=1,
    `den`=0.3
)
# Process per specificity
tmp.spcs <- names(ag.spc.cols)[names(ag.spc.cols) %in% plots.data[, consensus.pred]]
tmp.spcs <- c('All', tmp.spcs)
for(tmp.spc in tmp.spcs){
    # Set color variable
    tmp.data[, col.val:=clusters.tag]
    # Set scaling variable
    tmp.data[, scale.val:=ifelse(
        test=!is.na(consensus.pred) & consensus.pred==tmp.spc,
        yes='spc',
        no='no-spc'
    )]
    if(tmp.data[scale.val=='spc', .N<10]) next
    # Set row order, allowing for dots w/ predictions to come on top.
    setorderv(x=tmp.data, col='scale.val', order=-1)
    # Density data.
    plot.data <- if(tmp.spc!='All') tmp.data[!is.na(consensus.pred) & consensus.pred==tmp.spc] else tmp.data
    # Plot
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=UMAP_1, y=UMAP_2, color=col.val)) +
        geom_point(aes(alpha=scale.val), size=0.3) +
        stat_density_2d(
            data=plot.data,
            aes(fill=after_stat(level), alpha='den'),
            geom='polygon', linewidth=1,
            contour=TRUE,
            color='black', bins=5
        ) +
        scale_color_manual(values=cluster.cols) +
        scale_fill_gradient(low="#FFFFD4", high="#F6407F") +
        scale_alpha_manual(values=alpha.scale) +
        labs(x='UMAP 1', y='UMAP 2', col='Cluster +\nReactivity', size='Pred?', alpha='Pred?')
    tmp.lab <- paste0('ImpClonesOnUMAP_Spc-', tmp.spc)
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=dim.red.path, file.name=tmp.lab, type='tiff',
        blank.comp=blank.complement.1, do.legend=FALSE,
        height=10, width=10
    )
}


### -------------------------- Final report -------------------------- ###

# ---> Minimal report
# List barcodes and their predictions.
tmp.data.1 <- gex.meta[, .(barcode, donor.id.tag)]
tmp.data.2 <- cells.clons.info[, .(barcode, qry.clone.id)]
tmp.data.1 <- merge(x=tmp.data.1, y=tmp.data.2, by='barcode', all.x=TRUE, all.y=FALSE, sort=FALSE)
tmp.data.2 <- pred.summ[, .(qry.clone.id, donor.id.tag, consensus.pred)]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('qry.clone.id', 'donor.id.tag'), all.x=TRUE, all.y=FALSE, sort=FALSE)
tmp.data <- tmp.data[, .(barcode, clonotype.tag=qry.clone.id, consensus.pred)]
tmp.check <- tmp.data[, uniqueN(barcode)]==gex.meta[, uniqueN(barcode)] & tmp.data[, uniqueN(barcode)]==tmp.data[, .N]
if(!tmp.check) stop('Unexpected error w/ minimal report.\n')
# Randomly select a set of background cells (i.e., cells left without any specificity prediction)
if(tmp.data[is.na(consensus.pred), .N]<cell.bg.size) cell.bg.size <- tmp.data[is.na(consensus.pred), .N]
bg.cells <- sample(x=tmp.data[is.na(consensus.pred), barcode], size=cell.bg.size, replace=FALSE)
tmp.data[, ext.consensus.pred:=consensus.pred]
tmp.data[barcode %in% bg.cells, ext.consensus.pred:='Background']
# Output report.
tmp.file.name <- paste0(run.path, '/BarcodePredictions.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# ---> Full report.
final.output <- merge(
    x=clone.info, y=pred.summ,
    by.x=c('clonotype.tag', 'donor.id.tag'),
    by.y=c('qry.clone.id', 'donor.id.tag'),
    all=TRUE
)


# ---> Add HLA info.
# BROKEN. TO BE FIXED.
# @ Function to calculate match consensus on a reference donor basis.
# get.hla.match.consensus.1 <- function(ref.donor.info, final.output, ref.hla.info, qry.hla.info){
#     # Retrieve unique reference donors for each query clonotype.
#     ref.donor.info <- apply(X=ref.donor.info, MARGIN=1, FUN=function(x){
#         x <- x[!is.na(x)]
#         if(length(x)==0) return(NA)
#         x <- unlist(str_split(string=x, pattern='\\|'))
#         x <- unique(sort(x))
#         return(x)
#     })
#     # Compute HLA matches.
#     tmp.data <- unlist(lapply(X=1:length(ref.donor.info), FUN=function(i){
#         ref.donors <- ref.donor.info[[i]]
#         if(all(is.na(ref.donors))) return(NA)
#         qry.donor <- final.output[i, donor.id.tag]
#         ref.donors.hlas <- ref.hla.info[ref.donors]
#         int.hlas <- lapply(X=ref.donors.hlas, FUN=function(ref.donor.hlas){
#             tmp.data <- intersect(
#                 x=ref.donor.hlas,
#                 y=qry.hla.info[[qry.donor]]
#             )
#             if(length(tmp.data)==0) return(NA)
#             return(tmp.data)
#         })
#         if(all(is.na(int.hlas))) return('Mismatch')
#         # Attempt to make perfect matches.
#         int.hlas <- Reduce(x=int.hlas, f=intersect)
#         if(length(int.hlas)>0){
#             int.hlas <- paste0(int.hlas, collapse='|')
#             return(int.hlas)
#         }else{
#             return('Partial match')
#         }
#     }))
#     return(tmp.data)
# }
# # @ Function to calculate match consensus on a reference clonotype basis.
# get.hla.match.consensus.2 <- function(ref.hla.info, final.output, qry.hla.info){
#     # Retrieve unique reference HLAs for each query clonotype.
#     ref.hla.info <- apply(X=ref.hla.info, MARGIN=1, FUN=function(x){
#         x <- x[!is.na(x)]
#         if(length(x)==0) return(NA)
#         x <- unlist(str_split(string=x, pattern='\\|'))
#         x <- unique(sort(x))
#         return(x)
#     })
#     # Compute HLA matches.
#     tmp.data <- unlist(lapply(X=1:length(ref.hla.info), FUN=function(i){
#         ref.hlas <- ref.hla.info[[i]]
#         if(all(is.na(ref.hlas))) return(NA)
#         qry.donor <- final.output[i, donor.id.tag]
#         int.hlas <- intersect(
#             x=ref.hlas,
#             y=qry.hla.info[[qry.donor]]
#         )
#         if(length(int.hlas)==0) return('Mismatch')
#         int.hlas <- paste0(int.hlas, collapse='|')
#         return(int.hlas)
#     }))
#     return(tmp.data)
# }
# # @ Process on a case-by-case basis.
# if(!is.null(qry.hla.info)){
#     if(!is.null(ref.hla.info)){
#         # @ For predictions.
#         if(is.list(x=ref.hla.info)){
#             # ---> When HLA info is provided on a reference donor basis.
#             # Retrieve reference donor IDs for predictions.
#             pred.donor.cols <- paste(tool.names, ref.donor.id.tag, sep='|')
#             pred.donor.cols <- pred.donor.cols[pred.donor.cols %in% colnames(final.output)]
#             tmp.check <- length(pred.donor.cols)>0
#             if(!tmp.check) stop('Failed to find columns referring to donor IDs for predictions. (1)\n')
#             tmp.data <- final.output[, ..pred.donor.cols]
#             # Compute HLA matches.
#             tmp.data <- get.hla.match.consensus.1(
#                 ref.donor.info=tmp.data,
#                 final.output=final.output,
#                 ref.hla.info=ref.hla.info, qry.hla.info=qry.hla.info
#             )
#             set(x=final.output, j='pred.hla.consensus', value=tmp.data)
#             # Define general type of HLA match for predictions.
#             final.output[, gen.pred.hla:=pred.hla.consensus]
#             final.output[
#                 !is.na(gen.pred.hla) & gen.pred.hla!='Mismatch' & gen.pred.hla!='Partial match',
#                 gen.pred.hla:='Match'
#             ]
#         }else{
#             # ---> When HLA info is provided on a reference clonotype basis.
#             # Retrieve reference HLAs for predictions.
#             pred.hla.cols <- paste(tool.names, ref.hla.info, sep='|')
#             pred.hla.cols <- pred.hla.cols[pred.hla.cols %in% colnames(final.output)]
#             tmp.check <- length(pred.hla.cols)>0
#             if(!tmp.check) stop('Failed to find columns referring to HLA info for predictions. (1)\n')
#             tmp.data <- final.output[, ..pred.hla.cols]
#             # Compute HLA matches.
#             tmp.data <- get.hla.match.consensus.2(
#                 ref.hla.info=tmp.data,
#                 final.output=final.output,
#                 qry.hla.info=qry.hla.info
#             )
#             set(x=final.output, j='pred.hla.consensus', value=tmp.data)
#             # Define general type of HLA match for predictions.
#             final.output[, gen.pred.hla:=pred.hla.consensus]
#             final.output[
#                 !is.na(gen.pred.hla) & gen.pred.hla!='Mismatch' & gen.pred.hla!='Partial match',
#                 gen.pred.hla:='Match'
#             ]
#         }
#         # @ For CDR3 matches.
#         if(is.list(x=ref.hla.info)){
#             # ---> When HLA info is provided on a reference donor basis.
#             # Retrieve reference donor IDs for CDR3 matches.
#             cdr.donor.cols <- paste(c('MatchBeta', 'MatchAlpha'), ref.donor.id.tag, sep='|')
#             cdr.donor.cols <- cdr.donor.cols[cdr.donor.cols %in% colnames(final.output)]
#             tmp.check <- length(cdr.donor.cols)==2
#             if(!tmp.check) stop('Failed to find columns referring to donor IDs for CDR3 matches. (2)\n')
#             tmp.data <- final.output[, ..cdr.donor.cols]
#             # Compute HLA matches.
#             tmp.data <- get.hla.match.consensus.1(
#                 ref.donor.info=tmp.data,
#                 final.output=final.output,
#                 ref.hla.info=ref.hla.info, qry.hla.info=qry.hla.info
#             )
#             set(x=final.output, j='match.hla.consensus', value=tmp.data)
#             # Define general type of HLA match for predictions.
#             final.output[, gen.match.hla:=match.hla.consensus]
#             final.output[
#                 !is.na(gen.match.hla) & gen.match.hla!='Mismatch' & gen.match.hla!='Partial match',
#                 gen.match.hla:='Match'
#             ]
#         }else{
#             # ---> When HLA info is provided on a reference clonotype basis.
#             # Retrieve reference HLAs for CDR3 matches.
#             cdr.hla.cols <- paste(c('MatchBeta', 'MatchAlpha'), ref.hla.info, sep='|')
#             cdr.hla.cols <- cdr.hla.cols[cdr.hla.cols %in% colnames(final.output)]
#             tmp.check <- length(cdr.hla.cols)==2
#             if(!tmp.check) stop('Failed to find columns referring to HLA info for CDR3 matches. (1)\n')
#             tmp.data <- final.output[, ..cdr.hla.cols]
#             # Compute HLA matches.
#             tmp.data <- get.hla.match.consensus.2(
#                 ref.hla.info=tmp.data,
#                 final.output=final.output,
#                 qry.hla.info=qry.hla.info
#             )
#             set(x=final.output, j='match.hla.consensus', value=tmp.data)
#             # Define general type of HLA match for predictions.
#             final.output[, gen.match.hla:=match.hla.consensus]
#             final.output[
#                 !is.na(gen.match.hla) & gen.match.hla!='Mismatch' & gen.match.hla!='Partial match',
#                 gen.match.hla:='Match'
#             ]
#         }
#     }
# }


# ---> Final format.
# @ Set row order
setorderv(x=final.output, cols='rel.freq', order=-1)
# Save full report.
tmp.file.name <- paste0(run.path, '/IntersectedPredictionDetails.csv')
fwrite(file=tmp.file.name, x=final.output, na=NA, quote=TRUE)


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
