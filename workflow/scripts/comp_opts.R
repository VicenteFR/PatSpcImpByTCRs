############    -------   Predictions evaluation    -------    ############
cat('\n\n')
cat('############    -------   Predictions evaluation    -------    ############\n')

# By: Vicente Fajardo-Rosas

### -------------------------- Description -------------------------- ###


cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(optparse)


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
    make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Char, absolute path to directory to save results to.\n"),
    make_option(opt_str="--RunPath", type="character", default=NULL, dest="run.path", help="Char, absolute path to directory where clustering results have been saved to.\n"),
    make_option(opt_str="--OptsFile1", type="character", default=NULL, dest="yaml.file.1", help="Char, absolute path to yaml file for workflow; general options.\n"),
    make_option(opt_str="--OptsFile2", type="character", default=NULL, dest="yaml.file.2", help="Char, absolute path to yaml file for workflow; project-specific options.\n"),
    make_option(opt_str="--RefID", type="character", default=NULL, dest="ref.name", help="Char, indicates the reference ID for the results to be evaluated.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
reports.path <- opt$reports.path
run.path <- opt$run.path
yaml.file.1 <- opt$yaml.file.1
yaml.file.2 <- opt$yaml.file.2
ref.name <- opt$ref.name
# ---> Define constant vars.
set.seed(seed=1)

# ---> Options file.
# Load yaml file
if(!file.exists(yaml.file.1)) stop(paste0('Attempted definition of yaml file (general options) failed. Next file could not be found:\n', yaml.file.1, '\n'))
if(!file.exists(yaml.file.2)) stop(paste0('Attempted definition of yaml file (project-specific options) failed. Next file could not be found:\n', yaml.file.2, '\n'))
tmp.data.1 <- yaml::read_yaml(file=yaml.file.1)
tmp.data.2 <- yaml::read_yaml(file=yaml.file.2)
gen.wflow.opts <- c(tmp.data.2, tmp.data.1[!names(tmp.data.1) %in% names(tmp.data.2)])
# Define files and general options.
tholds.per.ref <- gen.wflow.opts$ref_tholds
ref.tholds <- tholds.per.ref[[ref.name]]
tmp.check <- is.numeric(ref.tholds)
if(!tmp.check) stop('Error defining reference-specific size thresholds.\n')
tholds.per.ref <- gen.wflow.opts$ucm_tholds
ucm.tholds <- tholds.per.ref[[ref.name]]
if(is.list(ucm.tholds)) ucm.tholds <- unlist(ucm.tholds)
tmp.check <- is.numeric(ucm.tholds)
if(!tmp.check) stop('Error defining reference-specific weight thresholds.\n')
match.opts.per.ref <- gen.wflow.opts$vdj_match_opts
vdj.match.opts <- match.opts.per.ref[[ref.name]]
if(is.list(vdj.match.opts)) vdj.match.opts <- unlist(vdj.match.opts)
vdj.match.opts <- as.logical(vdj.match.opts)
# @ Prediction output.
out.paths <- paste0(run.path, '/VDJMatch-', vdj.match.opts)
out.paths <- lapply(X=out.paths, FUN=function(tmp.path){
    out.paths <- paste0(tmp.path, '/RefThold-', ref.tholds)
    out.paths <- lapply(X=out.paths, FUN=function(tmp.path){
        tmp.paths <- paste0(tmp.path, '/UCMThold-', ucm.tholds, '/IntersectedPredictionDetails.csv')
        names(tmp.paths) <- paste0('Weight-', ucm.tholds)
        return(tmp.paths)
    })
    names(out.paths) <- paste0('Size-', ref.tholds)
    out.paths <- unlist(out.paths)
})
names(out.paths) <- paste0('VDJMatch-', vdj.match.opts)
out.paths <- unlist(out.paths)
tmp.check <- all(file.exists(out.paths))
if(!tmp.check){
    out.paths <- out.paths[!file.exists(out.paths)]
    tmp.err <- paste0('Error attempting to define prediction outputs. Next are the attempted and failed definitions:\n', paste0(out.paths, collapse='\n'))
    stop(tmp.err)
}


### --------------------------- Functions --------------------------- ###


cat('\n\n')
### --------------------------- Load data --------------------------- ###
### ----------------------- and preprocessing ----------------------- ###
cat('### --------------------------- Load data --------------------------- ###\n')
cat('### ----------------------- and preprocessing ----------------------- ###\n')

# ---> Prediction details.
pred.summ <- lapply(X=out.paths, FUN=fread)
pred.summ <- rbindlist(l=pred.summ, use.names=TRUE, fill=TRUE, idcol='opt')
tmp.cols <- c(
    'opt', 
    'clonotype.tag', 'donor.id.tag', 'abs.freq', 'rel.freq',
    'consensus.pred', 'consensus.approach',
    'DM.pred', 'DM.chain.class', 'DM.match.class', 'DMB|consensus.ag', 'DMA|consensus.ag',
    'RM.pred', 'RM.chain.class', 'RM.match.class',
    'UCM.pred', 'UCM.clust.class', 'UCM.ag.class', 'UCM.ratio', 'UCM.max.ags'
)
pred.summ <- pred.summ[, ..tmp.cols]
pred.summ[, vdj.opt:=str_extract(
        string=str_extract(string=opt, pattern='VDJMatch-[^\\.]+'),
        pattern='TRUE$|FALSE$'
    )
]
tmp.vals <- c(`TRUE`='W/ VDJ matches', `FALSE`='W/out VDJ matches')
pred.summ[,
    vdj.opt:=tmp.vals[vdj.opt]
]
pred.summ[, size.thold:=str_extract(
        string=str_extract(string=opt, pattern='Size-\\d+'),
        pattern='\\d+$'
    )
]
pred.summ[, weight.thold:=as.numeric(str_extract(
        string=str_extract(string=opt, pattern='Weight-\\d+(\\.\\d+)*'),
        pattern='\\d+(\\.\\d+)*$'
    ))
]
# "Mock" experimental prediction column added if necessary.
if(!'DM.pred' %in% colnames(pred.summ)) pred.summ[, DM.pred:=NA]


cat('\n\n')
### ------------------------- Main program -------------------------- ###
cat('### ------------------------- Main program -------------------------- ###\n')

cat('\n\n')
### ---------------- Individual assessment for each ----------------- ###
### ---------------------- prediction approach ---------------------- ###
cat('### ---------------- Individual assessment for each ----------------- ###\n')
cat('### ---------------------- prediction approach ---------------------- ###\n')

# ---> Number of predictions provided by each approach post consensus prediction criteria.
tmp.report <- pred.summ[
    !is.na(consensus.approach),
    .(
        clone.count=uniqueN(clonotype.tag),
        cell.count=sum(abs.freq)
    ),
    by=.(consensus.approach, opt, vdj.opt, size.thold, weight.thold)
]
up.bound <- tmp.report[, quantile(x=clone.count, probs=0.975)]
tmp.ggplot.1 <- ggplot(data=tmp.report, aes(x=consensus.approach, y=clone.count)) +
    geom_boxplot(width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0) +
    facet_wrap(facets=~vdj.opt+size.thold, nrow=tmp.report[, uniqueN(vdj.opt)]) +
    scale_y_continuous(expand=c(0, NA)) +
    coord_cartesian(ylim=c(NA, up.bound)) +
    labs(
        title='Clonotype count per consensus approach class',
        x='Consensus approach class',
        y='Number of clonotypes w/ final antigen specificity prediction',
        caption='Different quadrants represent different size threshold options.\nEach dot indicates a different size weight threshold value.\nZoom in on plot\'s y axis. Upper bound established by the 97.5 percentile.'
    ) +
    theme_bw() + theme(axis.text.x=element_text(angle=45))
up.bound <- tmp.report[, quantile(x=cell.count, probs=0.975)]
tmp.ggplot.2 <- ggplot(data=tmp.report, aes(x=consensus.approach, y=cell.count)) +
    geom_boxplot(width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0) +
    facet_wrap(facets=~vdj.opt+size.thold, nrow=tmp.report[, uniqueN(vdj.opt)]) +
    scale_y_continuous(expand=c(0, NA)) +
    coord_cartesian(ylim=c(NA, up.bound)) +
    labs(
        title='Cell count per consensus approach class',
        x='Consensus approach class',
        y='Number of cells w/ final antigen specificity prediction',
        caption='Different quadrants represent different size threshold options.\nEach dot indicates a different size weight threshold value.'
    ) +
    theme_bw() + theme(axis.text.x=element_text(angle=45))
tmp.width <- (2 * tmp.report[, uniqueN(size.thold)]) + 1 
tmp.height <- (2 * tmp.report[, uniqueN(vdj.opt)]) + 3
tmp.file.name <- paste0(reports.path, '/_1_ConsensusSummary.pdf')
pdf(file=tmp.file.name, width=tmp.width, height=tmp.height)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
dev.off()

# ---> Clustering by sequence similarity approach.
# @ Single over multiple ratio.
# Determine ratio, particularly for clustering results.
tmp.data <- pred.summ[!is.na(UCM.clust.class)]
tmp.report <- tmp.data[, 
    .(  
        single.clone.count=.SD[UCM.clust.class=='Single', uniqueN(clonotype.tag)]
    ),
    by=.(opt, vdj.opt, size.thold, weight.thold)
]
tmp.cols <- setdiff(x=colnames(tmp.report), y=c('opt', 'vdj.opt'))
tmp.report <- unique(tmp.report[, ..tmp.cols])
tmp.check <- tmp.report[, .N==uniqueN(weight.thold)*uniqueN(size.thold)]
if(!tmp.check) stop('Faulty reference match overlap results (1).\n')
tmp.ggplot.1 <- ggplot(data=tmp.report, aes(x=weight.thold, y=single.clone.count, group=size.thold)) +
    geom_line(linewidth=1.4, color='#0000b3') +
    geom_point(size=5) +
    facet_wrap(facets=~size.thold, nrow=1) +
    scale_y_continuous(limits=c(0, NA)) +
    labs(
        title='Weight threshold vs Single count',
        x='Weight threshold',
        y='Single clonotype count',
        caption='Different quadrants represent different size threshold options.'
    ) +
    theme_bw()
tmp.width <- (5 * tmp.report[, uniqueN(size.thold)]) + 1 
tmp.file.name <- paste0(reports.path, '/_2_ClusteringApproachAssessment.pdf')
pdf(file=tmp.file.name, width=tmp.width, height=6)
print(tmp.ggplot.1)
dev.off()


# ---> Reference match approach.
# Number of clonotypes provided by this approach for each size threshold.
tmp.data <- pred.summ[!is.na(RM.chain.class)]
tmp.report <- tmp.data[, 
    .(  
        clone.count=uniqueN(clonotype.tag),
        cell.count=sum(abs.freq)

    ),
    by=.(opt, vdj.opt, size.thold, weight.thold, gen.match.class=RM.chain.class)
]
tmp.report[,
    simp.match.class:=gen.match.class
]
tmp.vals <- c('Perfect beta match with alpha quasi match', 'Perfect alpha match with beta quasi match', 'Two-chain perfect match')
tmp.report[
    gen.match.class%in%tmp.vals,
    simp.match.class:='Both'
]
tmp.cols <- setdiff(x=colnames(tmp.report), y=c('opt', 'weight.thold'))
tmp.report <- unique(tmp.report[, ..tmp.cols])
# tmp.check <- tmp.report[, .N==uniqueN(vdj.opt)*uniqueN(size.thold)*uniqueN(gen.match.class)]
# if(!tmp.check) stop('Faulty reference match overlap results (1).\n') # NOTE: Found out this doesn't necessarily hold true, especially when dealing w/ small datasets.
tmp.lvls <- c(
    'Two-chain perfect match',
    'Perfect beta match with alpha quasi match', 'Perfect alpha match with beta quasi match',
    'Beta-only perfect match', 'Alpha-only perfect match'
)
tmp.report$gen.match.class <- factor(x=tmp.report$gen.match.class, levels=tmp.lvls)
tmp.lvls <- c( 'Both', 'Beta-only perfect match', 'Alpha-only perfect match')
tmp.report$simp.match.class <- factor(x=tmp.report$simp.match.class, levels=tmp.lvls)
# Plots.
tmp.ggplot.1 <- ggplot(data=tmp.report, aes(x=gen.match.class, y=clone.count)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
    facet_wrap(facets=~vdj.opt+size.thold, nrow=tmp.report[, uniqueN(vdj.opt)]) +
    scale_y_continuous(expand=c(0, NA)) +
    labs(
        title='Clonotype count per general match class',
        x='General class of match between query and reference clonotypes',
        y='Number of overlapping clonotypes between query and reference TCR sets',
        caption='Different quadrants represent different size threshold options.'
    ) +
    theme_bw() + theme(axis.text.x=element_text(angle=45))
tmp.ggplot.2 <- ggplot(data=tmp.report, aes(x=gen.match.class, y=cell.count)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
    facet_wrap(facets=~vdj.opt+size.thold, nrow=tmp.report[, uniqueN(vdj.opt)]) +
    scale_y_continuous(expand=c(0, NA)) +
    labs(
        title='Cell count per general match class',
        x='General class of match between query and reference clonotypes',
        y='Number of query cells w/ overlapping clonotypes between query and reference TCR sets',
        caption='Different quadrants represent different size threshold options.'
    ) +
    theme_bw() + theme(axis.text.x=element_text(angle=45))
tmp.ggplot.3 <- ggplot(data=tmp.report[gen.match.class!='Alpha'], aes(x=gen.match.class, y=clone.count)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
    facet_wrap(facets=~vdj.opt+size.thold, nrow=tmp.report[, uniqueN(vdj.opt)]) +
    scale_y_continuous(expand=c(0, NA)) +
    labs(
        title='Clonotype count per simplified match class',
        x='Simplified class of match between query and reference clonotypes',
        y='Number of overlapping clonotypes between query and reference TCR sets',
        caption='Different quadrants represent different size threshold options.'
    ) +
    theme_bw() + theme(axis.text.x=element_text(angle=45))
tmp.ggplot.4 <- ggplot(data=tmp.report[gen.match.class!='Alpha'], aes(x=gen.match.class, y=cell.count)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
    facet_wrap(facets=~vdj.opt+size.thold, nrow=tmp.report[, uniqueN(vdj.opt)]) +
    scale_y_continuous(expand=c(0, NA)) +
    labs(
        title='Cell count per simplified match class',
        x='Simplified class of match between query and reference clonotypes',
        y='Number of query cells w/ overlapping clonotypes between query and reference TCR sets',
        caption='Different quadrants represent different size threshold options.'
    ) +
    theme_bw() + theme(axis.text.x=element_text(angle=45))
tmp.width <- (1.5 * tmp.report[, uniqueN(size.thold)]) + 1 
tmp.height <- (2 * tmp.report[, uniqueN(vdj.opt)]) + 3
tmp.file.name <- paste0(reports.path, '/_3_MatchApproachAssessment_1.pdf')
pdf(file=tmp.file.name, width=tmp.width, height=tmp.height)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
dev.off()

# Comparison between VDJ options (if necessary)
if(tmp.data[, uniqueN(vdj.opt)]>1){
    tmp.report <- tmp.data[
        RM.chain.class!='Alpha-only perfect match',
        .(  
            clone.count=uniqueN(clonotype.tag),
            cell.count=sum(abs.freq)

        ),
        by=.(opt, vdj.opt, size.thold, weight.thold)
    ]
    tmp.cols <- setdiff(x=colnames(tmp.report), y=c('opt', 'weight.thold'))
    tmp.report <- unique(tmp.report[, ..tmp.cols])
    tmp.check <- tmp.report[, .N==uniqueN(vdj.opt)*uniqueN(size.thold)]
    if(!tmp.check) stop('Faulty reference match overlap results (1).\n')
    tmp.report.1 <- spread(data=tmp.report[, .(vdj.opt, size.thold, clone.count)], key=vdj.opt, value=clone.count)
    tmp.report.2 <- spread(data=tmp.report[, .(vdj.opt, size.thold, cell.count)], key=vdj.opt, value=cell.count)
    tmp.report <- merge(
        x=tmp.report.1,
        y=tmp.report.2,
        by='size.thold',
        suffixes=c('.clone', '.cell')
    )
    tmp.report$size.thold <- as.numeric(tmp.report$size.thold)
    up.bound <- max(c(
        tmp.report[, 'W/out VDJ matches.clone'], 
        tmp.report[, 'W/ VDJ matches.clone']
    ))
    tmp.ggplot.1 <- ggplot(data=tmp.report, aes(x=`W/out VDJ matches.clone`, y=`W/ VDJ matches.clone`)) +
        geom_point(aes(fill=size.thold), shape=21, size=4, stroke=2) +
        labs(
            title='Clonotype count comparison between VDJ gene-specific approached',
            x='Count disregarding VDJ gene matches',
            y='Count considering VDJ gene matches',
        ) +
        scale_x_continuous(limits=c(0, up.bound)) +
        scale_y_continuous(limits=c(0, up.bound)) +
        theme_bw()
    up.bound <- max(c(
        tmp.report[, 'W/out VDJ matches.cell'], 
        tmp.report[, 'W/ VDJ matches.cell']
    ))
    tmp.ggplot.2 <- ggplot(data=tmp.report, aes(x=`W/out VDJ matches.cell`, y=`W/ VDJ matches.cell`)) +
        geom_point(aes(fill=size.thold), shape=21, size=4, stroke=2) +
        scale_x_continuous(limits=c(0, up.bound)) +
        scale_y_continuous(limits=c(0, up.bound)) +
        labs(
            title='Cell count comparison between VDJ gene-specific approached',
            x='Cell disregarding VDJ gene matches',
            y='Cell considering VDJ gene matches',
        ) +
        theme_bw()
    tmp.file.name <- paste0(reports.path, '/_3_MatchApproachAssessment_2.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    dev.off()
}


# ---> Experimental validation approach.
# Apply if experimental data was available.
if(!pred.summ[, all(is.na(DM.pred))]){
    # @ Comparison between separate chains ("pairing" assessment) for clonotypes overlaping between lung and experimental data.
    #       For all query clonotypes with both chains matched to experimental data clonotypes
    tmp.data <- pred.summ[!is.na(DM.chain.class)]
    tmp.report <- tmp.data[
        DM.chain.class=='Two-chain perfect match',
        .(
            `Specificity match`=.SD[
                `DMB|consensus.ag`==`DMA|consensus.ag`,
                uniqueN(clonotype.tag)
            ],
            `Specificity mismatch`=.SD[
                `DMB|consensus.ag`!=`DMA|consensus.ag`,
                uniqueN(clonotype.tag)
            ]
        ),
        by=.(opt, size.thold, weight.thold)
    ]
    tmp.report <- unique(tmp.report[, .(`Specificity match`, `Specificity mismatch`)])
    tmp.check <- tmp.report[, .N==1]
    if(!tmp.check) stop('Faulty experimental overlap results (1).\n')
    tmp.report <- gather(data=tmp.report, key='key', value='value')
    tmp.ggplot.1 <- ggplot(data=tmp.report, aes(x=key, y=value)) +
        geom_bar(stat='identity', width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
        scale_y_continuous(expand=c(0, NA)) +
        labs(
            title='Specificty accordance between chains from experimental data',
            subtitle='All clonotypes',
            x='Specificity accordance between chains',
            y='Number of query clonotypes with both chains matched to experimental data clonotypes',
            caption='All query clonotypes with both chains matched to experimental data clonotypes were considered.\nSame across pipeline options.' 
        ) +
        theme_bw()
    #       For query clonotypes with both chains that are both 1) matched to experimental data clonotypes and 2) defined w/ a single reactivity for both matched chains.
    tmp.report <- tmp.data[
        DM.chain.class=='Two-chain perfect match' &
        !(str_detect(string=`DMB|consensus.ag`, pattern='\\|') | str_detect(string=`DMA|consensus.ag`, pattern='\\|')),
        .(
            `Specificity match`=.SD[
                `DMB|consensus.ag`==`DMA|consensus.ag`,
                uniqueN(clonotype.tag)
            ],
            `Specificity mismatch`=.SD[
                `DMB|consensus.ag`!=`DMA|consensus.ag`,
                uniqueN(clonotype.tag)
            ]
        ),
        by=.(opt, size.thold, weight.thold)
    ]
    tmp.report <- unique(tmp.report[, .(`Specificity match`, `Specificity mismatch`)])
    tmp.check <- tmp.report[, .N==1]
    if(!tmp.check) stop('Faulty experimental overlap results (1).\n')
    tmp.report <- gather(data=tmp.report, key='key', value='value')
    tmp.ggplot.2 <- ggplot(data=tmp.report, aes(x=key, y=value)) +
        geom_bar(stat='identity', width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
        scale_y_continuous(expand=c(0, NA)) +
        labs(
            title='Specificty accordance between chains from experimental data',
            subtitle='Clonotypes w/ single specificities',
            x='Specificity accordance between chains',
            y='Number of query clonotypes with both chains matched to experimental data clonotypes',
            caption='Only query clonotypes with both chains matched to experimental data clonotypes\nthat had single-reactivity annotates were considered.\nSame across pipeline options.' 
        ) +
        theme_bw()
    #       Final report.
    tmp.file.name <- paste0(reports.path, '/_4_ExperimentalApproachAssessment.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    dev.off()
}


# ---> Full report.
# Mock to indicate to the pipeline that this is the end of this step.
tmp.file.name <- paste0(reports.path, '/BestPredReport.csv')
tmp.cmmd <- paste0('touch ', tmp.file.name)
system(command=tmp.cmmd, intern=FALSE)


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
