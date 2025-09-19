############    -----   R handy code and functions    -----    ############


############    ------------   R handy code    ------------    ############

library(ggplot2)

# ---> Blank themes.
# Blank themes.
blank.complement.1 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_blank(), axis.line=element_blank()) # Default blank.
# ---> For standard plots.
blank.complement.3 <- theme(
    line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none',
    axis.line=element_line(size=2.5),
    axis.ticks=element_line(size=2.5), axis.ticks.length=unit(1, "cm")
)
# No x axis ticks.
blank.complement.3.1 <- theme(
    line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none',
    axis.line=element_line(size=2.5), axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(size=2.5), axis.ticks.length.y=unit(1, "cm")
)
# No axis ticks.
blank.complement.3.2 <- theme(
    line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none',
    axis.line=element_line(size=2.5),
    axis.ticks.x=element_blank(), axis.ticks.y=element_blank()
)
# For donut plots.
blank.complement.4 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none')
# Blank them for the repertoire barplots.
blank.complement.7 <- theme(
  text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none',
  axis.line.x=element_blank(), axis.line.y=element_line(size=1.5),
  axis.ticks.x=element_blank(), axis.ticks.y=element_line(size=1.5),
  axis.ticks.length.y=unit(0.4, units='cm'),
)
blank.complement.7.2 <- theme(
  text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none',
  axis.line.x=element_line(size=1.5), axis.line.y=element_line(size=1.5),
  axis.ticks.x=element_blank(), axis.ticks.y=element_line(size=1.5),
  axis.ticks.length.x=unit(0.4, units='cm'), axis.ticks.length.y=unit(0.4, units='cm'),
)



############    ----------   R handy functions    ---------    ############

# -------------------------------------------------------------------->
# Name: Subset seurat objects
# Dependencies:
#   Seurat and stringr
# Functions dependencies:
#   translate.ids, among others.
# Description ----------------------->
#   Given a seurat object, this program will subset it according to two kind of rules: 1) Rules per featurs and 2) rules per tags.
# Arguments ------------------------->
# seurat.obj - seurat object.
# tags.criteria - Data frame with criteria to discard info from tags in meta data.
# feats.criteria - Data frame with criteria to discard info from features described in meta data.
# Value ----------------------------->
# Subsetted seurat object.
# Function -------------------------->

subset.seurat.obj <- function(seurat.obj, tags.criteria, feats.criteria){
  # ---> Subset seurat object.
  # ---> Assess criteria per feature/tag
  # By tag criteria.
  if(!is.null(tags.criteria)){
    tags.evals <- sapply(X=colnames(tags.criteria), FUN=function(tmp.tag){
      tmp.criteria <- tags.criteria[, tmp.tag]
      to.do <- rownames(tags.criteria)[!is.na(tmp.criteria)]
      tmp.criteria <- tmp.criteria[!is.na(tmp.criteria)]
      tmp.criteria <- str_split(string=tmp.criteria, pattern=";")[[1]]
      # Check if NA value is part of the criteria.
      na.flag <- any(grepl(x=tmp.criteria, pattern='^[\'"]NA[\'"]$', perl=TRUE))
      tmp.criteria <- tmp.criteria[!grepl(x=tmp.criteria, pattern='^[\'"]NA[\'"]$', perl=TRUE)]
      # ---> Get evaluations per cell per tag value that need to be evaluated.
      # Check there's at least one value to evaluate.
      if(length(tmp.criteria)>0){ to.keep <- sapply(X=tmp.criteria, FUN=function(criteria.term) seurat.obj@meta.data[, tmp.tag]==criteria.term); normal.evals <- TRUE} else normal.evals <- FALSE
      # Evaluate NA value if required.
      if(na.flag){
        na.evals <- is.na(seurat.obj@meta.data[, tmp.tag])
        if(!normal.evals) to.keep <- na.evals else to.keep <- cbind(to.keep, na.evals)
      }
      if(to.do=='discard'){
        to.keep <- !to.keep
        # In case many tag values were evaluated.
        if(!is.null(dim(to.keep))) to.keep <- apply(X=to.keep, MARGIN=1, FUN=all, na.rm=TRUE)
      }else{
        # In case many tag values were evaluated.
        if(!is.null(dim(to.keep))) to.keep <- apply(X=to.keep, MARGIN=1, FUN=any, na.rm=TRUE)
      }
      return(to.keep)
    })
    if(!is.null(dim(tags.evals))) tags.evals <- apply(X=tags.evals, MARGIN=1, FUN=all)
  }else{
    tags.evals <- NULL
  }
  # By features criteria.
  if(!is.null(feats.criteria)){
    feats.evals <- sapply(X=colnames(feats.criteria), FUN=function(tmp.feat){
      upper.thold <- feats.criteria['upper', tmp.feat]
      lower.thold <- feats.criteria['lower', tmp.feat]
      if(is.na(upper.thold)) upper.thold <- max(seurat.obj@meta.data[, tmp.feat])
      if(is.na(lower.thold)) lower.thold <- min(seurat.obj@meta.data[, tmp.feat])
      to.keep <- seurat.obj@meta.data[, tmp.feat]>=lower.thold & seurat.obj@meta.data[, tmp.feat]<=upper.thold
      return(to.keep)
    })
    if(!is.null(dim(feats.evals))) feats.evals <- apply(X=feats.evals, MARGIN=1, FUN=all)
  }else{
    feats.evals <- NULL
  }
  # ---> Assess all evaluation over cells.
  if(is.null(feats.evals)|is.null(tags.evals)){
    if(is.null(feats.evals)) all.evals <- tags.evals else all.evals <- feats.evals
  }else{
    all.evals <- cbind(feats.evals, tags.evals)
    all.evals <- apply(X=all.evals, MARGIN=1, FUN=all)
  }
  names(all.evals) <- rownames(seurat.obj@meta.data)
  # ---> Subset seurat object accordingly.
  cells.to.keep <- names(all.evals)[all.evals]
  seurat.obj <- subset(x=seurat.obj, cells=cells.to.keep)
  # ---> Return subsetted seurat object.
  return(seurat.obj)
}


# -------------------------------------------------------------------->
# Name: Add new tags.
# Description ----------------------->
# This is the function that should be called straight from the main program.
#       For approach 1, based on the set of rules provided in the data frame 'new.tags.rules' and the execution of the function 'get.new.tag', this function obtains the new tags and adds them straight to the metadata of the input seurat object, returning such a modified object. This is achieved in three main steps:
#   1. Checking of the rules object to make sure it's a valid one.
#   2. Definition of the new tags (one at a time) through the execution of the function 'get.new.tag'.
#   3. Addition of the new information to the input object's metadata.
#       For approach 2, the function will merge the currently existing metadata with the new metadata provided (argument "new.meta.data") by barcode (subapproach 1) or groups of one or more already existing metadata columns (subapproach 2). To distinguish between one or the other subapproach, the function looks for the column "barcode" in the new metadata. If such a column is defined, then it'll proceed to merge based on the first subappraoch (after making sure that the rest of the new columns aren't redundant compared with what already exists in the current metadata). Otherwise, the function will attempt to find matching columns between the already existing metadata and the new metadata, and will merge according to such matching columns. If no matching columns can be found, an error will be output.
# Arguments ------------------------->
# @ seurat.obj - Input seurat object whose metadata should contain the tags to be used to create the new tags as defined in the rules object.
# @ new.tags.rules - Data frame, whose specificities are defined next. Leave as NULL when you want to proceed with approach 2. Otherwise, it is mandatory to proceed with approach 1.
#   The dataframe describes the rules to merge two or more different tags already defined in the input seurat object, thereby creating a new tag. The data frame should have next columns (in any order):
#     1) tag.names: For each row, describes the name of the new tag to be produced according to the specifications given in the rest of the columns. This is the only case where NAs are discarded before concatenating all values.
#     2) merge: For each row, sequence of already defined tags to be -somehow- combined, which should be seprated by and only by semicolon (e.g., 'virus.tag;donor.id.tag;batch.tag'). If rule equals 'join' or 'add' (see below), then the order given will be the one to use for the output.
#     3) rule: For each row, one of three options as follows:
#         a) join: Tag values are joined/concatenated (in the order provided) and separated by the pattern provided.
#         b) remove: Final tag values for the new tag will be produced by the appending of the columns specified in merge (separated by a dash and given in the order provided in that column), such that the tag values provided in 'pattern' will be discarded and, therefore, the final value for the cells having them will be set to NA. The tag values defined in 'pattern' and to be removed (i.e., to be set to NA) must be defined as valid values for any of the tags defined in 'merge', however, if any these 'pattern' values are defined for multiple tags defined in 'merge' (i.e., they're redundant), the program will fail outputting an error.
#           For example, if you provided as merge argument 'virus.tag;donor.id.tag;batch.tag' and you set pattern to 'Flu;P07', assuming that they're valid values for 'virus.tag' and 'donor.id.tag', respectively, then the final values for the new column will be composed by the appending of 'virus.tag-donor.id-batch.tag' (notice that '-' is 'joining' the classes) and any cell 1) with value 'Flu' and/or 'P07' for the tags 'virus.tag' and 'donor.id.tag', respectively or 2) with value set originally to NA for any of the original tags will be marked by NA in the final column/new tag.
#         c) keep: Same action as for 'remove', but insted of setting NA as final value to the cells being classified as the values in the 'pattern' argument, it will set NA to all values not being classified as such.
#         d) add.end: Concatenates 'pattern' at the end of the concatenation of the 'merge' columns provided (these ones always joined by '-'). That is, same as 'join', but will add the value given to pattern at the end of the result and always uses '-' to join everything.
#         e) add.beg: Same as add.end, but adding 'pattern' argument at the beggining instead.

# @ new.meta.data - Data frame, whose specificities are defined next. Leave as NULL when you want to proceed with approach 1. Otherwise, it is mandatory to proceed with approach 2.
# Value ----------------------------->
# Updated seurat object, that is, including the new tags gotten through the merging of the others based on the rules provided.
#   The dataframe lists new metadata to be merged with the already existing metadata at one of two possible levels/subapproaches: based on barcode or based on groups.
#       Subapproach 1, by barcode. For this approach to be taken, the column "barcode" must be defined in the dataframe and a value must be provided for > 50% of the cells in the currently existing metadata. Then, the current and new metadata dataframes are merged by barcode.
#       Subapproach 2, by groups. For this approach to be applied, at least one column must be matched between the current and new metadata dataframes. Note that for this one to work, no column named "barcode" must appear. Then, the dataframes are merged by these groups.
# Function -------------------------->

add.new.tags <- function(seurat.obj, new.tags.rules=NULL, new.meta.data=NULL, rem.tags=FALSE, check.min.fract=TRUE){

    ### -------------------------- Preflights --------------------------- ###

    # ---> Sanity check. Info for only one of the two possible approaches has to be provided.
    tmp.check <- is.null(new.tags.rules) & is.null(new.meta.data)
    if(tmp.check) stop('No rules or new metadata were indicated for the program to work properly. One and only one of the arguments "new.tags.rules" and "new.meta.data" has to be different thank NULL.\n')
    tmp.check <- !is.null(new.tags.rules) & !is.null(new.meta.data)
    if(tmp.check) stop('One and only one of the arguments "new.tags.rules" and "new.meta.data" has to be different thank NULL. Both of them were provided.\n')

    ### -------------------------- Approach 1 --------------------------- ###

    if(!is.null(new.tags.rules)){
        # ---> Preflights.
        if(!all(c('tag.name', 'merge', 'rule', 'pattern') %in% colnames(new.tags.rules))) stop(paste0('File describing rules to create new tags was provided (path below), but it does not have the format (necessary columns) required.\n'))
        tmp.tags <- unique(unlist(str_split(string=new.tags.rules$merge, pattern=';', simplify=FALSE)))
        if(!all(tmp.tags %in% colnames(seurat.obj@meta.data))){
            tmp.tags <- tmp.tags[!tmp.tags %in% colnames(seurat.obj@meta.data)]
            stop(paste0('File describing rules to create new tags was provided (path below), but any of the tags\' names provided in the column \'merge\' has not been previously defined in the seurat object as listed below:\n', paste0(tmp.tags, collapse='\n'), '\n'))
        }
        if(any(new.tags.rules$tag.name %in% colnames(seurat.obj@meta.data))){
            if(rem.tags){
                tmp.tags <- new.tags.rules$tag.name[new.tags.rules$tag.name %in% colnames(seurat.obj@meta.data)]
                for(tmp.tag in tmp.tags){
                    seurat.obj@meta.data[, tmp.tag] <- NULL
                }
            }else{
            stop('New tag names should not be previosly defined in the suerat object.\n')
            }
        }

        # ---> Get new tags.
        # New tag definition.
        new.tags <- lapply(X=1:nrow(new.tags.rules), FUN=function(idx){
            tag.name <- new.tags.rules[idx, 'tag.name']; merge <- new.tags.rules[idx, 'merge']; rule <- new.tags.rules[idx, 'rule']; pattern <- new.tags.rules[idx, 'pattern']
            new.tag <- get.new.tag(tag.name=tag.name, tags.to.merge=merge, rule=rule, pattern=pattern)
            return(new.tag)
        })
        names(new.tags) <- new.tags.rules$tag.name
        new.tags <- rbindlist(l=new.tags, use.names=TRUE, idcol='tag.name')
        new.tags <- spread(data=new.tags, key=tag.name, val=new.tag)
        new.tags <- as.data.frame(new.tags)
        row.names(new.tags) <- new.tags$cell; new.tags <- new.tags[, new.tags.rules$tag.name]
        new.tags <- new.tags[Cells(seurat.obj), ]
        # Add new tags to previous seurat object.
        meta.data <- merge(seurat.obj@meta.data, new.tags, by='row.names', sort=FALSE)
        row.names(meta.data) <- meta.data$Row.names
        meta.data$Row.names <- NULL
        if(!all(row.names(meta.data)==Cells(seurat.obj))) stop('Unexpected error 1. See program.\n')
        seurat.obj@meta.data <- meta.data
        rm(meta.data)
        rm(new.tags)
    }

    ### -------------------------- Approach 2 --------------------------- ###

    if(!is.null(new.meta.data)){
        # ---> Generalities.
        meta.data <- seurat.obj@meta.data
        meta.data$barcode <- Cells(seurat.obj)
        if('barcode' %in% colnames(new.meta.data)){
        # ---> Subapproach 1
            # Sanity checks. Unique barcodes must be provided in the new metadata and at least 50% of the current barcodes must be defined in the new metadata.
            tmp.check.1 <- nrow(new.meta.data) == length(unique(new.meta.data$barcode))
            tmp.check.2 <- sum(meta.data$barcode %in% new.meta.data$barcode) / nrow(meta.data)
            tmp.check.2 <- if(check.min.fract) tmp.check.2 > 0.5 else TRUE
            tmp.check <- tmp.check.1 & tmp.check.2
            if(!tmp.check) stop('Unexpected error while attempting to merge new metadata by barcode.\n')
            # Avoid redundant columns between the current and new metadata dataframes.
            tmp.cols <- setdiff(x=colnames(new.meta.data), y='barcode')
            tmp.cols <- tmp.cols[!tmp.cols %in% colnames(meta.data)]
            tmp.cols <- c('barcode', tmp.cols)
            new.meta.data <- new.meta.data[, tmp.cols]
            # Final merge.
            meta.data <- merge(x=meta.data, y=new.meta.data, by='barcode', all.x=TRUE, all.y=FALSE, sort=FALSE)
        }else{
        # ---> Subapproach 2
            # Identify matched groups between current and new metadata.
            tmp.cols <- intersect(x=colnames(meta.data), y=colnames(new.meta.data))
            if(length(tmp.cols)==0) stop('New metadata does not contain a "barcode" column nor a set of groups already defined in the current metadata to proceed with merging. Perhaps you meant to create new tags based on a set of rules? (see argument "new.tags.rules").\n')
            # Set columns as character for the input data.
            for(tmp.col in tmp.cols){
                new.meta.data[, tmp.col] <- as.character(new.meta.data[, tmp.col])
            }
            # Final merge
            meta.data <- merge(x=meta.data, y=new.meta.data, by=tmp.cols, all.x=TRUE, all.y=FALSE, sort=FALSE)
        }
        # Update metadata in seurat object.
        tmp.check <- all(meta.data$barcode %in% Cells(seurat.obj)) & all(Cells(seurat.obj) %in% meta.data$barcode)
        if(!tmp.check) stop('Unexpected error 2 while attempting to merge new metadata by approach 2.\n')
        row.names(meta.data) <- meta.data$barcode; meta.data$barcode <- NULL
        seurat.obj@meta.data <- meta.data[Cells(seurat.obj), ]
    }

    ### ---------------------------- Output ----------------------------- ###
    # Return modified seurat object.
    return(seurat.obj)
}


# -------------------------------------------------------------------->
# Name: Get new tag.
# Description ----------------------->
# For approach 1, by set of rules.
# This function is meant to work for every row for every row in the data.frame 'new.tags.rules'. It does retrieve the information from the original seurat object and based on the rules defined in a single row of the data frame (i.e., the rules meant to create one single new tag), it defines that new tag.
# Arguments ------------------------->
# All taken from the file provided as new.tags.file. Seurat object is taken from the main environment. The rest of the arguments are the ones defined in the data frame 'new.tags.rules' (or its equivalent) per row (see documentation above if needed).
# tag.name
# tags.to.merge - equvalent to 'merge'
# rule
# pattern
# Value ----------------------------->
# Data table listing the new tag's value for each cell.
# Function -------------------------->

get.new.tag <- function(tag.name, tags.to.merge, rule, pattern){

  ### ---------------------- Data preprocessing ----------------------- ###

  if(!any(rule %in% c('join', 'remove', 'keep', 'add.end', 'add.beg'))) stop(paste0('Invalid rule was provided in the file describing rules to create new tags.\nInvalid rule: ', rule, '\n\n'))
  # Get tags' names.
  tags.to.merge <- unlist(str_split(string=tags.to.merge, pattern=';'))
  if(!length(unique(tags.to.merge)) == length(tags.to.merge)) stop(paste0('For new tag named as ', tag.name, ', please provide unique tag names to merge.\n'))


  ### ------------------------- Main program -------------------------- ###

  # Create new tag accordingly.
  if(rule=='join'){
    # ---> Join.
    to.output <- tidyr::unite(data=seurat.obj[[tags.to.merge]], col=new.tag, sep=pattern, remove=TRUE)
  }else{
    if(rule %in% c('remove', 'keep')){
    # ---> Remove or keep.
      # Check tag values provided as pattern are valid in the sense that they are not redundant and are coming from a valid tag (notice that redundant tag values are not supported). Then, order pattern values according to tag origin.
      pattern <- unique(unlist(str_split(string=pattern, pattern=';')))
      tag.vals <- lapply(X=pattern, FUN=function(tmp.val){
        tmp.data <- seurat.obj@meta.data[, tags.to.merge]==tmp.val & !is.na(seurat.obj@meta.data[, tags.to.merge])
        tmp.data <- colSums(tmp.data)
        tmp.data <- names(tmp.data)[tmp.data>1]
        tmp.check <- length(tmp.data) > 1
        if(tmp.check){ tmp.err <- paste0('Tag value \'', tmp.val, '\' was found to be defined for multiple tags, mainly: ', paste0(tmp.data, collapse=', '), '. Such a redundancy cannot be resolved in the current version of the program. Please develope a new version to deal with it. Alternatively, try and rearrange the rules so that tag values are not redundant.\n'); stop(tmp.err) }
        tmp.check <- length(tmp.data) < 1
        if(tmp.check){ tmp.err <- paste0('Tag value \'', tmp.val, '\' is not valid for any of the tags defined for ', tag.name, ', the tag to be created. Please check it is a valid tag value.\n'); stop(tmp.err) }
        return(tmp.data)
      })
      tag.vals <- unlist(tag.vals)
      names(tag.vals) <- pattern
      tag.names <- unique(tag.vals); tag.vals <- lapply(X=tag.names, FUN=function(tmp.name) names(tag.vals)[tag.vals==tmp.name]); names(tag.vals) <- tag.names
      # Pattern values are then replaced by NA for any cell classified (for 'remove') or not classified (for 'keep') with those values.
      cells.to.modify <- lapply(X=names(tag.vals), FUN=function(tmp.tag){
        pattern.vals <- tag.vals[[tmp.tag]]
        to.output <- seurat.obj@meta.data[, tmp.tag] %in% pattern.vals
        if(rule=='keep') to.output <- !to.output
        to.output <- data.table(cell=Cells(seurat.obj), keep=to.output)
        return(to.output)
      }); names(cells.to.modify) <- names(tag.vals)
      cells.to.modify <- rbindlist(l=cells.to.modify, use.names=TRUE, idcol='tag')
      cells.to.modify <- spread(data=cells.to.modify, key=tag, value=keep)
      cells.to.modify <- apply(X=as.matrix(cells.to.modify, rownames='cell'), MARGIN=1, FUN=if(rule=='keep') any else all)
      # Proceed merging columns.
      to.output <- tidyr::unite(data=seurat.obj[[tags.to.merge]], col=new.tag, sep='-', remove=TRUE)
      # Set NAs according to the rules provided by 'remove' or 'keep'
      to.output[cells.to.modify[row.names(to.output)], 'new.tag'] <- NA
    }else{
    # ---> Add.
      to.output <- tidyr::unite(data=seurat.obj[[tags.to.merge]], col=new.tag, sep='-', remove=TRUE)
      if(rule=='add.end'){
        to.output$new.tag <- paste0(to.output$new.tag, pattern)
      }else{
        to.output$new.tag <- paste0(pattern, to.output$new.tag)
      }
    }
  }
  # Define in a established format and output.
  to.output$cell <- row.names(to.output)
  to.output <- as.data.table(to.output[, c('cell', 'new.tag')])
  # to.output[, tag.name:=tag.name]
  return(to.output)
}


# -------------------------------------------------------------------->
# Name: Tranform current tags.
# Dependencies:
#   Seurat, data.table and stringr
# Functions dependencies:
#   None.
# Description ----------------------->
#   The idea and main code body for this script comes from: /home/vfajardo/scripts/seurat_analysis/tag_transfer_analysis.1.0.R
#   This script is built to somehow transform the info defined for a given tag previously defined in a seurat object (input), this based on the rules in file stored to argument 'transform.rules.file'. This file must contain the next columns:
#   1) ref.tag - Name of a tag as defined in the input seurat object's metadata.
#   2) new.tag - Name of the new tag to store the original info in a modified way. If this value is the same as the original tag name, no action will be taken.
#   3) tag.vals - Tag values to be embedded into a single 'transformed' tag value (also called class). Values of a single class must be separated by a semicolon, while the values themselves must be separated by a '/'. For example, if we want two groups of values to be embedded into single classes, we may specify: '0/2/5;1/3/4/6/7' indicating that values 0, 2 and 5 will belong to the new class 1, while class 2 will be composed of values 1, 3, 4, 6 and 7. Notice that the names of the new classes must be provided in an independent column of this file (see below). Of note, all of the tag values defined in the original tag that are not specified for any class in this column won't be taken into account for any of the new classes and their values will be NA'd.
#   4) tag.classes - Names of the new classes to be defined post-aggregation of the original values. They should simply be separated by a semicolon and the number of names defined here must equal the number of groups defined in column 'tag.vals'. For example, if two new classes were defined in column 'tag.vals', their names may be defined in this column as 'class.1;class.2'.
# Arguments ------------------------->
#   seurat.obj - Seurat object.
#   transform.rules - Data frame with the characteristics described above.
# Value ----------------------------->
# Updated seurat object, that is, including the new tags gotten through the transformation of the others based on the rules provided.
# Function -------------------------->

transform.tags <- function(seurat.obj, transform.rules){

  ### ---------------------- Data preprocessing ----------------------- ###

  # ---> Objects' metadata.
  # Reference object.
  meta.data <- as.data.table(seurat.obj@meta.data)
  meta.data[, cell:=Cells(seurat.obj)]

  # ---> File listing rules to tranform incoming tags.
  # Confirm all necesary columns are defined.
  tmp.cols <- c('ref.tag', 'new.tag', 'tag.vals', 'tag.classes')
  tmp.cols <- tmp.cols[!tmp.cols %chin% colnames(transform.rules)]
  if(length(tmp.cols)>0){
    tmp.err <- paste0(tmp.cols, collapse='\n')
    tmp.err <- paste0('Next columns not appropriately defined in the rules file\n', tmp.err, '\n')
    stop(tmp.err)
  }
  # Confirm all reference tags are defined in seurat object's metadata.
  to.check <- all(transform.rules$ref.tag %chin% colnames(meta.data))
  if(!to.check){
    tmp.err <- paste0(unique(transform.rules$ref.tag[!transform.rules$ref.tag %chin% colnames(meta.data)]), collapse='\n')
    tmp.err <- paste0('Next tags listed in the rules file could be found in the reference seurat object\'s metadata:\n', tmp.err, '\n')
    stop(tmp.err)
  }
  # Confirm there are unique sets of rules for new tag values.
  transform.rules <- unique(x=transform.rules)
  if(length(unique(transform.rules$new.tag))!=nrow(transform.rules)){
    stop('Not unique sets of rules defined for new tags to be transferred on to the query seurat object. Check you provide unique names for the new tags.\n')
  }

  ### ------------------------- Main program -------------------------- ###

  # ---> Tranform tags to get new ones.

  # Get unique indexes and iterate over them to modify one tag at a time.
  idxs <- 1:nrow(transform.rules)
  for(tmp.idx in idxs){
    # Extract tag names as defined and to be defined.
    ref.tag <- transform.rules[tmp.idx, 'ref.tag']
    new.tag <- transform.rules[tmp.idx, 'new.tag']
    if(ref.tag==new.tag) next
    # Define old tag in the reference meta data with an intuitive name in a universal format.
    meta.data[, old.tag:=as.character(get(ref.tag))]
    meta.data[, new.tag:=NULL]
    # Extract tags values and names and confirm their sizes are the same for a set of rules.
    #   Tag values.
    tag.vals <- transform.rules[tmp.idx, 'tag.vals']
    tag.vals <- str_split(string=tag.vals, pattern=';', simplify=FALSE)[[1]]
    tag.vals <- lapply(X=tag.vals, FUN=function(x) str_split(string=x, pattern='/', simplify=FALSE)[[1]])
    #   Tag classes.
    tag.classes <- transform.rules[tmp.idx, 'tag.classes']
    tag.classes <- str_split(string=tag.classes, pattern=';', simplify=FALSE)[[1]]
    #   Confirm size equality.
    to.check <- length(tag.vals)==length(tag.classes)
    if(!to.check){
      tmp.err <- paste0('Not the same amount of classes defined as the amount of groups of values provided for set of rules to go from ', ref.tag, ' to ', new.tag, '.\n')
      stop(tmp.err)
    }else{
      names(tag.vals) <- tag.classes
    }
    # Transform reference tag.
    if(!any(is.na(tag.vals))){
      for(new.class in names(tag.vals)){
        old.classes <- tag.vals[[new.class]]
        meta.data[old.tag%chin%old.classes, new.tag:=new.class]
      }
    }else{
      meta.data[, new.tag:=old.tag]
    }
    # Rename new column.
    tmp.cols <- colnames(meta.data)
    tmp.cols[tmp.cols=='new.tag'] <- new.tag
    colnames(meta.data) <- tmp.cols
  }

  # ---> Tag transfer back to seurat object.
  # Tidy new metadata object and check cell order is kept.
  meta.data[, old.tag:=NULL]
  tmp.check <- all(meta.data[, cell]==Cells(seurat.obj))
  if(!to.check) stop('Unexpected error for "transform.tags". Something went wrong while tranforming current tags. Check function script.\n')
  meta.data <- as.data.frame(meta.data); row.names(meta.data) <- meta.data$cell; meta.data$cell <- NULL
  # Transfer metadata back to seurat object.
  seurat.obj@meta.data <- meta.data
  # Return updated seurat object.
  return(seurat.obj)
}



# -------------------------------------------------------------------->
# Name: Module scoring.
# Dependencies:
#   Seurat and stringr, ggplot2, stringr, data.table, R_handy_functions
# Functions dependencies:
#   translate.ids, among others.
# Description ----------------------->
# You can find an extended explanation on the AddModuleScore function from Seurat in the script associated to this module function.
# This program takes a seurat object and a directory where a set of files listing gene signatures should be and calculates module score for each signature.
# Arguments ------------------------->
# seurat.obj - seurat object.
# module.feats.dir - String for the absolute path to the directory storing one or more gene-signature lists. Each file should have a csv format (one column which must be named 'feature') with rows listing the set of features (usually gene names) that define the module and should be part of the seurat object; if they're not, they will be filtered out. Each file name must follow the next pattern '<Signature name -spaced substituted by underscores->_signature.csv
# is.ensembl - Logical indicating whether input seurat object has ENSEMBL IDs defined as gene IDs instead of common genen names.
# reports.path - Absolute path to where filtered gene signatures should be stored.
# Value ----------------------------->
# Seurat object with module scores calculated.
# Function -------------------------->

get.module.scores <- function(seurat.obj, module.feats.dir, is.ensembl=FALSE, reports.path=NULL){
  cat('\n\n')
  cat('############    -----------   Module Scoring    -----------    ############\n')
  cat('############    -Based on module-defining genes expression-    ############\n')

  ### --------------------------- Arguments --------------------------- ###
  this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # provided by Ciro.

  cat('\n\n')
  ### --------------------------- Load data --------------------------- ###
  cat('### --------------------------- Load data --------------------------- ###\n')
  # Module-defining genes file.
  module.feats.files <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
  modules.feats <- lapply(X=module.feats.files, FUN=function(tmp.file){
    tmp.feats <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
    if(!'feature' %in% colnames(tmp.feats)) stop(paste0('File ', tmp.file, ' not appropriately defined. Column listing the gene features should be named \'feature\'.\n\n'))
    tmp.feats <- tmp.feats$feature
    return(tmp.feats)
  })
  names(modules.feats) <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
  names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
  names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
  modules.names <- names(modules.feats)
  cat('Module-defining features loaded...\n')

  ### ----------------------- Module definition ----------------------- ###
  cat('### ----------------------- Module definition ----------------------- ###\n')
  data.path <- if(!is.null(reports.path)) paste0(reports.path, '/data') else NULL
  if(!is.null(data.path)) if(!dir.exists(data.path)) dir.create(data.path)

  # ---> Module-defining genes.
  # Make sure that all module-defining genes are part of the seurat object.
  cat('Filtering non-defined features...\n')
  modules.feats <- lapply(X=names(modules.feats), FUN=function(tmp.module){
    module.feats <- modules.feats[[tmp.module]]
    # Translate IDs if necessary.
    if(is.ensembl) module.feats <- translate.ids(ids=module.feats, ensembl=FALSE) else names(module.feats) <- module.feats
    # Define IDs as a data frame.
    module.feats <- module.feats[!is.na(names(module.feats))]
    module.feats <- data.frame(features=names(module.feats), original.features=module.feats, present=(names(module.feats) %in% rownames(seurat.obj)), stringsAsFactors=FALSE)
    module.feats <- module.feats[!is.na(module.feats$features), ]
    # Make sure there are any features present to continue.
    tmp.check <- all(!module.feats$present)
    if(tmp.check){
      tmp.warning <- paste0('There are not enough features from module ', tmp.module, '. It will be skipped.\n')
      warning(tmp.warning)
      return(NA)
    }
    # Then, output final list of module-defining genes.
    if(!is.null(data.path)){
      tmp.file.name <- paste0(data.path, '/', str_replace_all(string=toupper(tmp.module), pattern='\\.', replacement='_'), '_ModuleDefiningGenesInExpressionData.csv')
      write.csv(x=module.feats, file=tmp.file.name, quote=FALSE, row.names=FALSE)
    }
    # Check at least 90% of the original siganture features are present in the seurat object.
    if(sum(module.feats$present)/nrow(module.feats) < 0.85) warning(paste0('\nFor module/signature ', tmp.module, ', we couldn\'t recover for the seurat object more than 85% of genes in the original list. Something may be wrong either with the list or the seurat object. Please be careful about that when interpreting the results.\n'))
    # Output final features.
    return(module.feats$features[module.feats$present])
  })
  modules.names <- modules.names[!is.na(modules.feats)]
  modules.feats <- modules.feats[!is.na(modules.feats)]
  names(modules.feats) <- modules.names
  cat('List of module-defining features depicting presence in the gene expression data output.\n')

  cat('\n\n')
  ### -------------------- Module score inference --------------------- ###
  cat('### -------------------- Module score inference --------------------- ###\n')
  # ---> Module score.
  # Add module score across cells.
  seurat.obj <- AddModuleScore(object=seurat.obj, features=modules.feats, name=modules.names)
  # Fix module names given by seurat by removing the suffix number.
  modules.names <- paste0(modules.names, '.score')
  names(modules.names) <- paste0(str_replace(string=modules.names, pattern='\\.score', replacement=''), 1:length(modules.names))
  for(tmp.module in names(modules.names)){
    colnames(seurat.obj@meta.data) <- str_replace(string=colnames(seurat.obj@meta.data), pattern=tmp.module, replacement=modules.names[tmp.module])
  }
  names(modules.names) <- str_replace(string=names(modules.names), pattern='\\d+$', replacement='')
  names(modules.names) <- str_replace_all(string=names(modules.names), pattern='\\.', replacement='_')

  return(seurat.obj)
}


# ---------------------------------------->
# Name: Get publication figure.
# Description ----------------------->
# Provided a ggplot, the name of a file to save it to and any further pdf function arguments, this function will output the publication quality pdf of such a plot in both versions: 1) 'Complete' version: version with all labels; and 2) 'Blank' version: version to be used in the actual paper figures file. If requested, this function will also provide a file with the plot scale(s).
# Arguments ------------------------------>
# tmp.ggplot - ggplot object to be output.
# output.path - Character, absolute path to directory to save file to.
# file.name - Character, prefix to be given to all types of file names provided by this function.
# type - Character, indicates the extension of the file, either "pdf" or "tiff".
# do.legend - Logical, whether to output an extra pdf file showing only the legend scales (e.g., color and size). Note that legend is output as a pdf file regardless of whether the main plots are requested to be output as tiff files.
# legend.height, legend.width - Integer, height and width, respectively, for legend file. Default values have been fine-selected for a file where a single color scale is to be output.
# stat.cor - Logical, indicates whether a scatter plot should include a regression line and its corresponding statistics.
# cor.group - Character, relevant only when stat.cor is TRUE. Indicates the groups to apply regression to.
# repel.data, repel.label - TBD.
# Function -------------------------->

publish.plot <- function(
    tmp.ggplot, output.path, file.name, type='pdf',
    blank.comp, comp.comp=theme_bw(),
    do.rotate=FALSE,
    do.legend=FALSE,
    c.height=NULL, c.width=NULL,
    legend.height=1.5, legend.width=0.5,
    stat.cor=FALSE, cor.group=NULL, cor.method='pearson',
    repel.data=NULL, repel.label=NULL,
    ...
){
  # Blank version.
  tmp.file.name <- if(type=='pdf') paste0(output.path, '/', file.name, '.B.pdf') else paste0(output.path, '/', file.name, '.B.tiff')
  # if(type=='pdf') pdf(file=tmp.file.name, ...) else tiff(file=tmp.file.name, height=1400, width=1400)
  if(type=='pdf') pdf(file=tmp.file.name, ...) else tiff(file=tmp.file.name, units='in', res=600, compression='lzw', ...)
  print(tmp.ggplot + blank.comp)
  dev.off()
  # Complete version.
  tmp.ggplot <- if(stat.cor){
    if(!is.null(cor.group)) tmp.ggplot + stat_cor(aes_string(group=cor.group), method=cor.method) else tmp.ggplot + stat_cor(method=cor.method)
  }else{
    to.check <- !is.null(repel.data) & !is.null(repel.label)
    if(to.check){
      tmp.ggplot + ggrepel::geom_text_repel(data=repel.data, aes_string(label=repel.label), color='black')
    }else{
      tmp.ggplot
    }
  }
  tmp.file.name <- if(type=='pdf') paste0(output.path, '/', file.name, '.C.pdf') else paste0(output.path, '/', file.name, '.C.tiff')
  tmp.ggplot <- if(do.rotate) tmp.ggplot + comp.comp + theme(axis.text.x=element_text(angle=45)) else tmp.ggplot + comp.comp
  if(type=='pdf'){
    if(any(is.null(c(c.height, c.width)))){
      pdf(file=tmp.file.name, ...)
    }else{
      pdf(file=tmp.file.name, height=c.height, width=c.width)
    }
  }else{
    # tiff(file=tmp.file.name, height=1400, width=1400)
    tiff(file=tmp.file.name, units='in', res=600, compression='lzw', ...)
  }
  print(tmp.ggplot)
  dev.off()
  # Output legend independently.
  if(do.legend){
    tmp.ggplot <- get_legend(
        p=tmp.ggplot + theme(
            legend.background=element_blank(), legend.key=element_blank(), legend.text=element_blank(), legend.title=element_blank(),
            legend.position='right',
            legend.frame=element_rect(linewidth=0.5),
            legend.ticks=element_line(color='black', linewidth=0.5), legend.ticks.length=unit(0.5, "cm")
        )
    )
    tmp.file.name <- paste0(output.path, '/', file.name, '.L.pdf')
    pdf(file=tmp.file.name, heigh=legend.height, width=legend.width)
    tmp.ggplot <- as_ggplot(tmp.ggplot)
    print(tmp.ggplot)
    dev.off()
  }
}