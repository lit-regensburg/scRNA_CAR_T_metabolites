# Get root directories for pipeline

get_root_directories <- function(mode = "env") {
  root_directories <- list()

  # Use environment variables to get root directories
  if (mode == "env") {
    env_variables <- c(
      "analysis" = "ANALYSIS_DIR",
      "project" = "PROJECT_DIR",
      "data" = "TOPLEVEL_PROJECT_DIR"
    )

    for (root_dir in names(env_variables)) {
      root_directories[[root_dir]] <- Sys.getenv(env_variables[root_dir])

      if (root_directories[[root_dir]] == "") {
        stop(paste0(
          "The environment variable ", env_variables[root_dir],
          " is undefined, unable to get path for ", root_dir, " root directory"
        ))
      }
    }
  # Use rprojroot to get root directories
  } else if (mode == "rprojroot") {
    analysis_root_criterion <- ".analysisroot"
    project_root_criterion <- ".projectroot"
    data_root_criterion <- ".dataroot"

    project_root <- rprojroot::has_file(project_root_criterion)
    project_dir <- project_root$find_file()

    # Try to find analysis root. If there is none, use project root as analysis root

    analysis_dir <- tryCatch({
      analysis_root <- rprojroot::has_file(analysis_root_criterion)
      analysis_dir <- analysis_root$find_file()
      analysis_dir
    }, error = function(e) {
      analysis_dir <- project_dir
      message("Analysis root directory not explicitly set, using project root directory instead")
      analysis_dir
    })

    # Try to find data root. If there is none, use project root as data root

    data_dir <- tryCatch({
      data_root <- rprojroot::has_file(data_root_criterion)
      data_dir <- data_root$find_file()
      data_dir
    }, error = function(e) {
      data_dir <- project_dir
      message("Data root directory not explicitly set, using project root directory instead")
      data_dir
    })

    root_directories[["analysis"]] <- analysis_dir
    root_directories[["project"]] <- project_dir
    root_directories[["data"]] <- data_dir
  } else {
    stop(paste0("Invalid mode ", mode, " supplied"))
  }

  message('Analysis root directory: ', root_directories[["analysis"]])
  message('Project root directory: ', root_directories[["project"]])
  message('Data root directory: ', root_directories[["data"]])

  return(root_directories)
}

find_input_files = function(
  input_path, file_pattern = NULL, dir_pattern = NULL,
  cellranger_output = TRUE, ocm_output = FALSE, samples_to_exclude = NULL
) {
  if(is.null(file_pattern) & is.null(dir_pattern)){
    stop('At least one of `file_pattern` or `dir_pattern` must be specified')
  }
  #test if sample names are unique
  valid_names = function(snames){
    all(!anyDuplicated(snames) & grepl('^[A-Za-z0-9]+.*$', snames))
  }
  #remove "_library*"
  remove_pipe_additions = function(sample_names,added_pattern='_library_.*$'){
    sample_names_out = sub(added_pattern,"",sample_names)
    sample_names_out
  }

  #limit length of the sample names to `max_length` characters
  trim_names = function(sample_names, max_length = 20){
    maxl = max(sapply(sample_names, nchar))
    if(maxl <= max_length) return(sample_names)
    sample_names_trimmed = strtrim(sample_names, max_length)
    if(valid_names(sample_names_trimmed)) return(sample_names_trimmed)
    sample_names_trimmed = sapply(sample_names, function(x) intToUtf8(rev(head(rev(utf8ToInt(x)), max_length))))
    sample_names_trimmed
  }

  #replace underscores with dash (Seurat does not allow _, which will lead to problems in sample names)
  replace_underscore = function(sample_names) {
    gsub('_', '-', sample_names)
  }

  exclude_samples <- function(sample_names, samples_to_exclude) {
    for (samp in samples_to_exclude) {
        sample_match <- sample_names == samp
        
        if (!sum(sample_match)) {
          stop(paste0("Sample ", samp, ", which should be excluded, was not found in the input directory"))
        } else if (sum(sample_match) > 1) {
          stop(paste0("Sample ", samp, ", which should be excluded, was found multiple times in the input directory"))
        } else {
          sample_names <- sample_names[!sample_match]
        }
    }

    return(sample_names)
  }

  #get absolute paths to fragments.tsv.gz files
  if(!is.null(file_pattern) & !is.null(dir_pattern)){
    #search for files with `file_pattern`
    inputFiles_abs = gsub('//','/',list.files(input_path, pattern = file_pattern, full.names = T, recursive = T))
    #filter results for those files that additionally have `dir_pattern` in the file parent directory
    inputFiles_abs = inputFiles_abs[grepl(dir_pattern, basename(dirname(inputFiles_abs)))]

    if (length(inputFiles_abs) == 0) {
      stop(paste0(
        "No input files found for directory pattern '", dir_pattern, "' and file pattern '",
        file_pattern, "' in input directory ", input_path
      ))
    }
  }else if(!is.null(file_pattern)){
    inputFiles_abs = gsub('//','/',list.files(input_path, pattern = file_pattern, full.names = T, recursive = T))
    if(grepl("scTCR",inputFiles_abs[1])){
        inputFiles_abs = inputFiles_abs[grep('SC_MULTI_CORE', inputFiles_abs,invert=T)]}

    if (length(inputFiles_abs) == 0) {
      stop(paste0(
        "No input files found for file pattern '", file_pattern, "' in input directory ", input_path
      ))
    }
  }else if(!is.null(dir_pattern)){
    inputFiles_abs = gsub('//','/', list.dirs(input_path, full.names = T, recursive = T))
    inputFiles_abs = inputFiles_abs[grepl(dir_pattern, basename(inputFiles_abs))]

    if (length(inputFiles_abs) == 0) {
      stop(paste0(
        "No input files found for directory pattern '", dir_pattern, "' in input directory ", input_path
      ))
    }
  }

  #convert to relative paths

  inputFiles_rel = gsub(input_path, "", inputFiles_abs) #gsub('//','/',list.files(input_path, pattern = file_pattern, full.names = F, recursive = T))
  if(all(dirname(inputFiles_rel)=="/outs")){
      inputFiles_rel = gsub(dirname(input_path), "", inputFiles_abs)
  }
  
  #Split into parts
  input_split = strsplit(inputFiles_rel, '/')
  names_are_ok = FALSE
  #1: is file path standard cellranger output, where names are above `outs` folder?
  if(cellranger_output && ocm_output){
    if (
      all(sapply(input_split, function(x) 'outs' %in% x)) &&
      all(sapply(input_split, function(x) 'per_sample_outs' %in% x))
    ) {
      lib_names_raw = sapply(input_split, function(x) x[which('outs' == x) - 1])
      ocm_names_raw = sapply(input_split, function(x) x[which('per_sample_outs' == x) + 1])
      sample_names_raw = paste0(lib_names_raw, "__", ocm_names_raw)
      sample_names = remove_pipe_additions(sample_names_raw)
      sample_names = replace_underscore(sample_names)
      names_are_ok = valid_names(sample_names)
    }
  } else if(cellranger_output){
    if(all(sapply(input_split, function(x) 'outs' %in% x))){
      sample_names_raw = sapply(input_split, function(x) x[which('outs' == x) -1])
      sample_names = remove_pipe_additions(sample_names_raw)
      sample_names = replace_underscore(sample_names)
      sample_names = try(trim_names(sample_names))
      names_are_ok = valid_names(sample_names)
    }
  }
  if(!names_are_ok){
    #2: are sample names a substring of the fragments.tsv.gz files?
    if(is.null(file_pattern)) file_pattern = ''
    sample_names_raw = sapply(input_split, function(x) gsub(file_pattern, '', x[length(x)]))
    sample_names = replace_underscore(sample_names_raw)
    sample_names = try(trim_names(sample_names))
    names_are_ok = valid_names(sample_names)
  }
  if(!names_are_ok){
    #3: is there a valid name in any of the directory folders?
    maxl = max(sapply(input_split, length))
    i = 1
    while(i <= maxl){
      sample_names_raw = sapply(input_split, '[[', i)
      sample_names = replace_underscore(sample_names_raw)
      sample_names = try(trim_names(sample_names))
      names_are_ok = valid_names(sample_names)
      if(names_are_ok) break
      i = i + 1
    }
  }
  if(!names_are_ok){
    #4: if all is FALSE, create new sample names
    sample_names_raw = paste('Sample-', seq_along(inputFiles_abs))
    sample_names = sample_names_raw
  }

  names(inputFiles_abs) = sample_names

  if (!is.null(samples_to_exclude)) {
    names(sample_names) = sample_names_raw
    samples_to_include = exclude_samples(sample_names_raw, samples_to_exclude)
    sample_names = sample_names[samples_to_include]
    inputFiles_abs = inputFiles_abs[sample_names]

    if (length(inputFiles_abs) == 0) {
      stop(paste0("No samples left after excluding samples ", paste0("'", samples_to_exclude, "'", collapse = ", ")))
    }
  }

  inputFiles_abs
}

setIfNull = function(x, new_val){
  #set a value if x is null
  if(is.null(x)){
    x = new_val
  }
  x
}

has_valid_class = function(param, valid){
  #test if the `param` object has a valid class specified in `valid`
  if(is.null(param)){
    return(FALSE)
  }
  class(param) %in% valid
}

dircreatefun <- function(x, out, analysis_root) {
  #create a directory if it doesnt exist
  if(!grepl(out,x)) {
    relpath <- file.path(out,x)
  } else {
    relpath <- x
  }
  if (!dir.exists(file.path(analysis_root, relpath))) {
    dir.create(file.path(analysis_root, relpath))
  }
  return(relpath)
}

## Beta and alpha diversity:

clustercontrib_fun<-function(obj,clusterCol="seurat_clusters",spCol,cellnum=FALSE){
        # obj: Seurat Object
        # clusterCol: charcter vector specifing the cluster colname
        # spCol: character vector specifing the category to count for
        # https://stackoverflow.com/questions/59178124/problem-with-count-dplyr-within-function
        # solution if: With tidyverse, we can make use of {{}} if we pass an unquoted argument
        #          else if: we need to pass a quoted string as variable name, convert it to symbol with ensym and evaluate (!!)
        a=clusterCol
        if(is.na(a)){
            message(" clusterCol is NA")}
        b=spCol
        if(is.na(b)){
            message(" spCol is NA")}
        
        var1 <- rlang::ensym(clusterCol)
        var2 <- rlang::ensym(spCol)
        dat<-obj@meta.data[,c(clusterCol,spCol)]
        require(tidyverse)
        ## tidyverse code
          datt<-as_tibble(dat)
          df<-datt %>%
          dplyr::count(!!var1,!!var2)
          dfo<-as.data.frame(pivot_wider(df,names_from=!!var2,values_from=n))
        ##
        rownames(dfo)<-dfo[,1]
        dfo<-dfo[,-1]
        if(is.factor(dat[1,spCol])){
        ordvec<-levels(dat[,spCol])[which(levels(dat[,spCol])%in%colnames(dfo))]
        }else{
        ordvec<-unique(dat[,spCol])[which(unique(dat[,spCol])%in%colnames(dfo))]
        }
        dfo<-dfo[,ordvec]
        dfoprc<-(apply(dfo,2,function(x) round(x*100/rowSums(dfo,na.rm=T),2)))
        if(cellnum){
            dfout<-sapply(colnames(dfo),function(x) paste(dfo[,x],paste(dfoprc[,x],"%",sep=""),sep="\n"))
        }else{
            dfout<-dfoprc
            dfout[is.na(dfout)]<-0
        }
        return(dfout)
}



cluster_alphaDiv_fun<-function(obj,clusterCol="seurat_clusters",spCol,cellnum=FALSE){
        # similar to alpha Div in Microbioligy, see the composition of the sample in terms of cell types/clusters
        # obj: Seurat Object
        # 
        # clusterCol: charcter vector specifing the cluster colname
        # spCol: character vector specifing the category to count for
        # https://stackoverflow.com/questions/59178124/problem-with-count-dplyr-within-function
        # solution if: With tidyverse, we can make use of {{}} if we pass an unquoted argument
        #          else if: we need to pass a quoted string as variable name, convert it to symbol with ensym and evaluate (!!) 
        
        a=clusterCol
        if(is.na(a)){
            message(" clusterCol is NA")}
        b=spCol
        if(is.na(b)){
            message(" spCol is NA")}
        
        var1 <- rlang::ensym(clusterCol)
        var2 <- rlang::ensym(spCol)
        dat<-obj@meta.data[,c(clusterCol,spCol)]
        require(tidyverse)
        ## tidyverse code
          datt<-as_tibble(dat)
          df<-datt %>%
          dplyr::count(!!var1,!!var2)
          dfo<-as.data.frame(pivot_wider(df,names_from=!!var2,values_from=n))
        ##
        
        rownames(dfo)<-dfo[,1]
        dfo<-dfo[,-1]
        if(is.factor(dat[1,spCol])){
            ordvec<-levels(dat[,spCol])[which(levels(dat[,spCol])%in%colnames(dfo))]
        }else{
            ordvec<-unique(dat[,spCol])[which(unique(dat[,spCol])%in%colnames(dfo))]
        }    
        dfo<-dfo[,ordvec]
        dfoprc<-apply(dfo,2,function(x) round(x*100/sum(x,na.rm=T),2))
        if(cellnum){
            dfout<-dfo
        }else{
            dfout<-dfoprc
            dfout[is.na(dfout)]<-0
        }
        return(dfout)
}

## plot alphadiv

plottingDivalpha<-function(df){
    require(RColorBrewer)
    dat<-df%>%
        as.data.frame()%>%
        rownames_to_column(var = "cluster")%>%
        pivot_longer(!cluster,names_to="groupvar",values_to="value")
    # change order clusters to numeric order
    dat$cluster<-as.numeric(pull(dat,cluster))
    dat$cluster<-factor(pull(dat,cluster),levels=as.character(unique(pull(dat,cluster))))

    #define colours
    if(length(unique(dat$cluster))<9){
        pal<-c(brewer.pal(length(unique(dat$cluster)),"Pastel2"))
    }else if(length(unique(dat$cluster))<18){
        pal<-c(brewer.pal(8, "Pastel2"),
               brewer.pal((length(unique(dat$cluster))-8),"Pastel1"))
    }else{
        pal<-c(brewer.pal(8, "Pastel2"),
               brewer.pal(9,"Pastel1"),
               brewer.pal((length(unique(dat$cluster))-17),"Set3"))
    }
    # plot
    textsize=18
    pd<-ggplot(dat,aes(x=groupvar,y=value,fill=cluster))+
        geom_bar(stat="identity",position=position_stack(reverse=T))+
        theme_classic()+
        labs(y="cluster contribution [%] alpha div")+
        theme(text=element_text(size=textsize),
              axis.text.x=element_text(angle=90,vjust=0.5),
              axis.title.x=element_blank())+
        scale_fill_manual(values=pal)

     print(pd)

}



which_case = function(charvector){
  #determine the case of elements in a character vector (such as gene symbols)
  charvector = head(charvector, 21)
  cases = c('uppercase'=sum(grepl("^[A-Z0-9]+$", charvector)),
            'titlecase'=sum(grepl("^[A-Z0-9][a-z0-9]+$", charvector)),
            'lowercase'=sum(grepl("^[a-z0-9]+$", charvector)))
  names(cases)[cases==max(cases)]
}

archr_plot_markers = function(proj, marker_genes, impute_weights=T, reduction="UMAP"){
  #plot gene activities from ArchR project
  #archr_plot_markers(gmat, c('Cd8a','Cd8b'), impute_weights=getImputeWeights(proj), reduction=dr_df)
  #proj can be either archr project or gene expression matrix of n genes * m cells
  if(class(proj) == 'ArchRProject'){
    gene_mat = get_gene_mat(proj)
  }else if(is.matrix(as.matrix(proj[1:2,1:2]))){
    gene_mat = proj
  }else{
    stop('proj must be either archr project or matrix')
  }
  if(length(reduction) == 1){
    stopifnot(is.character(reduction))
    dr_df <- jj_get_reduction_coords(proj, reduction)
  }else{
    stopifnot(is.matrix(reduction) | is.data.frame(reduction))
    dr_df = reduction
  }
  
  marker_genes <- switch(which_case(rownames(gene_mat)),
                         "lowercase"=tolower(marker_genes),
                         "titlecase"=tools::toTitleCase(tolower(marker_genes)),
                         "uppercase"=toupper(marker_genes))
  marker_genes = marker_genes[ marker_genes %in% rownames(gene_mat)]
  
  if(length(impute_weights) == 1){
    if(impute_weights & class(proj) == 'ArchRProject'){
      #when magic is used:
      proj = addImputeWeights(proj)
      goi_mat = as.matrix(msGetGOIRNAMat(gene_mat, goi = marker_genes, normFold = F))
      goi_mat = imputeMatrix(goi_mat, imputeWeights = getImputeWeights(proj)) #use magic values instead of original gene activities
      dr_df = cbind(dr_df, t(goi_mat))
    }else{
      dr_df = jj_bind_features_with_dr_df(gene_mat, features = marker_genes, dr_df = dr_df)
    }
  }else{
    stopifnot(identical(names(impute_weights), c('Weights', 'Params'))) #should be output from getImputeWeights(proj)
    goi_mat = as.matrix(msGetGOIRNAMat(gene_mat, goi = marker_genes, normFold = F))
    goi_mat = imputeMatrix(goi_mat, imputeWeights = impute_weights) #use magic values instead of original gene activities
    dr_df = cbind(dr_df, t(goi_mat))
  }
  
  
  #marker_genes = gsub('-', '_', marker_genes)
  gg = jj_plot_features(seurat_obj = NULL, 
                        reduction = dr_df,
                        pt.size = 0.5,
                        meta_features = marker_genes, my_title = '', 
                        cont_or_disc = 'c', colorScale = 'viridis', cap_top = 'auto',
                        return_gg_object = T)
  
  #put the gene label to the title and remove gene name next to color scale
  gg = lapply(gg, function(x) x + labs(color='', title=x$labels$colour) + theme_minimal(base_size = 8))
  names(gg) = marker_genes
  return(gg)
}

quick_reduction_rna = function(seurat_obj, ndims = 10, nfeatures= 2000, umap_name="umap"){
  #perform dimensionality reduction with standard parameters in Seurat
  seurat_obj <- NormalizeData(seurat_obj,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = "vst",
                                     nfeatures = nfeatures)
  seurat_obj <- ScaleData(seurat_obj,
                          features = VariableFeatures(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj,
                       features = VariableFeatures(object = seurat_obj))
  seurat_obj <- RunUMAP(seurat_obj,
                        dims = 1:ndims,
                        reduction.name= umap_name)
  return(seurat_obj)
}

underscore = function(x, replace_pattern = '-| |\\.'){
  x = gsub(replace_pattern, '_', x)
  return(x)
} 

get_gene_mat = function(proj){
  gene_mat = ArchR::getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix')
  dr_df <- as.data.frame(proj@cellColData)
  rnames = rowData(gene_mat)
  gene_mat = assays(gene_mat)[[1]]
  rownames(gene_mat) = rnames$name
  gene_mat = gene_mat[, match(rownames(dr_df), colnames(gene_mat))]
  stopifnot(identical(rownames(dr_df), colnames(gene_mat)))
  #names_match(dr_df, gene_mat)
  return(gene_mat)
}

plot_dotplot = function(obj, features_use, group_vec, scale_data = FALSE, 
                        cluster_rows = TRUE,  cluster_columns = TRUE,
                        cap_top = NULL, cap_bottom = NULL, max_pt_size=5){
  
  get_fraction_expressing = function (rna_mat, features_use, group_vec){
    rna_mat = rna_mat[rownames(rna_mat) %in% features_use, , 
                      drop = FALSE] %>% as.matrix %>%  base::t()
    cell_exp_ct_mat <- rna_mat %>% as.data.frame %>% dplyr::mutate(cluster = group_vec) %>% 
      dplyr::group_by(cluster) %>% dplyr::summarise_all(list(function(x) sum(x >  0))) %>%
      as.data.frame %>% column_to_rownames("cluster") %>% 
      base::t(.) %>% as.data.frame %>% rownames_to_column("Gene") %>% 
      tidyr::pivot_longer(!Gene, names_to = "cluster", values_to = "cell_exp_ct")
    cell_ct_mat <- data.frame(cluster = group_vec) %>% dplyr::group_by(cluster) %>% 
      dplyr::summarise(cell_ct = n()) %>% as.data.frame
    cell_exp_ct_mat = cell_exp_ct_mat %>% dplyr::left_join(cell_ct_mat, 
                                                           by = "cluster")
    return(cell_exp_ct_mat)
  }
  
  get_mean_expressing = function (rna_mat, features_use, group_vec, scale_data = FALSE){
    rna_mat = rna_mat[rownames(rna_mat) %in% features_use, , 
                      drop = FALSE] %>% as.matrix %>%  base::t()
    
    expr_mat <- rna_mat %>% as.data.frame %>% dplyr::mutate(cluster = group_vec) %>% 
      dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean)) %>% 
      as.data.frame %>% column_to_rownames("cluster") 
    
    if (scale_data) {
      expr_mat = scale(expr_mat)
    }
    
    expr_mat = base::t(expr_mat) %>% 
      as.data.frame %>% 
      rownames_to_column("Gene") %>% 
      tidyr::pivot_longer(!Gene, names_to = "cluster", values_to = "count")
    return(expr_mat)
  }
  
  
  stopifnot((is(obj, "Matrix") | class(obj) == "matrix"))
  stopifnot(identical(ncol(obj), length(group_vec)))
  stopifnot(all(features_use %in% rownames(obj)))
  stopifnot(length(features_use) > 1)
  
  features_use = unique(features_use)
  cell_exp_ct_mat = get_fraction_expressing(obj, features_use, 
                                            group_vec)
  expr_mat = get_mean_expressing(obj, features_use, group_vec, 
                                 scale_data = scale_data)
  expr_mat = expr_mat %>% dplyr::left_join(cell_exp_ct_mat, 
                                           by = c("Gene", "cluster"))
  
  markers <- expr_mat$Gene %>% unique()
  mat <- expr_mat %>% dplyr::filter(Gene %in% markers) %>% 
    dplyr::select(Gene, cluster, count) %>% 
    pivot_wider(names_from = cluster, values_from = count) %>% as.data.frame()
  row.names(mat) <- mat$Gene
  mat <- mat[, -1]
  
  if(cluster_rows) {
    clust <- hclust(dist(mat %>% as.matrix()))
  }else {
    clust = data.frame(labels = features_use, order = 1:length(features_use))
  }
  
  if(cluster_columns){
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
  }else{
    v_clust = data.frame(labels = levels(as.factor(group_vec)), 
                         order = 1:length(unique(as.character(group_vec))))
  }
  
  expr_mat_use = expr_mat %>% dplyr::filter(Gene %in% markers) %>% 
    dplyr::mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 
                    100, Gene = factor(Gene, levels = clust$labels[clust$order]), 
                  cluster = factor(cluster, levels = v_clust$labels[v_clust$order]))
  
  if(!scale_data) {
    expr_mat_use = expr_mat_use %>% dplyr::filter(count > 0, `% Expressing` > 1)
  }else{
    expr_mat_use = expr_mat_use %>% dplyr::filter(`% Expressing` >  1)
  }
  expr_mat_use$count = jj::jj_cap_vals(expr_mat_use$count, cap_top = cap_top, cap_bottom = cap_bottom)
  
  dotplot <- ggplot(expr_mat_use, aes(x = cluster, y = Gene, 
                                      color = count, size = `% Expressing`)) + geom_point() + 
    theme_minimal() +
    theme(axis.line = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_blank()) + 
    viridis::scale_color_viridis() + scale_y_discrete(position = "right") + 
    scale_size_continuous(range = c(0,max_pt_size))
  if(scale_data){
    dotplot = dotplot + labs(colour='scaled\nmean expr', y='')
  }else{
    dotplot = dotplot + labs(colour='mean expr', y='')
  }
  dotplot
}


read_markers = function(filepath, skiprows = "#"){
  #read markers skipping rows starting with `skiprows`
  #marker file can be txt file with one gene per row
  #or csv file with column names 'group' and 'feature' specifying signature names and genes, respectively 
  #if the file is incorrect, the function returns a list of length 1 with zero entries
  make_empty_list = function(){
    l1_empty = list()
    l1_empty[[1]] = character()
    return(l1_empty)
  }
  
  if(!file.exists(filepath)){
    message('Marker file does not exist. Skipping plots.')
    return(make_empty_list())
  }
  
  filetype = rev(unlist(strsplit(basename(filepath),'\\.')))[1]
  if(filetype == 'txt'){
    signature_list = list()
    #get the name of the txt file
    filename = unlist(strsplit(basename(filepath),'\\.'))[-2]
    #add all marker genes from the file to a list using the filename as a signature name (skipping rows starting with #)
    skiprows = paste0("^", skiprows)
    signature_list[[filename]] = tryCatch({unique(grep(skiprows, readLines(filepath, skipNul = T), value=T, invert=T))}, error = function(x) {message('Marker txt file could not be read successfully. Skipping plots. Original error message:\n', x); return(character())})
  }else if(filetype == 'csv'){
    sigdf = tryCatch({read.csv(filepath, header = TRUE, stringsAsFactors = F, comment.char = '#', row.names = NULL)}, error = function(x) {message('Marker csv file is empty or another error occured while reading. Skipping plots. Original error message:\n', x); return(NULL)})
    if(is.null(sigdf)){
      return(make_empty_list())
    }    
    if(!all(c('group','feature') %in% colnames(sigdf))){
      message("Marker csv file should contain column names 'group' and 'feature'. Skipping plots.")
      return(make_empty_list())
    }
    sigdf = sigdf[complete.cases(sigdf), ]
    if(nrow(sigdf) == 0){
      message("Marker csv file should contains 0 entries after filtering for complete cases.")
      return(make_empty_list())
    }
    sigdf$group = as.character(sigdf$group)
    signature_list = split(sigdf, sigdf$group)
    signature_list = lapply(signature_list, function(x) pull(x, feature))
  }else{
    message('Marker file should be either a .txt or a .csv file. Skipping plots.')
    return(make_empty_list())
  }
  if(Hmisc::all.is.numeric(names(signature_list))){
    names(signature_list) = paste0('s', names(signature_list))
  }
  signature_list
}

get_point_size = function(n){
  ifelse(n < 5000, 1.25,
         ifelse(n < 10000, 1, 
                ifelse(n < 25000, 0.75,
                       ifelse(n < 50000, 0.5, 0.3))))
}

### Config file processing functions

# Most validation should be performed by direct schema validation, more complex
# validation will be handled in these functions (which should use the schema, e.g. to get
# default values)

get_clustering_dimred <- function(config) {
  cluster_on <- unlist(strsplit(config$scParam$cluster_on, '_'))
  stopifnot(!anyNA(as.integer(cluster_on[2:3])))

  return(cluster_on)
}

check_integration_config <- function(config, config_schema, analysis_root) {
  integrate <- config$integration$integrate

  if (integrate) {
    if (is.null(config$integration$methods)) {
      # Obtain possible values for integration methods from config schema

      config_schema_integration <- config_schema$properties$integration
      integration_methods_available <- unlist(config_schema_integration$properties$methods$items$enum)

      stop(
        "No integration methods supplied. Available methods are: ",
        paste(integration_methods_available, collapse = ", ")
      )
    }

    if (is.null(config$integration$batch_info_file)) {
      stop("No path to batch info file supplied")
    } else {
      batch_info_path <- file.path(analysis_root, config$integration$batch_info_file)
    }

    if (!file.exists(batch_info_path)) {
      stop(paste0("The batch info file ", batch_info_path, " does not exist"))
    }
  }
}

check_clustering_config <- function(config) {
  cluster_on <- get_clustering_dimred(config)

  integration_methods <- c()
  if (config$integration$integrate) {
    integration_methods <- config$integration$methods
  }
  stopifnot(cluster_on[1] %in% c('pca', integration_methods))
}

shorten_names_unique <- function(name_vec, nchar_start = 8) {
  nchar_shortened <- nchar_start
  name_vec_shortened <- substr(name_vec, 1, nchar_shortened)

  while(length(unique(name_vec_shortened)) < length(unique(name_vec))) {
    nchar_shortened <- nchar_shortened + 1
    name_vec_shortened <- substr(name_vec, 1, nchar_shortened)
  }

  return(name_vec_shortened)
}

firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

shorten_names <- function(input_names, max_length = 31, sep = NULL) {
    name_length <- sapply(input_names, nchar)
    names_too_long <- name_length > max_length

    if (any(names_too_long)) {
        if(!is.null(sep)) {
            split_names <- stringr::str_split(input_names, sep)
            max_length_half <- floor((max_length - nchar(sep))/2)
            shortened_names <- sapply(
                split_names,
                function (x) paste0(
                    substring(x[1], 1, max_length_half), "_vs_", substring(x[2], 1, max_length_half)
                )
            )
            names(shortened_names) <- input_names
        } else {
            shortened_names <- sapply(input_names, function(x) substring(x, 1, max_length))
        }
    }

    return(shortened_names)
}

convertLibNames <- function(libNames) {
    libNames <- str_replace_all(libNames, "-", "_")
    libNames <- str_remove_all(libNames, "^p[0-9]{2,4}_")
    libNames <- str_remove_all(libNames, "_run_[1-9][0-9]*$")

    return(libNames)
}

rearrangeFactorLevels <- function(factorLevels, ctrl) {
    if (!(ctrl %in% factorLevels)) {
        stop(paste0("Control level '", ctrl, "' not found in factor levels '", paste(factorLevels, collapse = "', '"), "'"))
    }

    levelsExCtrl <- factorLevels[factorLevels != ctrl]
    factorLevelsOrdered <- c(ctrl, levelsExCtrl)

    return(factorLevelsOrdered)
}

sample_metadata <- function(metadata, group_by_cols, id_col, size_col) {
    metadata_sampled <- metadata %>% dplyr::group_by(across(all_of(group_by_cols))) %>%
        dplyr::group_modify(
            function(df, keys) {
                cell_ids_sampled <- sample(df[[id_col]], size = df[[size_col]][1], replace = FALSE)
                df_sampled <- df[base::match(cell_ids_sampled, df[[id_col]]), ]
                
                return(df_sampled)
            }
        )

    return(metadata_sampled)
}

sample_replicates_metadata <- function(metadata, ident_col, replicate_col, cond_col = NULL, mode = "replicate") {
    if (is.null(mode)) {
        stop("The mode parameter must be supplied!")
    }

    metadata <- tibble::rownames_to_column(metadata, "cell_id")

    if (mode == "enforce-fractions") {
        if (is.null(cond_col)) {
            stop("In mode 'enforce-fractions' the parameter 'cond_col' must be supplied")
        }

        cell_counts_lib <- metadata %>% dplyr::group_by(across(all_of(c(ident_col, replicate_col)))) %>% dplyr::count()
        cell_counts_min <- cell_counts_lib %>% dplyr::group_by(across(all_of(ident_col))) %>%
            dplyr::slice(which.min(n)) %>% dplyr::rename(n_lib_min = n, lib_min = all_of(replicate_col))

        cell_counts <- metadata %>% dplyr::group_by(across(all_of(c(ident_col, replicate_col, cond_col)))) %>% dplyr::count()
        cell_counts <- cell_counts %>% dplyr::group_by(across(all_of(c(ident_col, replicate_col)))) %>%
            dplyr::mutate(prop_cond = prop.table(n))

        cell_counts <- dplyr::inner_join(
            cell_counts,
            cell_counts_min,
            by = ident_col
        )

        cell_counts <- cell_counts %>% dplyr::rename(n_lib_min_total = n_lib_min) %>%
            dplyr::mutate(n_lib_min = pmin(round(n_lib_min_total*prop_cond, 0), n))

        metadata <- dplyr::inner_join(
            metadata,
            cell_counts,
            by = c(ident_col, replicate_col, cond_col)
        )

        metadata_sampled <- sample_metadata(
            metadata, group_by_cols = c(ident_col, replicate_col, cond_col),
            id_col = "cell_id", size_col = "n_lib_min"
        )
    } else if (mode == "replicate-condition") {
        if (is.null(cond_col)) {
            stop("In mode 'replicate-condition' the parameter 'cond_col' must be supplied")
        }

        cell_counts <- metadata %>% dplyr::group_by(across(all_of(c(ident_col, replicate_col, cond_col)))) %>% dplyr::count()
        cell_counts_min <- cell_counts %>% dplyr::group_by(across(all_of(c(ident_col, cond_col)))) %>%
            dplyr::slice(which.min(n)) %>% dplyr::rename(n_lib_min = n, lib_min = all_of(replicate_col))

        metadata <- dplyr::inner_join(
            metadata,
            cell_counts_min,
            by = c(ident_col, cond_col)
        )

        metadata_sampled <- sample_metadata(
            metadata, group_by_cols = c(ident_col, replicate_col, cond_col),
            id_col = "cell_id", size_col = "n_lib_min"
        )
    } else if (mode == "replicate") {
        cell_counts <- metadata %>% dplyr::group_by(across(all_of(c(ident_col, replicate_col)))) %>% dplyr::count()
        cell_counts_min <- cell_counts %>% dplyr::group_by(across(all_of(ident_col))) %>%
            dplyr::slice(which.min(n)) %>% dplyr::rename(n_lib_min = n, lib_min = all_of(replicate_col))

        metadata <- dplyr::inner_join(
            metadata,
            cell_counts_min,
            by = ident_col
        )

        metadata_sampled <- sample_metadata(
            metadata, group_by_cols = c(ident_col, replicate_col),
            id_col = "cell_id", size_col = "n_lib_min"
        )
    } else {
        stop(paste0("Invalid mode supplied: ", mode))
    }

    metadata_sampled <- tibble::column_to_rownames(metadata_sampled, "cell_id")

    return(metadata_sampled)
}

plot_cell_counts <- function(
    cell_counts, ident_col, idents = NULL,
    cond_col, replicate_col, y_col, ncol = 3
) {
    if (is.null(idents)) {
        idents <- unique(as.character((sort(cell_counts[[ident_col]]))))
    }

    cell_count_plots <- list()

    for (ident in idents) {
        cell_count_plots[[ident]] <- ggplot2::ggplot(
            cell_counts %>% dplyr::filter(.data[[ident_col]] == ident),
            ggplot2::aes(x = .data[[cond_col]], y = .data[[y_col]], fill = .data[[replicate_col]])
        ) + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
            ggplot2::labs(title = ident)
    }

    cell_count_plots_combined <- patchwork::wrap_plots(cell_count_plots, ncol = ncol) + 
        patchwork::plot_layout(guides = "collect")

    return(list(list = cell_count_plots, patchwork = cell_count_plots_combined))
}

match_wildcard_string <- function(string, vec, value = FALSE) {
    string_regex <- stringr::str_replace_all(string, "\\.", "\\\\.")
    string_regex <- paste0("^", stringr::str_replace_all(string_regex, "\\*", ".*"), "$")
    wildcard_matches <- grep(string_regex, vec, value = value)

    return(wildcard_matches)
}

expand_ident_group <- function(group, identVec) {
    group_expanded <- c()

    for (ident in group) {
        if (grepl("\\*", ident)) {
            ident <- sort(unique(match_wildcard_string(
                ident, identVec, value = TRUE
            )))
        }
        group_expanded <- c(group_expanded, ident)
    }

    return(group_expanded)
}

split_number <- function(number, split_size, max_residual) {
    n_sets <- floor(number/split_size)
    set_counts <- rep(split_size, length = n_sets)
    residual_counts <- number %% split_size

    if (number < split_size) {
        set_counts <- residual_counts
    } else if (residual_counts > max_residual) {
        n_sets <- n_sets + 1
        set_counts[n_sets] <- residual_counts
    } else {
        set_counts[n_sets] <- set_counts[n_sets] + residual_counts
    }

    return(set_counts)
}

split_list <- function(x, split_size, max_residual) {
    set_counts <- split_number(
        length(x), split_size = split_size, max_residual = max_residual
    )

    x_subsets <- list()
    idx_start <- 1

    for (set_no in seq_along(set_counts)) {
        x_subsets[[set_no]] <- x[idx_start:(idx_start + set_counts[set_no] - 1)]
        idx_start <- idx_start + set_counts[set_no]
    }

    return(x_subsets)
}
