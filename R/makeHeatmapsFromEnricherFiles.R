makeHeatmapsFromEnricherFiles<-function(files_and_label_df, number_terms_per_file=20,
                                        apply_log10=F){
  
  enricher_data=mergeDataFromEnricherFiles(files_and_label_df , number_terms_per_file)
  makeHeatmapsFromEnricherData(enricher_data, apply_log10)
}

makeHeatmapsFromEnricherData<-function(merged_enricher_data, apply_log10=F){

    
    
    if (apply_log10){
      merged_enricher_data[ , 2:ncol(merged_enricher_data)]=-log10(merged_enricher_data[ ,
                                                            2:ncol(merged_enricher_data)])
      
      pheatmap(data.matrix(merged_enricher_data[ ,
                                                 2:ncol(merged_enricher_data)]),
               fontsize_row = 5,
               fontsize_col = 7,
               treeheight_row=0,
               treeheight_col=0 ,
               color = rev(heat.colors(100)),
               breaks=seq(0,8, length=101),
               cluster_cols = F,
               border_color = NA)
      
    } else {
  
      pheatmap(data.matrix(merged_enricher_data[ ,
                            2:ncol(merged_enricher_data)]),
               fontsize_row = 5,
               fontsize_col = 7,
               treeheight_row=0,
               treeheight_col=0 ,
               color = heat.colors(100),
               breaks=seq(0,0.1, length=101),
               cluster_cols = F,
               border_color = NA)
      }
}

mergeDataFromEnricherFiles<-function(files_and_label_df, number_terms_per_file=20){
  
  i=1
  blank_labels=c()
  full_dat=data.frame()
  cur_dat=data.frame()
  full_terms=c()
  cur_terms=c()
  
  cur_file=as.character(files_and_label_df[ i, 1])
  cur_label=as.character(files_and_label_df[ i, 2])
  
  if (!file.exists(cur_file)){
    blank_labels=c(blank_labels, cur_label)
  } else {
    cur_dat=read.table(cur_file, head=T, sep="\t")
    if (nrow(cur_dat)==0){
      blank_label=c(blank_labels, cur_coef)
      
    } else {
      cur_dat=cur_dat[ , c('Term','Adjusted.P.value' )]
      cur_terms=head(as.character(cur_dat$Term[ order(cur_dat$Adjusted.P.value)]), n=number_terms_per_file)
      
      colnames(cur_dat)[[2]]=cur_label
      
      full_dat=cur_dat
      full_terms=cur_terms
    }
    
  }
  
  for ( i in 2:nrow(files_and_label_df)){
    
    cur_file=as.character(files_and_label_df[ i, 1])
    cur_label=as.character(files_and_label_df[ i, 2])
    
    
    blank_labels=c(blank_labels, cur_label)
    if (!file.exists(cur_file)){
      blank_labels=c(blank_labels, cur_label)
    } else{
      cur_dat=read.table(cur_file, head=T, sep="\t")
      
      if (nrow(cur_dat)==0){
        blank_label=c(blank_labels, cur_coef)
        
      }else {
        cur_dat=cur_dat[ , c('Term','Adjusted.P.value' )]
        cur_terms=head(as.character(cur_dat$Term[ order(cur_dat$Adjusted.P.value)]), n=number_terms_per_file)
        full_terms=c(full_terms , cur_terms)
        
        colnames(cur_dat)[[2]]=cur_label
        
        full_dat=merge(full_dat, cur_dat, by='Term', all=T)
      }
    }
  }
  
  
  full_dat=subset(full_dat, Term %in% full_terms)
  if (length(blank_coefs)>0){
    for (b in blank_coefs){
      full_dat[, b]=1
    }
  }
  
  rownames(full_dat)=full_dat$Term
  full_dat[ is.na(full_dat)] =1
  
  return(full_dat)
}

