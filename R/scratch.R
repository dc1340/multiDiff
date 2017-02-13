p2_sexcomb_cd_diff=calculateMultiDiffMeth( p2_cd_meth, 
cd_design[ , c(2,3,4)], 'Is10ppm+IsCdHet+IsCdHet10ppm', num.cores=8)

p2_sexcomb_l_diff=calculateMultiDiffMeth( p2_l_meth, 
                                           l_design[ , c(2,3,4)], 'Is10ppm+IsLHet+IsLHet10ppm', num.cores=8)




cd_25_dmc=read.table('p2_output/processed_tables/Cd_25Diff_DMCs_Full_Summary.txt', head=T)

for (cur_file in list.files('p2_output/processed_tables',pattern='*.txt', full.names = T)){
  cur_lab=gsub('_DMCs_Full_Summary.txt', '', basename(cur_file))
  cur_dat=read.table(cur_file, head=T)
  cur_add_dat=cur_add_dat=data.frame(matrix(ncol=ncol(cur_dat), nrow=length(unique(cur_dat$Group))))
  
  colnames(cur_add_dat)=colnames(cur_dat)
  cur_add_dat$DMC_Type='All'
  i=1
  for (g in unique(cur_dat$Group)){
    cur_add_dat[ i, ncol(cur_add_dat)]=g
    cur_add_dat[ i, 2:(ncol(cur_add_dat)-1)]= colSums(data.matrix(subset(cur_dat, Group==g)))[ 2:(ncol(cur_add_dat)-1)]
    i=i+1
  }
  cur_dat=rbind(cur_dat, cur_add_dat)
  
  cur_dat=melt( cur_dat, id.vars =c( 'DMC_Type', 'Group'))

  cur_plot=ggplot(cur_dat, aes(fill=Group, y=value, x =variable ))+
    geom_bar(stat='identity', position='dodge')+
    theme(axis.text.x =
            element_text(angle = 90, hjust = 1, size = 10),
          legend.text = element_text(size = 10))+
    ggtitle(cur_lab)+
    xlab('Effect')+
    ylab('# DMCs')+ylim(0,210000)+facet_wrap(~DMC_Type)+
    scale_fill_brewer("Analysis", palette = 'RdYlGn')
    #scale_fill_discrete(guide = guide_legend(title = "Analysis"))
  
  pdf(paste0(dirname(cur_file), '/', cur_lab, '_Summary_DMC_Plot.pdf'), width=10, height=4)
    print(cur_plot)
  dev.off()
}



cur_plot=ggplot(cur_dat, aes(fill=Group, y=value, x =variable ))+
  geom_bar(stat='identity', position='dodge')+
  theme(axis.text.x =
          element_text(angle = 90, hjust = 1, size = 10),
        legend.text = element_text(size = 10))+
  ylab('# DMCs')+facet_wrap(~DMC_Type)+
  
  
#----------------------------------------------
##Colocation Analysis
  #Cd
  
  meth.thresh=25
  cur_dat=data.frame(DMCs=c(18015, 17102 , 3286,3188 ),
                     Analysis=final.labels)
  cur_het='IsCdHet'  
  write.table(cur_dat, file=paste0('p2_output/processed_tables/fa_colocation/P2_Cd_Is10ppm_',cur_het,'_',meth.thresh,'_DMC_Counts.txt')
        , sep="\t", quote=F , row=F)
  
  cur_plot=ggplot(cur_dat, aes(Analysis, DMCs, fill=Analysis))+
    geom_bar(stat='identity')+
    theme(axis.text.x =
            element_blank(),
          legend.text = element_text(size = 10))+
    ggtitle(paste('P2 Cd, FA and', cur_het,'Common DMCs at', meth.thresh, 'diff' ))+
    ylab('# DMCs')+ylim(0,65000)+
    scale_fill_brewer("Analysis", palette = 'RdYlGn')
  
  pdf(paste0('p2_output/processed_tables/fa_colocation/P2_Cd_Is10ppm_',cur_het,'_', meth.thresh, '_DMC_Counts.pdf'),
      width=7, height=4)
  print(cur_plot)
  dev.off()

  #L Analysis
  
  cur_dat=data.frame(DMCs=c(  14026 , 356, 5993, 134),
                     Analysis=final.labels)
  

  cur_het='IsLHet'  
  write.table(cur_dat, file=paste0('p2_output/processed_tables/fa_colocation/P2_L_Is10ppm_',cur_het,'_', meth.thresh, '_DMC_Counts.txt')
  , sep="\t", quote=F , row=F)
  
  
  cur_plot=ggplot(cur_dat, aes(Analysis, DMCs, fill=Analysis))+
    geom_bar(stat='identity')+
    theme(axis.text.x =element_blank(),
          legend.text = element_text(size = 10))+
    ggtitle(paste('P2 L, FA and', cur_het,'Common DMCs at ',meth.thresh,'diff' ))+
    ylab('# DMCs')+ylim(0,65000)+
    scale_fill_brewer("Analysis", palette = 'RdYlGn')
  
  pdf(paste0('p2_output/processed_tables/fa_colocation/P2_L_Is10ppm_',cur_het,'_', meth.thresh, '_DMC_Counts.pdf'),
      width=7, height=4)
    print(cur_plot)
  dev.off()
  
  #Cd -10meth .diff
  
  meth.thresh=10
  
  cur_dat=data.frame(DMCs=c(61302, 
                            59487 ,
                            45120 ,
                            44007 ),
                     Analysis=final.labels)
  cur_het='IsCdHet'  
  write.table(cur_dat, file=paste0('p2_output/processed_tables/fa_colocation/P2_Cd_Is10ppm_',cur_het,'_', meth.thresh,'_DMC_Counts.txt'),
              sep="\t", quote=F , row=F)
  
  cur_plot=ggplot(cur_dat, aes(Analysis, DMCs, fill=Analysis))+
    geom_bar(stat='identity')+
    theme(axis.text.x =
            element_blank(),
          legend.text = element_text(size = 10))+
    ggtitle(paste('P2 Cd, FA and', cur_het,'Common DMCs at', meth.thresh, 'diff' ))+
    ylab('# DMCs')+ylim(0,65000)+
    scale_fill_brewer("Analysis", palette = 'RdYlGn')
  
  pdf(paste0('p2_output/processed_tables/fa_colocation/P2_Cd_Is10ppm_',cur_het,'_', meth.thresh,'_DMC_Counts.pdf'),
      width=7, height=4)
  print(cur_plot)
  dev.off()
  
  #L Analysis
  
  cur_dat=data.frame(DMCs=c( 40482, 
                             742, 
                             43273, 
                             653) ,
                     Analysis=final.labels)
  
  
  cur_het='IsLHet'  
  write.table(cur_dat, file=paste0('p2_output/processed_tables/fa_colocation/P2_L_Is10ppm_',cur_het,'_', meth.thresh,'_DMC_Counts.txt')
              , sep="\t", quote=F , row=F)
  
  
  cur_plot=ggplot(cur_dat, aes(Analysis, DMCs, fill=Analysis))+
    geom_bar(stat='identity')+
    theme(axis.text.x =element_blank(),
          legend.text = element_text(size = 10))+
    ggtitle(paste('P2 L, FA and', cur_het,'Common DMCs at', meth.thresh,'diff' ))+
    ylab('# DMCs')+ylim(0,65000)+
    scale_fill_brewer("Analysis", palette = 'RdYlGn')
  
  pdf(paste0('p2_output/processed_tables/fa_colocation/P2_L_Is10ppm_',cur_het,'_', meth.thresh,'_DMC_Counts.pdf'),
      width=7, height=4)
  print(cur_plot)
  dev.off()
  
  

##Combine data
  cur_dat=data.frame(DMCs=c(61302, 
                            59487 ,
                            45120 ,
                            44007 ),
                     Analysis=final.labels, meth.thresh=10, BkGd='Cd')
  
  full_dat=cur_dat
  
  cur_dat=data.frame(DMCs=c( 40482, 
                             742, 
                             43273, 
                             653) ,
                     Analysis=final.labels, meth.thresh=10, BkGd='L')
  full_dat=rbind(full_dat, cur_dat)
  
  cur_dat=data.frame(DMCs=c(18015, 17102 , 3286,3188 ),
                     Analysis=final.labels, meth.thresh=25, BkGd='Cd')
  full_dat=rbind(full_dat, cur_dat)
  
  cur_dat=data.frame(DMCs=c(  14026 , 356, 5993, 134),
                     Analysis=final.labels, meth.thresh=25, BkGd='L')
  full_dat=rbind(full_dat, cur_dat)
  
  full_dat$meth.thresh=as.factor(full_dat$meth.thresh)
  full_dat$meth.thresh=relevel(full_dat$meth.thresh, c('25'))
  
  pdf('p2_output/processed_tables/fa_colocation/P2_Summary_FA_Lrp6_Common_DMCs.pdf',
      height=5, width = 7)
    ggplot(full_dat, aes(Analysis, DMCs, fill=Analysis))+
        geom_bar(stat='identity')+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())+
      facet_grid(meth.thresh~BkGd)+
      scale_fill_brewer("Analysis", palette = 'RdYlGn')+
      ylab('# DMCs')+ylim(0,65000)+
      ggtitle('P2 DMCs Effected by FA and Lrp6 Mutations')
  dev.off()
  
  #----Do requiring signature to be correct
  
  ##Combine data
  cur_dat=data.frame(DMCs=c(3581, 
                            3551, 
                            721, 
                            700 ),
                     Analysis=final.labels, meth.thresh=10, BkGd='Cd')
  
  full_dat=cur_dat
  
  cur_dat=data.frame(DMCs=c( 32371, 
                             645, 
                             40455, 
                             574 ) ,
                     Analysis=final.labels, meth.thresh=10, BkGd='L')
  full_dat=rbind(full_dat, cur_dat)
  
  cur_dat=data.frame(DMCs=c(608, 
                            595, 
                            5, 
                            5  ),
                     Analysis=final.labels, meth.thresh=25, BkGd='Cd')
  full_dat=rbind(full_dat, cur_dat)
  
  cur_dat=data.frame(DMCs=c( 12276, 
                             325, 
                             5941, 
                             129 ),
                     Analysis=final.labels, meth.thresh=25, BkGd='L')
  full_dat=rbind(full_dat, cur_dat)
  
  full_dat$meth.thresh=as.factor(full_dat$meth.thresh)
  full_dat$meth.thresh=relevel(full_dat$meth.thresh, c('25'))
  
  pdf('p2_output/processed_tables/fa_colocation/P2_Summary_FA_Lrp6_Common_DMCs_with_Expected_Signature.pdf',
      height=5, width = 7)
  ggplot(full_dat, aes(Analysis, DMCs, fill=Analysis))+
    geom_bar(stat='identity')+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    facet_grid(meth.thresh~BkGd)+
    scale_fill_brewer("Analysis", palette = 'RdYlGn')+
    ylab('# DMCs')+ylim(0,65000)+
    ggtitle('P2 DMCs With Expected Signature of FA and Lrp6 Mutations')
  dev.off()
  
  
  for (cur_coef in names(cur_multi[[3]])) {
    methDiff25p=getMethylDiff(getDiffFromMulti(
      cur_multi,coef=cur_coef),
      q=q.thresh, diff=meth.thresh)
    cat(cur_coef, min(abs(getData(methDiff25p)$meth.diff)),  nrow(methDiff25p), '\n')

  }
  
#----10ppm signaute analysis
  
  meth.thresh=25
  cur_coef='Is10ppm'
  cur_outdir='/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/plots/FA_Signature_Plots/'
  
  number_of_terms=10
  
  for (i in 1:nrow(selected_enricher_df)){
    cur_cat=selected_enricher_df[i , 'Category' ]
    cur_ann=selected_enricher_df[i , 'Annotation' ]
    
    cur_input_dir=paste0('p2_output/selected_enricher_3_dmc_prom/',cur_cat,'/')
    cur_files=paste0(cur_input_dir,
                     jobid_subset_df$job.id,
                     '_', meth.thresh,
                     '_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters_',
                     cur_ann,'.txt')
    
    
    
    cur_ann_df=data.frame(Files=cur_files, Labels=jobid_subset_df$Labels)
    
    cur_file=paste0(cur_dir,'P2_',cur_ann,'_',number_of_terms,'_Terms_',meth.thresh,'.pdf' )
    pdf(cur_file)
    makeHeatmapsFromEnricherFiles(cur_ann_df,
                                  num=number_of_terms)
    dev.off()
    
    cur_file=paste0(cur_dir,
                    'P2_',cur_ann,'_',number_of_terms,'_Terms_Logged_',meth.thresh,'.pdf' )
    pdf(cur_file)
    makeHeatmapsFromEnricherFiles(cur_ann_df,
                                  num=number_of_terms, app=T)
    dev.off()
  }
  
  #----10ppm NEGATIVE test signaute analysis
  meth.thresh=10
  cur_coef='Is10ppm'
  cur_outdir='/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/plots/FA_Signature_Plots/'
  
  number_of_terms=10
  
  for (i in 1:nrow(selected_enricher_df)){
    cur_cat=selected_enricher_df[i , 'Category' ]
    cur_ann=selected_enricher_df[i , 'Annotation' ]
    
    cur_input_dir=paste0('p2_output/selected_enricher_3_dmc_prom/',cur_cat,'/')
    cur_files=paste0(cur_input_dir,
                     jobid_subset_df$job.id,
                     '_Is10ppm_', meth.thresh,
                     '_Genes_with_gte_3_DMCs_in_Promoters_',
                     cur_ann,'.txt')

    
      
    cur_ann_df=data.frame(Files=cur_files, Labels=jobid_subset_df$Labels)
    
    cur_file=paste0(cur_dir,'P2_No_Sig_',cur_ann,'_',number_of_terms,'_Terms_',meth.thresh,'.pdf' )
    pdf(cur_file)
      makeHeatmapsFromEnricherFiles(cur_ann_df,
                                    num=number_of_terms)
    dev.off()
    
    cur_file=paste0(cur_dir,
                    'P2_No_Sig_',cur_ann,'_',number_of_terms,'_Terms_Logged_',meth.thresh,'.pdf' )
    pdf(cur_file)
      makeHeatmapsFromEnricherFiles(cur_ann_df,
                                    num=number_of_terms, app=T)
    dev.off()
  }
  

  P2_L_Sex_Combined_Sites_with_FA_IsLHet_Colocalization_and_Signature_Is10ppm_25_Genes_with_gte_3_DMCs_in_Promoters_Human_Phenotype_Ontology.txt
    
  tfiles=c('/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/Ontologies/P2_Cd_All_Interactions_10_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters_GO_Biological_Process_2015.txt',
  '/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/Ontologies/P2_Cd_All_Interactions_Filtered_Background_IsCd_10_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters_GO_Biological_Process_2015.txt',
  '/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/Ontologies/P2_Cd_Sex_Combined_10_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters_GO_Biological_Process_2015.txt',
  '/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/Ontologies/P2_Cd_Sex_Combined_Filtered_Wt_Combined_IsCd_10_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters_GO_Biological_Process_2015.txt')

  cur_dir='p2_output/selected_enricher_3_dmc_prom/'
  cur_prefix=P2_L_Sex_Combined_Sites_with_FA_IsLHet_Colocalization_and_Signature_',Is10ppm,'_25_Genes_with_gte_3_DMCs_in_Promoters_Human_Phenotype_Ontology.txt''
  P2_L_Sex_Combined_Sites_with_FA_IsLHet_Colocalization_and_Signature_Is10ppm_25_Genes_with_gte_3_DMCs_in_Promoters_Human_Phenotype_Ontology.txt
  P2_L_Sex_Combined_Sites_with_FA_IsLHet_Colocalization_and_Signature_Is10ppm_25_Genes_with_gte_3_DMCs_in_Promoters_Human_Phenotype_Ontology.txt

  meth.thresh=25
  number_of_terms=20
  for (i in 1:nrow(selected_enricher_df)){
    cur_cat=selected_enricher_df[i , 'Category' ]
    cur_ann=selected_enricher_df[i , 'Annotation' ]
    
  
    cur_files=paste0('/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/',
                     cur_cat,
                     '/P2_L_Sex_Combined_Sites_with_FA_IsLHet_Colocalization_and_Signature_',
                     cur_coefs,
                     '_25_Genes_with_gte_3_DMCs_in_Promoters_',
                     cur_ann,'.txt')
    
    cur_pdf=paste0('/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/plots/P2_L_Sex_Combined_Sites_with_FA_IsLHet_Colocalization_and_Signature_25_',cur_cat, '_',cur_ann,'.pdf')
    pdf(cur_pdf)
    makeHeatmapsFromEnricherFiles(data.frame(files=cur_files, label=names(cur_multi[[3]])),
                                  num=number_of_terms)
    dev.off()
    
    cur_pdf=paste0('/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/plots/P2_L_Sex_Combined_Sites_with_FA_IsLHet_Colocalization_and_Signature_25_',cur_cat, '_',cur_ann,'_Log.pdf')
    pdf(cur_pdf)
    makeHeatmapsFromEnricherFiles(data.frame(files=cur_files, label=names(cur_multi[[3]])),
                                  num=number_of_terms, apply_log10 = T)
    dev.off()
    
  }

#---FA Signature Gene Plots
cur_files=list.files('p2_output/', pattern='_25_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters.txt*', full.names = T)
cur_df=data.frame(Files=cur_files)
cur_df$Labels=gsub('_25_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters.txt', '', cur_files)
cur_df$Labels=gsub('P2_', '', cur_df$Labels)
cur_df$Labels=gsub('Filtered_Background_IsCd', 'Background_Filtered', cur_df$Labels)
cur_df$Labels=gsub('Filtered_Wt_Combined_IsCd', 'Background_Filtered', cur_df$Labels)
cur_df$Labels=gsub('p2_output//', '', cur_df$Labels)

i=1
cur_file=as.character(cur_df[ i, 'Files'])
job.id=as.character(cur_df[ i, 'Labels'])
if (grepl('^Cd', job.id)){
  cur_group="Cd"
  cur_het='IsCdHet'
  cur_sig=-1
}else{
  cur_group='L'
  cur_het='IsLHet'
  cur_sig=1
}

if (grepl('Filtered', job.id)){
  cur_filt="_Background_Filtered"
}else{
  cur_filt=''
}

if (grepl('Sex', job.id)){
  cur_sex="Sex_Combined"
}else{
  cur_sex='All_Interactions'
}

cur_dat=read.table(cur_file, head=F)
cur_dat=data.frame( Analysis=paste0(cur_sex, cur_filt),
                    Count=nrow(cur_dat),
                    BkGd=cur_group,
                    meth.diff=meth.thresh)
# colnames(cur_dat)='Gene'
# cur_dat$Analysis=cur_lab

full_fa_dat=cur_dat

for (i in 2:nrow(cur_df)){
  cur_file=as.character(cur_df[ i, 'Files'])
  job.id=as.character(cur_df[ i, 'Labels'])
  if (grepl('^Cd', job.id)){
    cur_group="Cd"
    cur_het='IsCdHet'
    cur_sig=-1
  }else{
    cur_group='L'
    cur_het='IsLHet'
    cur_sig=1
  }
  
  if (grepl('Filtered', job.id)){
    cur_filt="_Background_Filtered"
  }else{
    cur_filt=''
  }
  
  if (grepl('Sex', job.id)){
    cur_sex="Sex_Combined"
  }else{
    cur_sex='All_Interactions'
  }
  
  cur_dat=read.table(cur_file, head=F)
  cur_dat=data.frame( Analysis=paste0(cur_sex, cur_filt),
                      Count=nrow(cur_dat),
                      BkGd=cur_group,
                      meth.diff=meth.thresh)
  
  full_fa_dat=rbind(full_fa_dat, cur_dat)
}

#-----Combine Het10ppm Enricher

meth.thresh=25

cur_outdir='/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/plots/Het10ppm_Signature_Plots//'

number_of_terms=20

for (cur_group in c('L', 'Cd')){
  for (i in 1:nrow(selected_enricher_df)){
    cur_cat=selected_enricher_df[i , 'Category' ]
    cur_ann=selected_enricher_df[i , 'Annotation' ]
    
    cur_df=subset(jobid_subset_df, Group==cur_group)
    cur_input_dir=paste0('p2_output/selected_enricher_3_dmc_prom/',cur_cat,'/')
    cur_files=paste0(cur_input_dir,
                     cur_df$job.id,
                     '_Is', cur_group,'Het10ppm_', meth.thresh,
                     '_Genes_with_gte_3_DMCs_in_Promoters_',
                     cur_ann,'.txt')
    
    
    cur_ann_df=data.frame(Files=cur_files, Labels=cur_df$Labels)
    
    
    cur_file=paste0(cur_outdir,'P2_Is',cur_group, 'Het10ppm_',
                    cur_ann,'_',number_of_terms,'_Terms_',meth.thresh,'.pdf' )
    pdf(cur_file)
    makeHeatmapsFromEnricherFiles(cur_ann_df,
                                  num=number_of_terms)
    dev.off()
    
    cur_file=paste0(cur_outdir,'P2_Is',cur_group, 'Het10ppm_',
                    cur_ann,'_',number_of_terms,'_Terms_Logged_',meth.thresh,'.pdf' )
    pdf(cur_file)
    makeHeatmapsFromEnricherFiles(cur_ann_df,
                                  num=number_of_terms, app=T)
    dev.off()
  }
}

