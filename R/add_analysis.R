additional_analysis_list=c(
  'Add_Grad_Table',
  'Make_Grad_Plot',
  'Add_3_DMC_Promoter_Table',
  'Visualize_Enricher_5_dmc_Promoter',
  'Visualize_Enricher_3_dmc_Promoter',
  'Smart_Visualize_Enricher_3_dmc_Promoter',
  'Check_Gene_Intersection',
  'Make_Venn_Diagram',
  'Print_Count_of_DMCs_with_Colocal_Sig',
  '10ppm_Sites_Hitting_Wnt_Genes',
  'FA_Signature',
  'Mechanism_Signature_Table'
)

default_selection=F
selected_analyses=rep(default_selection, length(additional_analysis_list))
names(selected_analyses)=additional_analysis_list

#--------------------------------------
## SELECT ANALYSES HERE
#--------------------------------------

selected_analyses[ c('FA_Signature')
]=T

##Expect job.id, cur_multi, meth.thresh,
## and q.thresh to be set externally

outdir='p2_output'
out.string=paste0(outdir, '/', job.id, '_', meth.thresh)
out.string.nothresh=paste0(outdir, '/', job.id)

if (grepl('P2_Cd', job.id)){
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

cat(job.id, '\n')
# #----------------------------
# #Gradiated DMC Table
if (selected_analyses['Add_Grad_Table']==T){

  diff_sum=getGradiatedDiff(cur_multi, q=q.thresh)
  write.table(diff_sum, file=paste0(out.string.nothresh, '_Gradiated_DMC_Call_Table.txt'), quote=F, sep="\t", row.names = F)
}

# #----------------------------
# #Make Grad Plot
if (selected_analyses['Make_Grad_Plot']==T){
  
  cur_pdf=file=paste0(out.string.nothresh, '_Gradiated_DMC_Distr.pdf')
  pdf(cur_pdf)
    print(plotGradDiff(cur_multi, q.value = q.thresh, dist=T))
  dev.off()
}
 
# #----------------------------
# ##Visualizing Enricher results
if (selected_analyses['Visualize_Enricher_5_dmc_Promoter']==T){
  
 outdir='p2_output/selected_enricher/plots'
 job.string=paste0(job.id, '_', meth.thresh)
 out.string=paste0(outdir, '/', job.id, '_', meth.thresh)
 # out.string.nothresh=paste0(outdir, '/', job.id)

 selected_enricher_df=read.table('annotations/enricher/selected_analyses.txt', head=T, sep="\t")

 for (cur_cat in unique(selected_enricher_df$Category)){

   cur_dir=paste0('p2_output/selected_enricher/', cur_cat, '/')
   cur_ann_list=subset(selected_enricher_df, Category==cur_cat)$Annotation

   for (cur_ann in cur_ann_list){

     i=1
     blank_coefs=c()
     full_dat=data.frame()
     cur_dat=data.frame()
     full_terms=c()
     cur_terms=c()

     cur_coef=names(cur_multi[[3]])[[i]]
     cur_prefix=paste(job.id, cur_coef,meth.thresh,
                      'Genes_with_gte_5_DMCs_in_Promoters', cur_ann, sep='_')
     cur_file=paste0(cur_dir,cur_prefix, '.txt' )

     if (!file.exists(cur_file)){
       blank_coefs=c(blank_coefs, cur_coef)
     } else {
       cur_dat=read.table(cur_file, head=T, sep="\t")
       if (nrow(cur_dat)==0){
         blank_coefs=c(blank_coefs, cur_coef)

       } else {
         cur_dat=cur_dat[ , c('Term','Adjusted.P.value' )]
         cur_terms=head(as.character(cur_dat$Term[ order(cur_dat$Adjusted.P.value)]), n=20)

         colnames(cur_dat)[[2]]=cur_coef

         full_dat=cur_dat
         full_terms=cur_terms
       }

     }


 for ( i in 2:length(names(cur_multi[[3]]))){
       cur_coef=names(cur_multi[[3]])[[i]]
       cur_prefix=paste(job.id, cur_coef,meth.thresh,
                        'Genes_with_gte_5_DMCs_in_Promoters', cur_ann, sep='_')
       cur_file=paste0(cur_dir,cur_prefix, '.txt' )

       if (!file.exists(cur_file) ){
         blank_coefs=c(blank_coefs, cur_coef)
         next
       } else {
         cur_dat=read.table(cur_file, head=T, sep="\t")

         if (nrow(cur_dat)==0){
           blank_coefs=c(blank_coefs, cur_coef)

         } else {

           cur_dat=cur_dat[ , c('Term','Adjusted.P.value' )]

           cur_terms=head(as.character(cur_dat$Term[ order(cur_dat$Adjusted.P.value)]), n=20)
           full_terms=c(full_terms , cur_terms)

           colnames(cur_dat)[[2]]=cur_coef

           if (nrow(full_dat)==0){
             full_dat=cur_dat
           }else {
             full_dat=merge(full_dat, cur_dat, by='Term', all=T)
           }
         }


       }

     }


     #Filter data to most enriched results
     if (nrow(full_dat)==0) {
       cat("SKIPPING", job.id, cur_cat, cur_ann, 'with', nrow(full_dat), 'rows\n')
       next
     }

     full_dat=subset(full_dat, Term %in% full_terms)
     if (length(blank_coefs)>0){
       for (b in blank_coefs){
         full_dat[, b]=1
       }
     }

     rownames(full_dat)=full_dat$Term
     full_dat[ is.na(full_dat)] =1

     out_pdf=paste0(out.string,'_', cur_cat, '_', cur_ann,   '.pdf')
     pdf(out_pdf)
       pheatmap(data.matrix(full_dat[ , 2:ncol(full_dat)]),
                fontsize_row = 5,
                fontsize_col = 7,
                treeheight_row=0,
                treeheight_col=0 ,
                color = heat.colors(100),
                breaks=seq(0,1, length=101),
                cluster_cols = F,
                border_color = NA)
     dev.off()

     #sys.call(paste('pdftk  cat 2  ',out_pdf, 'out.pdf ; mv out.pdf' ,out_pdf))
     cat("Finished", job.id, cur_cat, cur_ann, 'with', nrow(full_dat), 'rows\n')

     ##Logged Plots
     # outdir='p2_output/selected_enricher/plots/logged'
     # job.string=paste0(job.id, '_', meth.thresh)
     # out.string=paste0(outdir, '/', job.id, '_', meth.thresh)

     sub.out.string=paste0(outdir, '/logged', job.id, '_', meth.thresh)
     out_pdf=paste0(sub.out.string,'_', cur_cat, '_', cur_ann,   '_Log.pdf')
     pdf(paste0(out.string,'_', cur_cat, '_', cur_ann,'_Log.pdf'))
       pheatmap(-log10(data.matrix(full_dat[ , 2:ncol(full_dat)])),
          fontsize_row = 5,
          fontsize_col = 7,
          treeheight_row=0,
          treeheight_col=0 ,
          color = rev(heat.colors(100)),
          breaks=seq(0,8, length=101),
          cluster_cols = F,
          border_color = NA)
     dev.off()
   }
 }

}

##----------------------------
##Add tables for 3 DMC promoters

if (selected_analyses['Add_3_DMC_Promoter_Table']==T){
 #cat(job.id, '\n')
 for (cur_coef in names(cur_multi[[3]])){
 
   #methDiff25p=get.methylDiff(getDiffFromMulti(cur_multi, coef=cur_coef), q=q.thresh, diff=meth.diff)
   methDiff25p=getMethylDiff(getDiffFromMulti(cur_multi, coef=cur_coef), q=q.thresh, diff=meth.thresh)
 
   job.string=paste0(job.id, '_', cur_coef)
   out.string=paste0(outdir, '/', job.string,'_', meth.thresh)
 
   #genic_meth_annotation=annotate.WithGenicParts( as(meth,"GRanges"),gene.obj)
   #genic_diff_annotation=annotate.WithGenicParts( methDiff25p,gene.obj)
   genic_diff_annotation=annotateWithGenicParts( methDiff25p,gene.obj)
 
   good_inds=which((1:nrow(methDiff25p)) %in% (genic_diff_annotation@dist.to.TSS)$target.row  )
 
   gene_table=genic_diff_annotation@members[ good_inds, ]
   gene_table=cbind(genic_diff_annotation@dist.to.TSS, gene_table)
   gene_table=merge(gene_table, mm10_refseq, by.x='feature.name', by.y='rsid')
   gene_table=gene_table[ , c('name', 'feature.name', 'dist.to.feature', 'prom', 'exon', 'intron') ]
 
 
 
 
   sub_table=subset(count(subset(gene_table, prom==1 )$name), freq>=3)
   write.table(sub_table$x, file=
                 paste0(out.string, '_Genes_with_gte_3_DMCs_in_Promoters.txt' ),
               quote=F, sep="\t", col=F, row=F)
 
 }
}
 

# #----------------------------
# # ## SMART Visualizing Enricher results for 3 DMC in Promoters
if (selected_analyses['Smart_Visualize_Enricher_3_dmc_Promoter']==T){
  

  job.string=paste0(job.id, '_', meth.thresh)
  
  
  selected_enricher_df=read.table('annotations/enricher/selected_analyses.txt', head=T, sep="\t")
  
  for (i in nrow(selected_enricher_df)){
    cur_cat=selected_enricher_df[ i , 'Category']
    
    cur_dir=paste0('p2_output/selected_enricher_3_dmc_prom/', cur_cat, '/')
    
  
    cur_ann=selected_enricher_df[ i , 'Annotation']
    sub.outdir=paste0('p2_output/selected_enricher_3_dmc_prom/plots/',cur_ann,'/')
    out.string=paste0(sub.outdir, '/', job.id, '_', meth.thresh)
    
    cur_coefs=names(cur_multi[[3]])
    cur_prefix=paste(job.id, cur_coefs,meth.thresh,
                     'Genes_with_gte_3_DMCs_in_Promoters', cur_ann, sep='_')
    cur_files=paste0(cur_dir,cur_prefix, '.txt' )
 
    out_pdf=paste0(out.string,'_', cur_cat, '_', cur_ann,   '.pdf')
    pdf(out_pdf)
      makeHeatmapsFromEnricherFiles(data.frame(Files=cur_files, Labels=names(cur_multi[[3]])))
    dev.off()
    
    ##Logged Plots
    
    sub.outdir=paste0('p2_output/selected_enricher_3_dmc_prom/plots/logged/',cur_ann,'/')
    sub.out.string=paste0(sub.outdir, job.id, '_', meth.thresh)
    out_pdf=paste0(sub.out.string,'_', cur_cat, '_', cur_ann,   '_Log.pdf')
    
    pdf(paste0(out.string,'_', cur_cat, '_', cur_ann,'_Log.pdf'))
    makeHeatmapsFromEnricherFiles(data.frame(Files=cur_files, Labels=names(cur_multi[[3]])),
                                  apply_log10 = T)
    dev.off()

    #cat("Finished", job.id, cur_cat, cur_ann, 'with', nrow(full_dat), 'rows\n')
  }
}


# #----------------------------
# # ##Visualizing Enricher results for 3 DMC in Promoters
if (selected_analyses['Visualize_Enricher_3_dmc_Promoter']==T){
 sub.outdir='p2_output/selected_enricher_3_dmc_prom/plots'
 job.string=paste0(job.id, '_', meth.thresh)
 out.string=paste0(sub.outdir, '/', job.id, '_', meth.thresh)
 
 
 selected_enricher_df=read.table('annotations/enricher/selected_analyses.txt', head=T, sep="\t")
 
 for (cur_cat in unique(selected_enricher_df$Category)){
 
   cur_dir=paste0('p2_output/selected_enricher_3_dmc_prom/', cur_cat, '/')
   cur_ann_list=subset(selected_enricher_df, Category==cur_cat)$Annotation
 
   for (cur_ann in cur_ann_list){
 
     i=1
     blank_coefs=c()
     full_dat=data.frame()
     cur_dat=data.frame()
     full_terms=c()
     cur_terms=c()
 
     cur_coef=names(cur_multi[[3]])[[i]]
     cur_prefix=paste(job.id, cur_coef,meth.thresh,
                      'Genes_with_gte_3_DMCs_in_Promoters', cur_ann, sep='_')
     cur_file=paste0(cur_dir,cur_prefix, '.txt' )
 
     if (!file.exists(cur_file)){
       blank_coefs=c(blank_coefs, cur_coef)
     } else {
       cur_dat=read.table(cur_file, head=T, sep="\t")
       if (nrow(cur_dat)==0){
         blank_coefs=c(blank_coefs, cur_coef)
 
       } else {
         cur_dat=cur_dat[ , c('Term','Adjusted.P.value' )]
         cur_terms=head(as.character(cur_dat$Term[ order(cur_dat$Adjusted.P.value)]), n=20)
 
         colnames(cur_dat)[[2]]=cur_coef
 
         full_dat=cur_dat
         full_terms=cur_terms
       }
 
     }
 
 
 for ( i in 2:length(names(cur_multi[[3]]))){
       cur_coef=names(cur_multi[[3]])[[i]]
       cur_prefix=paste(job.id, cur_coef,meth.thresh,
                        'Genes_with_gte_3_DMCs_in_Promoters', cur_ann, sep='_')
       cur_file=paste0(cur_dir,cur_prefix, '.txt' )
 
       if (!file.exists(cur_file) ){
         blank_coefs=c(blank_coefs, cur_coef)
         next
       } else {
         cur_dat=read.table(cur_file, head=T, sep="\t")
 
         if (nrow(cur_dat)==0){
           blank_coefs=c(blank_coefs, cur_coef)
 
         } else {
 
           cur_dat=cur_dat[ , c('Term','Adjusted.P.value' )]
 
           cur_terms=head(as.character(cur_dat$Term[ order(cur_dat$Adjusted.P.value)]), n=20)
           full_terms=c(full_terms , cur_terms)
 
           colnames(cur_dat)[[2]]=cur_coef
 
           if (nrow(full_dat)==0){
             full_dat=cur_dat
           }else {
             full_dat=merge(full_dat, cur_dat, by='Term', all=T)
           }
         }
 
 
       }
 
     }
 
}
 #     #Filter data to most enriched results
     if (nrow(full_dat)==0) {
       cat("SKIPPING", job.id, cur_cat, cur_ann, 'with', nrow(full_dat), 'rows\n')
       next
     }
 
     full_dat=subset(full_dat, Term %in% full_terms)
     if (length(blank_coefs)>0){
       for (b in blank_coefs){
         full_dat[, b]=1
       }
     }
 
     rownames(full_dat)=full_dat$Term
     full_dat[ is.na(full_dat)] =1
 
     out_pdf=paste0(out.string,'_', cur_cat, '_', cur_ann,   '.pdf')
     pdf(out_pdf)
       pheatmap(data.matrix(full_dat[ , 2:ncol(full_dat)]),
                fontsize_row = 5,
                fontsize_col = 7,
                treeheight_row=0,
                treeheight_col=0 ,
                color = heat.colors(100),
                breaks=seq(0,1, length=101),
                cluster_cols = F,
                border_color = NA)
     dev.off()
 
     #sys.call(paste('pdftk  cat 2  ',out_pdf, 'out.pdf ; mv out.pdf' ,out_pdf))
     cat("Finished", job.id, cur_cat, cur_ann, 'with', nrow(full_dat), 'rows\n')
 
      ##Logged Plots
     outdir='p2_output/selected_enricher_3_dmc_prom/plots'
     job.string=paste0(job.id, '_', meth.thresh)
     out.string=paste0(outdir, '/', job.id, '_', meth.thresh)
 
     sub.out.string=paste0(outdir, '/logged/', job.id, '_', meth.thresh)
     out_pdf=paste0(sub.out.string,'_', cur_cat, '_', cur_ann,   '_Log.pdf')
     pdf(paste0(out.string,'_', cur_cat, '_', cur_ann,'_Log.pdf'))
       pheatmap(-log10(data.matrix(full_dat[ , 2:ncol(full_dat)])),
          fontsize_row = 5,
          fontsize_col = 7,
          treeheight_row=0,
          treeheight_col=0 ,
          color = rev(heat.colors(100)),
          breaks=seq(0,8, length=101),
          cluster_cols = F,
          border_color = NA)
     dev.off()
   }
 }
 
# #---------------------------
# ## Checking intersection with gene lists

if (selected_analyses['Check_Gene_Intersection']==T){
 cur_dir='p2_output/'
 cur_ntd=data.frame( matrix( nrow=length(ntd_genes),
                             ncol=length(cur_multi[[3]])))
 rownames(cur_ntd)=ntd_genes
 colnames(cur_ntd)=names(cur_multi[[3]])

 cur_lrp6=data.frame( matrix( nrow=length(lrp6_genes),
                             ncol=length(cur_multi[[3]])))
 rownames(cur_lrp6)=lrp6_genes
 colnames(cur_lrp6)=names(cur_multi[[3]])

 cat(job.id, '\n')
 for (i in 1:length(cur_multi[[3]])){
   cur_coef=names(cur_multi[[3]])[[i]]
   cur_prefix=paste(job.id, cur_coef,meth.thresh,
                    'Genes_with_gte_3_DMCs_in_Promoters', sep='_')
   cur_file=paste0(cur_dir,cur_prefix, '.txt' )
   if (file.exists(cur_file) ){

     cur_dat=read.table(cur_file, head=F, sep="\t")$V1
     cur_ntd[ , i]=ntd_genes %in% cur_dat
     cur_lrp6[ , i]=lrp6_genes %in% cur_dat

   } else {
     cur_ntd[ , i]=F
     cur_lrp6[ , i]=F
   }

 }

 if (sum(data.matrix(cur_ntd))>1){
   cat('...has NTD Genes\n')
 }

 if (sum(data.matrix(cur_lrp6))>1){
   cat('...has Lrp6 Genes\n')
 }
}
 
# #----Make Venn Diagram

if (selected_analyses['Make_Venn_Diagram']==T){
  
  out.string=paste0(outdir, job.id, '_', meth.thresh)
  out_pdf=paste0(out.string, '_Venn_Diagram.pdf')
  pdf(out_pdf, width = 5, height =5)
    plot(makeVennDiagram(cur_multi, meth.thresh = meth.thresh,
                         q.thresh = q.thresh))
  dev.off()
}

if (selected_analyses['Print_Count_of_DMCs_with_Colocal_Sig']==T){

  if (grepl('P2_Cd', job.id)){
    cur_group="Cd"
    cur_het='IsCdHet'
    cur_sig=-1
  }else{
      cur_group='L'
      cur_het='IsLHet'
      cur_sig=1
  }
  cur_diff_call=getDiffCallMatrix(cur_multi, meth.thresh = meth.thresh,
                                  q.thresh = q.thresh)
  cur_diff=getDiffMatrix(cur_multi)
  cur_sig_diff=(sign(cur_diff[, 'Is10ppm'])*sign(cur_diff[, cur_het]))==cur_sig
  
  #cat(summary(cur_diff_call[  , 'Is10ppm'] &  cur_diff_call[ , cur_het])[ 'TRUE'], '\n')
  cat(summary(cur_sig_diff & cur_diff_call[  , 'Is10ppm'] &  cur_diff_call[ , cur_het])[ 'TRUE'], '\n')
}

#--------------------------
## Checking for Is10ppm genes hitting Wnt pathway genes
if (selected_analyses['10ppm_Sites_Hitting_Wnt_Genes']==T){
  
  cur_coef='Is10ppm'
  
  
  methDiff25p=get.methylDiff(getDiffFromMulti(cur_multi,
                               coef=cur_coef),
                            q=q.thresh, diff=meth.thresh)
  genic_diff_annotation=annotateWithGenicParts( methDiff25p,lrp6.obj)
  
  good_inds=which((1:nrow(methDiff25p)) %in% (genic_diff_annotation@dist.to.TSS)$target.row  )
  gene_table=genic_diff_annotation@members[ good_inds, ]
  gene_table=cbind(genic_diff_annotation@dist.to.TSS, gene_table)
  #gene_table=count(gene_table$feature.name)
  gene_table=count(gene_table[ , c('feature.name',
                                   'prom' ,
                                   'exon',
                                   'intron')])
  gene_table$Annotation='None'
  gene_table$Annotation[ gene_table$intron==1]='Intron'
  gene_table$Annotation[ gene_table$exon==1]='Exon'
  gene_table$Annotation[ gene_table$prom==1]='Promoter'
  
  gene_table=subset(gene_table, Annotation!='None')
  
  gene_table$BkGd=cur_group
  gene_table$Analysis=paste0(cur_sex, cur_filt)
  gene_table$meth.diff=meth.thresh
  
  colnames(gene_table)[[1]]='id'
  cur_cols=colnames(gene_table)
  write.table(gene_table,
              file=cur_table,
              append=T , quote=F, row=F, col.names = F, sep="\t")
  
  
  #gene_table=subset(gene_table, prom==1)
}

#----Looking for 10ppm signature sites
# Opposing interaction in Het10ppm condition should not be present

if (selected_analyses['FA_Signature']==T){
  
  sub.out.string=paste(job.string, 'FA_Signature', sep='_')

  
  cur_diff_call=getDiffCallMatrix(cur_multi, meth.thresh = meth.thresh, q.thresh = q.thresh)
  cur_diff=sign(getDiffMatrix(cur_multi))
  
  cur_10ppm_sig=cur_diff_call[ , 'Is10ppm']
  
  cur_strict_10ppm_sig=!( ( (cur_diff[ , 'Is10ppm'] +
                               cur_diff[ , paste0(cur_het, '10ppm')] )==0) &
                            (cur_diff_call[ , paste0(cur_het, '10ppm')]) )
  cur_strict_10ppm_sig=cur_10ppm_sig & cur_strict_10ppm_sig 
  
  cur_10ppm_sig=cur_strict_10ppm_sig
  
  cur_filt=selectMultiDiffPos(cur_multi, cur_10ppm_sig)
  genic_diff_annotation=annotateWithGenicParts( 
    getMethylDiff(getDiffFromMulti(cur_multi, 'Is10ppm'),
                   diff=meth.thresh,
                   q=q.thresh)
      ,gene.obj)
  
  good_inds=which((1:nrow(cur_multi[[1]])) %in% (genic_diff_annotation@dist.to.TSS)$target.row  )
  
  gene_table=genic_diff_annotation@members[ good_inds, ]
  gene_table=cbind(genic_diff_annotation@dist.to.TSS, gene_table)
  gene_table=merge(gene_table, mm10_refseq, by.x='feature.name', by.y='rsid')
  gene_table=gene_table[ , c('name', 'feature.name', 'dist.to.feature', 'prom', 'exon', 'intron') ]
  
  write.table(gene_table, file=paste0(sub.out.string, '_FA_Signature_Gene_Table.txt' ), quote=F, sep="\t", col=NA)
  sub_table=subset(count(subset(gene_table, prom==1 )$name), freq>=5)
  write.table(sub_table$x, file=
                paste0(sub.out.string, '_FA_Signature_Genes_with_gte_5_DMCs_in_Promoters.txt' ),
              quote=F, sep="\t", col=F, row=F)
  
  sub_table=subset(count(subset(gene_table, prom==1 )$name), freq>=3)
  write.table(sub_table$x, file=
                paste0(out.string, '_FA_Signature_Genes_with_gte_3_DMCs_in_Promoters.txt' ),
              quote=F, sep="\t", col=F, row=F)
  

}

#----------------
# Output Signature table
if (selected_analyses['Mechanism_Signature_Table']==T){
  
  # cur_table=paste0('p2_output/processed_tables/P2_Mechanism_Signature_DMC_Count_',meth.thresh,'.txt')
  # if (file.exists(cur_table)){
  #   file.remove(cur_table)
  # }
  
  cur_diff_call=getDiffCallMatrix(cur_multi, meth.thresh = meth.thresh, q.thresh = q.thresh)
  cur_diff=sign(getDiffMatrix(cur_multi))
  
  #Colocal sig
  ## Shoul occur in both FA and het, with background depdendent relative sign
  cur_colocal_sig=cur_diff_call[ , 'Is10ppm'] & cur_diff_call[ , cur_het]
    
  cur_strict_colocal_sig=cur_colocal_sig &
    ((cur_diff[, 'Is10ppm']*cur_diff[, cur_het])==cur_sig)
  
  cur_colocal_sig=sum(cur_colocal_sig)
  cur_strict_colocal_sig=sum(cur_strict_colocal_sig)
  
  #FA Signal
  ## Shoul occur in FA , and no contravening interaction
  cur_10ppm_sig=cur_diff_call[ , 'Is10ppm']
  
  cur_strict_10ppm_sig=!( ( (cur_diff[ , 'Is10ppm'] +
                          cur_diff[ , paste0(cur_het, '10ppm')] )==0) &
    (cur_diff_call[ , paste0(cur_het, '10ppm')]) )
  cur_strict_10ppm_sig=cur_10ppm_sig & cur_strict_10ppm_sig 
  
  cur_10ppm_sig=sum(cur_10ppm_sig)
  cur_strict_10ppm_sig=sum(cur_strict_10ppm_sig)
  
  
  
  #Interaction Signal
  cur_inter_sig=sum(cur_diff_call[ , paste0(cur_het, '10ppm')])
  
  #Strict Interaction Signal
  cur_strict_inter_sig=!(cur_diff_call[, 'Is10ppm']) & !(cur_diff_call[, cur_het])
  cur_strict_inter_sig=cur_diff_call[ , paste0(cur_het, '10ppm')] &
                    cur_strict_inter_sig
  
  cur_strict_inter_sig=sum(cur_strict_inter_sig)

  
  cur_df=data.frame(Colocal=cur_colocal_sig,
                    Strict_Colocal=cur_strict_colocal_sig,
                    FA=cur_10ppm_sig,
                    Strict_FA=cur_strict_10ppm_sig,
                    Interaction=cur_inter_sig,
                    Strict_Interaction=cur_strict_inter_sig,
                    Analysis=paste0(cur_sex, cur_filt),
                    Group=cur_group, 
                    meth.thresh=meth.thresh)
  cur_df=melt(cur_df, id=c('Analysis', 'Group', 'meth.thresh'))
  cur_df$Strict=grepl('Strict', cur_df$variable)
  write.table(cur_df, file=cur_table, append=T, quote=F, col=F, row=F)
}


