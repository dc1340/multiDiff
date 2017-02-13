#cur_multi=p2_sexcomb_l_diff
#job.string='P2_Combined_sex_L_JM_25'

# q.thresh=0.01
# meth.diff=25

cat(job.id, '\n')
for (cur_coef in names(cur_multi[[3]])){

  #methDiff25p=get.methylDiff(getDiffFromMulti(cur_multi, coef=cur_coef), q=q.thresh, meth.diff=meth.diff)
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
  
  write.table(gene_table, file=paste0(out.string, '_Gene_Table.txt' ), quote=F, sep="\t", col=NA)
  write.table(unique(gene_table$name), file=
                paste0(out.string, '_Genes_with_Any_DMCs.txt' ),
              quote=F, sep="\t", col=F, row=F)
  
  sub_table=subset(gene_table, prom==1 )
  write.table(unique(sub_table$name), file=
                paste0(out.string, '_Genes_with_DMCs_in_Promoters.txt' ),
              quote=F, sep="\t", col=F, row=F)
  
  sub_table=subset(count(subset(gene_table, prom==1 )$name), freq>=5)
  write.table(sub_table$x, file=
                paste0(out.string, '_Genes_with_gte_5_DMCs_in_Promoters.txt' ),
              quote=F, sep="\t", col=F, row=F)
  
  sub_table=subset(count(subset(gene_table, prom==1 )$name), freq>=3)
  write.table(sub_table$x, file=
                paste0(out.string, '_Genes_with_gte_3_DMCs_in_Promoters.txt' ),
              quote=F, sep="\t", col=F, row=F)
  
  sub_table=subset(gene_table, prom==1 | exon==1)
  write.table(unique(sub_table$name), file=paste0(out.string, '_Genes_with_DMCs_in_Promoters_or_Exons.txt' ), quote=F, sep="\t", col=F, row=F)
  
  sub_table=subset(gene_table, prom==1 | intron==1)
  write.table(unique(sub_table$name), file=paste0(out.string, '_Genes_with_DMCs_in_Promoters_or_Introns.txt' ), quote=F, sep="\t", col=F, row=F)
  
  
  #pdf(paste(job.string, "Genic_Annotation_of_All_Sites.pdf", sep="_"))
  #plotTargetAnnotation( genic_meth_annotation,precedence=TRUE, main="All Sites")
  #dev.off()
  
  pdf(paste0(out.string,  "_Genic_Annotation_of_Differential_Sites.pdf"))
    plotTargetAnnotation(genic_diff_annotation ,precedence=TRUE, main="Differentially Methylated Sites")
  dev.off()
  
#   cpgi_meth_annotation=annotate.WithFeature.Flank(as(meth,"GRanges"),cpg.obj$CpGi,
#                                                   cpg.obj$shores,
#                                                   feature.name="CpGi",
#                                                   flank.name="shores")
  
  # cpgi_diff_annotation=annotate.WithFeature.Flank( methDiff25p,cpg.obj$CpGi,
  #                                                  cpg.obj$shores,
  #                                                  feature.name="CpGi",
  #                                                  flank.name="shores")
  
  cpgi_diff_annotation=annotateWithFeatureFlank( methDiff25p,cpg.obj$CpGi,
                                                   cpg.obj$shores,
                                                   feature.name="CpGi",
                                                   flank.name="shores")
  
#   pdf(paste(job.id,"CpGi_Annotation_of_All_Sites.pdf", sep="_"))
#   plotTargetAnnotation(cpgi_meth_annotation,col=c("green","gray","white"), main="All Sites")
#   dev.off()
  
  pdf(paste0(out.string,  "_CpGi_Annotation_of_Differential_Sites.pdf"))
    plotTargetAnnotation(cpgi_diff_annotation ,col=c("green","gray","white"), main="Differentially Methylated Sites")
  dev.off()

}