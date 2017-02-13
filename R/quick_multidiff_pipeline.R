# cur_multi=tdiff
# job.id='test'
# outdir='test'
# 
# meth.thresh=25
# q.thresh=0.01
out.string=paste0(outdir, '/', job.id, '_', meth.thresh)
out.string.nothresh=paste0(outdir, '/', job.id)


#Summary Table
 diff_sum=getSummaryDiffCall(cur_multi, meth=meth.thresh, q=q.thresh)
 write.table(diff_sum, file=paste0(out.string, '_Summary_DMC_Call_Table.txt'), quote=F, sep="\t", col.names = NA)

 #Gradiated DMC Table
 diff_sum=getGradiatedDiff(p2_cd_jm_inter_diff, q=q.thresh)
 write.table(diff_sum, file=paste0(out.string.nothresh, '_Gradiated_DMC_Call_Table.txt'), quote=F, sep="\t", col.names = NA)
 
#Heatmap
png(paste0(out.string, '_Heatmap.png'), width=11, height=8.5, unit='in', res=200)
  makeHeatMap(cur_multi, meth.thresh = meth.thresh, q.thresh = q.thresh, cluster=F)
dev.off()

##Violin Plot
pdf(paste0(out.string, '_Violin_Plot.pdf'), width=11, height=8.5)
  makeViolinPlot(cur_multi, meth.thresh = meth.thresh, q.thresh = q.thresh)
dev.off()

##Venn Diagram
# cur_diff_call=getDiffCallMatrix(cur_multi)
# venn_els=c()
# venn_sets=c()
# venn_rows=rownames(cur_diff_call)
# 
# for (coef in colnames(cur_diff_call)) {
#   venn_els=c(venn_rows[ cur_diff_call[ , coef]==T ] , venn_els)
#   venn_sets=c( rep(coef, sum(cur_diff_call[ , coef]==T)  ) , venn_sets)
# }
# 
# pdf(paste0(outdir, '/', job.id,'_', meth.thresh, '_Venn_Diagram.pdf', sep='_'), width=5, height=5)
# plot(venneuler(data.frame(elements=venn_els, sets=venn_sets)))
# dev.off()
# 


