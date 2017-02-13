png(paste0(outdir, '/', job.id, '_', meth.thresh, '_Heatmap.png'), width=11, height=8.5, unit='in', res=200)
  makeHeatMap(cur_multi, meth.thresh = meth.thresh, q.thresh = q.thresh, cluster=F)
dev.off()