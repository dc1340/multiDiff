makeDendoDiff <-function (multiDiff, meth.thresh=15){
  require(reshape)
  pos_dat=multiDiff[[1]][ , c(1,2)]

  diff_dat=multiDiff[[2]][ , 1, ]

  
  colnames(diff_dat)=names(multiDiff[[3]])


  plot_dat=cbind(pos_dat, diff_dat)

  plot_dat=melt(plot_dat , id.vars = c('chr', 'start'))
  plot_dat=subset(plot_dat, abs(value)>=meth.thresh)

  ggplot(plot_dat, aes(start, value, col=variable , alpha=0.6))+geom_point()+
    facet_grid(chr~.)+
    theme(strip.text.y=element_text(size=8, angle=0))
}