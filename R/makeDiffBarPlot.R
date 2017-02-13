makeDiffBarPlot<- function(multiDiff, meth.thresh=0, make.plot=T, binwidth=NULL){ 
  
  
  
  cur_terms=names(multiDiff[[3]])
  
  
    cur_dat=Filter(function(x){ abs(x)>meth.thresh}, getData(getDiffFromMulti(multiDiff, cur_terms[[1]]))[ , 'meth.diff'])
    
    violin_dat=data.frame(coef=cur_terms[[1]], meth.diff=cur_dat)
    if (length(cur_terms)>1) {
        for (i in 2:length(cur_terms)){
          cur_dat=Filter(function(x){ abs(x)>meth.thresh}, getData(getDiffFromMulti( multiDiff, cur_terms[[i]]))[ , 'meth.diff'])
          violin_dat=rbind(violin_dat, data.frame(coef=cur_terms[[i]], meth.diff=cur_dat))
        }
    }
    
  if (make.plot){
   cur_plot=ggplot(violin_dat,  aes(meth.diff, fill=coef))+geom_bar(colour='black', binwidth=binwidth)+facet_grid(coef~.)
   print(cur_plot)
 }
  
}