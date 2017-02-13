makeVennDiagram<-function(multiDiff, meth.thresh=25, q.thresh = 0.01, return.dat=F){
  
  require(venneuler)
  cur_diff_call=getDiffCallMatrix(cur_multi,
                                  meth.thresh=meth.thresh,
                                  q.thresh = q.thresh)
  venn_els=c()
  venn_sets=c()
  venn_rows=rownames(cur_diff_call)

  for (coef in colnames(cur_diff_call)) {
    venn_els=c(venn_rows[ cur_diff_call[ , coef]==T ] , venn_els)
    venn_sets=c( rep(coef, sum(cur_diff_call[ , coef]==T)  ) , venn_sets)
  }
  
  venn_dat=data.frame(elements=venn_els,
                      sets=venn_sets)
  
  if (return.dat){
    return(venn_dat)
  } else{
    vennout=venneuler(venn_dat)
    
    return(vennout)
  }

  
}