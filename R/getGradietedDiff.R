getGradiatedDiff <-function (multiDiff, q.value=0.01, distr=F, remove.zero=F){
  q.calls=multiDiff[[2]][ ,3  ,  ]<q.value
  meth.vals=multiDiff[[2]][ , 1 ,  ]
  
  diff.out=data.frame(matrix(nrow=100, ncol=length(multiDiff[[3]])+1))
  colnames(diff.out)=c('meth.diff', names(multiDiff[[3]] ) )
  
  for (m in 1:100){
    cur_calls=q.calls & (meth.vals>=m)
    diff.out[ m, ]=c(m, colSums(cur_calls))
  }
  
  if (!distr){
    return(diff.out)
  }else{
    distr.out=data.frame(matrix(nrow=100, ncol=length(multiDiff[[3]])+1))
    colnames(distr.out)=c('meth.diff', names(multiDiff[[3]] ) )
    
    distr.out[ 1, ]=diff.out[ 1,  ]
    distr.out[ , 1]=0:99
    
    distr.out[ 2:100 , 2:ncol(distr.out)]=diff.out[ 1:99, 2:ncol(diff.out) ]-
      diff.out[ 2:100,2:ncol(diff.out)  ]
    
    if (remove.zero){
      distr.out=distr.out[ 2:100, ]
    }
    
    return(distr.out)
    
  }
  
}

plotGradDiff<-function(multiDiff, q.value=0.01, distr=F, remove.zero=F){
  cur_dat=getGradiatedDiff(multiDiff, q.value=q.value,
                           distr=distr,
                           remove.zero=remove.zero)
  cur_dat=melt(cur_dat, id=c('meth.diff'))
  cur_plot=ggplot(cur_dat, aes(meth.diff, value, col=variable))+geom_line()
  return(cur_plot)
}

