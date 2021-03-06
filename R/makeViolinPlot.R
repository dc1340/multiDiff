makeViolinPlot<- function(multiDiff, meth.thresh=0, q.thresh=0.01, make.violin=T, make.bar=T, make.plot=T, bw_plot=F){ 
  require(ggplot2)
  require(plyr)
  
  if (!make.violin & !make.bar){ return(NULL)}
  
  cur_terms=names(multiDiff[[3]])
 
  
  if (make.violin){
     violin_dat=data.frame(meth.diff=as.vector(multiDiff[[2]][ , 1 , ]),
                        coef=rep(cur_terms, each=nrow(multiDiff[[1]])),
                        q.value=as.vector(multiDiff[[2]][ , 3 , ])
                        )
  violin_dat=subset(violin_dat, abs(meth.diff)>=meth.thresh & q.value<=q.thresh)

    violin_dat$coef=factor(violin_dat$coef, levels=cur_terms)
    cur_violin=qplot(factor(coef), meth.diff, data = violin_dat, geom = "violin", fill=coef)+ylim(-100, 100)
  
 }
 
  if (make.bar){
    bar_dat=count(violin_dat$coef)
  
    colnames(bar_dat)=c('coef', 'site_count')
  
     
    bar_dat$coef=factor(bar_dat$coef, levels=cur_terms)
    #Changed to work with more recent versions of ggplot
    #cur_bar=qplot(x=factor(coef), y=site_count, data=bar_dat, geom='bar', fill=coef, stat='identity')+theme_set(theme_gray(base_size=18))
    
    cur_bar=qplot(data=violin_dat, x=(coef), geom='bar', fill=coef)+scale_y_continuous(label=scientific)
  }
  
  if (bw_plot){
    cur_violin=cur_violin+theme_bw()
    cur_bar=cur_bar+theme_bw()
    
  }
  if(!make.bar){
    cur_plot=cur_violin
  } else if (!make.violin){
    cur_plot=cur_bar
  } else{
    cur_plot=multiplot(cur_violin, cur_bar)
  }

  

  if (make.plot){ print(cur_plot)}
  


 
  return(cur_plot)
}

makeBarPlot<- function(multiDiff, meth.thresh=0, q.thresh=0.01,  bw_plot=F) { 
  require(ggplot2)
  require(plyr)
  require(scales)
  
  cur_terms=names(multiDiff[[3]])
  
  violin_dat=data.frame(meth.diff=as.vector(multiDiff[[2]][ , 1 , ]),
                        coef=rep(cur_terms, each=nrow(multiDiff[[1]])),
                        q.value=as.vector(multiDiff[[2]][ , 3 , ])
  )
  violin_dat=subset(violin_dat, abs(meth.diff)>=meth.thresh & q.value<=q.thresh)
  
  violin_dat$coef=factor(violin_dat$coef, levels=cur_terms)
  
  
  bar_dat=count(violin_dat$coef)
  
  
  colnames(bar_dat)=c('coef', 'site_count')
  
  
  bar_dat$coef=factor(bar_dat$coef, levels=cur_terms)
  
  #cur_bar=qplot(x=factor(coef), y=site_count, data=bar_dat, geom='bar', fill=coef)+theme_set(theme_gray(base_size=18))
  cur_bar=qplot(data=violin_dat, x=(coef), geom='bar', fill=coef)+scale_y_continuous(label=scientific)
  
  if (bw_plot){
    cur_bar=cur_bar+theme_bw()
    
  }
  
  
  return(cur_bar)
}

scientific_10 <- function(x) {
     parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

scientific <- function(x) {
  parse(text=scientific_format()(x))
}