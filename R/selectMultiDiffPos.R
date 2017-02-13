selectMultiDiffPosByChr <- function(multiDiff, chrs=NULL, start.g=0, end.g=1e9){
  pos_dat=multiDiff[[1]]
  
  if (is.null(chrs)){
    selected_pos=with(pos_dat, start >= start.g & end <= end.g )
  } else {
    selected_pos=with(pos_dat, chr %in% chrs &  start>- start.g & end <= end.g )
  }
  
  out_diff=multiDiff
  
  out_diff[[1]]=pos_dat[ selected_pos, ]
  out_diff[[2]]=out_diff[[2]][ selected_pos, , ]
  
  return(out_diff)
}

selectMultiDiffPos<- function(multiDiff, selected_pos=NULL){
  
  if (is.null(selected_pos)){
    return(multiDiff)
  } else{
    out_diff=multiDiff
    out_diff[[1]]=out_diff[[1]][ selected_pos,  ]
    out_diff[[2]]=out_diff[[2]][ selected_pos, , ]
    return(out_diff)
  }
  
}

