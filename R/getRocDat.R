getRocDat<-function( multiDiff, site_states,base_q.thresh=1, base_meth.thresh=0){
  cur_terms=names(multiDiff[[3]])
  num_sites=nrow(multiDiff[[1]])
  
  roc_dat=data.frame(matrix(nrow = 203* length(cur_terms), ncol=5) )
  colnames(roc_dat)=c('Term', 'type', 'thresh', 'FPR' , 'TPR' )
  diff_mat=abs(getDiffMatrix(multiDiff))
  q_mat=multiDiff[[2]][ , 3, ]
  cond_pos=colSums(site_states)
  
  base_q_mat=multiDiff[[2]][ , 3, ]<=base_q.thresh  
  
  i=1 
  
  for (thresh in 0:100){
      cur_diff_mat=diff_mat>=thresh
      
      cur_call=cur_diff_mat & base_q_mat
      for (j in 1:length(cur_terms)){
        term=cur_terms[[j]]
        cur_tpr=sum(site_states[ ,  j ]& cur_call[ , j ])/cond_pos[[j]]
        cur_fpr=sum((!site_states[ , j ]) & cur_call[ , j ])/(num_sites-cond_pos[[j]])
          roc_dat[ i,  ]=list(term , 'meth', thresh, cur_fpr, cur_tpr)
        i=i+1
      }
      
  }

  #base_meth_mat=diff_mat>=base_meth.thresh
  for (thresh in seq(-0.1,1, by=0.01)){
    cur_q_mat=q_mat<=thresh
    cur_call=cur_q_mat #& base_meth_mat
    for (j in 1:length(cur_terms)){
      term=cur_terms[[j]]
      cur_tpr=sum(site_states[ ,  j ]& cur_call[ , j ])/cond_pos[[j]]
      cur_fpr=sum((!site_states[ , j ]) & cur_call[ , j ])/(num_sites-cond_pos[[j]])
      roc_dat[ i,  ]=list(term , 'q.value', thresh, cur_fpr, cur_tpr)
      i=i+1
    }
  }

  
  
  return(roc_dat)
}


toBinState <-function( state_matrix){
  rowSums(matrix(as.numeric(state_matrix),
               nrow=nrow(state_matrix)) %*% diag(2^c(0:(ncol(state_matrix)-1))))+1
  
}