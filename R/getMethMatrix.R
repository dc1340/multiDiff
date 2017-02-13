getMethMatrix<-function(mbase){
  meth_matrix=matrix(nrow=nrow(mbase), ncol=length(mbase@sample.ids))
  mdata=getData(mbase)
  meth_matrix=mdata[ , mbase@numCs.index]/mdata[  , mbase@coverage.index]
  colnames(meth_matrix)=mbase@sample.ids
  return (meth_matrix)
}