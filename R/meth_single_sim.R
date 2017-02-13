require(plyr)
require(Matrix)
require(methylKit)

sim_simp_design=data.frame(Cov1=c(rep(1,8), rep(0,8)) , Cov2=c(rep(c(rep(1,4), rep(0,4)), 2)))

sim_simp_design$IsBoth=with(sim_simp_design, Cov1* Cov2)

sim_rev_design=-sim_simp_design+1
colnames(sim_rev_design)=c('NotCov1', 'NotCov2', 'Neither')
sim_rev_design$Neither=with(sim_rev_design, NotCov1 * NotCov2)

sim_possible_site_states=unique(sim_simp_design)
#sim_possible_site_states
sim_possible_site_states=rbind(sim_possible_site_states, c(1,1,0), c(1,0,1), c(0,1,1), c(0,0,1))

#Choose generative design, and the design used for differential
#sim_design=data.frame(Cov1=sim_simp_design$Cov1,  row.names=rownames(sim_simp_design))

##The generative model is specified in sim_design
sim_design=sim_simp_design

##The model used to calculate the differential 
## is specified with sim_diff_design
sim_diff_design=sim_design

sim_read_depth=100
sim_num_sites_per_cond=100

sim_base_meth=0.5

sim_mean_beta=1
sim_beta_multiplier=c(1,1.5,1.2) #Allows for covariate to have differing effects

sim_condition_noise=2
sim_biological_noise=0.05

sim_base_linear_meth=logit(sim_base_meth)

sim_num_sites=nrow(sim_possible_site_states)*sim_num_sites_per_cond

sim_num_samps=nrow(sim_design)


#Prepare data matrix
sim_dat=data.frame(matrix(nrow=sim_num_sites, ncol=4+3*sim_num_samps))


sim_cov_cols=4+seq(1, by=3, length.out = sim_num_samps)
sim_c_cols=4+seq(2, by=3, length.out = sim_num_samps)
sim_t_cols=4+seq(3, by=3, length.out = sim_num_samps)

base_cols=c('chr', 'start', 'end', 'strand')
colnames(sim_dat)[ 1:4]=base_cols
sim_dat$chr=rep(paste0('chr', 1:nrow(sim_possible_site_states)), each=sim_num_sites_per_cond)
sim_dat$start=1:sim_num_sites
sim_dat$end=1:sim_num_sites
sim_dat$strand='+'

sim_dat[ , sim_cov_cols]=sim_read_depth

colnames(sim_dat)[ sim_cov_cols]=paste0('coverage' , 1:sim_num_samps)
colnames(sim_dat)[ sim_c_cols]=paste0('numCs' , 1:sim_num_samps)
colnames(sim_dat)[ sim_t_cols]=paste0('numTs' , 1:sim_num_samps)

sim_site_states=Matrix( nrow=sim_num_sites, ncol=ncol(sim_design))

for ( i in 1:nrow(sim_possible_site_states)){
  sim_start=((i-1)*sim_num_sites_per_cond)+1
  sim_end=((i-1)*sim_num_sites_per_cond)+sim_num_sites_per_cond
  sim_cur_dat=Matrix(unlist(rep(sim_possible_site_states[ i , ],sim_num_sites_per_cond)), nrow=sim_num_sites_per_cond, ncol=ncol(sim_design), byrow = T)
  
  sim_site_states[ sim_start:sim_end , ]=sim_cur_dat
}
sim_beta_multipler_matrix=matrix((rep(sim_beta_multiplier, sim_num_sites)), nrow=sim_num_sites, ncol=ncol(sim_design), byrow=T)

sim_raw_betas=Matrix(rnorm(sim_num_sites_per_cond*ncol(sim_design), sim_mean_beta, sim_condition_noise ), nrow=sim_num_sites, ncol=ncol(sim_design))
sim_betas=sim_raw_betas *sim_beta_multipler_matrix* sim_site_states 
sim_linear_predictors=sim_base_linear_meth+(sim_betas %*% t(sim_design))+matrix(rnorm(sim_num_sites*sim_num_samps, 0, sim_biological_noise), nrow=sim_num_sites, ncol=sim_num_samps)

for (i in 1:sim_num_samps){
  
  sim_cur_c_count=aaply(ilogit(sim_linear_predictors[ , i]), 1, function (x) { sum(rbinom(n=sim_read_depth, prob = x, size=1))})
  sim_cur_t_count=sim_read_depth-sim_cur_c_count
  sim_dat[ , sim_c_cols[[i]]]=sim_cur_c_count
  sim_dat[ , sim_t_cols[[i]]]=sim_cur_t_count
}

sim_sample.ids=paste0('sample',1:sim_num_samps)

sim_meth_base=new("methylBase",sim_dat,
                  sample.ids=sim_sample.ids,assembly='mm10'
                  ,context='CpG',
                  treatment=rep(c(1,0), each=sim_num_samps/2),
                  coverage.index=sim_cov_cols, 
                  numCs.index=sim_c_cols,
                  numTs.index=sim_t_cols,
                  destranded=T,
                  resolution='base')

sim_diff=calculateMultiDiffMeth(sim_meth_base, sim_diff_design, paste(colnames(sim_diff_design), collapse = '+') , num.cores = 2);


makeHeatMap(sim_diff, cluster_cols = F)
sim_diff_matrix=getDiffMatrix(sim_diff)
