outdir='p2_output/'
q.thresh=0.01
meth.thresh=25


#cur_table=paste0('p2_output/processed_tables/P2_Wnt_Genes_Hit_by_FA_', meth.thresh,'.txt')
cur_table=paste0('p2_output/processed_tables/P2_Mechanism_Signature_DMC_Count_',meth.thresh,'.txt')
if (file.exists(cur_table)){
  file.remove(cur_table)
}

#Cd analysis
#Standard

cur_multi=p2_cd_jm_inter_diff
job.id='P2_Cd_All_Interactions'
source('R/add_analysis.R')

# job.id='test_fix'
# source('R/add_analysis.R')


##Filter against wildtype only results
wt_coef='IsCd'
cur_multi=filterMultiOnMulti(p2_cd_jm_inter_diff, p2_wt_diff, inv=T, coef=wt_coef)
job.id=paste('P2_Cd_All_Interactions_Filtered_Background', wt_coef, sep="_")
source('R/add_analysis.R')


## Sexcombined

cur_multi=p2_sexcomb_cd_diff
job.id='P2_Cd_Sex_Combined'
source('R/add_analysis.R')

# source('R/quick_genic_annotation.R')
# source('R/quick_multidiff_pipeline.R')
#source('R/quick_heatmap.R')


## Sexcombined, filtered against wildtype

wt_coef='IsCd'
cur_multi=filterMultiOnMulti(p2_sexcomb_cd_diff, p2_wt_diff, inv=T, coef=wt_coef)
job.id=paste('P2_Cd_Sex_Combined_Filtered_Wt_Combined', wt_coef, sep="_")
source('R/add_analysis.R')


#L analysis
#Standard

cur_multi=p2_l_jm_inter_diff
job.id='P2_L_All_Interactions'
source('R/add_analysis.R')
# source('R/quick_genic_annotation.R')
# source('R/quick_multidiff_pipeline.R')
#source('R/quick_heatmap.R')

##Filter against wildtype only results
wt_coef='IsCd'
cur_multi=filterMultiOnMulti(p2_l_jm_inter_diff, p2_wt_diff, coef=wt_coef)
job.id=paste('P2_L_All_Interactions_Filtered_Wt_Combined', wt_coef, sep="_")
source('R/add_analysis.R')


## Sexcombined

cur_multi=p2_sexcomb_l_diff
job.id='P2_L_Sex_Combined'
source('R/add_analysis.R')
# source('R/quick_genic_annotation.R')
# source('R/quick_multidiff_pipeline.R')
#source('R/quick_heatmap.R')

## Sexcombined, filtered against wildtype

wt_coef='IsCd'
cur_multi=filterMultiOnMulti(p2_sexcomb_l_diff, p2_wt_diff, coef=wt_coef)
job.id=paste('P2_L_Sex_Combined_Filtered_Background', wt_coef, sep="_")
source('R/add_analysis.R')


# cur_multi=p2_wt_diff
# job.id='P2_Wt_Combined'
# source('R/add_analysis.R')
# source('R/quick_genic_annotation.R')
# source('R/quick_multidiff_pipeline.R')


## Self intesecting on Is10ppm sites 
##Cd
# cur_multi=p2_cd_jm_inter_diff
# cur_coef='Is10ppm'
# cur_multi=filterMultiOnMulti(p2_cd_jm_inter_diff, p2_cd_jm_inter_diff, coef=cur_coef, meth.thresh = meth.thresh)
# cur_multi[[2]]=cur_multi[[2]][, , 1:3]
# cur_multi[[3]]=cur_multi[[3]][ 1:3]
# 
# job.id=paste('P2_Cd_All_Interactions_SelfIntersect' , cur_coef, sep="_")
# source('R/add_analysis.R')
# 
# ## Self intesecting on Is10ppm sites 
# ##L
# cur_multi=p2_l_jm_inter_diff
# cur_coef='Is10ppm'
# cur_multi=filterMultiOnMulti(p2_l_jm_inter_diff, p2_l_jm_inter_diff, coef=cur_coef, meth.thresh = meth.thresh)
# cur_multi[[2]]=cur_multi[[2]][, , 1:3]
# cur_multi[[3]]=cur_multi[[3]][ 1:3]
# 
# job.id=paste('P2_L_All_Interactions_SelfIntersect' , cur_coef, sep="_")
# source('R/add_analysis.R')


