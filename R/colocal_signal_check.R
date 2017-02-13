meth.thresh=10


cur_multi=p2_cd_jm_inter_diff
job.id='P2_Cd_All_Interactions'

cur_diff_call=getDiffCallMatrix(cur_multi, meth=meth.thresh)
cur_diff=abs(getDiffMatrix(cur_multi))
cur_diff=((cur_diff[, 'Is10ppm']+cur_diff[, 'IsCdHet'])==0)
cur_diff_call=cur_diff_call[ , 'Is10ppm'] & cur_diff_call[ , 'IsCdHet'] & cur_diff
(summary(cur_diff)) 
(summary(cur_diff_call))  


cur_multi=p2_l_jm_inter_diff
job.id='P2_L_All_Interactions'


cur_diff_call=getDiffCallMatrix(cur_multi, meth=meth.thresh)
cur_diff=abs(getDiffMatrix(cur_multi))
cur_diff=((cur_diff[, 'Is10ppm']-cur_diff[, 'IsLHet'])==0)
cur_diff_call=cur_diff_call[ , 'Is10ppm'] & cur_diff_call[ , 'IsLHet'] & cur_diff
(summary(cur_diff)) 
(summary(cur_diff_call))  



cur_multi=p2_sexcomb_cd_diff
job.id='P2_Cd_All_Interactions'

cur_diff_call=getDiffCallMatrix(cur_multi, meth=meth.thresh)
cur_diff=abs(getDiffMatrix(cur_multi))
cur_diff=((cur_diff[, 'Is10ppm']+cur_diff[, 'IsCdHet'])==0)
cur_diff_call=cur_diff_call[ , 'Is10ppm'] & cur_diff_call[ , 'IsCdHet'] & cur_diff
(summary(cur_diff)) 
(summary(cur_diff_call))  


cur_multi=p2_sexcomb_l_diff
job.id='P2_L_All_Interactions'

cur_diff_call=getDiffCallMatrix(cur_multi, meth=meth.thresh)
cur_diff=abs(getDiffMatrix(cur_multi))
cur_diff=((cur_diff[, 'Is10ppm']-cur_diff[, 'IsLHet'])==0)
cur_diff_call=cur_diff_call[ , 'Is10ppm'] & cur_diff_call[ , 'IsLHet'] & cur_diff
(summary(cur_diff)) 
(summary(cur_diff_call))  