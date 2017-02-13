require(methylKit)
cpg.obj=read.feature.flank('/Users/dhruvachandramohan/GitHub/multiDiff/annotations/mouse.cpgi.bed', feature.flank.name = c('CpGi', 'shores') )
gene.obj=read.transcript.features('/Users/dhruvachandramohan/GitHub/multiDiff/annotations/mouse.refseq.bed.txt')
mm10_refseq=read.table('/Users/dhruvachandramohan/GitHub/multiDiff/annotations/refseq_id_and_name.mm10.txt', head=F)
colnames(mm10_refseq)=c('rsid', 'name')

