refseq=read.table('mouse.refseq.bed.txt', head=F)
refnames=read.table('refseq_id_and_name.mm10.txt')
colnames(refnames)=c('id', 'name')

ntd_genes=read.table('ntd_gene_symbols.txt')$V1
lrp6_genes=read.table('lrp6_genes.txt')$V1


cur_tab=(subset( refseq, V4  %in% refnames$id [ refnames$name %in% ntd_genes ] ))
cur_tab=merge(cur_tab, refnames [ refnames$name %in% ntd_genes,  ], by.x='V4', by.y='id')
cur_tab$V4=cur_tab$name
cur_tab=cur_tab[ , c(2:4, 1,5:13 )]

write.table(cur_tab[ , 1:12] , file='refseq.ntd.mm10.txt', quote=F, col=F, sep="\t" , row=F)

cur_tab=(subset( refseq, V4  %in% refnames$id [ refnames$name %in% lrp6_genes ] ))
cur_tab=merge(cur_tab, refnames [ refnames$name %in% lrp6_genes,  ], by.x='V4', by.y='id')
cur_tab$V4=cur_tab$name
cur_tab=cur_tab[ , c(2:4, 1,5:13 )]
write.table(cur_tab[ , 1:12] , file='refseq.lrp.mm10.txt', quote=F, col=F, sep="\t" , row=F)