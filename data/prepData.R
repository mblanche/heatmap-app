library(biomaRt)

mart <- useMart("ensembl","dmelanogaster_gene_ensembl")

atts <- c('ensembl_gene_id',
          'external_gene_id',
          'description',
          'gene_biotype')

gene2name <-  getBM(attributes=atts,
                    mart=mart)

gene2name$description <- sub("\\s\\[Source.+","",gene2name$description)

saveRDS(gene2name,"gene2name.rds")
