
#' ---
#' title: Real data KRAS
#' subtitle: Definizione dataset KRAS per verificare normalita'
#' author:
#' date:
#' output:
#'    html_document:
#'      toc: true
#' ---

knitr::opts_chunk$set(echo = F)
knitr::opts_chunk$set(warnings=F)
knitr::opts_chunk$set(message =F)


#---- Libraries -----------------------------------------------------

library(parallel)
library(Matrix)
library(graphite)
library(SourceSet2)
library(ggplot2)
library(tidyr)
library(gRbase)
library(reshape2)
library(scry)
library(ggsci)
library(ggpubr)
library(stringr)
library(EDASeq)
library(sctransform)
library(SingleCellExperiment)
library(BDgraph)
theme_set(theme_bw())



path = 'data/Kras'



#---- DATA ---------------------------------------------------------------------

# Decidiamo di tenere i cluster 10, 3, 5, 4, 1 di tutte e tre le linee cellulari.
cell.annotation = read.csv(paste0(path,'/GSE137912_cell.annotation.csv.gz'))

# Totale di cellule nel cluster
sum(cell.annotation$Cluster==1)
sum(cell.annotation$Cluster==3)
sum(cell.annotation$Cluster==4)
sum(cell.annotation$Cluster==5)
sum(cell.annotation$Cluster==10)

# Barcode associati a queste cellule
# 1
string = cell.annotation[which(cell.annotation$Cluster==1), ]$X
idcell_1 = read.table(text=string, sep='_')[,-2]
# 3
string = cell.annotation[which(cell.annotation$Cluster==3), ]$X
idcell_3 = read.table(text=string, sep='_')[,-2]
# 4
string = cell.annotation[which(cell.annotation$Cluster==4), ]$X
idcell_4 = read.table(text=string, sep='_')[,-2]
# 5
string = cell.annotation[which(cell.annotation$Cluster==5), ]$X
idcell_5 = read.table(text=string, sep='_')[,-2]
# 10
string = cell.annotation[which(cell.annotation$Cluster==10), ]$X
idcell_10 = read.table(text=string, sep='_')[,-2]

# all e aggiungo var cluster
idcell_keep = cbind(rbind(idcell_1, idcell_3, idcell_4, idcell_5, idcell_10),
                    'cluster' = c(rep(1, nrow(idcell_1)),
                                  rep(3, nrow(idcell_3)),
                                  rep(4, nrow(idcell_4)),
                                  rep(5, nrow(idcell_5)),
                                  rep(10,nrow(idcell_10))))
rm(string, cell.annotation)


## Carico i dati di tutte le linee cellulari

## H0 
# H358
H358 = as.matrix(readMM(file = paste0(path,'/GSM4094251_H358.1.matrix.mtx.gz')))
colnames(H358) = read.table(paste0(path,'/GSM4094251_H358.1.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358) = read.table(paste0(path,'/GSM4094251_H358.1.genes.tsv.gz'), header=F)[,1]
# H2122
H2122 = as.matrix(readMM(file = paste0(path,'/GSM4094255_H2122.0.matrix.mtx.gz')))
colnames(H2122) = read.table(paste0(path,'/GSM4094255_H2122.0.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122) = read.table(paste0(path,'/GSM4094255_H2122.0.genes.tsv.gz'), header=F)[,1]
# SW1573
SW1573 = as.matrix(readMM(file = paste0(path,'/GSM4094259_SW1573.0.matrix.mtx.gz')))
colnames(SW1573) = read.table(paste0(path,'/GSM4094259_SW1573.0.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573) = read.table(paste0(path,'/GSM4094259_SW1573.0.genes.tsv.gz'), header=F)[,1]

## Tengo solo i samples dei cluster di interesse
# H358
idcol = read.table(text = colnames(H358), sep='-')[,1]
H358 = H358[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H358']]
# H2122
idcol = read.table(text = colnames(H2122), sep='-')[,1]
H2122 = H2122[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H2122']]
# SW1573
idcol = read.table(text = colnames(SW1573), sep='-')[,1]
SW1573 = SW1573[, idcol %in% idcell_keep$V3[idcell_keep$V1=='SW1573']]


## H4
# H358
H358_1 = as.matrix(readMM(file = paste0(path,'/GSM4094252_H358.2.matrix.mtx.gz')))
colnames(H358_1) = read.table(paste0(path,'/GSM4094252_H358.2.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358_1) = read.table(paste0(path,'/GSM4094252_H358.2.genes.tsv.gz'), header=F)[,1]
# H2122
H2122_1 = as.matrix(readMM(file = paste0(path,'/GSM4094256_H2122.4.matrix.mtx.gz')))
colnames(H2122_1) = read.table(paste0(path,'/GSM4094256_H2122.4.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122_1) = read.table(paste0(path,'/GSM4094256_H2122.4.genes.tsv.gz'), header=F)[,1]
# SW1573
SW1573_1 = as.matrix(readMM(file = paste0(path,'/GSM4094260_SW1573.4.matrix.mtx.gz')))
colnames(SW1573_1) = read.table(paste0(path,'/GSM4094260_SW1573.4.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573_1) = read.table(paste0(path,'/GSM4094260_SW1573.4.genes.tsv.gz'), header=F)[,1]

## Tengo solo i samples di interesse
# H358
idcol = read.table(text = colnames(H358_1), sep='-')[,1]
H358_1 = H358_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H358']]
# H2122
idcol = read.table(text = colnames(H2122_1), sep='-')[,1]
H2122_1 = H2122_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H2122']]
# SW1573
idcol = read.table(text = colnames(SW1573_1), sep='-')[,1]
SW1573_1 = SW1573_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='SW1573']]

# table(rownames(H358) == rownames(H358_1))
# table(rownames(H2122) == rownames(H2122_1))
# table(rownames(SW1573) == rownames(SW1573_1))

H358  = cbind(H358, H358_1); rm(H358_1)
H2122 = cbind(H2122, H2122_1); rm(H2122_1)
SW1573= cbind(SW1573, SW1573_1); rm(SW1573_1)


## H24
# H358
H358_1 = as.matrix(readMM(file = paste0(path,'/GSM4094253_H358.3.matrix.mtx.gz')))
colnames(H358_1) = read.table(paste0(path,'/GSM4094253_H358.3.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358_1) = read.table(paste0(path,'/GSM4094253_H358.3.genes.tsv.gz'), header=F)[,1]
# H2122
H2122_1 = as.matrix(readMM(file = paste0(path,'/GSM4094257_H2122.24.matrix.mtx.gz')))
colnames(H2122_1) = read.table(paste0(path,'/GSM4094257_H2122.24.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122_1) = read.table(paste0(path,'/GSM4094257_H2122.24.genes.tsv.gz'), header=F)[,1]
# SW1573
SW1573_1 = as.matrix(readMM(file = paste0(path,'/GSM4094261_SW1573.24.matrix.mtx.gz')))
colnames(SW1573_1) = read.table(paste0(path,'/GSM4094261_SW1573.24.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573_1) = read.table(paste0(path,'/GSM4094261_SW1573.24.genes.tsv.gz'), header=F)[,1]

## Tengo solo i samples di interesse
# H358
idcol = read.table(text = colnames(H358_1), sep='-')[,1]
H358_1 = H358_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H358']]
# H2122 -  nessun sample di interesse
# idcol = read.table(text = colnames(H2122_1), sep='-')[,1]
# H2122_1 = H2122_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H2122']]
# SW1573
idcol = read.table(text = colnames(SW1573_1), sep='-')[,1]
SW1573_1 = SW1573_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='SW1573']]

# table(rownames(H358) == rownames(H358_1))
# table(rownames(H2122) == rownames(H2122_1))
# table(rownames(SW1573) == rownames(SW1573_1))

H358  = cbind(H358, H358_1); rm(H358_1)
# H2122 = cbind(H2122, H2122_1); rm(H2122_1)
SW1573= cbind(SW1573, SW1573_1); rm(SW1573_1)


## H72
# H358
H358_1 = as.matrix(readMM(file = paste0(path,'/GSM4094254_H358.5.matrix.mtx.gz')))
colnames(H358_1) = read.table(paste0(path,'/GSM4094254_H358.5.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358_1) = read.table(paste0(path,'/GSM4094254_H358.5.genes.tsv.gz'), header=F)[,1]
# H2122
H2122_1 = as.matrix(readMM(file = paste0(path,'/GSM4094258_H2122.72.matrix.mtx.gz')))
colnames(H2122_1) = read.table(paste0(path,'/GSM4094258_H2122.72.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122_1) = read.table(paste0(path,'/GSM4094258_H2122.72.genes.tsv.gz'), header=F)[,1]
# SW1573
SW1573_1 = as.matrix(readMM(file = paste0(path,'/GSM4094262_SW1573.72.matrix.mtx.gz')))
colnames(SW1573_1) = read.table(paste0(path,'/GSM4094262_SW1573.72.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573_1) = read.table(paste0(path,'/GSM4094262_SW1573.72.genes.tsv.gz'), header=F)[,1]

## Tengo solo i samples di interesse
# H358
idcol = read.table(text = colnames(H358_1), sep='-')[,1]
H358_1 = H358_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H358']]
# H2122
idcol = read.table(text = colnames(H2122_1), sep='-')[,1]
H2122_1 = H2122_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='H2122']]
# SW1573
idcol = read.table(text = colnames(SW1573_1), sep='-')[,1]
SW1573_1 = SW1573_1[, idcol %in% idcell_keep$V3[idcell_keep$V1=='SW1573']]

# table(rownames(H358) == rownames(H358_1))
# table(rownames(H2122) == rownames(H2122_1))
# table(rownames(SW1573) == rownames(SW1573_1))

H358  = cbind(H358, H358_1); rm(H358_1)
H2122 = cbind(H2122, H2122_1); rm(H2122_1)
SW1573= cbind(SW1573, SW1573_1); rm(SW1573_1)


## Unisco i dataset delle tre linee cellulari
cell_line = as.factor(c(rep('H358', ncol(H358)), rep('H2122', ncol(H2122)), rep('SW1573', ncol(SW1573))))
# controllo se ci sono barcode duplicati nelle tre linee cellulari: 32
# quindi cambio nome ai barcode
# table(duplicated(c(colnames(H358), colnames(SW1573), colnames(H2122))))
colnames(H358) = paste0('H358_', 1:ncol(H358))
colnames(H2122) = paste0('H2122_', 1:ncol(H2122))
colnames(SW1573) = paste0('SW1573_', 1:ncol(SW1573))

## Geni in comune (ne perdo 1506)
#table(rownames(H358) %in% rownames(SW1573))
rows_keep = intersect(rownames(H358), rownames(SW1573))
rows_keep = intersect(rows_keep, rownames(H2122))
SW1573 = SW1573[rows_keep, ]; H358 = H358[rows_keep, ]; H2122 = H2122[rows_keep, ]

#table(rownames(H358) ==  rownames(SW1573))
#table(rownames(H358) ==  rownames(H2122))
#table(rownames(SW1573) ==  rownames(H2122))

counts = cbind(H358, SW1573, H2122)
rm(rows_keep, H358, SW1573, H2122)
#dim(counts) # 31232 geni e 5505 samples

# Ci sono geni con il nome duplicato? no
table(duplicated(rownames(counts)))

## Togliamo tutti i geni che hanno tutte le conte=0 (9410)
table(rowSums(counts>0)>0)
counts = counts[which(rowSums(counts>0)>0), ]
dim(counts)

# elimino anche quelli che hanno meno di 10 conte diverse da zero
# table(rowSums(counts>0)>10)
# counts = counts[which(rowSums(counts>0)>10), ]
# dim(counts)


#---- ENTREZID ----------------------------------------------------------

#' I geni nei pathway sono codificati con ENTREZID, quindi trasformiamo i nomi nel nostro dataset.
#' Scarichiamo da biomRt la tabella con le corrispondnze.
#' (una sessantina di geni hanno doppia corrispondenza con entrezid, tolgo la seconda)

## Carico dataset con nome geni
names_converter = read.csv(paste0(path,'/names_converter_KIM.csv'), he=T)[,-1]

# Elimino questi perche' non hanno corrispondenza in mapped
names_converter$entrezgene_id[names_converter$entrezgene_id=='107986084'] = NA
names_converter$entrezgene_id[names_converter$entrezgene_id=='723788'] = NA
names_converter$entrezgene_id[names_converter$entrezgene_id=='84953'] = NA

# geni che non hanno la corrispondenza ensemble-entrezid e sono NA, li tolgo
names_converter = names_converter[-which(is.na(names_converter$entrezgene_id)), ] 

# geni che hanno doppia corrispondenza ensemble-entrezid, tengo la prima entry
table(duplicated(names_converter$ensembl_gene_id))
# names_converter[which(duplicated(names_converter$ensembl_gene_id)), ]
names_converter = names_converter[-which(duplicated(names_converter$ensembl_gene_id)), ]

# library(org.Hs.eg.db)
# mapped.genes.symbol <- as.list(org.Hs.egSYMBOL[as.character(names_converter$entrezgene_id)])
# 
# geneREC = c('HBEGF', 'EGFR', 'SHP2', 'GRB2',
#             'SOS1',  'KRAS', 'PTEN', 'RAF',
#             'PI3K',  'P21',  'RB1',  'G1', 'S')
# 
# mapped.genes.symbol[sapply(mapped.genes.symbol, function(x) x %in% geneREC)]


# geni che hanno doppia corrispondenza entrez-ensemble, tolgo dopo il meno variabile
table(duplicated(names_converter$entrezgene_id))
# names_converter[duplicated(names_converter$entrezgene_id),]

# geni che hanno una corrispondenza
table(rownames(counts) %in% names_converter$ensembl_gene_id)

# Nel dataset counts tengo solo le righe che hanno poi una corrispondenza
counts = counts[which(rownames(counts) %in% names_converter$ensembl_gene_id),]

# Ordino i nomi secondo l'ordine del dataset
names_converter = names_converter[order(match(names_converter$ensembl_gene_id, rownames(counts))), ]
# table(names_converter$ensembl_gene_id == rownames(counts))

# Assegno nuovo nome
rownames(counts) = names_converter$entrezgene_id

# Elimino i duplicati
counts = counts[-which(duplicated(rownames(counts))), ]

dim(counts)


#---- PATHWAYS -----------------------------------------------------------------------------------
#' #### Pathways

# Ora scarico i pathways su cui poi applichero' sourceset.
pathKEGG = pathways("hsapiens", "kegg")[] # specie e database
pathWIKI = pathways("hsapiens", "wikipathways")[] # specie e database

graphsKEGG = lapply(pathKEGG, function(p) pathwayGraph(p)) # lista con tutti i pathways (graphs) 321
graphsWIKI = lapply(pathWIKI, function(p) pathwayGraph(p)) # lista con tutti i pathways (graphs) 601

for(i in 1:length(graphsKEGG)) graph::nodes(graphsKEGG[[i]]) = gsub("ENTREZID:","", graph::nodes(graphsKEGG[[i]]))
for(i in 1:length(graphsWIKI)) graph::nodes(graphsWIKI[[i]]) = gsub("ENTREZID:","", graph::nodes(graphsWIKI[[i]]))

# Eliminiamo dai grafi i nodi non presenti nei nostri dataset
new_nodesK = lapply(graphsKEGG, function(x,y){nodes(x)[which(nodes(x) %in% rownames(y))]}, y=counts)
new_nodesW = lapply(graphsWIKI, function(x,y){nodes(x)[which(nodes(x) %in% rownames(y))]}, y=counts)

graphsKEGG = mapply(graphsKEGG, new_nodesK, FUN = function(x,y){ x = graph::subGraph(y,x)})
graphsWIKI = mapply(graphsWIKI, new_nodesW, FUN = function(x,y){ x = graph::subGraph(y,x)})


# Elimino i grafi con zero nodi
nodiK = sapply(graphsKEGG, function(g) length(nodes(g))!=0)
graphsKEGG = graphsKEGG[nodiK]

nodiW = sapply(graphsWIKI, function(g) length(nodes(g))!=0)
graphsWIKI = graphsWIKI[nodiW]

#' Abbiamo un totale di 316 pathways KEGG e 284 WikiPAth.

#' Tengo nel dataset solo i geni presenti nei pathways
new_nodesK = unique(unlist(sapply(graphsKEGG, function(x){nodes(x)}), use.names = FALSE))
new_nodesW = unique(unlist(sapply(graphsWIKI, function(x){nodes(x)}), use.names = FALSE))

# sum(new_nodesK %in% new_nodesW)
# length(new_nodesK); length(new_nodesW)

new_nodes = unique(c(new_nodesK, new_nodesW))

counts = counts[rownames(counts) %in% new_nodes, ]
rm(new_nodes, pathways, i)
dim(counts)


#' Seleziono i pathways con i geni che ci interessano

library(org.Hs.eg.db)
mapped.genes.symbol <- as.list(org.Hs.egSYMBOL[rownames(counts)])
mapped.genes = sapply(mapped.genes.symbol, function(x) x)

# list of genes of interest
gene.name = c('HBEGF', 'EGFR', 'SHP2', 'GRB2',
            'SOS1', 'KRAS', 'PTEN', 'RAF1', 'PIK3CA',
            'PI3K', 'p21', 'RB1', 'PTPN11', 'CDKN1A')

map.genes = sapply(mapped.genes.symbol[sapply(mapped.genes.symbol, function(x) x %in% gene.name)], function(y) y)
map.genes = data.frame('g.name' = map.genes, 'g.entr'=names(map.genes))
rownames(map.genes) = NULL
map.genes


# Seleziono i grafi KEGG con questi dentro
new_nodesK = lapply(graphsKEGG, function(x,y){nodes(x)[which(nodes(x) %in% map.genes$g.entr)]}, y=counts)
graphs_keep = unname(sapply(new_nodesK, function(x) length(x)))
max(graphs_keep)

# Questi grafi hanno almeno 6 degli 11 geni
graphs_int_KEGG = graphsKEGG[graphs_keep>=6]
names(graphs_int_KEGG)

# geni contenuti nei pathways: tutti ma mai insieme
nodesG = lapply(graphs_int_KEGG, function(x,y){nodes(x)[which(nodes(x) %in% map.genes$g.entr)]}, y=counts)
unique(unname(unlist(nodesG)))


# Seleziono i grafi WIKI con questi dentro
new_nodesW = lapply(graphsWIKI, function(x,y){nodes(x)[which(nodes(x) %in% map.genes$g.entr)]}, y=counts)
graphs_keep = unname(sapply(new_nodesW, function(x) length(x)))
max(graphs_keep)

# Questi grafi hanno almeno 6 degli 11 geni
graphs_int_WIKI = graphsWIKI[graphs_keep>=6]
names(graphs_int_WIKI)

# geni contenuti nei pathways: tutti tranne rb1
nodesG = lapply(graphs_int_WIKI, function(x,y){nodes(x)[which(nodes(x) %in% map.genes$g.entr)]}, y=counts)
map.genes$g.name[!(map.genes$g.entr %in% unique(unname(unlist(nodesG))))]


# Graphs with KRAS

#' ## Pathways with KRAS (3845)
# selector = function(x){ if(3845 %in% nodes(x)){ x } }
# graphsKRAS = lapply(graphs, selector)
# graphsKRAS = graphsKRAS[-which(sapply(graphsKRAS, is.null))]
# names(graphsKRAS)[order(names(graphsKRAS))]

# Seleziono solo i geni che ci sono anche nel dataset TUTTI
#unname(sapply(graphsKRAS, function(x) sum(nodes(x) %in% rownames(counts))/length(nodes(x))))


rm(list=ls()[!ls() %in% c('counts','graphsKEGG', 'graphsWIKI','graphs_int_KEGG','graphs_int_WIKI',
                         'idcell_keep','map.genes','gene.name','path')])

save.image(paste0(path,'/Kras_dataset.RData'))

