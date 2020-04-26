# requirement 
library(sSeq)
library(edgeR)
library(org.Hs.eg.db)
library(gplots)
library(biomaRt)
library(enrichR)
library(pheatmap)
library(ggplot2)
library(M3C)

# Extraction of data ------------------------------------------------------

seq <-  read.delim('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct')
rownames(seq) <- seq$Name
seq <- subset(seq, select = -c(Name))


# IDs from obtained data --------------------------------------------------

pancreas <- c('GTEX.11TT1.0326.SM.5LUAY','GTEX.ZTPG.1026.SM.5DUWP','GTEX.PX3G.1026.SM.48TZW',
              'GTEX.O5YT.1026.SM.3MJGF','GTEX.QESD.0226.SM.447BH','GTEX.15ER7.1326.SM.6LPK9',
              'GTEX.144GM.0726.SM.79OJQ','GTEX.UPIC.0726.SM.3GADW','GTEX.11EQ9.1026.SM.5H134',
              'GTEX.1QPFJ.2226.SM.EXOJW','GTEX.ZYT6.1326.SM.5E453','GTEX.1399R.0426.SM.5IJE3',
              'GTEX.S33H.1226.SM.4AD69','GTEX.U8XE.2026.SM.3DB8S','GTEX.WFON.0626.SM.4LVLX',
              'GTEX.17MF6.1426.SM.7IGOX','GTEX.1GF9X.1126.SM.7MKHC','GTEX.1CAMQ.1626.SM.7MGXG',
              'GTEX.1H11D.0926.SM.9OSXO','GTEX.Q734.0426.SM.48TZX','GTEX.1B932.1626.SM.731EA',
              'GTEX.WFON.0626.SM.5S2SX','GTEX.15RJ7.0826.SM.6LLI8','GTEX.PWOO.0626.SM.48TZH',
              'GTEX.WYVS.0926.SM.4SOJV','GTEX.ZAB4.1726.SM.5HL8C','GTEX.YEC4.1326.SM.5IFHG',
              'GTEX.Y5LM.0526.SM.4V6G3','GTEX.ZPU1.0226.SM.4WWC9','GTEX.YB5E.0526.SM.4VDSD',
              'GTEX.ZPU1.0226.SM.4WWC9','GTEX.X3Y1.0726.SM.3P5YU','GTEX.WYVS.0926.SM.4SOJV',
              'GTEX.14AS3.0326.SM.5Q5DB','GTEX.ZPU1.0226.SM.4WWC9','GTEX.13FTW.0526.SM.5IFIP',
              'GTEX.ZVP2.0726.SM.59HKY','GTEX.ZAB5.0826.SM.5P9FU','GTEX.QDVN.0926.SM.2I5GL',
              'GTEX.QEL4.1326.SM.447AD','GTEX.ZF2S.0426.SM.4WKGP','GTEX.ZEX8.1026.SM.4WKHE',
              'GTEX.WY7C.1026.SM.4OND3','GTEX.QV44.0426.SM.4R1KF','GTEX.R53T.0426.SM.48FEM',
              'GTEX.S32W.0826.SM.5S2SI','GTEX.S32W.0826.SM.4AD5Z','GTEX.TKQ2.0426.SM.4DXUO',
              'GTEX.VUSG.0326.SM.3GIJ7','GTEX.XXEK.1726.SM.5S2SR','GTEX.VUSG.0326.SM.5S2ST',
              'GTEX.XXEK.1726.SM.4BRVB','GTEX.QV44.0426.SM.4R1KF','GTEX.RM2N.0326.SM.48FD8',
              'GTEX.1K9T9.1726.SM.CXZK1','GTEX.13FLV.0626.SM.5IFEY','GTEX.13SLX.1326.SM.5S2QS',
              'GTEX.17HGU.2026.SM.79OKS','GTEX.17HHE.0526.SM.7DUGR','GTEX.18A7A.1726.SM.7LT93',
              'GTEX.1A8FM.1826.SM.7MGXO','GTEX.13N11.0226.SM.5KM3C','GTEX.1AMFI.0726.SM.731D9',
              'GTEX.14PHX.2026.SM.6872C','GTEX.1B933.1726.SM.731FC','GTEX.139YR.1526.SM.5IFJ1',
              'GTEX.146FH.1826.SM.5QGQ7','GTEX.1C6VQ.1726.SM.7IGLQ','GTEX.1GN1U.1126.SM.9WPPZ',
              'GTEX.132AR.1826.SM.5EGHR', 'GTEX.131XE.1926.SM.5IFER','GTEX.12WSL.0426.SM.5GCNX',
              'GTEX.12WSG.1026.SM.5EGII', 'GTEX.ZF29.1126.SM.4WKGO','GTEX.13O1R.1826.SM.5KM3B',
              'GTEX.13PVR.0726.SM.5S2PX','GTEX.12WSD.1626.SM.5GCNR','GTEX.Y5V5.1026.SM.5LUAH',
              'GTEX.11GSP.0426.SM.5A5KX', 'GTEX.145MO.2126.SM.5Q5CZ','GTEX.ZF29.1126.SM.4WKGO',
              'GTEX.OIZF.1026.SM.2AXUK','GTEX.WQUQ.2126.SM.4OOSO','GTEX.1497J.0426.SM.5Q5CO',
              'GTEX.OOBJ.1026.SM.3NB2L','GTEX.OIZF.1026.SM.7P8QV','GTEX.1J8EW.1926.SM.CGQGR',
              'GTEX.1HCU6.1926.SM.A8NA6','GTEX.1H23P.2026.SM.9OSXY', 'GTEX.1CB4F.1226.SM.7DHKU',
              'GTEX.148VI.0826.SM.5TDDI','GTEX.1B996.1326.SM.731EO','GTEX.1AX8Z.1726.SM.731DE',
              'GTEX.S4Z8.1126.SM.4AD5C','GTEX.18A6Q.1726.SM.7LT9A','GTEX.18A66.2126.SM.7189D',
              'GTEX.15DYW.1626.SM.6LLI1','GTEX.OIZF.1026.SM.2HML5','GTEX.14PJO.1826.SM.69LPR',
              'GTEX.1B98T.2126.SM.79OOG', 'GTEX.ZYY3.0826.SM.5E44R','GTEX.1ICG6.1126.SM.ADEIF')

liver <- c('GTEX.11TT1.1726.SM.5EQLJ','GTEX.ZTPG.1426.SM.51MT3','GTEX.PX3G.0826.SM.48TZS',
           'GTEX.O5YT.0826.SM.3TW8N','GTEX.QESD.2026.SM.447BI','GTEX.15ER7.1826.SM.6LLI7',
           'GTEX.144GM.1326.SM.5LU5E','GTEX.UPIC.0926.SM.4IHLV','GTEX.11EQ9.0526.SM.5A5JZ',
           'GTEX.1QPFJ.1126.SM.E9U5V','GTEX.ZYT6.0626.SM.5E45V','GTEX.1399R.1226.SM.5P9GF',
           'GTEX.S33H.1626.SM.4AD68','GTEX.U8XE.1526.SM.4E3HT','GTEX.WFON.1726.SM.4LVMQ',
           'GTEX.17MF6.1126.SM.7DUFS','GTEX.1GF9X.0426.SM.7MKHN','GTEX.1CAMQ.1426.SM.7MKFL',
           'GTEX.1H11D.0826.SM.9OSWB','GTEX.Q734.0326.SM.48U15','GTEX.1B932.1426.SM.793AN',
           'GTEX.WFON.1726.SM.4LVMQ','GTEX.15RJ7.2026.SM.6LPJ4','GTEX.PWOO.0826.SM.48TCL',
           'GTEX.WYVS.1926.SM.4PQZ2','GTEX.ZAB4.0826.SM.5LU9D','GTEX.YEC4.0826.SM.5P9FV',
           'GTEX.Y5LM.0426.SM.4VBRO','GTEX.ZPU1.0826.SM.57WG2','GTEX.YB5E.0326.SM.5IFHU',
           'GTEX.ZPU1.0826.SM.7DHMN','GTEX.X3Y1.2726.SM.4PQZH','GTEX.WYVS.1926.SM.4SOIC',
           'GTEX.14AS3.0126.SM.5Q5F4','GTEX.ZPU1.0826.SM.73KXI','GTEX.13FTW.1126.SM.5J2NV',
           'GTEX.ZVP2.0626.SM.51MSO','GTEX.ZAB5.0426.SM.5CVMI','GTEX.QDVN.0826.SM.48TZ2',
           'GTEX.QEL4.1226.SM.447A4','GTEX.ZF2S.3026.SM.4WWCH','GTEX.ZEX8.0826.SM.4WKHK',
           'GTEX.WY7C.0726.SM.4ONCB','GTEX.QV44.0326.SM.4R1KD','GTEX.R53T.0326.SM.48FEC',
           'GTEX.S32W.1926.SM.4AD63','GTEX.S32W.1926.SM.4AD63','GTEX.TKQ2.1726.SM.4DXUP',
           'GTEX.VUSG.0126.SM.4KL1X','GTEX.XXEK.1126.SM.4BRUX','GTEX.VUSG.0126.SM.4KL1X',
           'GTEX.XXEK.1126.SM.4BRUX','GTEX.QV44.0326.SM.C1YQX','GTEX.RM2N.1926.SM.48FCU',
           'GTEX.1K9T9.0726.SM.CXZJZ','GTEX.13FLV.0326.SM.5N9DJ','GTEX.13SLX.1226.SM.5S2Q6',
           'GTEX.17HGU.1826.SM.7IGQM','GTEX.17HHE.0126.SM.79398','GTEX.18A7A.1526.SM.72D69',
           'GTEX.1A8FM.1226.SM.7IGMD','GTEX.13N11.0926.SM.5IJG2','GTEX.1AMFI.0826.SM.731DV',
           'GTEX.14PHX.0526.SM.664NW','GTEX.1B933.1226.SM.731DW','GTEX.139YR.0226.SM.5IFEM',
           'GTEX.146FH.1526.SM.5NQBU','GTEX.1C6VQ.1626.SM.79OO1','GTEX.1GN1U.0926.SM.9WPPY',
           'GTEX.132AR.0426.SM.5IFH8','GTEX.131XE.0326.SM.5LZVO','GTEX.12WSL.0226.SM.5CVMJ',
           'GTEX.12WSG.0626.SM.5FQTQ','GTEX.ZF29.2026.SM.DNZYW','GTEX.13O1R.2026.SM.5KM3N',
           'GTEX.13PVR.0126.SM.5S2PY','GTEX.12WSD.1426.SM.5GCN9','GTEX.Y5V5.0926.SM.4VBPZ',
           'GTEX.11GSP.0626.SM.5986T','GTEX.145MO.2326.SM.5NQ9K','GTEX.ZF29.2026.SM.4WWB7',
           'GTEX.OIZF.0826.SM.3MJGO','GTEX.WQUQ.1926.SM.4OOSA','GTEX.1497J.0726.SM.5Q5D1',
           'GTEX.OOBJ.0826.SM.3NB2K','GTEX.OIZF.0826.SM.3MJGO','GTEX.1J8EW.1626.SM.C1YQ8',
           'GTEX.1HCU6.2026.SM.A8NA7','GTEX.1H23P.1826.SM.A8NA5','GTEX.1CB4F.1726.SM.7MGXH',
           'GTEX.148VI.0626.SM.5TDDH','GTEX.1B996.1426.SM.7EPIA','GTEX.1AX8Z.1026.SM.7189R',
           'GTEX.S4Z8.0526.SM.4AD4T','GTEX.18A6Q.2026.SM.718AQ','GTEX.18A66.2026.SM.7189C',
           'GTEX.15DYW.1326.SM.6LPIV','GTEX.OIZF.0826.SM.3MJGO','GTEX.14PJO.1726.SM.68719',
           'GTEX.1B98T.1126.SM.7EPHO','GTEX.ZYY3.0626.SM.5NQ6W','GTEX.1ICG6.1026.SM.B2LY6')

# it appears that not all of indexes could be found in seq
pancreas <- unique(colnames(seq)[match(pancreas, colnames(seq))])
liver <- unique(colnames(seq)[match(liver, colnames(seq))])

# subsetting of dataset
seq_pancreas <- seq[,c('Description',pancreas[!is.na(pancreas)])]
seq_liver <- seq[,c('Description',liver[!is.na(liver)])]


# For pancreas ------------------------------------------------------------

# compute standart edgeR format data
yp <-  DGEList(counts=seq_pancreas[,2:83], genes=rownames(seq_pancreas),
              remove.zeros = T)

# take only genes with HGNC IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=sub("[.][0-9]*","",yp$genes$genes),mart= mart)
m <- match(sub("[.][0-9]*","",yp$genes$genes), G_list$ensembl_gene_id)
yp$genes$EntrezGene <- G_list$hgnc_symbol[m]

# reorder data
new_order <-  order(rowSums(yp$counts), decreasing=TRUE)
yp <- yp[new_order,]

# remove duplicated genes
str(yp$genes$EntrezGene)
d <-  duplicated(yp$genes$EntrezGene)
table(d)
yp <- yp[!d,]

# change counts
yp$samples$lib.size <- colSums(yp$counts)
yp <- calcNormFactors(yp, method = "TMM")

# plot correlation for noise estimation
cor.mtx <- cor(yp$counts, method = "spearman")
heatmap.2(cor.mtx, trace = "none", density.info = "none")
pheatmap(cor.mtx)


# Clustering part ---------------------------------------------------------

d <- dist(t(yp$counts), method = "euclidean")
hc5 <- hclust(d, method = "ward.D2" )
sub_grp <- cutree(hc5, k = 2)

plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 2, border = 2:5)


# Dimensionality reduction ------------------------------------------------

# PCA
pca.res <- prcomp(yp$counts, scale = TRUE,rank. = 3)
pca.res <- data.frame(pca.res$rotation)
ggplot(data = pca.res, aes(x = PC1, y = PC2, label = rownames(pca.res))) +
  geom_point()

# tSNE
m <- match(colnames(seq_pancreas[,2:83]), pancreas_bins$X)
labels_bin <- pancreas_bins$binned[m]
tsne(yp$counts, perplex = 25,seed = 42, labels = labels_bin)


# DEGs analysis -----------------------------------------------------------

# building linear model
fat_percentage <-  factor(labels_bin)
anno <-  data.frame(Sample=colnames(yp),fat_percentage)
design <-  model.matrix(~fat_percentage)
rownames(design) <-  colnames(yp)

yp <- estimateGLMCommonDisp(yp, design, verbose=TRUE)
yp <-  estimateGLMTrendedDisp(yp, design)
yp <-  estimateGLMTagwiseDisp(yp, design)

fit = glmFit(yp, design)
lrt = glmLRT(fit)

# estimating DEGs
tt<-topTags(lrt,n=100000)
tt<-tt$table
DEG<-tt[abs(tt$logFC)>0.5 & tt$FDR<0.1,]
DEG<-data.frame(DEG,Symbol=DEG$EntrezGene,
                stringsAsFactors = F)
DEG<-na.omit(DEG)

DEG.up<-DEG[DEG$logFC>0,8]
DEG.down<-DEG[DEG$logFC<0,8]
DEG.up.down<-DEG$Symbol

# gene set enrichment analysis according to different db
db<-c("KEGG_2019_Human")
enriched <- enrichr(DEG.up.down, db)
DEP<-enriched$KEGG.2019.Human
DEP<-DEP[DEP$Adjusted.P.value<0.1,]
DEP <- DEP[order(DEP$Adjusted.P.value),]
row.names(DEP) <- rownames(DEP[order(DEP$P.value),])

# plot gene set enrichment
ggplot(DEP, aes(x= reorder(Term, -log10(Adjusted.P.value)) , y=-log10(Adjusted.P.value))) + 
  geom_bar(stat='identity', aes(fill=rownames(DEP)), width=.5,
           position="dodge", show.legend = F) +
  coord_flip() +
  ylab("-log10(Adjusted P-value)") + xlab("Enriched pathways") +
  labs(title = "Pancreas fat % for 0-10 vs. 10-80")+
  theme(axis.text.x = element_text(color = "grey20", size = 14),
        axis.text.y = element_text(color = "grey20", size = 14),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))


# For liver ---------------------------------------------------------------

# compute standart edgeR format data
y <-  DGEList(counts=seq_liver[,2:87], genes=rownames(seq_liver),
              remove.zeros = T)

# take only genes with HGNC IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=sub("[.][0-9]*","",y$genes$genes),mart= mart)
m <- match(sub("[.][0-9]*","",y$genes$genes), G_list$ensembl_gene_id)
y$genes$EntrezGene <- G_list$hgnc_symbol[m]

# reorder data
new_order <-  order(rowSums(y$counts), decreasing=TRUE)
y <- y[new_order,]

# remove duplicated genes
str(y$genes$EntrezGene)
d <-  duplicated(y$genes$EntrezGene)
y <- y[!d,]

# change counts
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method = "TMM")

# plot correlation for noise estimation
cor.mtx <- cor(y$counts, method = "spearman")
heatmap.2(cor.mtx, trace = "none", density.info = "none")

# Clustering part ---------------------------------------------------------

d <- dist(t(y$counts), method = "euclidean")
hc5 <- hclust(d, method = "ward.D2" )
sub_grp <- cutree(hc5, k = 2)

plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 2, border = 2:5)

# Dimensionality reduction ------------------------------------------------

# PCA
pca.res <- prcomp(y$counts, scale = TRUE,rank. = 3)
pca.res <- data.frame(pca.res$rotation)
ggplot(data = pca.res, aes(x = PC1, y = PC2, label = rownames(pca.res))) +
  geom_point()

# tSNE
m <- match(colnames(seq_liver[,2:87]), liver_bin$SAMPID_liver)
labels_bins <- liver_bin$binned[m]
tsne(y$counts, perplex = 25,seed = 42, labels = labels_bins)

# DEGs analysis -----------------------------------------------------------

# building linear model
fat_percentage_liv <-  factor(labels_bins)
anno <-  data.frame(Sample=colnames(y),fat_percentage_liv)
design <-  model.matrix(~fat_percentage_liv)
rownames(design) <-  colnames(y)

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <-  estimateGLMTrendedDisp(y, design)
y <-  estimateGLMTagwiseDisp(y, design)

fit = glmFit(y, design)
lrt = glmLRT(fit)
topTags(lrt)

# estimating DEGs
tt<-topTags(lrt,n=100000)
tt<-tt$table
DEG<-tt[abs(tt$logFC)>0.5 & tt$FDR<0.1,]

DEG<-data.frame(DEG,Symbol=DEG$EntrezGene, stringsAsFactors = F)
DEG<-na.omit(DEG)

DEG.up<-DEG[DEG$logFC>0,8]
DEG.down<-DEG[DEG$logFC<0,8]
DEG.up.down<-DEG$Symbol

# gene set enrichment analysis according to different db
db<-c("KEGG_2019_Human")
enriched <- enrichr(DEG.up.down, db)
DEP<-enriched$KEGG_2019_Human
DEP<-DEP[DEP$Adjusted.P.value<0.1,]

# plot gene set enrichment
ggplot(DEP, aes(x= reorder(Term, -log10(Adjusted.P.value)) , y=-log10(Adjusted.P.value))) + 
  geom_bar(stat='identity', aes(fill=rownames(DEP)), width=.5,
           position="dodge", show.legend = F) +
  coord_flip() +
  ylab("-log10(Adjusted P-value)") + xlab("Enriched pathways") +
  labs(title = "Liver fat % for 0-25 vs. 25-80")+
  theme(axis.text.x = element_text(color = "grey20", size = 14),
        axis.text.y = element_text(color = "grey20", size = 14),
        axis.title.x = element_text(color = "grey20", size = 16),
        axis.title.y = element_text(color = "grey20", size = 16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))

# analysis was performed by polinashpudeiko
