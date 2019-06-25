
####### PERMANOVA

# Load libraries
library(microbiome)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/Generate_microbiome_compositional_data')

load(file = 'mice.data.20181031.RData') #phy.genus: a phyloseq object
ind = (sample_data(phy.genus)$Week==13 & 
         (sample_data(phy.genus)$Treatment != 'STAT'))
pseq = subset_samples(phy.genus,ind)
pseq.rel = microbiome::transform(pseq, "compositional")
otu = abundances(pseq.rel)
meta = meta(pseq.rel)

dim(otu) # 37 68; 37 kinds of microbiome, 68 subjects
dim(meta) # 68, 53; 68 subjects and 53 features

permanova <- adonis(t(otu) ~ Treatment,
                    data = meta, permutations=99, method = "bray")

png('2.png')
p = plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Treatment", size = 3)
print(p)
dev.off()

# Check that variance homogeneity assumptions hold 
# (to ensure the reliability of the results)
dist = vegdist(t(otu))
anova(betadisper(dist, meta$Treatment))

# Show coefficients for the top taxa separating the groups
coef = coefficients(permanova)["Treatment1",]
top.coef = coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 15.5, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")


rownames(otu)[rev(order(abs(coef)))[1:10]] = '1'
rownames(otu)[rownames(otu)!=1] = '2'
colnames(otu) = meta$Treatment

X = as.matrix(otu)
otu[1:3,1:3]
save(X, file = 'Xmicrobiom.RData')

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation1/result 0415')

col_types = colnames(X)
row_types = rownames(X)
col_ty = as.numeric(factor(col_types))
row_ty = as.numeric(factor(row_types))

cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)
heatmap(X, col=hmcols,Colv=NA, Rowv=NA,
        ColSideCol=cols[col_ty],
        RowSideColors = cols[row_ty],
        main = 'Fig 1: Raw data (control & PAT)',scale ='none')

heatmap(X, Colv=NA, Rowv=NA, scale='none')

library(gplots)
png('raw heatmap.png')
heatmap.2(X, Rowv=TRUE, Colv=TRUE,trace='none')
dev.off()


png('raw heatmap2.png')
heatmap.2(X, trace='none',
          ColSideCol=cols[col_ty],col=hmcols,
          RowSideColors = cols[row_ty])
dev.off()

png('raw eric heatmap.png')
heatmap.2(sol$U[[ix]], trace='none',col=hmcols,
          ColSideCol=cols[col_ty],
          RowSideColors = cols[row_ty])
dev.off()

png('raw eric smooth heatmap.png')
heatmap.2(M, trace='none',col=hmcols,
          ColSideCol=cols[col_ty],
          RowSideColors = cols[row_ty])
dev.off()


png('admm.png')
heatmap.2(res_admm1$A, trace='none',col=hmcols,
          ColSideCol=cols[col_ty],
          RowSideColors = cols[row_ty])
dev.off()


png('admm2.png')
heatmap.2(MM$M, trace='none',col=hmcols,
          ColSideCol=cols[col_ty],
          RowSideColors = cols[row_ty])
dev.off()


x = rnorm(100, 0,1)
y = x * 2 + rnorm(100, 0, 5)

a = lm(y~x)
summary(a)
