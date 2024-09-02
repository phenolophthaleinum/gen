#### DEMO bulk-segregant-analysis ("BSA-seq", "QTL-seq")


library(ggplot2)
library(magrittr)

### Fenotypy
phen <- read.table("out/phen_cov", header=F)
head(phen)
str(phen)
summary(phen)

colnames(phen) <- c("ID","Subpopulacja","Fenotyp","GT")

plot(phen$GT, phen$Fenotyp)  ## modulowany przez "penetrance" w skrypcie + error

lm(Fenotyp~GT, phen) %>% summary  ## R^2 - odziedziczalność (heritability)

### Rozkład fenotypu
ggplot(phen, aes(x=Fenotyp)) + geom_density() + theme_minimal()

#high_bulk <- phen[order(phen$Fenotyp),]$ID %>% tail(30)
#low_bulk <- phen[order(phen$Fenotyp),]$ID %>% head(30)


### Wczytajmy warianty (output z bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\t[%AD\t]\n' variants.vcf )
### Warto sprawdzić, czy kolejność próbek zgodna z headerem vcf...

var_data <- read.table("out/var_data", header=F)
colnames(var_data) <- c("chrom","pos","ref","alt","gt_p2","gt_p1","gt_l","gt_h","ad_p2","ad_p1","ad_l","ad_h" )
head(var_data)

### Odfiltrowanie wariantów multiallelicznych (oddzielone przecinkiem)
var_data <- var_data[sapply(var_data$alt, function(x){!grepl(",",x,fixed=T)}),]

### Odfiltrowanie wariantów homomorficznych (u rodziców) - uwaga, jeżeli segregują w populacji - prawdopodpbinie błąd, imputacja
var_data <- var_data[var_data$gt_p1 != var_data$gt_p2,]

### Wektor określający haplotyp rodzicielski względem referencji (dla p2)
p2_phase <- sapply(var_data$gt_p2, function(x){
  unlist(strsplit(x, "/"))[1] %>% as.numeric + 1
})

### Tabelka z pozycją, głębokością alleli (Ref,Alt) i wektorem haplotypu
indices <- cbind(var_data$pos , var_data$ad_l, var_data$ad_h, p2_phase) %>% data.frame
colnames(indices) <- c("pos", "low","high","phase")
rownames(indices) <- rownames(var_data)
indices %>% head

### SNP index dla próby z niską wartością fenotypu
indices$low_index <- apply(indices, 1, function(x){
  counts <- unlist(strsplit(x[["low"]], ",")) %>% as.numeric
  counts[as.numeric(x[["phase"]])]/sum(counts)
})

indices

### SNP index dla próby z wysoką wartością fenotypu

indices$high_index <- apply(indices, 1, function(x){
  counts <- unlist(strsplit(x[["high"]], ",")) %>% as.numeric
  counts[as.numeric(x[["phase"]])]/sum(counts)
})


indices

### delta SNP index 

indices$d_index <- indices$low_index - indices$high_index

### Plot 

ggplot(indices, aes(x=as.numeric(pos), y=d_index)) + geom_line() + theme_minimal()

#### CI (np. permutacja, G-score)

