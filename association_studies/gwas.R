#### DEMO GWAS - bardzo ogólne wprowadzenie

library(ggplot2)
library(magrittr)

### Fenotypy

phen <- read.table("out/phen_cov", header=F)
colnames(phen) <- c("ID","Subpopulacja","Fenotyp","GT")
str(phen)

phen$Subpopulacja <- as.factor(phen$Subpopulacja)


### Rozkład fenotypu - normalizacja?
ggplot(phen, aes(x=Fenotyp)) + geom_density() + theme_minimal()


### Genotypy

gen <- read.table("out/snp_table", header=T)

### Dataset:
ncol(gen)
nrow(phen)

### W klasycznym gwas mamy z reguły do czynienia z p>>n
### Tutaj dodatkowo niezbalansowany design


### Tabela zliczeń genotypów
sapply(gen, table) 


######## Błędy w genotypach: N's, Singletony, Filtrowanie, imputacja
######## Błędy w fenotypach: Np heterozygotyczność, N's


### Odfiltrowanie multiallelicznych wariantów
biallelic_mask <- sapply(gen, function(x){
  gts <- table(x) %>% names %>% paste0(., collapse = "") %>% strsplit(.,"") %>% table %>% length
  if(gts == 2) TRUE
  else FALSE
  })

gen_biallelic <- gen[,biallelic_mask]

### Odfiltrowanie wariantów z MAF < 5% - Moc statystyczna
maf_mask <- sapply(gen_biallelic, function(x){
  sortab <-  x %>% unlist %>% paste0(collapse ="") %>% strsplit(.,"") %>% table %>% sort
  if((sortab[1]*100)/(sortab[1]+sortab[2]) > 5) TRUE
  else FALSE
})

gen_maf <- gen_biallelic[,maf_mask]


### Odfiltrowanie wariantów dla których nie obserwujemy wszystkich trzech genotypów
### Jedynie w celu uproszczenia przykładu
### (Normalnie np. filtrowanie po HWE)

all_gts_mask <- sapply(gen_maf, function(x){
  tab <- table(x) %>% names %>% sort   ## het zawsze w srodku
  if(length(tab) == 3) TRUE
  else FALSE
})

gen_all_gt <- gen_maf[,all_gts_mask]



### Numeryczne zakodowanie genotypów. Zakładamy addytywny model, więc możemy wykorzystać 0,1,2.

gen_numcoded <- sapply(gen_all_gt, function(x){
    tab <- table(x) %>% names %>% sort   ## het zawsze w srodku (bez phasingu/haplotypowania)
    tr_table <- c(0:2)
    names(tr_table) <- as.vector(tab)
    sapply(x, function(y){
      tr_table[y] %>% as.numeric()
    })
})  %>% data.frame


rownames(gen_numcoded) <- rownames(gen_maf)


### Odfiltrowanie w pełni skorelowanych markerów (można zbinować)

cor_matrix <- cor(gen_numcoded)
cor_matrix[upper.tri(cor_matrix, diag = T)] <- 0

gen_numcoded <- gen_numcoded[,!apply(cor_matrix,2, function(x) any(round(abs(x), digits = 4)==1))] 


### Regresja liniowa dla markera nr 1 
### Pytanie - co jeżeli mamy cechę jakościową?

marker_num = 1
m1 <- lm(Fenotyp ~ gen_numcoded[,marker_num], data=phen)
m1
summary(m1)
plot(m1)



nm <- colnames(gen_numcoded)[marker_num] 
nm
alleles <- gen_all_gt[nm] %>% table %>% sort %>% names

to_plot <- cbind(gen_numcoded[,marker_num], phen$Fenotyp) %>% data.frame
ggplot(to_plot, aes(x=X1, y = X2)) + geom_point() + theme_minimal() + geom_smooth(method="lm", se = F) +
  ylab("Fenotyp") + xlab("Genotyp") + scale_x_continuous(breaks=c(0,1,2), labels=c(alleles[1],alleles[2],alleles[3]))



### N niezależnych modelów liniowych dla każdego z markerów:


out <- sapply(1:ncol(gen_numcoded), function(x){
  df <- cbind(phen$Fenotyp, gen_numcoded[,x] ) %>% data.frame
  colnames(df) <- c("phen","gen")
  m <- lm(phen ~ gen, data=df)
  -log10(summary(m)$coefficients[,4][2])
})


### Manhattan plot

to_plot <- cbind(out, as.numeric(gsub("X","",colnames(gen_numcoded),fixed=T)))

to_plot = as.data.frame(to_plot)

ggplot(to_plot, aes(x=V2, y=out)) + geom_point() + theme_minimal()

### Testowanie Zbioru Hipotez (Multiple Hypothesis Testing) - korekcja Bonferroniego (konserwatywne)

thresh <- -log10(0.01/ncol(gen_numcoded)) 

ggplot(to_plot, aes(x=V2, y=out)) + geom_point() + theme_minimal() + geom_hline(yintercept = thresh, col = "orange")


### Struktura populacji - jeden z głównych czynników zakłucającyh, innymi problemami są np.: efekt poligenowy, odstępstwa od addytywności
### https://pubmed.ncbi.nlm.nih.gov/10673763/ Beware the chopsticks gene


### PCA

gt_pca <- prcomp(gen_numcoded)
autoplot(gt_pca, data=phen, col = 'Subpopulacja') + theme_minimal()


#K = 1-dist(gen_numcoded, method="euclidean") %>% as.matrix 
#heatmap(K)

### Pokrewieństwo (Kinship)

`GAPIT.kinship.VanRaden` <-
  function(snps,hasInbred=TRUE) {
    # Object: To calculate the kinship matrix using the method of VanRaden (2009, J. Dairy Sci. 91:4414???C4423)
    # Input: snps is n individual rows by m snps columns
    # Output: n by n relationship matrix
    # Authors: Zhwiu Zhang
    # Last update: March 2, 2016 
    ############################################################################################## 
    print("Calculating kinship with VanRaden method...")
    #Remove invariants
    fa=colSums(snps)/(2*nrow(snps))
    index.non=fa>=1| fa<=0
    snps=snps[,!index.non]
    
    nSNP=ncol(snps)
    nInd=nrow(snps)
    n=nInd 
    
    ##allele frequency of second allele
    p=colSums(snps)/(2*nInd)
    P=2*(p-.5) #Difference from .5, multiple by 2
    snps=snps-1 #Change from 0/1/2 coding to -1/0/1 coding
    
    print("substracting P...")
    Z=t(snps)-P#operation on matrix and vector goes in direction of column
    print("Getting X'X...")
    #K=tcrossprod((snps), (snps))
    K=crossprod((Z), (Z)) #Thanks to Peng Zheng, Meng Huang and Jiafa Chen for finding the problem
    
    print("Adjusting...")
    adj=2*sum(p*(1-p))
    K=K/adj
    
    print("Calculating kinship with VanRaden method: done")
    
    return(K)
  }
#=============================================================================================



K = GAPIT.kinship.VanRaden(gen_numcoded)
heatmap(K)



### Efekt populacji

to_plot <- cbind(gen_numcoded[,marker_num], phen$Fenotyp, phen$Subpopulacja) %>% data.frame
ggplot(to_plot, aes(x=X1, y = X2)) + 
  geom_point(aes(col=factor(X3)), position=position_dodge(width=0.3)) + 
  geom_smooth(method="lm", aes(col=factor(X3)), se=F) +
  labs(col="Subpopulacja") + 
  ylab("Fenotyp") + xlab("Genotyp") + scale_x_continuous(breaks=c(0,1,2), labels=c(alleles[1],alleles[2],alleles[3])) +
  theme_minimal()


### Pierwsze trzy PC jako współzmienne (kowariaty) (efekt stały, "fixed effect")
### Częściej włącza się macierz pokrewieństwa jako efekt losowy ("random effect") w modelu mieszanym

out <- sapply(1:ncol(gen_numcoded), function(x){
  df <- cbind(phen$Fenotyp, gen_numcoded[,x],  gt_pca$x[,1],  gt_pca$x[,2],  gt_pca$x[,3] ) %>% data.frame
  colnames(df) <- c("phen","gen","P1","P2","P3")
  m <- lm(phen ~ gen + P1 + P2 + P3, data=df)
  -log10(summary(m)$coefficients[,4][2])
})




to_plot <- cbind(out, as.numeric(gsub("X","",colnames(gen_numcoded),fixed=T)))


thresh <- -log10(0.01/ncol(gen_numcoded)) 

ggplot(to_plot, aes(x=V2, y=out)) + geom_point() + theme_minimal() + geom_hline(yintercept = thresh, col = "orange")


### QQ plot

qqplot(-log10(ppoints(ncol(gen_numcoded))), out , 
       xlab = "theoretical", ylab = "observed")
abline(0, 1, col="red")

### Multi-marker association
### W naszym przypadku niemal całkowicie kontroluje strukturę. 

m3 <- lm(phen$Fenotyp~as.matrix(gen_numcoded))
summary(m3)                         ## R2 -  0.81  --> Odziedziczalność, heritability
out <- -log10(summary(m3)$coefficients[,4])

to_plot <- cbind(out[2:length(out)], as.numeric(gsub("X","",colnames(gen_numcoded),fixed=T))) %>% data.frame

thresh <- -log10(0.01) 

ggplot(to_plot, aes(x=X2, y=X1)) + geom_point() + theme_minimal() + geom_hline(yintercept = thresh, col = "orange")


qqplot(-log10(ppoints(ncol(gen_numcoded))), -log10(summary(fullmodel)$coefficients[,4] ), 
       xlab = "theoretical", ylab = "observed")
abline(0, 1, col="red")




nm <- which(out == max(out))
nm
nm2 <- sub("as.matrix(gen_numcoded)","",nm,fixed=TRUE)
alleles <- gen_all_gt[nm] %>% table %>% sort %>% names

to_plot <- cbind(gen_numcoded[,nm], phen$Fenotyp) %>% data.frame
ggplot(to_plot, aes(x=X1, y = X2)) + geom_point() + theme_minimal() + geom_smooth(method="lm", se = F) +
  ylab("Fenotyp") + xlab("Genotyp") + scale_x_continuous(breaks=c(0,1,2), labels=c(alleles[1],alleles[2],alleles[3]))


to_plot <- cbind(gen_numcoded[,nm], phen$Fenotyp, phen$Subpopulacja) %>% data.frame
ggplot(to_plot, aes(x=X1, y = X2)) + 
  geom_point(aes(col=factor(X3)), position=position_dodge(width=0.3)) + 
  geom_smooth(method="lm", aes(col=factor(X3)), se=F) +
  labs(col="Subpopulacja") + 
  ylab("Fenotyp") + xlab("Genotyp") + scale_x_continuous(breaks=c(0,1,2), labels=c(alleles[1],alleles[2],alleles[3])) +
  theme_minimal()

