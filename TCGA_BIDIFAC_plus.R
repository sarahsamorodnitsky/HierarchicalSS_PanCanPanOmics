#Load pan-cancer expression, proteomic, miRNA, and methylation files
#First, download the files below from https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin
pan.exp <- read.table('~/PanTCGA/PanCanExp.tsv', sep='\t', header=TRUE)
pan.prot <- read.table('~/PanTCGA/TCGA-RPPA-pancan-clean.txt', header=TRUE)
pan.mirna <- read.csv('~/PanTCGA/PanCan_miRNA.csv', header=TRUE)
pan.meth <-  read.table('~/PanTCGA/PANCAN_meth_merged.tsv', sep='\t', header=TRUE)

#save feature labels
genes <- pan.exp[,1]
cpgs <- pan.meth[,1]
mirnas <- pan.mirna[,1]
#find sample IDs common to all datasets and match across datasets
ids.exp <- names(pan.exp)
ids.prot <- as.character(pan.prot$SampleID)
ids.mirna <-  names(pan.mirna)
ids.meth <- names(pan.meth)
ids.exp <- gsub('\\.','-',ids.exp)
ids.meth <- gsub('\\.','-',ids.meth)
ids.mirna <- gsub('\\.','-',ids.mirna)
ids.exp <- substr(ids.exp,1,12)
ids.prot <- substr(ids.prot,1,12)
ids.mirna <- substr(ids.mirna,1,12)
ids.meth <- substr(ids.meth,1,12)
ids.int <- intersect(ids.exp,ids.prot)
ids.int <- intersect(ids.int,ids.mirna)
ids.int <- intersect(ids.meth,ids.int)
match.exp <- match(ids.int,ids.exp)
exp.int <- pan.exp[,match.exp]
match.meth <- match(ids.int,ids.meth)
meth.int <- pan.meth[,match.meth]
match.mirna <- match(ids.int,ids.mirna)
mirna.int <- pan.mirna[,match.mirna]
tum.type <- pan.prot$TumorType  ##use proteomic IDs to get tumor types
match.prot <- match(ids.int,ids.prot)
tum.type <- tum.type[match.prot]
prot.mat <- t(as.matrix(pan.prot[,3:200]))
colnames(prot.mat) = ids.prot
prot.int <- prot.mat[,match.prot]
rownames(exp.int) <- genes
rownames(meth.int) <- cpgs
rownames(mirna.int) <- mirnas
exp.int <- as.matrix(exp.int)
meth.int <- as.matrix(meth.int)
mirna.int <- as.matrix(mirna.int)

#exp.int, prot.int, meth.int and mirna.int give final files

tum.type <- as.factor(as.character(tum.type)) 
#log-normalize expression and miRNA
exp.int <- log(1+exp.int)
mirna.int <- log(1+mirna.int)
n.source <- 4
n.type <- 30 #will reduce to 29 below
#filter genes and miRNAs
gene.sd <- rep(0,dim(exp.int)[1])
for(i in 1:dim(exp.int)[1]){
  gene.sd[i] <- sd(lm(exp.int[i,]~tum.type)$resid)
}
thresh <- sort(gene.sd,decreasing=TRUE)[1001]
genestokeep <- c(1:dim(exp.int)[1])[gene.sd>thresh]
genes.filt <- genes[genestokeep]

meth.sd <- rep(0,dim(meth.int)[1])
for(i in 1:dim(meth.int)[1]){
  meth.sd[i] <- sd(lm(meth.int[i,]~tum.type)$resid)
}
thresh <- sort(meth.sd,decreasing=TRUE)[1001]
methtokeep <- c(1:dim(meth.int)[1])[meth.sd>thresh]
cpgs.filt <- cpgs[methtokeep]

exp.int <- exp.int[genestokeep,]
meth.int <- meth.int[methtokeep,]

#center and scale each data block
X <- matrix(list(),n.source,n.type)

for(j in 1:n.type){
  X[[1,j]] <- exp.int[,as.numeric(tum.type)==j]
  X[[2,j]] <- mirna.int[,as.numeric(tum.type)==j]
  X[[3,j]] <- meth.int[,as.numeric(tum.type)==j]
  X[[4,j]] <- prot.int[,as.numeric(tum.type)==j]
}

for(i in 1:n.source){ for(j in 1:n.type){
  X[[i,j]] <- t(scale(t(X[[i,j]]), scale=FALSE))
  X[[i,j]][is.na(X[[i,j]])] <- 0
  X[[i,j]] <- X[[i,j]]/sd(X[[i,j]])
}}

#get indices
n.vec <- c()
p.vec <- c()
for(i in 1:n.source){p.vec[i] <- dim(X[[i,1]])[1]}
for(j in 1:n.type){n.vec[j] <- dim(X[[1,j]])[2]}

N <- sum(n.vec)
P <- sum(p.vec)

n.ind <- list()
p.ind <- list()
p.ind[[1]] <- c(1:1000)
p.ind[[2]] <- c(1001:(p.vec[1]+p.vec[2]))
p.ind[[3]] <- c((sum(p.vec[1:2])+1):sum(p.vec[1:3]))
p.ind[[4]] <- c((sum(p.vec[1:3])+1):sum(p.vec[1:4]))

for(j in 1:n.type){n.ind[[j]] <- c(1:N)[as.numeric(tum.type)==j]}

#build full datasert X0
X0 <- matrix(nrow=P,ncol=N)
for(i in 1:n.source){ for(j in 1:n.type){
  X0[p.ind[[i]], n.ind[[j]]] <- X[[i,j]]
}}
#remove last cancer type (UVM) due to lack of data (only 12 samples)
X0 <- X0[,-c(n.ind[[30]])]
n.ind[[30]] <- NULL
ids.int = ids.int[as.numeric(tum.type)!=30]
tum.type =tum.type[as.numeric(tum.type)!=30]
n.type=29
N=length(tum.type)
for(j in 1:n.type){n.ind[[j]] <- c(1:N)[as.numeric(tum.type)==j]}

#run bidifac+ for initial iterations (this can take substantial time, 1-2 days)
res=bidifac.plus(X0,p.ind,n.ind,max.comb=50,num.comp=25,max.iter=20,temp.iter=10, conv.thresh=1)
#run bidifac+ for more iterations using given row and column sets (again, this can take several hours)
res2=bidifac.plus.given(X0,p.ind,n.ind,p.ind.list=res$p.ind.list,n.ind.list=res$n.ind.list,S=res$S,pen=res$pen,given.inits=TRUE,max.iter=80,conv.thresh=1)
#check convergence
plot(c(res$obj.vec,res2$obj.vec))
S=res2$S
Sums=res2$Sums
n.ind.list = res$n.ind.list
p.ind.list = res$p.ind.list

#Collect results and list modules
wtype <- list()
for(i in 1:50){ wtype[[i]] <- levels(tum.type)[1:29][colSums(Sums[i,,])>0]}

wsource <- list()
for(i in 1:50){ wsource[[i]] <- c('mRNA','miRNA','Methylation','Protein')[rowSums(Sums[i,,])>0]}

#order by variance explained
wvar <-c()
for(i in 1:50){ wvar[i] = sum(Sums[i,,])}

Order <- order(wvar,decreasing=TRUE)
S.list = S[Order]
p.ind.list.order = p.ind.list[Order]
n.ind.list.order = n.ind.list[Order]
#S.list=S
tum.type.list <- list()
sample.id.list <- list()
for(i in 1:50){
  S.list[[i]] <- S.list[[i]][p.ind.list.order[[i]],n.ind.list.order[[i]]]
  tum.type.list[[i]] <- as.character(tum.type)[n.ind.list.order[[i]]]
  sample.id.list[[i]] <- ids.int[n.ind.list.order[[i]]]
}

type.list = wtype[Order]
source.list = wsource[Order]
save(source.list,type.list,S.list,tum.type.list,sample.id.list,file='pan.fac.results_v2.rda')


