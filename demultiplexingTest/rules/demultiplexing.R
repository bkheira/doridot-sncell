# This script handles the demultiplexing of one snRNAseq run

#Date: 29/03/2022
### version modified by LD
### cleanded, updated and adapted to snakemake by Alix Silvert

##TODO remove Ludivine's note once I'm done updating the script

#Group: Gp2B (encapsulation du 8/12/2021)
#4 HashTag (HT1--> souris 2-2 placentas G3G5 (Female), HT2--> souris 2-3 placentas G1G6 (Male),
#HT3--> souris 5-2 placenta G4 (Female), HT5--> souris 5-5 placentas G2G3 (Male)

###	modifié pour prendre en compte le background des HT
### partir de la matrice raw qui contient les goutelettes vide (ie très peu de gènes détectés)


library(Matrix)
library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

#input parameters
path_dir = args TODO
sample_group= args TODO
sample_dir= args TODO
samples = args TODO
output_dir= args TODO

filterthreshold= args TODO


matrix_dir = paste0(path_dir, sample_dir, "outs/filtered_feature_bc_matrix/")
#prefix=paste0(format(Sys.Date(), format="%y%m%d"),"_",sample_group,"_")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
feature.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

#Loading data
mat <- readMM(file = matrix.path)
feature.names = read.delim(feature.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
dim(mat)
#[1] 32290  7273
tail(mat)
head(mat)

### load unfiltered matrix to extract empty droplet in order to assess hashtag background
unfilteredmatrix_dir=paste0(path_dir, sample_dir, "outs/raw_feature_bc_matrix/")
unfilteredbarcode.path <- paste0(unfilteredmatrix_dir, "barcodes.tsv.gz")
unfilteredfeature.path <- paste0(unfilteredmatrix_dir, "features.tsv.gz")
unfilteredmatrix.path <- paste0(unfilteredmatrix_dir, "matrix.mtx.gz")
unfilteredmat=readMM(file = unfilteredmatrix.path)
unfilteredfeature.names = read.delim(unfilteredfeature.path,
                                     header = FALSE,
                                     stringsAsFactors = FALSE)
unfilteredbarcode.names = read.delim(unfilteredbarcode.path,
                                     header = FALSE,
                                     stringsAsFactors = FALSE)
colnames(unfilteredmat) = unfilteredbarcode.names$V1
rownames(unfilteredmat) = unfilteredfeature.names$V1
dim(unfilteredmat)
#[1]

## subset of the unfiltered matrix to only keep empty droplet
# Keep cells with at least 1 HT detected
##empty droplet if <10 genes and <100 RNA (>1RNA)

SumAboveOne <- function(x){
  if(sum(x)>=1){
    result = TRUE
  } else {
    result = FALSE
  }
  return(result)
}

unfilteredmatAC=unfilteredmat[(nrow(unfilteredmat)-3):nrow(unfilteredmat),]
unfilteredmatRNA=unfilteredmat[1:(nrow(unfilteredmat)-4),]
TokeepAC=apply(unfilteredmatAC,2,function(x) SumAboveOne(x))
sum(TokeepAC)
matAC1=unfilteredmatAC[,TokeepAC]
matRNA_AC1=unfilteredmatRNA[,TokeepAC]
str(matRNA_AC1)
SumRNA=colSums(matRNA_AC1)
Ngenes=colSums(matRNA_AC1 != 0)
df=data.frame(SumRNA,Ngenes)
#ggplot(df, aes(x=SumRNA, y=Ngenes)) + geom_point()
df_empty=df[(df$SumRNA<=filterthreshold[1])&(df$Ngenes<=filterthreshold[2])&(df$SumRNA>0),]
dim(df_empty)
matAC1_emptydroplet=matAC1[,rownames(df_empty)]
Tokeep=apply(matAC1_emptydroplet,1,function(x) max(x)<1000)
#### outlier with 1079 hashtag 2-2 detected...
matAC1_emptydroplet2=matAC1_emptydroplet[,Tokeep]
backgroundHT=apply(matAC1_emptydroplet2,1,max)


#Remove background signal from matAC
# background defined as the max in all empty droplets
# empty droplets= from unfiltered matrix with more than 1 UMI for at least one ADT but less than 100 UMI for RNA/10genes

whichMaxPerso <- function(x){
  result = which.max(x)
  if(sum(x == x[result])>1){
    result = 0
  }
  return(result)
}

matAC=mat[(nrow(mat)-3):nrow(mat),] #matrix subset with only the antibody values
matAC=apply(matAC, 2, function(x) x-backgroundHT)
signalHT=apply(matAC,1,mean)

M <-apply(matAC,2,function(x) whichMaxPerso(x))
matAC[]

#Rajouter la ligne qui donne les infos sur l'AC max
mat <- rbind(mat,M)
mat


#Extraire les matrices avec les cellules pour chaque AC
i=1;
for(s in samples) {
  submat=mat[1:(nrow(mat)-5),mat[nrow(mat),]==i]
  assign(s, submat)
  writeMM(submat,file=paste0(output_dir, 'mat', s, '.mtx'))
  i=i+1
}

#Creéer un dossier par echantillon dans le groupe avec les fichiers matrix, fetaures and barcodes correspondant (les fichiers par sample),
#et gzipper ce qui doit l'être
#git / gitlab ou github ???

features_dt = fread(feature.path)
features_without_ab_dt = features_dt[1:(.N-4)]

for (s in samples){
  new_folder = paste(output_dir, 'mat', s, "/", sep="")
  dir.create(new_folder)
  file.copy(paste(output_dir,'mat', s,".mtx",sep=""),paste(new_folder, "matrix.mtx", sep=""))
  system(paste("gzip ", new_folder , "matrix.mtx", sep=""))
  barcodes_temp = data.frame(V1 = colnames(get(s)))
  write.table(barcodes_temp,
              file=paste0(new_folder, "barcodes.tsv"),
              quote=FALSE,
              row.names=FALSE,
              col.names=FALSE)
  system(paste("gzip ", new_folder , "barcodes.tsv", sep=""))
  write.table(features_without_ab_dt, paste0(new_folder, "features.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
  system(paste("gzip ", new_folder , "features.tsv", sep=""))
}
