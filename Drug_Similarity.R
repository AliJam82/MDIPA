library(ChemmineR)
library(ChemmineOB)
sdfset <- read.SDFstr("DATA.sdf") 
valid <- validSDF(sdfset); sdfset <- sdfset[valid]
apset <- sdf2ap(sdfset)
fpset <- desc2fp(apset, descnames=1024, type="FPset")
Tanimoto <- matrix(0:0, nrow = 2287, ncol = 2287)
for (i in 1:2287) {
  for (j in i:2287) {
    Tanimoto[i,j] <- fpSim(fpset[i], fpset[j], method="Tanimoto")
  }
}
write.csv(Tanimoto, "Tanimoto2287.csv")


#To copy Transpose of matrix below the diagonal
for (i in 1:2287){
  for (j in 1:i) {
    Tanimoto[i, j] <- Tanimoto[j,i]
  }
}
write.csv(Tanimoto, "Tanimoto2287-2.csv")

#To creat the labels
for (i in 1:2287){
  label[i] <- sdfset@SDF[[i]]@datablock[["SYNONYMS"]]
  }
write.csv(label, "Label.csv")
