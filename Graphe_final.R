
library(igraph)
library(Hmisc)

# Colors

couleurs1 <- adjustcolor( c("skyblue", "red"))
couleurs2 <- adjustcolor( c("grey", "yellowgreen", "lightgreen", "pink","lightblue"), alpha=1)

# Data reading

data_raw <- read.csv2("D:/REORG-Neuropsy/Data_Neuropsy_full.csv", header=T, as.is=T,sep=';')
data_raw[1:28,] <- transform(data_raw[1:28,], PAT = sprintf('L%s', PAT))
data_raw[29:51,] <- transform(data_raw[29:51,], PAT = sprintf('R%s', PAT))

data <- na.omit(data_raw[,c(1,3,4,6,8,12,16,19,21,23,25,26,60,65,2,27:58)])

data[,15:47] <- lapply(data[,15:47] , as.numeric)
data[,1:14] <- lapply(data[,1:14] , as.factor)

data[,15:47] <- scale(data[,15:47])


scores <- names(data)

# Binary variables

gen <- as.numeric(as.factor(data$GEN))-1
hds <- as.numeric(as.factor(data$HDS))-1
hem <- as.numeric(as.factor(data$HEM))-1
hs <- as.numeric(as.factor(data$HS))-1
dur <- as.numeric(as.factor(data$DUR_bin))-1
sev <- as.numeric(as.factor(data$SEV_bin))-1
bdi <- 2-as.numeric(as.factor(data$BDI_bin))
sta <- 2-as.numeric(as.factor(data$STA_bin))
stb <- 2-as.numeric(as.factor(data$STB_bin))
thy <- 2-as.numeric(as.factor(data$THY_split))
eng <- 2-as.numeric(as.factor(data$ENG_bin))
nam <- 2-as.numeric(as.factor(data$DIFF_NAM_1))



################################################################################
## Fonctions
# Simple Matching Coefficient

binary_table = function(x, y) {
  result = matrix(0, 2, 2)
  result[1, 1] = sum((x == 0) * (y == 0))
  result[1, 2] = sum((x == 0) * (y == 1))
  result[2, 1] = sum((x == 1) * (y == 0))
  result[2, 2] = sum((x == 1) * (y == 1))
  colnames(result) = c("y0", "y1")
  rownames(result) = c("x0", "x1")
  return(result)
}

SMC = function(x, y) {
  bt = binary_table(x, y)
  return((bt[1, 1] + bt[2, 2])/sum(bt))
}


################################################################################
# PRI, VMI, VCI, AMI, NAM, SFL, PFL, TMT, STR


distance_1 <- as.matrix(dist(data[c("PRI","VMI","VCI","AMI","NAM","SFL","PFL","TMT","STR")], diag = TRUE))
proximite_1 <- (max(distance_1) - distance_1)/(max(distance_1)-min(distance_1))
proximite_1 <- distance_1
colnames(proximite_1)<-data$PAT
rownames(proximite_1)<-data$PAT

graphe <- graph_from_adjacency_matrix(proximite_1, weighted=T, mode="undirected", diag=F)
layouts <- layout_with_fr(graphe)
titre = "Sous-jeu"


# Graphe original
plot(graphe, edge.arrow.mode=0,vertex.color=couleurs1[as.factor(data$HEM)],vertex.label.color="black",vertex.label.cex=.6, layout=layouts)
title(main = paste("Graphe de patients :",titre), font.main= 1)
legend(x=-1.5, y=0, c("Left", "Right"), pch=21, col="#777777", pt.bg=couleurs1, pt.cex=2, cex=.8, bty="n", ncol=1)

# Weighted Community detection based on Multilevel
louvain <- cluster_louvain(graphe)
com <- 1-(as.numeric(membership(louvain))-1)
SMCdur <- 1-round(SMC(com,dur),digits=2)
SMChs <- 1-round(SMC(com,1-hs),digits=2)
SMCsev <- 1-round(SMC(com,1-sev),digits=2)
SMCgen <- 1-round(SMC(com,gen),digits=2)
SMChem <- 1-round(SMC(com,hem),digits=2)
SMChds <- 1-round(SMC(com,hds),digits=2)

SMCbdi <- 1-round(SMC(com,bdi),digits=2)
SMCsta <- 1-round(SMC(com,sta),digits=2)
SMCstb <- 1-round(SMC(com,stb),digits=2)
SMCthy <- round(SMC(com,thy),digits=2)
SMCeng <- round(SMC(com,1-eng),digits=2)
SMCnam <- round(SMC(com,1-nam),digits=2)

plot(louvain, graphe ,edge.arrow.mode=0, vertex.label.cex=.6, layout=layouts)
title(main = paste("Multilevel/Louvain :",titre), 
      sub = paste("SMC Durée :",SMCdur," SMC Sclérose :",SMChs," SMC Sévérité :",SMCsev,
                  "\n SMC Genre :",SMCgen," SMC Hémisphère :",SMChem," SMC Latéralité :",SMChds,
                  "\n SMC BDI :",SMCbdi," SMC STA :",SMCsta," SMC STB :",SMCstb,
                  "\n SMC Thy :",SMCthy," SMC Engel :",SMCeng," SMC Diff_NAM :",SMCnam
      ))


## Enregistrement du graphe
# attributs
graphe <- set_vertex_attr(graph = graphe, name = "DUR_bin", index=V(graphe), value = data$DUR_bin)
graphe <- set_vertex_attr(graph = graphe, name = "HS", index=V(graphe), value = data$HS)
graphe <- set_vertex_attr(graph = graphe, name = "SEV_bin", index=V(graphe), value = data$SEV_bin)
graphe <- set_vertex_attr(graph = graphe, name = "GEN", index=V(graphe), value = data$GEN)
graphe <- set_vertex_attr(graph = graphe, name = "HEM", index=V(graphe), value = data$HEM)
graphe <- set_vertex_attr(graph = graphe, name = "HDS", index=V(graphe), value = data$HDS)
graphe <- set_vertex_attr(graph = graphe, name = "BDI_bin", index=V(graphe), value = data$BDI_bin)
graphe <- set_vertex_attr(graph = graphe, name = "STA_bin", index=V(graphe), value = data$STA_bin)
graphe <- set_vertex_attr(graph = graphe, name = "STB_bin", index=V(graphe), value = data$STB_bin)
graphe <- set_vertex_attr(graph = graphe, name = "THY_bin", index=V(graphe), value = data$THY_bin)
graphe <- set_vertex_attr(graph = graphe, name = "ENG_bin", index=V(graphe), value = data$ENG_bin)
graphe <- set_vertex_attr(graph = graphe, name = "DIFF_NAM_1", index=V(graphe), value = data$DIFF_NAM_1)
write_graph(graphe, "./graphe_Patients_ENG.graphml" , format = "graphml")
##



par(mfrow=c(1,2))

plot(louvain, graphe ,edge.arrow.mode=0, vertex.label.cex=.6, layout=layouts)
title(main = paste("Multilevel/Louvain :",titre), 
      sub = paste("SMC Durée :",SMCdur," SMC Sclérose :",SMChs," SMC Sévérité :",SMCsev,
                  "\n SMC Genre :",SMCgen," SMC Hémisphère :",SMChem," SMC Latéralité :",SMChds,
                  "\n SMC BDI :",SMCbdi," SMC STA :",SMCsta," SMC STB :",SMCstb,
                  "\n SMC Thy :",SMCthy," SMC Engel :",SMCeng," SMC Diff_NAM :",SMCnam
      ))

# Graphe spectral embedding
c2 = cluster_leading_eigen(graphe)
plot(c2, graphe ,edge.arrow.mode=0, vertex.label.cex=.6, layout=layouts)
com <- 1-(as.numeric(membership(c2))-1)
SMCdur <- round(SMC(com,dur),digits=2)
SMChs <- round(SMC(com,1-hs),digits=2)
SMCsev <- round(SMC(com,1-sev),digits=2)
SMCgen <- round(SMC(com,gen),digits=2)
SMChem <- round(SMC(com,hem),digits=2)
SMChds <- round(SMC(com,hds),digits=2)

SMCbdi <- round(SMC(com,bdi),digits=2)
SMCsta <- round(SMC(com,sta),digits=2)
SMCstb <- round(SMC(com,stb),digits=2)
SMCthy <- round(SMC(com,thy),digits=2)
SMCeng <- 1-round(SMC(com,1-eng),digits=2)
SMCnam <- round(SMC(com,1-nam),digits=2)

title(main = paste("Spectral clustering :",titre),
      sub = paste("SMC Durée :",SMCdur," SMC Sclérose :",SMChs," SMC Sévérité :",SMCsev,
                  "\n SMC Genre :",SMCgen," SMC Hémisphère :",SMChem," SMC Latéralité :",SMChds,
                  "\n SMC BDI :",SMCbdi," SMC STA :",SMCsta," SMC STB :",SMCstb,
                  "\n SMC Thy :",SMCthy," SMC Engel :",SMCeng," SMC Diff_NAM :",SMCnam
      ))

par(mfrow=c(1,1))
plot_dendrogram(c2)
title(main = paste("Spectral clustering :",titre))


# Spinglass

c3 = cluster_spinglass(graphe)
plot(c3, graphe ,edge.arrow.mode=0, vertex.label.cex=.6, layout=layouts)
com <- 1-(as.numeric(membership(c3))-1)
SMCdur <- 1-round(SMC(com,dur),digits=2)
SMChs <- round(SMC(com,1-hs),digits=2)
SMCsev <- round(SMC(com,1-sev),digits=2)
SMCgen <- round(SMC(com,gen),digits=2)
SMChem <- round(SMC(com,hem),digits=2)
SMChds <- round(SMC(com,hds),digits=2)

SMCbdi <- round(SMC(com,bdi),digits=2)
SMCsta <- round(SMC(com,sta),digits=2)
SMCstb <- round(SMC(com,stb),digits=2)
SMCthy <- round(SMC(com,thy),digits=2)
SMCeng <- 1-round(SMC(com,1-eng),digits=2)
SMCnam <- 1-round(SMC(com,1-nam),digits=2)

title(main = paste("Spinglass :",titre),
      sub = paste("SMC Durée :",SMCdur," SMC Sclérose :",SMChs," SMC Sévérité :",SMCsev,
                  "\n SMC Genre :",SMCgen," SMC Hémisphère :",SMChem," SMC Latéralité :",SMChds,
                  "\n SMC BDI :",SMCbdi," SMC STA :",SMCsta," SMC STB :",SMCstb,
                  "\n SMC Thy :",SMCthy," SMC Engel :",SMCeng," SMC Diff_NAM :",SMCnam
      ))

################################################################################
# Sous-jeu VMI,NAM,TMT,AMI

eng_1 = subset(data,data$ENG_bin==1)[c("VMI","AMI","NAM","TMT")]
eng_2 = subset(data,data$ENG_bin==2)[c("VMI","AMI","NAM","TMT")]
#

# distance_1 <- as.matrix(dist(data[c("VMI","AMI","NAM","TMT")], diag = TRUE))
distance_1 <- as.matrix(dist(eng_2, diag = TRUE))
proximite_1 <- (max(distance_1) - distance_1)/(max(distance_1)-min(distance_1))
colnames(proximite_1)<-subset(data,data$ENG_bin==2)$PAT
rownames(proximite_1)<-subset(data,data$ENG_bin==2)$PAT
colnames(proximite_1)<-data$PAT
rownames(proximite_1)<-data$PAT

graphe <- graph_from_adjacency_matrix(proximite_1, weighted=T, mode="undirected", diag=F)


# attributs
graphe <- set_vertex_attr(graph = graphe, name = "DUR_bin", index=V(graphe), value = data$DUR_bin)
graphe <- set_vertex_attr(graph = graphe, name = "HS", index=V(graphe), value = data$HS)
graphe <- set_vertex_attr(graph = graphe, name = "SEV_bin", index=V(graphe), value = data$SEV_bin)
graphe <- set_vertex_attr(graph = graphe, name = "GEN", index=V(graphe), value = data$GEN)
graphe <- set_vertex_attr(graph = graphe, name = "HEM", index=V(graphe), value = data$HEM)
graphe <- set_vertex_attr(graph = graphe, name = "HDS", index=V(graphe), value = data$HDS)
graphe <- set_vertex_attr(graph = graphe, name = "BDI_bin", index=V(graphe), value = data$BDI_bin)
graphe <- set_vertex_attr(graph = graphe, name = "STA_bin", index=V(graphe), value = data$STA_bin)
graphe <- set_vertex_attr(graph = graphe, name = "STB_bin", index=V(graphe), value = data$STB_bin)
graphe <- set_vertex_attr(graph = graphe, name = "THY_bin", index=V(graphe), value = data$THY_bin)
graphe <- set_vertex_attr(graph = graphe, name = "ENG_bin", index=V(graphe), value = data$ENG_bin)
graphe <- set_vertex_attr(graph = graphe, name = "DIFF_NAM_1", index=V(graphe), value = data$DIFF_NAM_1)
write_graph(graphe, "./graphe_Patients_NAM2.graphml" , format = "graphml")
##


layouts <- layout_with_fr(graphe)
titre = "Sous-jeu"

par(mfrow=c(1,2))

# Graphe original
plot(graphe, edge.arrow.mode=0,vertex.color=couleurs1[as.factor(subset(data,data$ENG_bin==2)$HEM)],vertex.label.color="black",vertex.label.cex=.6, layout=layouts)
title(main = paste("Graphe de patients :",titre), font.main= 1)
legend(x=-1.5, y=0, c("Left", "Right"), pch=21, col="#777777", pt.bg=couleurs1, pt.cex=2, cex=.8, bty="n", ncol=1)


gen <- as.numeric(as.factor(subset(data,data$ENG_bin==2)$GEN))-1
hds <- as.numeric(as.factor(subset(data,data$ENG_bin==2)$HDS))-1
hem <- as.numeric(as.factor(subset(data,data$ENG_bin==2)$HEM))-1
hs <- as.numeric(as.factor(subset(data,data$ENG_bin==2)$HS))-1
dur <- as.numeric(as.factor(subset(data,data$ENG_bin==2)$DUR_bin))-1
sev <- as.numeric(as.factor(subset(data,data$ENG_bin==2)$SEV_bin))-1
bdi <- 2-as.numeric(as.factor(subset(data,data$ENG_bin==2)$BDI_bin))
sta <- 2-as.numeric(as.factor(subset(data,data$ENG_bin==2)$STA_bin))
stb <- 2-as.numeric(as.factor(subset(data,data$ENG_bin==2)$STB_bin))
thy <- 2-as.numeric(as.factor(subset(data,data$ENG_bin==2)$THY_split))
eng <- 2-as.numeric(as.factor(subset(data,data$ENG_bin==2)$ENG_bin))
nam <- 2-as.numeric(as.factor(subset(data,data$ENG_bin==2)$DIFF_NAM_1))

# Weighted Community detection based on Multilevel
louvain <- cluster_louvain(graphe)
com <- 1-(as.numeric(membership(louvain))-1)
SMCdur <- round(SMC(com,dur),digits=2)
SMChs <- round(SMC(com,1-hs),digits=2)
SMCsev <- 1-round(SMC(com,1-sev),digits=2)
SMCgen <- 1-round(SMC(com,gen),digits=2)
SMChem <- 1-round(SMC(com,hem),digits=2)
SMChds <- round(SMC(com,hds),digits=2)

SMCbdi <- 1-round(SMC(com,bdi),digits=2)
SMCsta <- round(SMC(com,sta),digits=2)
SMCstb <- round(SMC(com,stb),digits=2)
SMCthy <- round(SMC(com,thy),digits=2)
SMCeng <- round(SMC(com,1-eng),digits=2)
SMCnam <- round(SMC(com,1-nam),digits=2)

plot(louvain, graphe ,edge.arrow.mode=0, vertex.label.cex=.6, layout=layouts)
title(main = paste("Multilevel/Louvain :",titre), 
      sub = paste("SMC Durée :",SMCdur," SMC Sclérose :",SMChs," SMC Sévérité :",SMCsev,
                  "\n SMC Genre :",SMCgen," SMC Hémisphère :",SMChem," SMC Latéralité :",SMChds,
                  "\n SMC BDI :",SMCbdi," SMC STA :",SMCsta," SMC STB :",SMCstb,
                  "\n SMC Thy :",SMCthy," SMC Engel :",SMCeng," SMC Diff_NAM :",SMCnam
      ))



################################################################################
# Sous-jeu : VCI,NAM,SFL,AMI

distance_1 <- as.matrix(dist(data[c("VCI","AMI","NAM","SFL")], diag = TRUE))
proximite_1 <- (max(distance_1) - distance_1)/(max(distance_1)-min(distance_1))
colnames(proximite_1)<-data$PAT
rownames(proximite_1)<-data$PAT


graphe <- graph_from_adjacency_matrix(proximite_1, weighted=T, mode="undirected", diag=F)
layouts <- layout_with_fr(graphe)
titre = "Sous-jeu 5.1 réduit dif_nam"

# png(paste("Graphe_",scores,".png"), width=10, height=6, units="in", res=150)
par(mfrow=c(1,2))

# Graphe original
plot(graphe, edge.arrow.mode=0,vertex.color=couleurs1[as.factor(data$HEM)],vertex.label.color="black",vertex.label.cex=.6, layout=layouts)
title(main = paste("Graphe de patients :",titre), font.main= 1)
legend(x=-1.5, y=0, c("Left", "Right"), pch=21, col="#777777", pt.bg=couleurs1, pt.cex=2, cex=.8, bty="n", ncol=1)

# Weighted Community detection based on Multilevel
louvain <- cluster_louvain(graphe)
com <- 1-(as.numeric(membership(louvain))-1)
SMCdur <- round(SMC(com,dur),digits=2)
SMChs <- round(SMC(com,1-hs),digits=2)
SMCsev <- 1-round(SMC(com,1-sev),digits=2)
SMCgen <- 1-round(SMC(com,gen),digits=2)
SMChem <- 1-round(SMC(com,hem),digits=2)
SMChds <- 1-round(SMC(com,hds),digits=2)

SMCbdi <- 1-round(SMC(com,bdi),digits=2)
SMCsta <- round(SMC(com,sta),digits=2)
SMCstb <- 1-round(SMC(com,stb),digits=2)
SMCthy <- round(SMC(com,thy),digits=2)
SMCeng <- round(SMC(com,1-eng),digits=2)
SMCnam <- round(SMC(com,1-nam),digits=2)

plot(louvain, graphe ,edge.arrow.mode=0, vertex.label.cex=.6, layout=layouts)
title(main = paste("Multilevel/Louvain :",titre), 
      sub = paste("SMC Durée :",SMCdur," SMC Sclérose :",SMChs," SMC Sévérité :",SMCsev,
                  "\n SMC Genre :",SMCgen," SMC Hémisphère :",SMChem," SMC Latéralité :",SMChds,
                  "\n SMC BDI :",SMCbdi," SMC STA :",SMCsta," SMC STB :",SMCstb,
                  "\n SMC Thy :",SMCthy," SMC Engel :",SMCeng," SMC Diff_NAM :",SMCnam
      ))





################################################################################
# graphe de similarité basé sur la corrélation des tests neuropsy
################################################################################

########### Labels des groupes thématiques de tests ###########

# Verbal / Non-verbal
Groupe1=c(1,2,1,2,1,1,1,2,1,1,1,1,2,2,2,1,1,1,1,2,2,2,2,2,2,1,1,1,1,2,1,1)
Couleurs1 <- adjustcolor( c( "chartreuse2", "deepskyblue1"))
Legend1 = c('Verbal', 'No Verbal')
Label1 = Legend1[Groupe1]
Color1 = Couleurs1[Groupe1]

# Memoire, Langage, etc.
Groupe2=c(1,1,2,2,3,3,3,5,5,3,3,3,4,4,4,2,2,2,2,2,2,2,2,5,5,2,5,3,3,4,3,3)
Couleurs2 <- c("#F7F7F7","#D01C8B","#F1B6DA","#B8E186","#4DAC26")
Legend2 = c('IQ','Memory','Langage','Visuo-spatial','Executive')
Label2 = Legend2[Groupe2]
Color2 = Couleurs2[Groupe2]

# Memoire + Langage
Groupe2=c(1,1,2,2,3,3,3,5,5,3,3,3,4,4,4,2,2,2,2,2,2,2,2,5,5,2,5,3,3,4,3,3)
Couleurs3 <- c("#F1B6DA","#B8E186")
Legend3 = c(1,2,2,1,1)
Label3 = Legend3[Groupe2]
Color3 = Couleurs3[Groupe2]

# Memoire
Groupe2=c(1,1,2,2,3,3,3,5,5,3,3,3,4,4,4,2,2,2,2,2,2,2,2,5,5,2,5,3,3,4,3,3)
Legend4 = c(1,2,1,1,1)
Label4 = Legend4[Groupe2]

# Langage
Groupe2=c(1,1,2,2,3,3,3,5,5,3,3,3,4,4,4,2,2,2,2,2,2,2,2,5,5,2,5,3,3,4,3,3)
Legend5 = c(1,1,2,1,1)
Label5 = Legend5[Groupe2]

# Mémoire, Langage, les autres réunis
Groupe2=c(1,1,2,2,3,3,3,5,5,3,3,3,4,4,4,2,2,2,2,2,2,2,2,5,5,2,5,3,3,4,3,3)
Legend6 = c(1,2,3,1,1)
Label6 = Legend6[Groupe2]

################################################################################

eng_1 = subset(data[,16:47],data$ENG_bin==1) # 30 patients (rmin = 0,49)
eng_2 = subset(data[,16:47],data$ENG_bin==2) # 17 patients (rmin = 0,64)
nam_1 = subset(data[,16:47],data$DIFF_NAM_1==1) # 28 patients (rmin = 0,51)
nam_2 = subset(data[,16:47],data$DIFF_NAM_1==2) # 19 patients (rmin = 0,6)

df <- nam_2 # to change

####

library(psych)
corr_data_adjusted <- corr.test(as.matrix(df),method="pearson",adjust='holm')
corr_matrix_adjusted <- corr_data_adjusted$r
p_vals_adjusted <- corr_data_adjusted$p
corr_matrix_adjusted[p_vals_adjusted >= 0.005] <- 0
corr_matrix_adjusted[corr_matrix_adjusted < 0] <- 0
###

proximite <- corr_matrix_adjusted

colnames(proximite)<-names(data[,16:47])
rownames(proximite)<-names(data[,16:47])

graphe <- graph_from_adjacency_matrix(proximite, weighted=T, mode="undirected", diag=F)


################################################################################
# construction des graphes

graphe <- set_vertex_attr(graph = graphe, name = "z-score", index=V(graphe), value = colMeans(df))

graphe <- set_vertex_attr(graph = graphe, name = "Label1", index=V(graphe), value = Label1)
graphe <- set_vertex_attr(graph = graphe, name = "Label2", index=V(graphe), value = Label2)
graphe <- set_vertex_attr(graph = graphe, name = "Label3", index=V(graphe), value = Label3)
graphe <- set_vertex_attr(graph = graphe, name = "Label4", index=V(graphe), value = Label4)
graphe <- set_vertex_attr(graph = graphe, name = "Label5", index=V(graphe), value = Label5)
graphe <- set_vertex_attr(graph = graphe, name = "Label6", index=V(graphe), value = Label6)
graphe <- set_vertex_attr(graph = graphe, name = "Color1", index=V(graphe), value = Color1)
graphe <- set_vertex_attr(graph = graphe, name = "Color2", index=V(graphe), value = Color2)
graphe <- set_vertex_attr(graph = graphe, name = "Color3", index=V(graphe), value = Color3)
write_graph(graphe, "./graphe_NAM2_Holm_corrected.graphml" , format = "graphml")
Isolated = which(degree(graphe)==0)
graphe_compact = delete.vertices(graphe, Isolated)


################################################################################
# affichage des graphes

library(eigenmodel)

l <- layout_with_fr(graphe)

name = "Réseau Neuropsy ENG_1"

for (i in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
  
  net <- delete.edges(graphe, which(E(graphe)$weight < i))

  par(mfrow=c(1,2))
  
  sub_text = paste("Pearson, p-value:0.005 Holm corrected")
  # Graphe Color1
  plot(net, edge.arrow.mode=0,vertex.label=V(net)$name,vertex.color=V(net)$Color1,
       vertex.label.color="black",vertex.label.cex=.8, layout=l)
  title(main = paste(name), font.main= 1, sub = sub_text)
  legend(x=-1.5, y=0, legend = Legend1, pch=21,
         col="#777777", pt.bg=Couleurs1, pt.cex=2, cex=.8, bty="n", ncol=1)
  
  # Graphe Color2
  plot(net, edge.arrow.mode=0,vertex.label=V(net)$name,vertex.color=V(net)$Color2,
       vertex.label.color="black",vertex.label.cex=.8, layout=l)
  title(main = paste(name), font.main= 1, sub = sub_text)
  legend(x=-1.5, y=0, legend = Legend2, pch=21,
         col="#777777", pt.bg=Couleurs2, pt.cex=2, cex=.8, bty="n", ncol=1)
  
  # Détection de communautés
  
  # Weighted Community detection based on Multilevel
  louv_G <- cluster_louvain(net)
  c_m <- membership(louv_G)
  plot(louv_G, net,edge.arrow.mode=0, layout=l)
  title(main = paste("Multilevel modularité : ",round(modularity(louv_G),2)))
  
  
  # Modèle latent
  A <- get.adjacency(net, sparse = FALSE)
  leig.fit1 <- eigenmodel_mcmc(A, R=2, S=11000,burn=10000)
  lat.sp.1 <- eigen(leig.fit1$ULU_postmean)$vec[, 1:2]
  coords <- format(round(apply(leig.fit1$L_postsamp, 2, mean), 2), nsmall = 2) 
  
  plot(net, edge.arrow.mode=0,vertex.label=V(net)$Label2,vertex.color=V(net)$Color1,
       vertex.label.color="black",vertex.label.cex=.8, layout=lat.sp.1)
  title(main = paste("Eigenmodel ", name), font.main= 1, sub = paste("composantes :",coords[1]," ",coords[2]))
  legend(x=-1.5, y=0, legend = Legend1, pch=21,
         col="#777777", pt.bg=Couleurs1, pt.cex=2, cex=.8, bty="n", ncol=1)
  
}

#### Mesures globales sur graphe neuropsy seuillé

hist_modularity <- array()
hist_density <- array()
hist_average_path_length <- array()
k <- 0

graphe <- graphe_eng2
for (i in c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8)){
  net <- delete.edges(graphe, which(E(graphe)$weight < i))
  k<-k+1
  
  # Weighted Community detection based on Multilevel
  louv_G <- cluster_louvain(net)
  hist_modularity[k] <- round(modularity(louv_G),2)
  hist_density[k] <- edge_density(net)
  hist_average_path_length[k] <- average.path.length(net)
  
}


plot(y=hist_modularity,x=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8), main='Graph Louvain Modularity', ylab='Modularity', xlab='Correlation Threshold', type='b')
plot(y=hist_density,x=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8), main='Graph Density', ylab='Density', xlab='Correlation Threshold', type='b')
plot(y=hist_average_path_length,x=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8), main='Graph Average Path Length', ylab='Average Path Length', xlab='Correlation Threshold', type='b')

# multi lines plot
# Create a first line
plot(y=modularity_nam1,x=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8),
      main='Louvain/Multilevel Modularity', ylab='Modularity', xlab='Correlation Threshold',
      type='b', frame = FALSE, pch = 19,
      col = couleur_wes[1])
# Add a second line
lines(y=modularity_nam2,x=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8),
      pch = 18, col = couleur_wes[2], type = "b", lty = 2)
# Add a third line
lines(y=modularity_eng1,x=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8),
      pch = 18, col = couleur_wes[3], type = "b", lty = 2)
# Add a fourth line
lines(y=modularity_eng2,x=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8),
      pch = 18, col = couleur_wes[4], type = "b", lty = 2)
# Add a legend to the plot
legend("left", legend=c("Naming 1", "Naming 2", "Engel 1", "Engel 2"),
       col=c(couleur_wes[1], couleur_wes[2],couleur_wes[3],couleur_wes[4]),
       lty = 1:2, cex=0.8, inset = 0.02)

################################################################################
# Boostrap sur les caractéristiques du graphe de tests
################################################################################

# Fonctions

library(psych)

# Closeness
library(CINNA)
node_harmonic <- function(x,d) {
  data = x[d,]
  # corr_data <- rcorr(as.matrix(data),type="pearson")
  corr_data <- corr.test(as.matrix(data),method="pearson",adjust='holm')
  corr_matrix <- corr_data$r
  p_vals <- corr_data$P
  corr_matrix[p_vals >= 0.005] <- 0
  corr_matrix[corr_matrix < 0] <- 0
  network=graph.adjacency(corr_matrix, weighted=T, mode="undirected", diag=F)
  result = harmonic_centrality(network)
  return(result)
}


# Weighted clustering Coefficient
trans_local <- function(x,d) {
  data = x[d,]
  # corr_data <- rcorr(as.matrix(data),type="pearson")
  corr_data <- corr.test(as.matrix(data),method="pearson",adjust='holm')
  corr_matrix <- corr_data$r
  p_vals <- corr_data$P
  corr_matrix[p_vals >= 0.005] <- 0
  corr_matrix[corr_matrix < 0] <- 0
  network=graph.adjacency(corr_matrix, weighted=T, mode="undirected", diag=F)
  result = transitivity(network, type="weighted",weights = E(network)$weight)
  return(result)
}


# Node Betweenness
node_betweenness <- function(x,d) {
  data = x[d,]
  # corr_data <- rcorr(as.matrix(data),type="pearson")
  corr_data <- corr.test(as.matrix(data),method="pearson",adjust='holm')
  corr_matrix <- corr_data$r
  p_vals <- corr_data$P
  corr_matrix[p_vals >= 0.005] <- 0
  corr_matrix[corr_matrix < 0] <- 0
  network=graph.adjacency(corr_matrix, weighted=T, mode="undirected", diag=F)
  result = betweenness(network, weights = E(network)$Weight)
  return(result)
}

# Strength
node_strength <- function(x,d) {
  data = x[d,]
  # corr_data <- rcorr(as.matrix(data),type="pearson")
  corr_data <- corr.test(as.matrix(data),method="pearson",adjust='holm')
  corr_matrix <- corr_data$r
  p_vals <- corr_data$P
  corr_matrix[p_vals >= 0.005] <- 0
  corr_matrix[corr_matrix < 0] <- 0
  network=graph.adjacency(corr_matrix, weighted=T, mode="undirected", diag=F)
  result = strength(network, weights = E(network)$Weight)
  return(result)
}

# Degree
node_degree <- function(x,d) {
  data = x[d,]
  # corr_data <- rcorr(as.matrix(data),type="pearson")
  corr_data <- corr.test(as.matrix(data),method="pearson",adjust='holm')
  corr_matrix <- corr_data$r
  p_vals <- corr_data$P
  corr_matrix[p_vals >= 0.005] <- 0
  corr_matrix[corr_matrix < 0] <- 0
  network=graph.adjacency(corr_matrix, weighted=T, mode="undirected", diag=F)
  result = degree(network)
  return(result)
}


library(wesanderson)
couleur3=wes_palette("GrandBudapest2", 4)

###############

library(boot)
library(reshape2)
library(ggpubr)
library(rstatix)

###

Clustering_coeff = boot(nam_1, trans_local, 1000)

# boxplot(Clustering_coeff$t[,1:9],col = Couleurs3[1], names = names(data[,16:24]), las=2, ylim = c(0, 1))
# title(main = "NAM_1 Clustering Coefficient pondéré")

Clustering_coeff_2 = boot(nam_2, trans_local, 1000)

# boxplot(Clustering_coeff_2$t[,1:9],col = Couleurs3[1], names = names(data[,16:24]), las=2, ylim = c(0, 1))
# title(main = "NAM_2 Clustering Coefficient pondéré")


# pvalue <- array()
# for (i in 1:32){  
#   t <- t.test(Clustering_coeff$t[,i],Clustering_coeff_2$t[,i])
#   pvalue[i] <-  format(round(t$p.value, 3), nsmall = 3) 
# }
# pvalue_Clustering_coeff <- rbind(names(data[,16:47]), pvalue)
# write.csv(pvalue_Clustering_coeff,"pvalue_Clustering_ENG.csv")


df1 <- Clustering_coeff$t[,1:9]
colnames(df1) <- names(data[,16:24])

df2 <- Clustering_coeff_2$t[,1:9]
colnames(df2) <- names(data[,16:24])

# soustraction pour mesure d'écart entre sélectionnés et non-sélectionnés en ML
delta=abs(df1-df2)
# write.csv(delta,"Clustering_NAM1-NAM2.csv")
delta_ML=subset(delta, select= c("VMI","AMI","NAM","TMT"))
# delta_ML=subset(delta, select= c("VCI","AMI","NAM","SFL"))
delta_noML=subset(delta, select= c("VCI","PRI","PFL","SFL","STR"))
# delta_noML=subset(delta, select= c("VMI","PRI","PFL","TMT","STR"))
delta_ML_mean = rowMeans(delta_ML)
delta_noML_mean = rowMeans(delta_noML)
t.test(delta_ML_mean,delta_noML_mean,alternative="less")
boxplot(cbind(delta_ML_mean,delta_noML_mean),ylab="Mean Clustering absolute difference",xlab="ENG - Nodes group",ylim=c(0,0.25))

# boxplot(delta,ylab="Mean Clustering absolute difference NAM1-NAM2",xlab="NAM - Nodes", col=couleur3[c(2,1,2,1,2,1,2,1,1)], ylim=c(0,0.35))
boxplot(delta,ylab="Mean Clustering absolute difference ENG1-ENG2",xlab="ENG - Nodes", col=couleur3[c(1,1,2,2,2,1,1,2,1)], ylim=c(0,0.4))
###

df1 <- melt(df1)
df2 <- melt(df2)

data_Clustering_coeff <- rbind(df1,df2)
data_Clustering_coeff <- cbind(data_Clustering_coeff,rep(1:2, each=9000))
colnames(data_Clustering_coeff) <- c("Index","Node","Clustering_coeff","NAM_groupe")
# data_Clustering_coeff$ENG_groupe <- as.factor(data_Clustering_coeff$ENG_groupe)
data_Clustering_coeff$NAM_groupe <- as.factor(data_Clustering_coeff$NAM_groupe)



stat.test <- data_Clustering_coeff %>%
  group_by(Node) %>%
  t_test(Clustering_coeff ~ NAM_groupe) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")
stat.test

# Créer un box plot
bxp <- ggboxplot(
  data_Clustering_coeff, x = "Node", y = "Clustering_coeff", 
  color = "NAM_groupe", palette = c("#00AFBB", "#E7B800")
)

# Ajoutez des p-values sur les graphiques en box plot
stat.test <- stat.test %>%
  add_xy_position(x = "Node", dodge = 0.8)

bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}{p.adj.signif}", 
  tip.length = 0.02, step.increase = 0.05, hide.ns = FALSE
)



###############

###############

node_betweenness_1 = boot(eng_1, node_betweenness, 1000)

# boxplot(node_betweenness_1$t[,1:9],col = Couleurs3[2], names = names(data[,16:24]), las=2, ylim = c(0, 100))
# title(main = "NAM_1 Node Betweenness")

node_betweenness_2 = boot(eng_2, node_betweenness, 1000)

# boxplot(node_betweenness_2$t[,1:9],col = Couleurs3[2], names = names(data[,16:24]), las=2, ylim = c(0, 100))
# title(main = "NAM_2 Node Betweenness")



# pvalue <- array()
# for (i in 1:32){  
#   t <- t.test(node_betweenness_1$t[,i],node_betweenness_2$t[,i])
#   pvalue[i] <-  format(round(t$p.value, 3), nsmall = 3) 
# }
# pvalue_node_betweenness <- rbind(names(data[,16:47]), pvalue)
# write.csv(pvalue_node_betweenness,"pvalue_node_betweenness_ENG.csv")


df1 <- node_betweenness_1$t[,1:9]
colnames(df1) <- names(data[,16:24])

df2 <- node_betweenness_2$t[,1:9]
colnames(df2) <- names(data[,16:24])

# soustraction pour mesure d'écart entre sélectionnés et non-sélectionnés en ML
delta=abs(df1-df2)
# write.csv(delta,"Node_Betweenness_NAM1-NAM2.csv")
delta_ML=subset(delta, select= c("VMI","AMI","NAM","TMT"))
# delta_ML=subset(delta, select= c("VCI","AMI","NAM","SFL"))
delta_noML=subset(delta, select= c("VCI","PRI","PFL","SFL","STR"))
# delta_noML=subset(delta, select= c("VMI","PRI","PFL","TMT","STR"))
delta_ML_mean = rowMeans(delta_ML)
delta_noML_mean = rowMeans(delta_noML)
t.test(delta_ML_mean,delta_noML_mean,alternative="greater")
boxplot(cbind(delta_ML_mean,delta_noML_mean),ylab="Mean Node Betweenness absolute difference",xlab="ENG Nodes group",ylim=c(0,110))

# boxplot(delta,ylab="Mean Node Betweenness absolute difference NAM1-NAM2",xlab="NAM - Nodes", col=couleur3[c(2,1,2,1,2,1,2,1,1)], ylim=c(0,180))
boxplot(delta,ylab="Mean Node Betweenness absolute difference ENG1-ENG2",xlab="ENG - Nodes", col=couleur3[c(1,1,2,2,2,1,1,2,1)], ylim=c(0,230))
###

df1 <- melt(df1)
df2 <- melt(df2)

data_node_betweenness <- rbind(df1,df2)
data_node_betweenness <- cbind(data_node_betweenness,rep(1:2, each=9000))
colnames(data_node_betweenness) <- c("Index","Node","node_betweenness","ENG_groupe")
data_node_betweenness$ENG_groupe <- as.factor(data_node_betweenness$ENG_groupe)

stat.test <- data_node_betweenness %>%
  group_by(Node) %>%
  t_test(node_betweenness ~ ENG_groupe) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")
stat.test

# Créer un box plot
bxp <- ggboxplot(
  data_node_betweenness, x = "Node", y = "node_betweenness", 
  color = "ENG_groupe", palette = c("#00AFBB", "#E7B800")
)

# Ajoutez des p-values sur les graphiques en box plot
stat.test <- stat.test %>%
  add_xy_position(x = "Node", dodge = 0.8)

bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}{p.adj.signif}", 
  tip.length = 0.02, step.increase = 0.05, hide.ns = FALSE
)


###############


###############

node_strength_1 = boot(nam_1, node_strength, 1000)


# boxplot(node_strength_1$t[,1:9],col = Couleurs3[3], names = names(data[,16:24]), las=2, ylim=c(0,15))
# title(main = "ENG_1 Node Strength")

node_strength_2 = boot(nam_2, node_strength, 1000)


# boxplot(node_strength_2$t[,1:9],col = Couleurs3[3], names = names(data[,16:24]), las=2, ylim=c(0,15))
# title(main = "ENG_2 Node Strength")


# pvalue <- array()
# for (i in 1:32){  
#   t <- t.test(node_strength_1$t[,i],node_strength_2$t[,i])
#   pvalue[i] <-  format(round(t$p.value, 3), nsmall = 3) 
# }
# pvalue_node_strength <- rbind(names(data[,16:47]), pvalue)
# write.csv(pvalue_node_strength,"pvalue_node_strength_ENG.csv")


df1 <- node_strength_1$t[,1:9]
colnames(df1) <- names(data[,16:24])

df2 <- node_strength_2$t[,1:9]
colnames(df2) <- names(data[,16:24])

# soustraction pour mesure d'écart entre sélectionnés et non-sélectionnés en ML
delta=abs(df1-df2)
write.csv(delta,"Strength_NAM1-NAM2.csv")
# delta_ML=subset(delta, select= c("VMI","AMI","NAM","TMT"))
delta_ML=subset(delta, select= c("VCI","AMI","NAM","SFL"))
# delta_noML=subset(delta, select= c("VCI","PRI","PFL","SFL","STR"))
delta_noML=subset(delta, select= c("VMI","PRI","PFL","TMT","STR"))
delta_ML_mean = rowMeans(delta_ML)
delta_noML_mean = rowMeans(delta_noML)
t.test(delta_ML_mean,delta_noML_mean,alternative="two.sided")
boxplot(cbind(delta_ML_mean,delta_noML_mean),ylab="Mean Strength absolute difference",xlab="NAM Nodes group",ylim=c(0,5.5))

# boxplot(delta,ylab="Mean Strength absolute difference ENG1-ENG2",xlab="ENG - Nodes", col=couleur3[c(1,1,2,2,2,1,1,2,1)], ylim=c(0,13))
boxplot(delta,ylab="Mean Strength absolute difference NAM1-NAM2",xlab="NAM - Nodes", col=couleur3[c(2,1,2,1,2,1,2,1,1)], ylim=c(0,10))
###

df1 <- melt(df1)
df2 <- melt(df2)

data_node_strength <- rbind(df1,df2)
data_node_strength <- cbind(data_node_strength,rep(1:2, each=9000))
colnames(data_node_strength) <- c("Index","Node","Strength","NAM_groupe")
data_node_strength$NAM_groupe <- as.factor(data_node_strength$NAM_groupe)


stat.test <- data_node_strength %>%
  group_by(Node) %>%
  t_test(Strength ~ NAM_groupe) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")
stat.test

# Créer un box plot
bxp <- ggboxplot(
  data_node_strength, x = "Node", y = "Strength", 
  color = "NAM_groupe", palette = c("#00AFBB", "#E7B800")
)

# Ajoutez des p-values sur les graphiques en box plot
stat.test <- stat.test %>%
  add_xy_position(x = "Node", dodge = 0.8)

bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}{p.adj.signif}", 
  tip.length = 0.02, step.increase = 0.05, hide.ns = FALSE
)




###############

###############

node_harmonic_1 = boot(eng_1, node_harmonic, 1000)


# boxplot(node_harmonic_1$t[,1:9],col = Couleurs3[3], names = names(data[,16:24]), las=2, ylim=c(0,30))
# title(main = "eng_1 Node harmonic")

node_harmonic_2 = boot(eng_2, node_harmonic, 1000)


# boxplot(node_harmonic_2$t[,1:9],col = Couleurs3[3], names = names(data[,16:24]), las=2, ylim=c(0,30))
# title(main = "ENG_2 Node harmonic")


# pvalue <- array()
# for (i in 1:32){  
#   t <- t.test(node_harmonic_1$t[,i],node_harmonic_2$t[,i])
#   pvalue[i] <-  format(round(t$p.value, 3), nsmall = 3) 
# }
# pvalue_node_harmonic <- rbind(names(data[,16:47]), pvalue)
# write.csv(pvalue_node_harmonic,"pvalue_node_harmonic_ENG.csv")


df1 <- node_harmonic_1$t[,1:9]
colnames(df1) <- names(data[,16:24])

df2 <- node_harmonic_2$t[,1:9]
colnames(df2) <- names(data[,16:24])


# soustraction pour mesure d'écart entre sélectionnés et non-sélectionnés en ML
delta=abs(df1-df2)
# write.csv(delta,"Strength_NAM1-NAM2.csv")
delta_ML=subset(delta, select= c("VMI","AMI","NAM","TMT"))
# delta_ML=subset(delta, select= c("VCI","AMI","NAM","SFL"))
delta_noML=subset(delta, select= c("VCI","PRI","PFL","SFL","STR"))
# delta_noML=subset(delta, select= c("VMI","PRI","PFL","TMT","STR"))
delta_ML_mean = rowMeans(delta_ML)
delta_noML_mean = rowMeans(delta_noML)
t.test(delta_ML_mean,delta_noML_mean,alternative="greater")
boxplot(cbind(delta_ML_mean,delta_noML_mean),ylab="Mean Node Harmonic absolute difference",xlab="ENG Nodes group",ylim=c(0,900))

boxplot(delta,ylab="Mean Node Harmonic absolute difference ENG1-ENG2",xlab="ENG - Nodes", col=couleur3[c(1,1,2,2,2,1,1,2,1)], ylim=c(0,900))
# boxplot(delta,ylab="Mean Node Harmonic absolute difference NAM1-NAM2",xlab="NAM - Nodes", col=couleur3[c(2,1,2,1,2,1,2,1,1)], ylim=c(0,1100))
###

df1 <- melt(df1)
df2 <- melt(df2)

data_node_harmonic <- rbind(df1,df2)
data_node_harmonic <- cbind(data_node_harmonic,rep(1:2, each=9000))
colnames(data_node_harmonic) <- c("Index","Node","harmonic","NAM_groupe")
data_node_harmonic$NAM_groupe <- as.factor(data_node_harmonic$NAM_groupe)


stat.test <- data_node_harmonic %>%
  group_by(Node) %>%
  t_test(harmonic ~ NAM_groupe) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")
stat.test

# Créer un box plot
bxp <- ggboxplot(
  data_node_harmonic, x = "Node", y = "harmonic", 
  color = "NAM_groupe", palette = c("#00AFBB", "#E7B800")
)

# Ajoutez des p-values sur les graphiques en box plot
stat.test <- stat.test %>%
  add_xy_position(x = "Node", dodge = 0.8)

bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}{p.adj.signif}", 
  tip.length = 0.02, step.increase = 0.05, hide.ns = FALSE
)

###############
###############

node_degree_1 = boot(eng_1, node_degree, 1000)


# boxplot(node_degree_1$t[,1:9],col = Couleurs3[3], names = names(data[,16:24]), las=2, ylim=c(0,15))
# title(main = "ENG_1 Node degree")

node_degree_2 = boot(eng_2, node_degree, 1000)


# boxplot(node_degree_2$t[,1:9],col = Couleurs3[3], names = names(data[,16:24]), las=2, ylim=c(0,15))
# title(main = "ENG_2 Node degree")


# pvalue <- array()
# for (i in 1:32){  
#   t <- t.test(node_degree_1$t[,i],node_degree_2$t[,i])
#   pvalue[i] <-  format(round(t$p.value, 3), nsmall = 3) 
# }
# pvalue_node_degree <- rbind(names(data[,16:47]), pvalue)
# write.csv(pvalue_node_degree,"pvalue_node_degree_ENG.csv")


df1 <- node_degree_1$t[,1:9]
colnames(df1) <- names(data[,16:24])

df2 <- node_degree_2$t[,1:9]
colnames(df2) <- names(data[,16:24])

# soustraction pour mesure d'écart entre sélectionnés et non-sélectionnés en ML
delta=abs(df1-df2)
write.csv(delta,"degree_NAM1-NAM2.csv")
delta_ML=subset(delta, select= c("VMI","AMI","NAM","TMT"))
# delta_ML=subset(delta, select= c("VCI","AMI","NAM","SFL"))
delta_noML=subset(delta, select= c("VCI","PRI","PFL","SFL","STR"))
# delta_noML=subset(delta, select= c("VMI","PRI","PFL","TMT","STR"))
delta_ML_mean = rowMeans(delta_ML)
delta_noML_mean = rowMeans(delta_noML)
t.test(delta_ML_mean,delta_noML_mean,alternative="greater")
boxplot(cbind(delta_ML_mean,delta_noML_mean),ylab="Mean degree absolute difference",xlab="ENG Nodes group",ylim=c(0,9))

boxplot(delta,ylab="Mean degree absolute difference ENG1-ENG2",xlab="ENG - Nodes", col=couleur3[c(1,1,2,2,2,1,1,2,1)], ylim=c(0,20))
# boxplot(delta,ylab="Mean degree absolute difference NAM1-NAM2",xlab="NAM - Nodes", col=couleur3[c(2,1,2,1,2,1,2,1,1)], ylim=c(0,20))
###

df1 <- melt(df1)
df2 <- melt(df2)

data_node_degree <- rbind(df1,df2)
data_node_degree <- cbind(data_node_degree,rep(1:2, each=9000))
colnames(data_node_degree) <- c("Index","Node","degree","NAM_groupe")
data_node_degree$NAM_groupe <- as.factor(data_node_degree$NAM_groupe)


stat.test <- data_node_degree %>%
  group_by(Node) %>%
  t_test(degree ~ NAM_groupe) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")
stat.test

# Créer un box plot
bxp <- ggboxplot(
  data_node_degree, x = "Node", y = "degree", 
  color = "NAM_groupe", palette = c("#00AFBB", "#E7B800")
)

# Ajoutez des p-values sur les graphiques en box plot
stat.test <- stat.test %>%
  add_xy_position(x = "Node", dodge = 0.8)

bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}{p.adj.signif}", 
  tip.length = 0.02, step.increase = 0.05, hide.ns = FALSE
)




###############

###############

# ABAQUE

######################################################################################

data_eng <- subset(data, select = c("VMI","AMI","NAM","TMT","ENG_bin"))

library("rms")

t.data <- datadist(data_eng)

options(datadist = 't.data')

fit <- lrm(formula = ENG_bin ~ VMI + AMI + NAM + TMT , data = data_eng, maxit=1000, x=TRUE, y=TRUE)

summary(fit)


plot(anova(fit), what='proportion chisq') # relative importance
plot(Predict(fit, fun=plogis)) # predicted values
rms::validate(fit, method="boot", B=1000) # bootstrapped validation
my.calib <- rms::calibrate(fit, method="boot", B=1000) # model calibration
plot(my.calib, las=1)

plot(nomogram(fit, fun = function(x)plogis(x)))

val <- validate(fit, method="boot", B=1000)
(c_opt_corr <- 0.5 * (val[1, 5] + 1)) # AUC ?


######

data_nam <- subset(data, select = c("VCI","AMI","NAM","SFL","DIFF_NAM_1"))


library("rms")

t.data2 <- datadist(data_nam)

options(datadist = 't.data2')

fit <- lrm(formula = DIFF_NAM_1 ~ VCI + AMI + NAM + SFL , data = data_nam, maxit=1000, x=TRUE, y=TRUE)

summary(fit)


plot(anova(fit), what='proportion chisq') # relative importance
plot(Predict(fit, fun=plogis)) # predicted values
rms::validate(fit, method="boot", B=1000) # bootstrapped validation
my.calib <- rms::calibrate(fit, method="boot", B=1000) # model calibration
plot(my.calib, las=1)

plot(nomogram(fit, fun = function(x)plogis(x)))

residuals(fit,type="gof")

val <- validate(fit, method="boot", B=1000)
(c_opt_corr <- 0.5 * (val[1, 5] + 1)) # AUC ?
## doc : https://www.nicholas-ollberding.com/post/an-introduction-to-the-harrell-verse-predictive-modeling-using-the-hmisc-and-rms-packages/

library("rms.gof")

######################################################################################
