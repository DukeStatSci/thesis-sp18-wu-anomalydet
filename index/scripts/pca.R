setwd("~/Desktop/Stats Thesis/Dataset/")
library(plyr)
library(data.table)
library(factoextra)
library(ggbiplot)
argus = read.csv("argus-anon-20170201.csv")#[1:50000,] 
summary(argus)

sapply(argus, class)
argus = transform(argus, 
                  Sport = as.factor(Sport),
                  Dport = as.factor(Dport))
argus = subset(argus, select = c("Flgs", "SrcAddr", "Sport", "DstAddr", "Dport", 
                                 "SrcPkts", "DstPkts", "SrcBytes", "DstBytes", "State"))
attach(argus)
categorical = c("Flgs", "SrcAddr", "Sport", "DstAddr", "Dport", "State")
continuous = c("SrcPkts", "DstPkts", "SrcBytes", "DstBytes")

#set counts
s_num = 25
d_num = 25
combo_num = 25
#get freqs
Sport_table = table(Sport)
Sport_table = as.data.frame(Sport_table)
Sport_table = Sport_table[order(-Sport_table$Freq),]
top_Sport = (head(Sport_table$Sport, s_num))

#get freqs
Dport_table = table(Dport)
Dport_table = as.data.frame(Dport_table)
Dport_table = Dport_table[order(-Dport_table$Freq),]
top_Dport = (head(Dport_table$Dport, d_num))

#subset data
argus_maxes = argus[is.element(Sport, top_Sport) & is.element(Dport, top_Dport), ]
argus_maxes = transform(argus_maxes, 
                        Sport = as.numeric(as.character(Sport)),
                        Dport = as.numeric(as.character(Dport)))
max_combinations = as.data.frame(table(argus_maxes$Sport, argus_maxes$Dport))

top_combinations = head(max_combinations[order(-max_combinations$Freq),], combo_num)
top_combinations$Sport = top_combinations$Var1
top_combinations$Dport = top_combinations$Var2
top_combinations$Var1 = NULL
top_combinations$Var2 = NULL
top_combinations = transform(top_combinations, 
                             Sport = as.numeric(as.character(Sport)),
                             Dport = as.numeric(as.character(Dport)))

extract_intersection = function(sport, dport){
  argus_subset = argus[Sport == sport & Dport == dport,]
  return (argus_subset)
}

generate_combinations_matrix = function(top_combinations){
  n = dim(top_combinations)[1]
  combinations = c()
  for (i in 1:n){
    sport = as.numeric(top_combinations[i,]$Sport)
    dport = as.numeric(top_combinations[i,]$Dport)
    combo = extract_intersection(sport, dport)
    combinations = c(combinations, list(combo))
  }
  return (combinations)
}
combinations = generate_combinations_matrix(top_combinations)
saveRDS(combinations, "combinations.rds")

#Normal Transformation
# normal_transform = function(y){
#   return (qnorm(rank(y, na.last = "keep")/(sum(!is.na(y)) + 1)))
# }
nscore <- function(x) {
  # Takes a vector of values x and calculates their normal scores. Returns 
  # a list with the scores and an ordered table of original values and
  # scores, which is useful as a back-transform table. See backtr().
  nscore = qqnorm(x, plot.it = FALSE)$x  # normal score 
  trn.table = data.frame(x=sort(x),nscore=sort(nscore))
  return (list(nscore=nscore, trn.table=trn.table))
}

#PCA
pca_analysis = function(SrcBytes, SrcPkts, DstBytes, DstPkts, Sport){
  pca_cont_vars = cbind(SrcBytes, SrcPkts, DstBytes, DstPkts)
  pca = prcomp(pca_cont_vars, center = TRUE, scale. = TRUE)
  print(pca$rotation)
  print((summary(pca)))
  screeplot(pca, type="lines",col=3)
  g = ggbiplot(pca, obs.scale = 1, var.scale = 1,
                ellipse = TRUE,
                circle = TRUE)
  g = g + scale_color_discrete(name = '')
  g = g + theme(legend.direction = 'horizontal',
                 legend.position = 'top')
  print(g)
  return(pca$rotation)
}

combo_table = combinations[1]
combo_table = transform(combo_table,
                        SrcBytes = as.numeric(SrcBytes),
                        SrcPkts = as.numeric(SrcPkts),
                        DstBytes = as.numeric(DstBytes),
                        DstPkts = as.numeric(DstPkts))
cat("Sport:", combo_table$Sport[1],"\t")
cat("Dport:", combo_table$Dport[1],"\n")
SrcBytes_norm = nscore(SrcBytes)$nscore
SrcPkts_norm = nscore(SrcPkts)$nscore
DstBytes_norm = nscore(DstBytes)$nscore
DstPkts_norm = nscore(DstPkts)$nscore
pca_analysis(SrcBytes_norm, SrcPkts_norm, DstBytes_norm, DstPkts_norm)

#Investigating Combinations
for (i in 1:combo_num){
  combo_table = combinations[i]
  combo_table = transform(combo_table, 
                          SrcBytes = as.numeric(SrcBytes),
                          SrcPkts = as.numeric(SrcPkts),
                          DstBytes = as.numeric(DstBytes),
                          DstPkts = as.numeric(DstPkts))
  cat("Sport:", combo_table$Sport[1],"\t")
  cat("Dport:", combo_table$Dport[1],"\n")
  SrcBytes_norm = normal_transform(SrcBytes)
  SrcPkts_norm = normal_transform(SrcPkts)
  DstBytes_norm = normal_transform(DstBytes)
  DstPkts_norm = normal_transform(DstPkts)
  pca_analysis(SrcBytes_norm, SrcPkts_norm, DstBytes_norm, DstPkts_norm, combo_table$Sport)
}
