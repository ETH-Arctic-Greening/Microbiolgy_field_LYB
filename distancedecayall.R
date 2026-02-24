
#############################################################
#ok new try!!
ps.propr <- ps.prorra
ps.propr.a <- ps.propr %>%
  subset_samples(site1 %in% c("ADT"))
ps.propr.b <- ps.propr %>%
  subset_samples(site1 %in% c("BJB"))
ps.propr.l <- ps.propr %>%
  subset_samples(site1 %in% c("LOD"))
###############
#ADT
abrel_bray.a <- phyloseq::distance(ps.propr.a, method = "bray")
abrel_bray.a <- as.matrix(abrel_bray.a)
#bray-curtis distance matrix
upper_triangular.a <- upper.tri(abrel_bray.a)
indices.a <- which(upper_triangular.a, arr.ind = TRUE)
# Create a list of comparisons
comparisons.a <- data.frame(
  Point1 = rownames(abrel_bray.a)[indices.a[,1]],
  Point2 = rownames(abrel_bray.a)[indices.a[,2]],
  Distance_braya = abrel_bray.a[upper_triangular.a]
)

#calculate distance between all plots ADT
coorda <- read.csv("coordinatesADT.csv")
library(geosphere)
distance.a <- distm(coorda[, c("Long","Lat")])
rownames(distance.a) <- coorda$Name
colnames(distance.a) <- coorda$Name

upper_triangular <- upper.tri(distance.a)

# Get the row and column indices of the upper triangular elements
indices.a <- which(upper_triangular, arr.ind = TRUE)
# Create a list of comparisons
comparisonsa <- data.frame(
  Point1 = rownames(distance.a)[indices.a[,1]],
  Point2 = rownames(distance.a)[indices.a[,2]],
  Distance = distance.a[upper_triangular]
)
distdeca <- cbind(comparisons.a,comparisonsa$Distance)
theme_set(theme_bw())
pa <- ggplot(distdeca, aes(x=comparisonsa$Distance, y=1-Distance_braya)) + 
  geom_point()+
  geom_smooth(method=lm)+
  ylim(0,1)+
  xlim(0,180)+
  labs(title="Distance Decay Adventdalen", x="Distance between plots [m]", y="1-Bray-Curtis Dissimilarity ASV level")
pa
model_linear <- lm((1-distdeca$Distance_braya) ~ comparisonsa$Distance)
summary(model_linear)
###############
#BJB
abrel_bray.b <- phyloseq::distance(ps.propr.b, method = "bray")
abrel_bray.b <- as.matrix(abrel_bray.b)
#bray-curtis distance matrix
upper_triangular.b <- upper.tri(abrel_bray.b)
indices.b <- which(upper_triangular.b, arr.ind = TRUE)
# Create a list of comparisons
comparisons.b <- data.frame(
  Point1 = rownames(abrel_bray.b)[indices.b[,1]],
  Point2 = rownames(abrel_bray.b)[indices.b[,2]],
  Distance_brayb = abrel_bray.b[upper_triangular.b]
)

#calculate distance between all plots ADT
coordb <- read.csv("coordinatesBJB.csv")
library(geosphere)
distance.b <- distm(coordb[, c("Long","Lat")])
rownames(distance.b) <- coordb$Name
colnames(distance.b) <- coordb$Name

upper_triangular <- upper.tri(distance.b)

# Get the row and column indices of the upper triangular elements
indices.b <- which(upper_triangular, arr.ind = TRUE)
# Create a list of comparisons
comparisonsb <- data.frame(
  Point1 = rownames(distance.b)[indices.b[,1]],
  Point2 = rownames(distance.b)[indices.b[,2]],
  Distance = distance.b[upper_triangular]
)
distdecb <- cbind(comparisons.b,comparisonsb$Distance)
theme_set(theme_bw())
pb <- ggplot(distdecb, aes(x=comparisonsb$Distance, y=1-Distance_brayb)) + 
  geom_point()+
  geom_smooth(method=lm)+
  ylim(0,1)+
  xlim(0,180)+
  labs(title="Distance Decay Bjørndalen", x="Distance between plots [m]", y="1-Bray-Curtis Dissimilarity ASV level")
###############
#LOD
abrel_bray.l <- phyloseq::distance(ps.propr.l, method = "bray")
abrel_bray.l <- as.matrix(abrel_bray.l)
#bray-curtis distance matrix
upper_triangular.l <- upper.tri(abrel_bray.l)
indices.l <- which(upper_triangular.l, arr.ind = TRUE)
# Create a list of comparisons
comparisons.l <- data.frame(
  Point1 = rownames(abrel_bray.l)[indices.l[,1]],
  Point2 = rownames(abrel_bray.l)[indices.l[,2]],
  Distance_brayl = abrel_bray.l[upper_triangular.l]
)

#calculate distance between all plots LOD
coordl <- read.csv("coordinatesLOD.csv")
library(geosphere)
distance.l <- distm(coordl[, c("Long","Lat")])
rownames(distance.l) <- coordl$Name
colnames(distance.l) <- coordl$Name

upper_triangular <- upper.tri(distance.l)

# Get the row and column indices of the upper triangular elements
indices.l <- which(upper_triangular, arr.ind = TRUE)
# Create a list of comparisons
comparisonsl <- data.frame(
  Point1 = rownames(distance.l)[indices.l[,1]],
  Point2 = rownames(distance.l)[indices.l[,2]],
  Distance = distance.l[upper_triangular]
)
distdecl <- cbind(comparisons.l,comparisonsl$Distance)
theme_set(theme_bw())
pl <- ggplot(distdecl, aes(x=comparisonsl$Distance, y=1-Distance_brayl)) + 
  geom_point()+
  geom_smooth(method=lm)+
  ylim(0,1)+
  xlim(0,180)+
  labs(title="Distance Decay Gamlebyen", x="Distance between plots [m]", y="1-Bray-Curtis Dissimilarity ASV level")
library(gridExtra)
grid.arrange(pa, pb, pl, ncol = 3)

#############################################
#######################
##### now do this on class level!
#start from "glomrclass"
saveRDS(glomr, "glomrclass.rds")
ps.propr <- transform_sample_counts(glomrclass, function(x) x/sum(x))

ps.propr.a <- ps.propr %>%
  subset_samples(site1 %in% c("ADT"))
ps.propr.b <- ps.propr %>%
  subset_samples(site1 %in% c("BJB"))
ps.propr.l <- ps.propr %>%
  subset_samples(site1 %in% c("LOD"))
###############
#ADT
abrel_bray.a <- phyloseq::distance(ps.propr.a, method = "bray")
abrel_bray.a <- as.matrix(abrel_bray.a)
#bray-curtis distance matrix
upper_triangular.a <- upper.tri(abrel_bray.a)
indices.a <- which(upper_triangular.a, arr.ind = TRUE)
# Create a list of comparisons
comparisons.a <- data.frame(
  Point1 = rownames(abrel_bray.a)[indices.a[,1]],
  Point2 = rownames(abrel_bray.a)[indices.a[,2]],
  Distance_braya = abrel_bray.a[upper_triangular.a]
)

#calculate distance between all plots ADT
coorda <- read.csv("coordinatesADT.csv")
library(geosphere)
distance.a <- distm(coorda[, c("Long","Lat")])
rownames(distance.a) <- coorda$Name
colnames(distance.a) <- coorda$Name

upper_triangular <- upper.tri(distance.a)

# Get the row and column indices of the upper triangular elements
indices.a <- which(upper_triangular, arr.ind = TRUE)
# Create a list of comparisons
comparisonsa <- data.frame(
  Point1 = rownames(distance.a)[indices.a[,1]],
  Point2 = rownames(distance.a)[indices.a[,2]],
  Distance = distance.a[upper_triangular]
)
distdeca <- cbind(comparisons.a,comparisonsa$Distance)
theme_set(theme_bw())
pac <- ggplot(distdeca, aes(x=comparisonsa$Distance, y=1-Distance_braya)) + 
  geom_point()+
  geom_smooth(method=lm)+
  ylim(0,1)+
  xlim(0,180)+
  labs(title="Distance Decay Adventdalen", x="Distance between plots [m]", y="1-Bray-Curtis Dissimilarity Class level")
pac
###############
#BJB
abrel_bray.b <- phyloseq::distance(ps.propr.b, method = "bray")
abrel_bray.b <- as.matrix(abrel_bray.b)
#bray-curtis distance matrix
upper_triangular.b <- upper.tri(abrel_bray.b)
indices.b <- which(upper_triangular.b, arr.ind = TRUE)
# Create a list of comparisons
comparisons.b <- data.frame(
  Point1 = rownames(abrel_bray.b)[indices.b[,1]],
  Point2 = rownames(abrel_bray.b)[indices.b[,2]],
  Distance_brayb = abrel_bray.b[upper_triangular.b]
)

#calculate distance between all plots ADT
coordb <- read.csv("coordinatesBJB.csv")
library(geosphere)
distance.b <- distm(coordb[, c("Long","Lat")])
rownames(distance.b) <- coordb$Name
colnames(distance.b) <- coordb$Name

upper_triangular <- upper.tri(distance.b)

# Get the row and column indices of the upper triangular elements
indices.b <- which(upper_triangular, arr.ind = TRUE)
# Create a list of comparisons
comparisonsb <- data.frame(
  Point1 = rownames(distance.b)[indices.b[,1]],
  Point2 = rownames(distance.b)[indices.b[,2]],
  Distance = distance.b[upper_triangular]
)
distdecb <- cbind(comparisons.b,comparisonsb$Distance)
theme_set(theme_bw())
pbc <- ggplot(distdecb, aes(x=comparisonsb$Distance, y=1-Distance_brayb)) + 
  geom_point()+
  geom_smooth(method=lm)+
  ylim(0,1)+
  xlim(0,180)+
  labs(title="Distance Decay Bjørndalen", x="Distance between plots [m]", y="1-Bray-Curtis Dissimilarity Class level")
###############
#LOD
abrel_bray.l <- phyloseq::distance(ps.propr.l, method = "bray")
abrel_bray.l <- as.matrix(abrel_bray.l)
#bray-curtis distance matrix
upper_triangular.l <- upper.tri(abrel_bray.l)
indices.l <- which(upper_triangular.l, arr.ind = TRUE)
# Create a list of comparisons
comparisons.l <- data.frame(
  Point1 = rownames(abrel_bray.l)[indices.l[,1]],
  Point2 = rownames(abrel_bray.l)[indices.l[,2]],
  Distance_brayl = abrel_bray.l[upper_triangular.l]
)

#calculate distance between all plots LOD
coordl <- read.csv("coordinatesLOD.csv")
library(geosphere)
distance.l <- distm(coordl[, c("Long","Lat")])
rownames(distance.l) <- coordl$Name
colnames(distance.l) <- coordl$Name

upper_triangular <- upper.tri(distance.l)

# Get the row and column indices of the upper triangular elements
indices.l <- which(upper_triangular, arr.ind = TRUE)
# Create a list of comparisons
comparisonsl <- data.frame(
  Point1 = rownames(distance.l)[indices.l[,1]],
  Point2 = rownames(distance.l)[indices.l[,2]],
  Distance = distance.l[upper_triangular]
)
distdecl <- cbind(comparisons.l,comparisonsl$Distance)
theme_set(theme_bw())
plc <- ggplot(distdecl, aes(x=comparisonsl$Distance, y=1-Distance_brayl)) + 
  geom_point()+
  geom_smooth(method=lm)+
  ylim(0,1)+
  xlim(0,180)+
  labs(title="Distance Decay Gamlebyen", x="Distance between plots [m]", y="1-Bray-Curtis Dissimilarity Class level")
plc
library(gridExtra)
grid.arrange(pa, pb, pl, pac, pbc, plc, ncol = 3)

################################################
##### do it with nutrients calculate difference between plots!
# Step 1: Load the data
###ADT
data <- read.csv('predADT.csv')
#pH
# Step 2: Set up the difference matrix
# Extract values and sample names
samples <- data$sample
values <- data$pH.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularph <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesph.a <- which(upper_triangularph, arr.ind = TRUE)
# Create a list of comparisons
comparisonsaph <- data.frame(
  Point1 = rownames(difference_matrix)[indicesph.a[,1]],
  Point2 = rownames(difference_matrix)[indicesph.a[,2]],
  Distance = difference_matrix[upper_triangularph]
)

# Step 4: Save to a CSV or display

write.csv(comparisonsaph, 'distancesapH.csv')  # Save to CSV
distdecaph <- cbind(comparisons.a,comparisonsaph$Distance)
library(dplyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())
paph <- ggplot(distdecaph, aes(x=comparisonsaph$Distance, y=1-Distance_braya)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Native Tundra pH", x="difference in pH", y="1-Bray-Curtis Dissimilarity")
paph


##### TC
# Extract values and sample names
samples <- data$sample
values <- data$TC.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulartc <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicestc.a <- which(upper_triangulartc, arr.ind = TRUE)
# Create a list of comparisons
comparisonsatc <- data.frame(
  Point1 = rownames(difference_matrix)[indicesph.a[,1]],
  Point2 = rownames(difference_matrix)[indicesph.a[,2]],
  Distance = difference_matrix[upper_triangulartc]
)

write.csv(comparisonsatc, 'distancesaTC.csv')  # Save to CSV
distdecatc <- cbind(comparisons.a,comparisonsatc$Distance)

theme_set(theme_bw())
patc <- ggplot(distdecatc, aes(x=comparisonsatc$Distance, y=1-Distance_braya)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Native Tundra TC", x="difference in TC", y="1-Bray-Curtis Dissimilarity")
patc

##### TN
samples <- data$sample
values <- data$TN.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulartn <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicestn.a <- which(upper_triangulartn, arr.ind = TRUE)
# Create a list of comparisons
comparisonsatn <- data.frame(
  Point1 = rownames(difference_matrix)[indicestn.a[,1]],
  Point2 = rownames(difference_matrix)[indicestn.a[,2]],
  Distance = difference_matrix[upper_triangulartn]
)

write.csv(comparisonsatn, 'distancesaTN.csv')  # Save to CSV
distdecatn <- cbind(comparisons.a,comparisonsatn$Distance)

theme_set(theme_bw())
patn <- ggplot(distdecatn, aes(x=comparisonsatn$Distance, y=1-Distance_braya)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Native Tundra TN", x="difference in TN", y="1-Bray-Curtis Dissimilarity")
patn

#### CEC
samples <- data$sample
values <- data$CEC.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularcec <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicescec.a <- which(upper_triangularcec, arr.ind = TRUE)
# Create a list of comparisons
comparisonsacec <- data.frame(
  Point1 = rownames(difference_matrix)[indicescec.a[,1]],
  Point2 = rownames(difference_matrix)[indicescec.a[,2]],
  Distance = difference_matrix[upper_triangularcec]
)

write.csv(comparisonsacec, 'distancesaCEC.csv')  # Save to CSV
distdecacec <- cbind(comparisons.a,comparisonsacec$Distance)

theme_set(theme_bw())
pacec <- ggplot(distdecacec, aes(x=comparisonsacec$Distance, y=1-Distance_braya)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Native Tundra CEC", x="difference in CEC", y="1-Bray-Curtis Dissimilarity")
pacec
library(gridExtra)
grid.arrange(paph, patn, patc, pacec, ncol = 4)


##### do it with nutrients calculate difference between plots!
# Step 1: Load the data
#BJB
data <- read.csv('predBJB.csv')
#pH
# Step 2: Set up the difference matrix
# Extract values and sample names
samples <- data$Name
values <- data$pH.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularph <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesph.b <- which(upper_triangularph, arr.ind = TRUE)
# Create a list of comparisons
comparisonsbph <- data.frame(
  Point1 = rownames(difference_matrix)[indicesph.b[,1]],
  Point2 = rownames(difference_matrix)[indicesph.b[,2]],
  Distance = difference_matrix[upper_triangularph]
)

# Step 4: Save to a CSV or display

write.csv(comparisonsbph, 'distancesbpH.csv')  # Save to CSV
distdecbph <- cbind(comparisons.b,comparisonsbph$Distance)
library(dplyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())
pbph <- ggplot(distdecbph, aes(x=comparisonsbph$Distance, y=1-Distance_brayb)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Bird Cliff pH", x="difference in pH", y="1-Bray-Curtis Dissimilarity")
pbph

##### TC
# Extract values and sample names

values <- data$TC.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulartc <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicestc.b <- which(upper_triangulartc, arr.ind = TRUE)
# Create a list of comparisons
comparisonsbtc <- data.frame(
  Point1 = rownames(difference_matrix)[indicesph.b[,1]],
  Point2 = rownames(difference_matrix)[indicesph.b[,2]],
  Distance = difference_matrix[upper_triangulartc]
)

write.csv(comparisonsbtc, 'distancesbTC.csv')  # Save to CSV
distdecbtc <- cbind(comparisons.b,comparisonsbtc$Distance)

theme_set(theme_bw())
pbtc <- ggplot(distdecbtc, aes(x=comparisonsbtc$Distance, y=1-Distance_brayb)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Bird Cliff TC", x="difference in TC", y="1-Bray-Curtis Dissimilarity")
pbtc

##### TN

values <- data$TN.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulartn <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicestn.b <- which(upper_triangulartn, arr.ind = TRUE)
# Create a list of comparisons
comparisonsbtn <- data.frame(
  Point1 = rownames(difference_matrix)[indicestn.b[,1]],
  Point2 = rownames(difference_matrix)[indicestn.b[,2]],
  Distance = difference_matrix[upper_triangulartn]
)

write.csv(comparisonsbtn, 'distancesbTN.csv')  # Save to CSV
distdecbtn <- cbind(comparisons.b,comparisonsbtn$Distance)

theme_set(theme_bw())
pbtn <- ggplot(distdecbtn, aes(x=comparisonsbtn$Distance, y=1-Distance_brayb)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Bird Cliff TN", x="difference in TN", y="1-Bray-Curtis Dissimilarity")
pbtn

#### CEC

values <- data$CEC.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularcec <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicescec.b <- which(upper_triangularcec, arr.ind = TRUE)
# Create a list of comparisons
comparisonsbcec <- data.frame(
  Point1 = rownames(difference_matrix)[indicescec.b[,1]],
  Point2 = rownames(difference_matrix)[indicescec.b[,2]],
  Distance = difference_matrix[upper_triangularcec]
)

write.csv(comparisonsbcec, 'distancesbCEC.csv')  # Save to CSV
distdecbcec <- cbind(comparisons.b,comparisonsbcec$Distance)

theme_set(theme_bw())
pbcec <- ggplot(distdecbcec, aes(x=comparisonsbcec$Distance, y=1-Distance_brayb)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Bird Cliff CEC", x="difference in CEC", y="1-Bray-Curtis Dissimilarity")
pbcec

grid.arrange(pbph, pbtn, pbtc, pbcec, ncol = 4)

############
#LOD

data <- read.csv('predLOD.csv')
#pH
# Step 2: Set up the difference matrix
# Extract values and sample names
samples <- data$Name
values <- data$pH.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularph <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesph.l <- which(upper_triangularph, arr.ind = TRUE)
# Create a list of comparisons
comparisonslph <- data.frame(
  Point1 = rownames(difference_matrix)[indicesph.l[,1]],
  Point2 = rownames(difference_matrix)[indicesph.l[,2]],
  Distance = difference_matrix[upper_triangularph]
)

# Step 4: Save to a CSV or display

write.csv(comparisonslph, 'distanceslpH.csv')  # Save to CSV
distdeclph <- cbind(comparisons.l,comparisonslph$Distance)
library(dplyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())
plph <- ggplot(distdeclph, aes(x=comparisonslph$Distance, y=1-Distance_brayl)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Disturbed pH", x="difference in pH", y="1-Bray-Curtis Dissimilarity")
plph

##### TC
# Extract values and sample names

values <- data$TC.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulartc <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicestc.l <- which(upper_triangulartc, arr.ind = TRUE)
# Create a list of comparisons
comparisonsltc <- data.frame(
  Point1 = rownames(difference_matrix)[indicesph.l[,1]],
  Point2 = rownames(difference_matrix)[indicesph.l[,2]],
  Distance = difference_matrix[upper_triangulartc]
)

write.csv(comparisonsltc, 'distanceslTC.csv')  # Save to CSV
distdecltc <- cbind(comparisons.l,comparisonsltc$Distance)

theme_set(theme_bw())
pltc <- ggplot(distdecltc, aes(x=comparisonsltc$Distance, y=1-Distance_brayl)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Disturbed TC", x="difference in TC", y="1-Bray-Curtis Dissimilarity")
pltc

##### TN

values <- data$TN.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulartn <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicestn.l <- which(upper_triangulartn, arr.ind = TRUE)
# Create a list of comparisons
comparisonsltn <- data.frame(
  Point1 = rownames(difference_matrix)[indicestn.l[,1]],
  Point2 = rownames(difference_matrix)[indicestn.l[,2]],
  Distance = difference_matrix[upper_triangulartn]
)

write.csv(comparisonsltn, 'distanceslTN.csv')  # Save to CSV
distdecltn <- cbind(comparisons.l,comparisonsltn$Distance)

theme_set(theme_bw())
pltn <- ggplot(distdecltn, aes(x=comparisonsltn$Distance, y=1-Distance_brayl)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Disturbed TN", x="difference in TN", y="1-Bray-Curtis Dissimilarity")
pltn

#### CEC

values <- data$CEC.prediction

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularcec <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicescec.l <- which(upper_triangularcec, arr.ind = TRUE)
# Create a list of comparisons
comparisonslcec <- data.frame(
  Point1 = rownames(difference_matrix)[indicescec.l[,1]],
  Point2 = rownames(difference_matrix)[indicescec.l[,2]],
  Distance = difference_matrix[upper_triangularcec]
)

write.csv(comparisonslcec, 'distanceslCEC.csv')  # Save to CSV
distdeclcec <- cbind(comparisons.l,comparisonslcec$Distance)

theme_set(theme_bw())
plcec <- ggplot(distdeclcec, aes(x=comparisonslcec$Distance, y=1-Distance_brayl)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Disturbed CEC", x="difference in CEC", y="1-Bray-Curtis Dissimilarity")
plcec

grid.arrange(plph, pltn, pltc, plcec, ncol = 4)
################
#now with Graminoids and Dwarf shrubs
#LOD

data <- read.csv('predLOD.csv')
#graminoids
# Step 2: Set up the difference matrix
# Extract values and sample names
samples <- data$Name
values <- data$graminoid

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularg <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesg.l <- which(upper_triangularg, arr.ind = TRUE)
# Create a list of comparisons
comparisonslg <- data.frame(
  Point1 = rownames(difference_matrix)[indicesg.l[,1]],
  Point2 = rownames(difference_matrix)[indicesg.l[,2]],
  Distance = difference_matrix[upper_triangularg]
)

# Step 4: Save to a CSV or display

write.csv(comparisonslg, 'distanceslg.csv')  # Save to CSV
distdeclg <- cbind(comparisons.l,comparisonslg$Distance)
library(dplyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())
plg <- ggplot(distdeclg, aes(x=comparisonslg$Distance, y=1-Distance_brayl)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Disturbed Graminoids", x="difference in graminoids", y="1-Bray-Curtis Dissimilarity")
plg

#dwarf shrubs
# Step 2: Set up the difference matrix
# Extract values and sample names
samples <- data$Name
values <- data$dwarf

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulard <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesd.l <- which(upper_triangulard, arr.ind = TRUE)
# Create a list of comparisons
comparisonsld <- data.frame(
  Point1 = rownames(difference_matrix)[indicesd.l[,1]],
  Point2 = rownames(difference_matrix)[indicesd.l[,2]],
  Distance = difference_matrix[upper_triangulard]
)

# Step 4: Save to a CSV or display

write.csv(comparisonsld, 'distancesld.csv')  # Save to CSV
distdecld <- cbind(comparisons.l,comparisonsld$Distance)
library(ggpmisc)
library(devtools)
pld <- ggplot(distdecld, aes(x=comparisonsld$Distance, y=1-Distance_brayl)) + 
  geom_point()+ 
  geom_smooth(method=lm)+
  labs(title="Distance Decay Disturbed Dwarf Shrubs", x="difference in dwarf shrubs", y="1-Bray-Curtis Dissimilarity")
pld
################
#now with Graminoids and Dwarf shrubs
#ADT

data <- read.csv('predADT.csv')
#graminoids
# Step 2: Set up the difference matrix
# Extract values and sample names
samples <- data$sample
values <- data$graminoid

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularg <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesg.a <- which(upper_triangularg, arr.ind = TRUE)
# Create a list of comparisons
comparisonsag <- data.frame(
  Point1 = rownames(difference_matrix)[indicesg.a[,1]],
  Point2 = rownames(difference_matrix)[indicesg.a[,2]],
  Distance = difference_matrix[upper_triangularg]
)

# Step 4: Save to a CSV or display

write.csv(comparisonsag, 'distancesag.csv')  # Save to CSV
distdecag <- cbind(comparisons.a,comparisonsag$Distance)
library(dplyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())
pag <- ggplot(distdecag, aes(x=comparisonsag$Distance, y=1-Distance_braya)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Tundra Graminoids", x="difference in graminoids", y="1-Bray-Curtis Dissimilarity")
pag

#dwarf shrubs
# Step 2: Set up the difference matrix
# Extract values and sample names
#samples <- data$Name
values <- data$dwarf

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulard <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesd.a <- which(upper_triangulard, arr.ind = TRUE)
# Create a list of comparisons
comparisonsad <- data.frame(
  Point1 = rownames(difference_matrix)[indicesd.a[,1]],
  Point2 = rownames(difference_matrix)[indicesd.a[,2]],
  Distance = difference_matrix[upper_triangulard]
)

# Step 4: Save to a CSV or display

write.csv(comparisonsad, 'distancesad.csv')  # Save to CSV
distdecad <- cbind(comparisons.a,comparisonsad$Distance)
library(ggpmisc)
library(devtools)
pad <- ggplot(distdecad, aes(x=comparisonsad$Distance, y=1-Distance_braya)) + 
  geom_point()+ 
  geom_smooth(method=lm)+
  labs(title="Distance Decay Tundra Dwarf Shrubs", x="difference in dwarf shrubs", y="1-Bray-Curtis Dissimilarity")
pad

##### BJB
data <- read.csv('predBJB.csv')
#graminoids

# Step 2: Set up the difference matrix
# Extract values and sample names

samples <- data$Name
values <- data$graminoid

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangularg <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesg.b <- which(upper_triangularg, arr.ind = TRUE)
# Create a list of comparisons
comparisonsbg <- data.frame(
  Point1 = rownames(difference_matrix)[indicesg.b[,1]],
  Point2 = rownames(difference_matrix)[indicesg.b[,2]],
  Distance = difference_matrix[upper_triangularg]
)

# Step 4: Save to a CSV or display

write.csv(comparisonsbg, 'distancesbg.csv')  # Save to CSV
distdecbg <- cbind(comparisons.b,comparisonsbg$Distance)
library(dplyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())
pbg <- ggplot(distdecbg, aes(x=comparisonsbg$Distance, y=1-Distance_brayb)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(title="Distance Decay Bird Cliff Graminoids", x="difference in graminoids", y="1-Bray-Curtis Dissimilarity")
pbg

#dwarf shrubs
# Step 2: Set up the difference matrix
# Extract values and sample names
#samples <- data$Name
values <- data$dwarf

# Create a matrix of absolute differences
difference_matrix <- outer(values, values, FUN = function(x, y) abs(x - y))

# Step 3: Set row and column names
rownames(difference_matrix) <- samples
colnames(difference_matrix) <- samples

upper_triangulard <- upper.tri(difference_matrix)
# Get the row and column indices of the upper triangular elements
indicesd.b <- which(upper_triangulard, arr.ind = TRUE)
# Create a list of comparisons
comparisonsbd <- data.frame(
  Point1 = rownames(difference_matrix)[indicesd.b[,1]],
  Point2 = rownames(difference_matrix)[indicesd.b[,2]],
  Distance = difference_matrix[upper_triangulard]
)

# Step 4: Save to a CSV or display

write.csv(comparisonsbd, 'distancesbd.csv')  # Save to CSV
distdecbd <- cbind(comparisons.b,comparisonsbd$Distance)
library(ggpmisc)
library(devtools)
pbd <- ggplot(distdecbd, aes(x=comparisonsbd$Distance, y=1-Distance_brayb)) + 
  geom_point()+ 
  geom_smooth(method=lm)+
  labs(title="Distance Decay Bird Cliff Dwarf Shrubs", x="difference in dwarf shrubs", y="1-Bray-Curtis Dissimilarity")
pbd

grid.arrange(pag, pad, pbg, pbd, plg, pld, ncol = 2)
