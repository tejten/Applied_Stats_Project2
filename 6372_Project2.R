# Include Libraries
library(tidyverse)
library(caret)
library(ggcorrplot)
library(kernlab)
library(ggplot2)

bc<-read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",header=F,sep=",")
names(bc)<- c('id_number', 'diagnosis', 'radius_mean', 
              'texture_mean', 'perimeter_mean', 'area_mean', 
              'smoothness_mean', 'compactness_mean', 
              'concavity_mean','concave_points_mean', 
              'symmetry_mean', 'fractal_dimension_mean',
              'radius_se', 'texture_se', 'perimeter_se', 
              'area_se', 'smoothness_se', 'compactness_se', 
              'concavity_se', 'concave_points_se', 
              'symmetry_se', 'fractal_dimension_se', 
              'radius_worst', 'texture_worst', 
              'perimeter_worst', 'area_worst', 
              'smoothness_worst', 'compactness_worst', 
              'concavity_worst', 'concave_points_worst', 
              'symmetry_worst', 'fractal_dimension_worst')

# Data Summary
summary(bc)

# Normalize Data
bc.clean <- bc[,-c(1)]
normalize <- function(x){
  return (( x - min(x))/(max(x) -min(x)))
}  
bc.clean.normalized <- as.data.frame(
  lapply(bc.clean[,2:31],normalize)
)  
bc.clean.normalized <- cbind(
  bc.clean[,1],
  bc.clean.normalized
)
names(bc.clean.normalized)[1] <- "diagnosis"

summary(bc.clean.normalized)

## Highly Correlated Features
#Highly correlated features sometime shows feature duplication. Dropping those feature would not affect end result.For this analysis I am using correlation cutoff 0.85. If any correlation is higher than cutoff then I will simply drop that feature(s) from analysis.
#```{r}
correlationMatrix <- cor(bc.clean.normalized[,2:31])
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix,cutoff = 0.85)
print(highlyCorrelated)
names(bc.clean.normalized[,2:31])
post.cor.test.cancer <- cbind(
  bc.clean.normalized[,2:31][,-c(highlyCorrelated)],
  bc.clean.normalized[,1]
)
names(post.cor.test.cancer)[18] <-  "diagnosis"

names(bc.clean.normalized[,2:31])[c(highlyCorrelated)]

#Getting a look at the distribution
table(bc$diagnosis)

# Malignant and Benign Distribution
m_and_b <- bc.clean %>% 
  group_by(diagnosis) %>%
  summarise(n = n()) %>%
  mutate(percentage = signif((100 * n/sum(n)),2))

ggplot(data = m_and_b) +
  geom_bar(
    mapping = aes(x = "",y = percentage, fill = diagnosis), 
    stat = "identity", 
    width = 1) +
  geom_text(
    mapping = aes(x = c(1,1), y = c(69,18), 
                  label = paste(percentage,"%")), 
    size = 3) +
  coord_polar("y")

#Scatter plots color coded by response for just the first few variables
#NOTE: updating this to include features 3 - 11. Strong correlation with Area, Perimeter, Radius
pairs(bc[,3:11],col=bc$diagnosis)

#Seeing has how a lot of the attributes are strongly correlated, we will use PCA to convert attributes into a set of uncorrelated components.

#Conduct PCA
pc.bc<-prcomp(bc.clean.normalized[,-c(1)],scale.=TRUE)
pc.bc.scores<-pc.bc$x

#Adding the response column to the PC's data frame
pc.bc.scores<-data.frame(pc.bc.scores)
pc.bc.scores$Diagnosis<-bc$diagnosis


#Use ggplot2 to plot the first few pc's
library(ggplot2)
ggplot(data = pc.bc.scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(col=Diagnosis), size=1)+
  ggtitle("PCA of Breast Cancer Tumor Biopsies")

ggplot(data = pc.bc.scores, aes(x = PC2, y = PC3)) +
  geom_point(aes(col=Diagnosis), size=1)+
  ggtitle("PCA of Breast Cancer Tumor Biopsies")

summary(pc.bc)
library(factoextra)
fviz_eig(pc.bc, addlabels = TRUE, ylim = c(0,100), barfill = "steelblue1", line="navy") + 
  theme_classic() +
  labs(x = "Principal Components", y = "% of Explained Variance", title = "WDBC - Principal Components")

# We see that 44.3% of the variance is explained by the first principal component.

## Model Creation

# Models will be created using 5-fold cross-validation, given the relatively small sample size of the dataset. Setting parameters below:
  
# Setting up 5-fold cross-validation
ctrl <- trainControl(method = "cv",
                     number = 5)

# Function for plotting confusion matrices
cm_plot <- function(ml, title) {
  confusionMatrix(ml)$table %>%
    round(1) %>%
    fourfoldplot(
      color = c("#CC6666", "#99CC99"),
      main=title, 
      conf.level=0, 
      margin=1
    )
}

### Logistic Regression

#Conduct PCA
pca_wdbc <- princomp(bc.clean.normalized[,-c(1)]) # PCA on attributes
pc_wdbc <- pca_wdbc$scores # PCA scores
pc_wdbc_c <- bc$diagnosis # WDBC class attribute
full_wdbc <- data.frame(pc_wdbc,pc_wdbc_c) # Combining PC with class attribute

summary(pca_wdbc)
fviz_eig(pca_wdbc, addlabels = TRUE, ylim = c(0,100), barfill = "steelblue1", line="navy") + 
  theme_classic() +
  labs(x = "Principal Components", y = "% of Explained Variance", title = "WDBC - Principal Components")

# We see that 44.3% of the variance is explained by the first principal component.

logit.ml <- train(pc_wdbc_c~., full_wdbc, method = "glm", family = "binomial", trControl =ctrl)
logit.cm <- confusionMatrix(logit.ml)
cm_plot(logit.ml, "Logistic Regression")
logit.metrics <- data.frame (
  "Model" = "Logistic Regression",
  "Accuracy" = (logit.cm$table[1,1] + logit.cm$table[2,2])/100,
  "Recall" = logit.cm$table[2,2] / (logit.cm$table[2,2] + logit.cm$table[1,2]),
  "Precision" = logit.cm$table[2,2] / (logit.cm$table[2,1] + logit.cm$table[2,2]),
  "FNR" = (logit.cm$table[1,2] / (logit.cm$table[2,2] + logit.cm$table[1,2])),
  "Fscore" = (2 * logit.cm$table[2,2]) / (2 * logit.cm$table[2,2] + logit.cm$table[1,2] + logit.cm$table[2,1])
)
logit.metrics

### k-Nearest Neighbours

knn.ml <- train(pc_wdbc_c~., full_wdbc, method = "knn", trControl =ctrl)
knn.cm <- confusionMatrix(knn.ml)
cm_plot(knn.ml, "kNN")
knn.metrics <- data.frame (
  "Model" = "k-NN",
  "Accuracy" = (knn.cm$table[1,1] + knn.cm$table[2,2])/100,
  "Recall" = knn.cm$table[2,2] / (knn.cm$table[2,2] + knn.cm$table[1,2]),
  "Precision" = knn.cm$table[2,2] / (knn.cm$table[2,1] + knn.cm$table[2,2]),
  "FNR" = (knn.cm$table[1,2] / (knn.cm$table[2,2] + knn.cm$table[1,2])),
  "Fscore" = (2 * knn.cm$table[2,2]) / (2 * knn.cm$table[2,2] + knn.cm$table[1,2] + knn.cm$table[2,1])
)
knn.metrics

### Bagging - Random Forest

rf.ml <- train(pc_wdbc_c~., full_wdbc, method = "rf", trControl =ctrl)
rf.cm <- confusionMatrix(rf.ml)
cm_plot(rf.ml, "Random Forest")
rf.metrics <- data.frame (
  "Model" = "Random Forest",
  "Accuracy" = (rf.cm$table[1,1] + rf.cm$table[2,2])/100,
  "Recall" = rf.cm$table[2,2] / (rf.cm$table[2,2] + rf.cm$table[1,2]),
  "Precision" = rf.cm$table[2,2] / (rf.cm$table[2,1] + rf.cm$table[2,2]),
  "FNR" = (rf.cm$table[1,2] / (rf.cm$table[2,2] + rf.cm$table[1,2])),
  "Fscore" = (2 * rf.cm$table[2,2]) / (2 * rf.cm$table[2,2] + rf.cm$table[1,2] + rf.cm$table[2,1])
)
rf.metrics

## Model Performance

metrics1 <- rbind(knn.metrics,logit.metrics, rf.metrics)
metrics1 # Taking a look at everything together



ggplot(metrics1, aes(Model, Accuracy)) + geom_bar(stat="identity", aes(fill=Model)) + coord_cartesian(ylim=c(0.9,1)) + ggtitle("Accuracy")
ggplot(metrics1, aes(Model, Recall)) + geom_bar(stat="identity", aes(fill=Model)) + coord_cartesian(ylim=c(0.9,1)) + ggtitle("Recall")
ggplot(metrics1, aes(Model, Precision)) + geom_bar(stat="identity", aes(fill=Model)) + coord_cartesian(ylim=c(0.8,1)) + ggtitle("Precision")
ggplot(metrics1, aes(Model, FNR)) + geom_bar(stat="identity", aes(fill=Model)) + coord_cartesian(ylim=c(0,0.05)) + ggtitle("False Negative Rate")
ggplot(metrics1, aes(Model, Fscore)) + geom_bar(stat="identity", aes(fill=Model)) + coord_cartesian(ylim=c(0.9,1)) + ggtitle("F score")

#Taking a look at all confusion matrices:
  
par(mfrow=c(3,3))
cm_plot(knn.ml, "k-NN")
cm_plot(logit.ml, "Logistic Regression")
cm_plot(rf.ml, "Random Forest")

