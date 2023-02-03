## -----------------------------------------------------------------
##
## Script name: alzheimer_insilico.R
##
## Purpose of script: All the execution of insilico final project
##
## Author: Alba Vallejo Urruchi
##
## Date Created: 2023-12-01
##
## Email: alvau@alumni.uv.es
##
## -----------------------------------------------------------------


install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pheatmap")

library(pheatmap)
library(GEOquery)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(limma)

## GSE138214 
## GSE138261 

## Extracción de datos

id_gse <- "GSE138261"
gse138261_raw <- getGEO(id_gse)

length(gse138261_raw)

gse138261_raw <- gse138261_raw[[1]]

pData(gse138261_raw) ## información de variables
fData(gse138261_raw) ## anotaciones genes
exprs(gse138261_raw) ## datos de expresión

summary(exprs(gse138261_raw))

## Boxplot
# Antes se aplica el logaritmo en base 2

exprs(gse138261_raw) <- log2(exprs(gse138261_raw))
boxplot(exprs(gse138261_raw),outline=FALSE)

## Miramos qué variables nos interesan

sampleInfo <- pData(gse138261_raw)
sampleInfo

## Las tres columnas (variables) de interés son:
# "age:ch1","disease state:ch1","gender:ch1"

sampleInfo <- select(sampleInfo, "age:ch1","disease state:ch1","gender:ch1")

## Renombramos las columnas con nombres más fáciles:
# age, disease_state, gender
sampleInfo <- rename(sampleInfo,age = "age:ch1", disease_state="disease state:ch1",gender="gender:ch1")

#Clustering y PCA
# Generamos la matriz de correlación
corMatrix <- cor(exprs(gse138261_raw),use="c")
pheatmap(corMatrix)   

rownames(sampleInfo)
colnames(corMatrix)

## Ahora hacemos que las filas también tenga la información de gender, age, disease_state

rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)    

# PCA

pca <- prcomp(t(exprs(gse138261_raw)))

cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=disease_state,label=paste("Age", age))) + geom_point() + geom_text_repel()

## Expresion diferencial

# Normalización

gse138261 <- gse138261_raw
exprs(gse138261) = normalizeBetweenArrays(exprs(gse138261_raw))

## Nivel medio de expresión

probenames = rownames(gse138261)
pData(gse138261)
length(colnames(gse138261))
sampleId <- c("AD_1","AD_2","AD_3","AD_4","AD_5","AD_6","AD_7","AD_8","AD_9","AD_10","AD_11","AD_12","AD_13","AD_14","AD_15","AD_16","AD_17",
              "control_1", "control_2", "control_3", "control_4", "control_5", "control_6", "control_7", "control_8", "control_9", "control_10",
              "control_11","control_12","control_13","control_14","control_15","control_16", "control_17","control_18","control_19")
colnames(gse138261) <- sampleId
head(gse138261)

## Calculamos los valores medios de expresión para cada genotipo/condición 
## sumando las correspondientes columnas y dividiendo por el número de réplicas
expression.level = exprs(gse138261)
control = (expression.level[,"control_1"] + expression.level[,"control_1"]+ expression.level[,"control_2"]+ expression.level[,"control_3"]
           + expression.level[,"control_4"]+ expression.level[,"control_5"]+ expression.level[,"control_6"]
           + expression.level[,"control_7"]+ expression.level[,"control_8"]+ expression.level[,"control_9"]
           + expression.level[,"control_10"]+ expression.level[,"control_11"]+ expression.level[,"control_12"]
           + expression.level[,"control_13"]+ expression.level[,"control_14"]+ expression.level[,"control_15"]
           + expression.level[,"control_16"]+ expression.level[, "control_17"]+ expression.level[,"control_18"]
           + expression.level[,"control_19"])/19

AD = (expression.level[,"AD_1"] + expression.level[,"AD_1"]+ expression.level[,"AD_2"]+ expression.level[,"AD_3"]
      + expression.level[,"AD_4"]+ expression.level[,"AD_5"]+ expression.level[,"AD_6"]
      + expression.level[,"AD_7"]+ expression.level[,"AD_8"]+ expression.level[,"AD_9"]
      + expression.level[,"AD_10"]+ expression.level[,"AD_11"]+ expression.level[,"AD_12"]
      + expression.level[,"AD_13"]+ expression.level[,"AD_14"]+ expression.level[,"AD_15"]
      + expression.level[,"AD_16"]+ expression.level[, "AD_17"])/17

## Creamos una matriz que contenga por columna la expresión media para cada 
## condición o genotipo. Nombramos las filas con el nombre de las sondas (genes) 
## y la columnas con la condición o genotipo. 
mean.expression <- matrix(c(control,AD),ncol=2)
conditions.id <- c("control","AD")
rownames(mean.expression) <- names(control)
colnames(mean.expression) <- conditions.id
head(mean.expression)

## Scatterplots o gráficos de dispersión para la comparción de distintos 
## genotipos/condiciones. Este tipo de gráficos nos permite obtener una visión 
## global de la comparación entre genotipos/condiciones.
plot(control,AD,xlab="control",ylab="AD")

# Generamos el model.matrix (control vs AD)
design <- model.matrix(~0+sampleInfo$disease_state)
design

#Cambiamos el nombre de columnas
colnames(design) <- c("AD","Control")

summary(exprs(gse138261))

## Nivel medio de expresión
cutoff <- median(exprs(gse138261))

## Miramos si está por encima o por debajo del nivel medio de expresión
is_expressed <- exprs(gse138261) > cutoff

## Nos quedamos con los genes que están en más de 2 muestras

keep <- rowSums(is_expressed) > 2

## En esta tabla se ve cuáles se van a eliminar
table(keep)

## Subconjunto de los elegidos
gse138261 <- gse138261[keep,]

fit <- lmFit(exprs(gse138261), design)
head(fit$coefficients)

#Especificamos un contraste
contrasts <- makeContrasts(AD - Control, levels=design)

fit_1 <- contrasts.fit(fit, contrasts)
fit_1 <- eBayes(fit_1)

topTable(fit_1,coef = 1,sort.by="logFC")

decideTests(fit_1)
table(decideTests(fit_1))

length(fit_1[fit_1$p.value < 1E-5,]$p.value)

select_fit <- fit_1[fit_1$p.value < 1E-5,]$p.value

## Valores atípicos
# arrayWeights asigna una puntuación a cada muestra. 1 mismo peso, < 1 a la baja y > 1 por arriba
aw <- arrayWeights(exprs(gse138261),design)
head(aw)

fit <- lmFit(exprs(gse138261), design,
             weights = aw)
contrasts <- makeContrasts(AD - Control, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
table(decideTests(fit2))

## Enriquecimiento funcional

#Se seleccionan los que tenga p-valor < 1E-5
select_fit <- fit_1[fit_1$p.value < 1E-5,]$p.value

#Columnas de anotacion seleccionadas
final <- fData(eset)[rownames(select_fit),c("ID","GENE_NAME","GENE_SYMBOL","DESCRIPTION","GO_ID","ENSEMBL_ID")]

datos <- final[final$GENE_SYMBOL!= "",]
#Extraemos p-valor
p_value <- select_fit[datos$ID,]
datos$p_value <- p_value
#Ordenamos por p-valor
datos <- datos[order(datos$p_value),]
dim(datos)

