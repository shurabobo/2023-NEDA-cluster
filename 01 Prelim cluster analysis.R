library(tidyverse)
library(ggokabeito)
library(skimr)

# 1 Load data ----
NuCecho_clean <- read_rds("./data/NedaUComet_cleanEcho_v73.rds")

# 2 Select variables ---
## 2.1 Select variables ====
data.frame(colnames(NuCecho_clean))

## 2.2 Choose variables ====
skim(NuCecho_clean)
varChosen <- colnames(NuCecho_clean)[c(5, 9, 45, 46, 52, 53, 56, 62:64, 
                                       85, 96, 124, 109, 130, 126,
                                       22, 32, 13, 19, 40, 41, 35)]
varChosen

## 2.3 Create new dataframe for cluster analysis ====
dfCluster1 <- NuCecho_clean %>% 
  select(varChosen)
skim(dfCluster1)

data.frame(colnames(dfCluster1))

dfCluster2 <- dfCluster1 %>% 
  select(-c(ICULos:Status))

dfCluster2f <- dfCluster2 %>% 
  mutate_if(is.character, as.factor)


# 3 Cluster analysis ------------------------------------------------------

## 3.1 Calculate gower distance (as opposed to Euclidean distance for all numerical variables) ====
library(cluster)
gower_df <- daisy(dfCluster2f,
                  metric = "gower" , stand = T)  # log-transform variable #2
summary(gower_df)

## 3.2 Find the optimal number of clusters using PAM (partitioning aroudn medioids) ====
# Use silhouette width to select the optimal number of clusters

silhouette <- c()
silhouette = c(silhouette, NA)
for(i in 2:10){
  pam_clusters = pam(as.matrix(gower_df),
                     diss = TRUE,
                     k = i)
  silhouette = c(silhouette ,pam_clusters$silinfo$avg.width)
}
plot(1:10, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width")
lines(1:10, silhouette)

## 3.3 perform pam ====
ncluster = 4
pam_neda = pam(gower_df, diss = TRUE, k = ncluster)
pam_neda

# detailed summary results by cluster====
pam_summary <- dfCluster2f %>%
  mutate(cluster = pam_neda$clustering) %>%
  group_by(cluster) %>%
  do(cluster_summary = summary(.))
pam_summary$cluster_summary[[2]]

## 1.5 Visualize results by t-distribute stochastic neighbor embedding (t-SNE) technique ====
# https://learnopencv.com/t-sne-t-distributed-stochastic-neighbor-embedding-explained/
library(Rtsne)  #https://www.r-bloggers.com/2019/05/quick-and-easy-t-sne-analysis-in-r/ ;  https://github.com/jkrijthe/Rtsne
library(ggplot2)
tsne_object <- Rtsne(gower_df, is_distance = TRUE)

tsne_df <- tsne_object$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_neda$clustering))

ggplot(aes(x = X, y = Y), data = tsne_df) +
  geom_point(aes(color = cluster)) +
  scale_colour_okabe_ito()

# 
set.seed(2020)
gower.pam <- pam(gower_df, diss = T, k = ncluster)

gower.mds <- as.data.frame(cmdscale(gower_df, ncluster))
gower.mds$gower_cluster <- as.factor(gower.pam$clustering)

ggplot(gower.mds,aes(x=V1,y=V2,color=gower_cluster)) + 
  geom_point(alpha = 0.6) + theme_ipsum() +
  labs(title="MDS plot",
       subtitle="Colored by PAM cluster") +
  scale_color_okabe_ito()


# compare clusters

clusterTable <- dfCluster1 %>% 
  mutate(cluster = gower.mds$gower_cluster) %>% 
  mutate_if(is.character, as.factor)

table <- tbl_summary(clusterTable, 
                      by = cluster,
                      missing = "ifany") %>% 
  add_n() %>% 
  add_p()
table

# Survival by cluster
library(survival)
library(survminer)
library(ggfortify)

kmCluster <- survfit(Surv(ICULos, ICU_status) ~ cluster, data = clusterTable)


pSurv2 <- ggsurvplot(kmCluster,
                     conf.int = T,
                     #legend.labs = c("Cardiovascular", "Neurological", "Respiratory", "Sepsis"),
                     pval = TRUE,
                     pval.size = 4,
                     pval.method = TRUE,
                     pval.method.size = 4,
                     log.rank.weights = "1")
pSurv2 +  labs(x = "Time (days)")

# ****ML clustering -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cluster)
library(factoextra)
library(gmodels)
library(gtsummary)


load("./data/00_masterLV.rda")

dfDay1 <- masterLV %>% 
  filter(EchoNo == 1) %>% 
  dplyr::select(LVEF, LVEDV, HR, MAP, ACP.imputed, IcuStatus)


df <- dfDay1 %>% 
  drop_na()

dfscale <- scale(dfCluster2f[,-ncol(dfCluster2f)])


distance <- get_dist(dfscale)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


# Compute k-means cluster -------------------------------------------------
## Find optimal number of clusters ====

fviz_nbclust(dfscale, kmeans, method = "wss")
fviz_nbclust(dfscale, kmeans, method = "silhouette")

set.seed(123)
gap_stat <- clusGap(dfscale, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

## 2 clusters ====

k2 <- kmeans(dfscale, centers = 2, nstart = 25)
fviz_cluster(k2, geom = "point", data = dfscale)

dfWcluster1 <- df %>% 
  mutate(cluster = k2$cluster)


### Test mortality within each cluster ====

table1 <- tbl_summary(dfWcluster1, 
                      by = cluster,
                      missing = "ifany") %>% 
  add_n() %>% 
  add_p()
table1

### Plot distribution ====
df %>%
  mutate(cluster = k2$cluster) %>%
  ggplot(aes(cluster, IcuStatus, color = factor(IcuStatus), label = IcuStatus)) +
  geom_jitter()


## 3 clusters -------------------------------------------------------------------------

k3 <- kmeans(dfscale, centers = 3, nstart = 25)
fviz_cluster(k3, geom = "point", data = dfscale)

dfWcluster2 <- df %>% 
  mutate(cluster = k3$cluster)


### Test mortality within each cluster ====

table2 <- tbl_summary(dfWcluster2, 
                      by = cluster,
                      missing = "ifany") %>% 
  add_n() %>% 
  add_p()
table2

### Plot distribution: 3 clusters ====

df %>%
  mutate(cluster = k3$cluster) %>%
  ggplot(aes(cluster, IcuStatus, color = factor(IcuStatus), label = IcuStatus)) +
  geom_jitter()

### check mortality: effects of LVEDV and HR ====

modglm1 <- glm(IcuStatus ~ HR + LVEF + LVEDV, data = df, family = "binomial")
summary(modglm1)
exp(coef(modglm1))

### Pt characteristics according to the clusters

dfListdf <- dfDay1 %>% 
  select(record_id, LVEF, LVEDV, HR, IcuStatus) %>% 
  drop_na()

PtClusterList <- dfListdf$record_id

dfCluster <- masterLV %>% 
  filter(EchoNo == 1) %>%
  filter(record_id %in% PtClusterList) %>% 
  mutate(cluster = k3$cluster)

data.frame(colnames(dfCLuster))

dfClusterFocus <- dfCluster %>% 
  select(record_id, Sex, Age, BMI, LOS_icu, IcuStatus, HosStatus, 
         LVEF, LV_function_vis, LVEDV, LV_size_vis, 
         TAPSE, RV_function_vis, RVAtoLVA_ratio, PSM, 
         IVC_diameter, CVP, highVP,
         ACP.imputed, RVfailure, RVFtapse, cluster)

tbl_summary(dfClusterFocus, 
            by = cluster,
            missing = "ifany") %>% 
  add_n() %>% 
  add_p()

