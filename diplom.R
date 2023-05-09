install.packages("readr")
install.packages("readxl")
install.packages("tidyr")
install.packages("tibble")
install.packages("stringr")
install.packages("dplyr")
install.packages("rstatix")
install.packages("tidygraph")
install.packages("ggraph")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("cluster")
install.packages("factoextra")
install.packages("dendextend")
install.packages("tibble")
library("stringr")
library("dplyr")
library("tidyr")
library("rstatix")
library("tidygraph")
library("ggraph")
library("ggplot2")
library("tidyverse")  # data manipulation
library("cluster")    # clustering algorithms
library("factoextra") # clustering visualization
library("dendextend") # for comparing two dendrograms
library("tibble")

Tible_1 <- read.csv("data_AE_Anastasiia_and_Natalia.csv")
clinical_data <- read.delim("E-MTAB-3732.sdrf.txt")

#Таблица с экспрессией без Настиного гена
expression <- Tible_1[ -c(6,7,8,63), ]

#Таблица по нормальным образцам
clinical_data1 <- clinical_data[!clinical_data$Characteristics.organism.part. != 	
                                  'pancreas', ]
normal <- clinical_data1[!clinical_data1$Characteristics.disease. != 'normal', ]

#из таблицы нормы мы достаем номера семплов как вектор, отдельно семпл норма

samples_normal <- select(normal, 1) %>%
  unlist()

samples_normal <- as.numeric(str_extract(samples_normal, "(\\d)+"))

#Таблица по опухолевым образцам
cancer1 <- clinical_data1 %>% filter(Characteristics.disease. %in% c('pancreatic cancer', 'pancreatic adenocarcinoma', 
                                                                     'pancreatic ductal adenocarcinoma', 'pancreatic carcinoma'))

#из таблицы рака мы достаем номера семплов как вектор, отдельно семплы рак

samples_cancer1 <- select(cancer1, 1) %>%
  unlist()

samples_cancer1 <- as.numeric(str_extract(samples_cancer1, "(\\d)+"))

#перевернем таблицу экспрессии чтобы было удобнее с ней работать 

expr_updated <- expression %>%
  select(-1) %>%
  pivot_longer(starts_with("Sample"), names_to = "Samples", values_to = "Expression") 

#найдем медиану экспрессии по поворяющимся генам по всем семплам

new_results <- expr_updated %>%
  group_by(HGNC.symbol, Samples) %>%
  summarise(Expression = median(Expression, na.rm = T)) %>%
  pivot_wider(names_from = HGNC.symbol, values_from = Expression)

#убираем из столбца семлов слово семпл, чтобы отфильтровать по нашим векторам

new_results_num_samples <- new_results %>%
  mutate(Samples = str_remove(Samples, "Sample.")) %>%
  mutate(Samples = as.numeric(Samples))

#затем фильтруем экспресию по векторам, в результате две таблички (экспр по опухол и норм) с экспрессией только по нужным семплам

expr_normal <- new_results_num_samples %>% filter(Samples %in% samples_normal) #таблица с нормальной экспрессией

expr_cancer <- new_results_num_samples %>% filter(Samples %in% samples_cancer1) #таблица с опухолевой экспрессией

#удаляем первый столбик с номерами семплов в обеих таблицах

expr_normal_no_sampl <- expr_normal %>% select(-1)

expr_cancer_no_sampl <- expr_cancer %>% select(-1)

#Подсчет корреляции между всеми генами по двум выборкам отдельно - должны получится 2 корреляционные матрицы (пакет rstatix, cor_mat, Корреляция Спирмана)

cor_normal <- expr_normal_no_sampl %>% cor_mat(method = "spearman")

cor_cancer <- expr_cancer_no_sampl %>% cor_mat(method = "spearman")

#Построить 2 коррелограммы по двум выборкам (поменять цвет)

cor_normal %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot()

cor_cancer %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot()

#Сделать 2 таблицы с информацией о корреляции по парам генов (пакет rstatix, cor_test, Корреляция Спирмана) 

cor_normal_all <-  expr_normal_no_sampl %>% cor_test(method = "spearman")

cor_cancer_all <-  expr_cancer_no_sampl %>% cor_test(method = "spearman")

#отобрать только значимые корреляции p < 0,05

cor_normal_all_1 <- filter(cor_normal_all, p < 0.05)

cor_cancer_all_1 <- filter(cor_cancer_all, p < 0.05)

#убрать р, метод и статитику

cor_normal_all_2 <- cor_normal_all_1  %>% select(-4 & -5 & -6)

cor_cancer_all_2 <- cor_cancer_all_1 %>% select(-4 & -5 & -6)

#строим сеть и визулизируем

normal.graph <- as_tbl_graph(cor_normal_all_2, directed = FALSE)

network_normal <- ggraph(normal.graph) + 
  geom_edge_link(aes(color = cor, width = cor)) + 
  geom_node_point(size = 6) +
  geom_node_text(aes(label = name), size = 6, repel = TRUE) +
  theme_graph() +
  scale_edge_color_gradient2(low = "blue", high = "red", mid = "white")

ggsave(filename = "normal.png", 
       plot = network_normal, scale = 1.8)


cancer.graph <- as_tbl_graph(cor_cancer_all_2, directed = FALSE)

network_cancer <- ggraph(cancer.graph) + 
  geom_edge_link(aes(color = cor, width = cor)) + 
  geom_node_point(size = 6) +
  geom_node_text(aes(label = name), size = 6, repel = TRUE) +
  theme_graph() +
  scale_edge_color_gradient2(low = "blue", high = "red", mid = "white")

ggsave(filename = "cancer.png", 
       plot = network_cancer, scale = 1.8)


#--------------------------------Анализ сетей ко-экспрессии--------------------------

# по норме

normal.graph <- normal.graph %>% activate(edges) %>% mutate (cor_abs = abs(cor))

normal.graph <- normal.graph %>% activate(nodes) %>%
  mutate(hub=centrality_hub(weights = cor_abs, scale = TRUE, options = igraph::arpack_defaults)) %>%
  mutate (betweenness = centrality_betweenness(
    weights = cor_abs ,
    directed = FALSE,
    cutoff = NULL,
    normalized = FALSE
  ))  %>% mutate (group = group_edge_betweenness(weights = cor_abs, directed = FALSE, n_groups = NULL)) 

hub_and_between_normal <- normal.graph %>% activate(nodes) %>% as_tibble()

# по раку

cancer.graph <- cancer.graph %>% activate(edges) %>% mutate (cor_abs = abs(cor))

cancer.graph <- cancer.graph %>% activate(nodes) %>%
  mutate(hub=centrality_hub(weights = cor_abs, scale = TRUE, options = igraph::arpack_defaults)) %>%
  mutate (betweenness = centrality_betweenness(
    weights = cor_abs ,
    directed = FALSE,
    cutoff = NULL,
    normalized = FALSE
  ))  %>% mutate (group = group_edge_betweenness(weights = cor_abs, directed = FALSE, n_groups = NULL)) 

hub_and_between_cancer <- cancer.graph %>% activate(nodes) %>% as_tibble()

#отбираем верхние 5%
#норма
per95_hub_normal = quantile(hub_and_between_normal$hub, 0.95)
per95_between_normal = quantile(hub_and_between_normal$betweenness, 0.95)

table_hub_normal <- filter(hub_and_between_normal, hub >= 0.962361)
table_between_normal <- filter(hub_and_between_normal, betweenness >= 44.575)
table_hubandbetween_normal <- full_join(table_hub_normal, table_between_normal)

#рак
per95_hub_cancer = quantile(hub_and_between_cancer$hub, 0.95)
per95_between_cancer = quantile(hub_and_between_cancer$betweenness, 0.95)

table_hub_cancer <- filter(hub_and_between_cancer, hub >= 0.9922763)
table_between_cancer <- filter(hub_and_between_cancer, betweenness >= 45.48405)
table_hubandbetween_cancer <- full_join(table_hub_cancer, table_between_cancer)

#визуализация сетей
#норма
ko_network_n <- ggraph(normal.graph) +
  geom_edge_link(aes(color = cor, width = cor))  +
                   scale_edge_width(range=c(1,3)) +
  geom_node_point(aes(size =  log(hub)))  +
  geom_node_text(aes(label=ifelse(group == 1, name, NA), size = betweenness*10), repel = TRUE) +
  theme_graph() +
  scale_edge_color_gradient2(low = "blue", high = "red", mid = "white")

ggsave(filename = "normal2.png", 
       plot = ko_network_n, scale = 2)
#рак
ko_network_c <- ggraph(cancer.graph) +
  geom_edge_link(aes(color = cor, width = cor))  +
  scale_edge_width(range=c(1,3)) +
  geom_node_point(aes(size =  log(hub)))  +
  geom_node_text(aes(label=ifelse(group == 1, name, NA), size = betweenness*10), repel = TRUE) +
  theme_graph() +
  scale_edge_color_gradient2(low = "blue", high = "red", mid = "white")

ggsave(filename = "cancer2.png", 
       plot = ko_network_c, scale = 2)

# Построение дендрограмм

cor_normal_1 <- column_to_rownames(cor_normal, var = "rowname") #перевод первого столбца в названия строк, для отображения генов на дендрограмме
cor_cancer_1 <- column_to_rownames(cor_cancer, var = "rowname")

d_normal <- dist(cor_normal_1, method = "euclidean")
hc1_normal <- hclust(d_normal, method = "complete" )
plot(hc1_normal, cex = 0.6, hang = -1)

d_cancer <- dist(cor_cancer_1, method = "euclidean")
hc1_cancer <- hclust(d_cancer, method = "complete" )
plot(hc1_cancer, cex = 0.6, hang = -1)

denrogramma <- tanglegram(hc1_normal, hc1_cancer) #объединили дендограммы

png("дендр.png", width = 12, height = 15, res = 300, units = "cm")
plot(tanglegram ( hc1_normal, hc1_cancer, lab.cex = 0.7))
dev.off()

cor_cophenetic(hc1_normal, hc1_cancer) #сравнили дендрограммы 0.1492133

# Установка пакета BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Установка пакетов для анализа
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("biomaRt")

library(ggplot2)
library(stringr)
library(readr)

library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(biomaRt)

# Здесь должен быть список ваших генов
x <- c(str_split("BRCA1, BRCA2, ATM, ATR, CDK12, CHEK1, CHEK2, EMC2, FANCA, BAP1, BARD1, 
                 BRIP1, AKT1, CTNNB1, ERCC4, ABRAXAS1, FANCD2, FANCE, FANCI, FANCL, KRAS, 
                 MLH1, MRE11, MSH2, MSH6, MUTYH, NBN, PALB2, RAD50, RAD51, RAD51B, RAD51C, 
                 RAD51D, RAD52, RAD54B, RAD54L, PARP1, TP53, TP53BP1, XRCC2, XRCC3, PIK3CA, 
                 PPP2R2A, PTEN", ", ", simplify = T))

mart=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",mart = useMart("ensembl"))

# Нахождение entrezgene_id по списку генов
ids <- unique(getBM(attributes = c("hgnc_symbol", "entrezgene_id"),    
                    filters = "hgnc_symbol",
                    values = x,
                    mart = ensembl))

# Две части скрипта по базе KEGG и GO
# По каждой части будут сохранены таблица и график
# Функция enrichKEGG() может не сработать, а enrichGO() точно должна

edo1 <- enrichKEGG(ids$entrezgene_id, keyType = "kegg", pvalueCutoff = 0.05)
plot1 <- cnetplot(edo1, circular = T, color.params = list(edge = 0.05))
edo_result1 <- edo1@result
write_csv2(edo_result1, "KEGG_result.csv")
ggsave(plot = plot1, filename = "go_concept_plot_KEGG.png", dpi = 600, scale = 2)

edo2 <- enrichGO(ids$entrezgene_id, 'org.Hs.eg.db', keyType = "ENTREZID", pvalueCutoff = 0.05, readable=T)
plot2 <- cnetplot(edo2, circular = T, colorEdge = T)
edo_result2 <- edo2@result
write_csv2(edo_result2, "GO_result.csv")
ggsave(plot = plot2, filename = "go_concept_plot_GO.png", dpi = 600, scale = 2)
