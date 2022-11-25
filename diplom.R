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
library("stringr")
library("dplyr")
library("tidyr")
library("rstatix")
library("tidygraph")
library("ggraph")
library("ggplot2")



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
       plot = network_normal, scale = 1.8)

