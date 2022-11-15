install.packages("readr")
install.packages("readxl")
install.packages("tidyr")
install.packages("tibble")
install.packages("stringr")
install.packages("dplyr")
library("stringr")
library("dplyr")


Tible_1 <- read.csv("data_AE_Anastasiia_and_Natalia.csv")
clinical_data <- read.delim("E-MTAB-3732.sdrf.txt")

#Таблица с экспрессией без Настиного гена
expression <- Tible_1[ -c(6,7,8,63), ]


#Таблица по нормальным образцам
clinical_data1 <- clinical_data[!clinical_data$Characteristics.organism.part. != 	
                                  'pancreas', ]
normal <- clinical_data1[!clinical_data1$Characteristics.disease. != 'normal', ]

samples_normal <- select(normal, 1) %>%
  unlist()

samples_normal <- as.numeric(str_extract(samples_normal, "(\\d)+"))

#Таблица по опухолевым образцам
cancer1 <- clinical_data1 %>% filter(Characteristics.disease. %in% c('pancreatic cancer', 'pancreatic adenocarcinoma', 
                                                                     'pancreatic ductal adenocarcinoma', 'pancreatic carcinoma'))


#из рака и из нормы мы достаем номера семплов как вектор, отдельно семплы рак и отдельно семпл норма
#затем фильтруем экспресию по векторам, в результате две таблички с экспрессией только по нужным семплам

