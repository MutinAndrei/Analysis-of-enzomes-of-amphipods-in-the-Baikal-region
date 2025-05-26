#Функция которая на вход получает 3 переменные 
x = EConlynMATRIX4numbers_EV1forvoolcano <- read.csv("30-01-2024_EConlynMATRIX4numbers_EV1.csv")
y = (Eve; Gla; Ecy)
z = b <- read.csv("T1_5andB12A.csv", row.names=1)
nameConrol = #Названия из матрицы и нижнее подчеркивание EveB12_ которое будет выступать контролем
nameExperiment = #Названия из матрицы и нижнее подчеркивание EveЕ1.5_ которое будет выступать опытом

#функция для подсчета дифф экспрессии 
dif_express_function <- function(DEC, x, y, z, nameConrol, nameExperiment) {
  if (DEC == "DECRE") {
    EConlynMATRIX4numbers_EV1forvoolcano <- read.csv(x)
    EConlynMATRIX4numbers_EV1forvoolcano <- EConlynMATRIX4numbers_EV1forvoolcano[ , c(-1, -3)]
    colnames(EConlynMATRIX4numbers_EV1forvoolcano)[1] <- "gen"
    library(dplyr)
    a2 <- EConlynMATRIX4numbers_EV1forvoolcano
    ##удаление NANfull
    a2[is.na(a2)] <- 0
    percisoformfull_NANRMa <- a2
    percisoformfull_NANRMaEXPRESSION <- percisoformfull_NANRMa
    #Преобразует таблицу в матрицу для dds
    library(tidyverse)
    afull <- percisoformfull_NANRMaEXPRESSION %>% 
      group_by(Temperature, gen) %>% 
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = Expression)
    write.csv(afull, "Expression.csv")
    afull <- afull[ , -3]
    afull2 <- afull %>% select(contains(c("gen", "EC.number", nameConrol, nameExperiment)))
    afull4 <- na.omit(afull2)
    afull4[,-c(1:2)] <- round(afull4[,-c(1:2)])
    #Сохраняем готовую матрицу
    write.csv(afull4, "afull4_1_kB12Ae.csv")
    #Подгружаем готовую матрицу 
    afull4 <- read.csv("afull4_1_kB12Ae.csv")
    afull5 <- afull4[ , -1]
    #пакеты 
    library(dplyr)
    library(tidyverse)
    #Меняем название
    afull5 <- afull5[!duplicated(afull5), ]
    afull6 <- afull5 %>% remove_rownames %>% column_to_rownames(var="gen")
    afull6 <- afull6[, -1]
    ##нумеруем
    afull4 <- afull4[ , -1]
    afull4[, -c(1:2)] <- lapply(afull4[, -c(1:2)], as.numeric)
    #Подгружаем столбцы с названием и распределением
    b <- read.csv(z, row.names=1)
    afull61 <- na.omit(afull6)
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = afull61,
                                  colData = b,
                                  design = ~ dex)
    
    dds <- DESeq(dds)
    res <- results(dds)
    head(results(dds, tidy=TRUE))
    summary(res)
    res <- res[order(res$padj),]
    head(res)
    par(mfrow=c(1,1))
    ddsdata <- (results(dds, tidy=TRUE))
    sum(is.na(ddsdata))
    ddsdata <- na.omit(ddsdata)
  } else if (DEC == "DECPE") {
    EConlynMATRIX4numbers_EV1forvoolcano <- read.csv(x)
    EConlynMATRIX4numbers_EV1forvoolcano <- EConlynMATRIX4numbers_EV1forvoolcano[ , c(-1, -3)]
    colnames(EConlynMATRIX4numbers_EV1forvoolcano)[1] <- "gen"
    library(dplyr)
    a1 <- EConlynMATRIX4numbers_EV1forvoolcano
    ##Процентный состав для всей выборки (каждая повторность)
    a2 <- a1 %>%
      group_by(Temperature, EC.number) %>%
      mutate(percent = Expression/sum(Expression) *100)
    ##запись файла с процентным вкладом
    write.csv(a2, "percent_table.csv")
    ##удаление NANfull
    a2[is.na(a2)] <- 0
    percisoformfull_NANRMa <- a2
    percisoformfull_NANRMaEXPRESSION <- percisoformfull_NANRMa[ , -4]
    #Преобразует таблицу в матрицу для dds
    library(tidyverse)
    afull <- percisoformfull_NANRMaEXPRESSION %>% 
      group_by(Temperature, gen) %>% 
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = percent)
    write.csv(afull, "Expression.csv")
    afull <- afull[ , -3]
    afull2 <- afull %>% select(contains(c("gen", "EC.number", nameConrol, nameExperiment)))
    afull4 <- na.omit(afull2)
    afull4[,-c(1:2)] <- round(afull4[,-c(1:2)])
    #Сохраняем готовую матрицу
    write.csv(afull4, "afull4_1_kB12Ae.csv")
    #Подгружаем готовую матрицу 
    afull4 <- read.csv("afull4_1_kB12Ae.csv")
    afull5 <- afull4[ , -1]
    #пакеты 
    library(dplyr)
    library(tidyverse)
    #Меняем название
    afull5 <- afull5[!duplicated(afull5), ]
    afull6 <- afull5 %>% remove_rownames %>% column_to_rownames(var="gen")
    afull6 <- afull6[, -1]
    ##нумеруем
    afull4 <- afull4[ , -1]
    afull4[, -c(1:2)] <- lapply(afull4[, -c(1:2)], as.numeric)
    #Подгружаем столбцы с названием и распределением
    b <- read.csv(z, row.names=1)
    afull61 <- na.omit(afull6)
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = afull61,
                                  colData = b,
                                  design = ~ dex)
    
    dds <- DESeq(dds)
    res <- results(dds)
    head(results(dds, tidy=TRUE))
    summary(res)
    res <- res[order(res$padj),]
    head(res)
    par(mfrow=c(1,1))
    ddsdata <- (results(dds, tidy=TRUE))
    sum(is.na(ddsdata))
    ddsdata <- na.omit(ddsdata)
  } else if (DEC == "DECSE") {
    EConlynMATRIX4numbers_EV1forvoolcano <- read.csv(x)
    EConlynMATRIX4numbers_EV1forvoolcano <- EConlynMATRIX4numbers_EV1forvoolcano[ , c(-1, -3)]
    colnames(EConlynMATRIX4numbers_EV1forvoolcano)[1] <- "gen"
    library(dplyr)
    a1 <- EConlynMATRIX4numbers_EV1forvoolcano
    a1 <- a1[ , -1]
    ##Процентный состав для всей выборки (каждая повторность)
    a2 <- a1 %>%
      group_by(Temperature, EC.number) %>%
      mutate(sumECC = sum(Expression))
    ##удаление NANfull
    a2[is.na(a2)] <- 0
    percisoformfull_NANRMa <- a2
    percisoformfull_NANRMaEXPRESSION <- percisoformfull_NANRMa[ ,-3]
    #Преобразует таблицу в матрицу для dds
    library(tidyverse)
    afull <- percisoformfull_NANRMaEXPRESSION %>% 
      group_by(Temperature, EC.number) %>% 
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = sumECC)
    write.csv(afull, "Expression.csv")
    afull <- subset(afull, select = -c(id))
    afull2 <- afull %>% select(contains(c("EC.number", nameConrol, nameExperiment)))
    afull4 <- na.omit(afull2)
    afull4[,-c(1)] <- round(afull4[,-c(1)])
    #Сохраняем готовую матрицу
    write.csv(afull4, "afull4_1_kB12Ae.csv")
    #Подгружаем готовую матрицу 
    afull4 <- read.csv("afull4_1_kB12Ae.csv")
    afull5 <- afull4[ , -1]
    #пакеты 
    library(dplyr)
    library(tidyverse)
    #Меняем название
    afull5 <- afull5[!duplicated(afull5), ]
    afull6 <- afull5 %>% remove_rownames %>% column_to_rownames(var="EC.number")
    ##нумеруем
    afull4 <- afull4[ , -1]
    afull4[, -c(1)] <- lapply(afull4[, -c(1)], as.numeric)
    #Подгружаем столбцы с названием и распределением
    b <- read.csv(z, row.names=1)
    b[, 1] <- as.factor(b[ , 1])
    afull61 <- na.omit(afull6)
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = afull61,
                                  colData = b,
                                  design = ~ dex)
    
    dds <- DESeq(dds)
    res <- results(dds)
    head(results(dds, tidy=TRUE))
    summary(res)
    res <- res[order(res$padj),]
    head(res)
    par(mfrow=c(1,1))
    ddsdata <- (results(dds, tidy=TRUE))
    sum(is.na(ddsdata))
    ddsdata <- na.omit(ddsdata)
  } else {print("ancorected data (DEC)")}
}

#Функция для фильтрации таблицы по logFoldChange и p-Value

ddsdata = #Таблица после deseq2
  logFoldchange = #число для логфолдчанже
  gen_type = #Down or Up
  pvalue = #P-value, entered in digits (example: pvalue = 0.001) Already adjusted for multiple comparisons 
  table_for_transcript = #Table with expression data with EC number "30-01-2024_EConlynMATRIX4numbers_EV1.csv"
  
  
  
filterlfpv <- function(ddsdata, logFoldchange, gen_type, pvalue, table_for_transcript, type_analisis) {
    ddsdata <- ddsdata
    if (type_analisis == "DECRE" | type_analisis == "DECPE") {
      DownUp <- if (gen_type == "Down" | gen_type == "DOWN" | gen_type == "DoWn" | gen_type == "down") {
        ddsdata[ddsdata$padj <= 0.001 & ddsdata$log2FoldChange <= -1, ]
      } else {ddsdata[ddsdata$padj <= 0.001 & ddsdata$log2FoldChange >= 1, ]}
      EConlynMATRIX4numbers_EV11 <- read.csv(table_for_transcript)
      EConlynMATRIX4numbers_EV11 <- EConlynMATRIX4numbers_EV11[ , -c(1)]
      colnames(EConlynMATRIX4numbers_EV11)[1] <- 'row'
      EConlynMATRIX4numbers_EV11 <- EConlynMATRIX4numbers_EV11[ , -c(4, 5)]
      EConlynMATRIX4numbers_EV11 <- EConlynMATRIX4numbers_EV11[!duplicated(EConlynMATRIX4numbers_EV11), ]
      Up_or_Down <- merge(EConlynMATRIX4numbers_EV11, DownUp)
    } else if (type_analisis == "DECSE") {
      DownUp <- if (gen_type == "Down" | gen_type == "DOWN" | gen_type == "DoWn" | gen_type == "down") {
        ddsdata[ddsdata$padj <= 0.001 & ddsdata$log2FoldChange <= -1, ]
      } else {ddsdata[ddsdata$padj <= 0.001 & ddsdata$log2FoldChange >= 1, ]}
      Up_or_Down <- DownUp
    }
  }



###функция считающая уникальные EC номера в табличках####
#df - это таблица созданная из всех значимо понижающихся или повышающихся генов
#gen_type - тип гена "Down" или "Up"

UniqECnumber <- function(df, gen_type, type_analisis) {
  if(type_analisis == "DECRE" | type_analisis == "DECPE") {
    ifelse(gen_type == "Down", df$gens_type <- "Down", df$gens_type <- "Up")
    alltable <- df[ , c("row","Prot.name.Uniprot.","EC.number", "gens_type")]
    alltable1 <- alltable[!duplicated(alltable), ]
    uniqECnumbers <- unique(alltable1$EC.number)
  } else if (type_analisis == "DECSE") {
    ifelse(gen_type == "Down", df$gens_type <- "Down", df$gens_type <- "Up")
    alltable <- df[ , c("row", "gens_type")]
    alltable1 <- alltable[!duplicated(alltable), ]
    uniqECnumbers <- unique(alltable1$row)
  } else {
    print("ancorected type analisis")
  }
}


##функциия для создания фасты коррелирующих пар3###

fastacreate_FUN <- function(Data_table, patern_column, type_analisis, name_fasta) {
  if (type_analisis == "Expression" | type_analisis == "E") {
    Evpersent <- read.csv(Data_table)
    Evpersent <- Evpersent[ , -c(1, 6) ]
    
    
    Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
    library(tidyverse)
    Evpersent2 <- Evpersent %>%
      group_by(Temperature, gen) %>%
      mutate(meanExpression = median(Expression, na.rm = TRUE))
    
    
    Evpersent2 <- Evpersent2[ , -4]
    Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]
    
    Evpersent4 <- Evpersent3 %>% 
      group_by(Temperature, gen) %>%
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = meanExpression)
    Evpersent4 <- Evpersent4[ , -3]
    Evpersent4 <- Evpersent4[ , -c(6, 4, 5)]
    
    Evpersent4 <- Evpersent4[ , patern_column]
    
    Evpersent4 <- Evpersent4[duplicated(Evpersent4[,2]) | duplicated(Evpersent4[,2], fromLast=TRUE),]
    Evpersent4$gen <- sub("\\|$", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\gb|", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\|", "", Evpersent4$gen)
    Evpersent5 <- Evpersent4
    
    Evpersent6 <- split(Evpersent5, Evpersent5$EC.number)
    
    
    a <- lapply(Evpersent6[1:length(Evpersent6)], function(x) data.frame(t(x[-c(1, 2)])))
    
    test <- lapply(Evpersent6[c(TRUE)], `[[`, 1)
    colnames(a[[1]]) <- test[[1]]
    
    for (i in 1:length(a)) {
      colnames(a[[i]]) <- test[[i]]
    }
    
    
    
    corr <- lapply(a[1:length(a)], function(x) cor(x, use = "complete.obs"))
    diagNA <- lapply(corr[1:length(corr)], function(x) `diag<-`(x, NA))
    miiiin <- lapply(diagNA[1:length(diagNA)], function(x) min(x, na.rm = TRUE))
    maaaax <- lapply(diagNA[1:length(diagNA)], function(x) max(x, na.rm = TRUE))
    sumsCOUNT <- lapply(diagNA[1:length(diagNA)], function(x) rowSums(x <= -0.9, na.rm = TRUE))
    sumsCOUNT_2 <- lapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x))
    
    
    
    dfmin <- as.data.frame(miiiin)
    dfmax <- as.data.frame(maaaax)
    sumsCOUNT_df <- as.data.frame(sumsCOUNT_2)
    dfmin_1 <- data.frame(t(dfmin[-1]))
    dfmax_1 <- data.frame(t(dfmax[-1]))
    sumsCOUNT_df <- data.frame(t(sumsCOUNT_df[-1]))
    write.csv(dfmin_1, "dfmin_1.csv")
    write.csv(dfmax_1, "dfmax_1.csv")
    write.csv(sumsCOUNT_df, "sumsCOUNT_df.csv")
    df_min_cor <- read.csv("dfmin_1.csv")
    df_max_cor <- read.csv("dfmin_1.csv")
    sumsCOUNT_df <- read.csv("sumsCOUNT_df.csv")
    colnames(df_min_cor)[2] <- "min_cor"
    colnames(df_max_cor)[2] <- "max_cor"
    colnames(sumsCOUNT_df)[2] <- "sum_cor"
    
    sumsCOUNT_df_2 <- sumsCOUNT_df[apply(sumsCOUNT_df!=0, 1, all),]
    sumsCOUNT_df_2$X <- gsub("X","",as.character(sumsCOUNT_df_2$X))
    sumsCOUNT_df_2$sum_cor <- sumsCOUNT_df_2$sum_cor / 2
    colnames(sumsCOUNT_df_2)[1] <- "Prot.name(Uniprot)"
    colnames(sumsCOUNT_df_2)[1] <- "EC.number"
    write.csv(sumsCOUNT_df_2, "sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- read.csv("sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- sumsCOUNT_df_2[ , -1]
    library(brendaDb)
    df <- ReadBrenda("brenda_2023_1.txt")
    ID <- df[df$field == "RECOMMENDED_NAME",]
    ID <- ID[ , -c(2)]
    colnames(ID)[1] <- "EC.number"
    sum_coun_w_ID <- merge(sumsCOUNT_df_2, ID)
    sum(sum_coun_w_ID$sum_cor == 1)
    write.csv(sum_coun_w_ID, "prot_name_w_count_corr_EV.csv")
    df222 <- subset(sum_coun_w_ID, sum_coun_w_ID$sum_cor == 1)
    write.csv(df222, "Table_of_only_par_for_EV.csv")
    
    
    test_221 <- sumsCOUNT[sapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x) == 2)]
    test_222 <- lapply(1:length(test_221), function(x) names(test_221[[x]][test_221[[x]] == 1]))
    test_223 <- unlist(test_222)
    test_224 <- as.data.frame(test_223)
    gen_EC <- Evpersent[ , -c(3,4)]
    gen_EC_2 <- gen_EC[!duplicated(gen_EC), ]
    colnames(test_224)[1] <- "gen"
    gen_EC_and_par <- merge(test_224, gen_EC_2)
    gen_EC_and_par$gen <- gsub('gb', '', gen_EC_and_par$gen)
    write.csv(gen_EC_and_par, "gen_EC_and_par_EV_exp_median.csv")
    gen_EC_and_par <- read.csv("gen_EC_and_par_EV_exp_median.csv")
    gen_EC_and_par <- gen_EC_and_par[ , -1]
    library(reutils)
    f <- efetch(uid = test_224$gen, db = "nucleotide", rettype = "fasta", retmode = "text")
    write(content(f), name_fasta)
  } else if (type_analisis == "Percent" | type_analisis == "P") {
    
    Evpersent <- read.csv(Data_table)
    Evpersent <- Evpersent[ , -c(1, 5) ]
    
    Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
    library(tidyverse)
    Evpersent2 <- Evpersent %>%
      group_by(Temperature, gen) %>%
      mutate(medianPercent = median(percent, na.rm = TRUE))
    
    
    Evpersent2 <- Evpersent2[ , -4]
    Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]
    
    Evpersent4 <- Evpersent3 %>% 
      group_by(Temperature, gen) %>%
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = medianPercent)
    Evpersent4 <- Evpersent4[ , -3]
    Evpersent4 <- Evpersent4[ , -c(6, 4, 5)]
    
    Evpersent4 <- Evpersent4[ , patern_column]
    
    Evpersent4 <- Evpersent4[duplicated(Evpersent4[,2]) | duplicated(Evpersent4[,2], fromLast=TRUE),]
    Evpersent4$gen <- sub("\\|$", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\gb|", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\|", "", Evpersent4$gen)
    Evpersent5 <- Evpersent4
    
    Evpersent6 <- split(Evpersent5, Evpersent5$EC.number)
    
    
    a <- lapply(Evpersent6[1:length(Evpersent6)], function(x) data.frame(t(x[-c(1, 2)])))
    
    test <- lapply(Evpersent6[c(TRUE)], `[[`, 1)
    colnames(a[[1]]) <- test[[1]]
    
    for (i in 1:length(a)) {
      colnames(a[[i]]) <- test[[i]]
    }
    
    
    
    corr <- lapply(a[1:length(a)], function(x) cor(x))
    diagNA <- lapply(corr[1:length(corr)], function(x) `diag<-`(x, NA))
    miiiin <- lapply(diagNA[1:length(diagNA)], function(x) min(x, na.rm = TRUE))
    maaaax <- lapply(diagNA[1:length(diagNA)], function(x) max(x, na.rm = TRUE))
    sumsCOUNT <- lapply(diagNA[1:length(diagNA)], function(x) rowSums(x <= -0.9, na.rm = TRUE))
    sumsCOUNT_2 <- lapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x))
    
    
    
    dfmin <- as.data.frame(miiiin)
    dfmax <- as.data.frame(maaaax)
    sumsCOUNT_df <- as.data.frame(sumsCOUNT_2)
    dfmin_1 <- data.frame(t(dfmin[-1]))
    dfmax_1 <- data.frame(t(dfmax[-1]))
    sumsCOUNT_df <- data.frame(t(sumsCOUNT_df[-1]))
    write.csv(dfmin_1, "dfmin_1.csv")
    write.csv(dfmax_1, "dfmax_1.csv")
    write.csv(sumsCOUNT_df, "sumsCOUNT_df.csv")
    df_min_cor <- read.csv("dfmin_1.csv")
    df_max_cor <- read.csv("dfmin_1.csv")
    sumsCOUNT_df <- read.csv("sumsCOUNT_df.csv")
    colnames(df_min_cor)[2] <- "min_cor"
    colnames(df_max_cor)[2] <- "max_cor"
    colnames(sumsCOUNT_df)[2] <- "sum_cor"
    
    sumsCOUNT_df_2 <- sumsCOUNT_df[apply(sumsCOUNT_df!=0, 1, all),]
    sumsCOUNT_df_2$X <- gsub("X","",as.character(sumsCOUNT_df_2$X))
    sumsCOUNT_df_2$sum_cor <- sumsCOUNT_df_2$sum_cor / 2
    colnames(sumsCOUNT_df_2)[1] <- "Prot.name(Uniprot)"
    colnames(sumsCOUNT_df_2)[1] <- "EC.number"
    write.csv(sumsCOUNT_df_2, "sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- read.csv("sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- sumsCOUNT_df_2[ , -1]
    library(brendaDb)
    df <- ReadBrenda("brenda_2023_1.txt")
    ID <- df[df$field == "RECOMMENDED_NAME",]
    ID <- ID[ , -c(2)]
    colnames(ID)[1] <- "EC.number"
    sum_coun_w_ID <- merge(sumsCOUNT_df_2, ID)
    sum(sum_coun_w_ID$sum_cor == 1)
    write.csv(sum_coun_w_ID, "prot_name_w_count_corr_EV.csv")
    df222 <- subset(sum_coun_w_ID, sum_coun_w_ID$sum_cor == 1)
    write.csv(df222, "Table_of_only_par_for_EV.csv")
    
    
    test_221 <- sumsCOUNT[sapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x) == 2)]
    test_222 <- lapply(1:length(test_221), function(x) names(test_221[[x]][test_221[[x]] == 1]))
    test_223 <- unlist(test_222)
    test_224 <- as.data.frame(test_223)
    gen_EC <- Evpersent[ , -c(3,4)]
    gen_EC_2 <- gen_EC[!duplicated(gen_EC), ]
    colnames(test_224)[1] <- "gen"
    gen_EC_and_par <- merge(test_224, gen_EC_2)
    gen_EC_and_par$gen <- gsub('gb', '', gen_EC_and_par$gen)
    write.csv(gen_EC_and_par, "gen_EC_and_par_EV_per.csv")
    gen_EC_and_par <- read.csv("gen_EC_and_par_EV_per.csv")
    gen_EC_and_par <- gen_EC_and_par[ , -1]
    library(reutils)
    f1 <- efetch(uid = test_224[1:500 ,], db = "nucleotide", rettype = "fasta", retmode = "text")
    f2 <- efetch(uid = test_224[501:nrow(test_224), ], db = "nucleotide", rettype = "fasta", retmode = "text")
    f1_1 <- c(content(f1),  content(f2))
    write(f1_1, name_fasta)
  }
  }


pTT_EV <- c("gen", "EC.number", "EveT1.5", "EveB12", "EveB12A", "EveT12", "EveT18", "EveT24")
pTT_EC <- c("gen", "EC.number", "EcyT1.5", "EcyB12", "EcyB12A", "EcyT12", "EcyT18", "EcyT24")
pTT_GL <- c("gen", "EC.number", "GlaT1.5", "GlaB12", "GlaB12A", "GlaT12", "GlaT18", "GlaT24")

fastacreate_FUN(Data_table = "percent_tableEV.csv", patern_column = pTT_EV, type_analisis = "E", name_fasta = "EV_exp_17_01_2025.fa")
fastacreate_FUN(Data_table = "percent_tableEV.csv", patern_column = pTT_EV, type_analisis = "P", name_fasta = "EV_per_17_01_2025.fa")
fastacreate_FUN(Data_table = "percent_tableEC.csv", patern_column = pTT_EC, type_analisis = "E", name_fasta = "EC_exp_17_01_2025.fa")
fastacreate_FUN(Data_table = "percent_tableEC.csv", patern_column = pTT_EC, type_analisis = "P", name_fasta = "EC_per_17_01_2025.fa")
fastacreate_FUN(Data_table = "percent_tableGL.csv", patern_column = pTT_GL, type_analisis = "E", name_fasta = "GL_exp_17_01_2025.fa")
fastacreate_FUN(Data_table = "percent_tableGL.csv", patern_column = pTT_GL, type_analisis = "P", name_fasta = "GL_per_17_01_2025.fa")


#### Полученные файлы .fasta мы обрабатывали с помощью transdecoder in BASH
#EV_exp
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t EV_exp_median.fa --complete_orfs_only
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t EV_exp_median.fa --single_best_only
#EC_exp
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t EC_exp_median.fa --complete_orfs_only
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t EC_exp_median.fa --single_best_only
#GL_exp
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t GL_exp_median.fa --complete_orfs_only
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t GL_exp_median.fa --single_best_only
#EV_per
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t EV_per_median.fa --complete_orfs_only
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t EV_per_median.fa --single_best_only
#EC_per
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t EC_per_median.fa --complete_orfs_only
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t EC_per_median.fa --single_best_only
#GL_per
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t GL_per_median.fa --complete_orfs_only
#/media/secondary/apps/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t GL_per_median.fa --single_best_only




####функция для подсчета корреляций###

corrFUN <- function(table_sp, pattern_Sp, trancdecoder_sp, analis, Pairs_or_not, cor_lvl) {
  if (analis == "E") {
    Evpersent <- read.csv(table_sp)
    Evpersent <- Evpersent[ , -c(1, 6) ]
    
    Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
    library(tidyverse)
    
    
    Evpersent2 <- Evpersent %>%
      group_by(Temperature, gen) %>%
      mutate(meanExpression = median(Expression, na.rm = TRUE))
    
    
    Evpersent2 <- Evpersent2[ , -4]
    Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]
    
    Evpersent4 <- Evpersent3 %>% 
      group_by(Temperature, gen) %>%
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = meanExpression)
    Evpersent4 <- Evpersent4[ , -3]
    Evpersent4 <- Evpersent4[ , -c(6, 4, 5)]
    
    Evpersent4 <- Evpersent4[ , pattern_Sp]
    
    Evpersent5 <- Evpersent4[duplicated(Evpersent4[,2]) | duplicated(Evpersent4[,2], fromLast=TRUE),]
    Evpersent5$gen <- sub("\\|$", "", Evpersent5$gen)
    Evpersent5$gen <- sub("\\gb|", "", Evpersent5$gen)
    Evpersent5$gen <- sub("\\|", "", Evpersent5$gen)
    Evpersent5 <- Evpersent5
    
    Evpersent6 <- split(Evpersent5, Evpersent5$EC.number)
    
    
    a <- lapply(Evpersent6[1:length(Evpersent6)], function(x) data.frame(t(x[-c(1, 2)])))
    
    test <- lapply(Evpersent6[c(TRUE)], `[[`, 1)
    
    colnames(a[[1]]) <- test[[1]]
    
    for (i in 1:length(a)) {
      colnames(a[[i]]) <- test[[i]]
    }
    
    
    
    corr <- lapply(a[1:length(a)], function(x) cor(x, use = "complete.obs"))
    diagNA <- lapply(corr[1:length(a)], function(x) `diag<-`(x, NA))
    miiiin <- lapply(diagNA[1:length(a)], function(x) min(x, na.rm = TRUE))
    maaaax <- lapply(diagNA[1:length(a)], function(x) max(x, na.rm = TRUE))
    sumsCOUNT <- lapply(diagNA[1:length(a)], function(x) rowSums(x <= cor_lvl, na.rm = TRUE))
    sumsCOUNT_2 <- lapply(sumsCOUNT[1:length(a)], function(x) sum(x))
    
    
    
    dfmin <- as.data.frame(miiiin)
    dfmax <- as.data.frame(maaaax)
    sumsCOUNT_df <- as.data.frame(sumsCOUNT_2)
    dfmin_1 <- data.frame(t(dfmin[-1]))
    dfmax_1 <- data.frame(t(dfmax[-1]))
    sumsCOUNT_df <- data.frame(t(sumsCOUNT_df[-1]))
    write.csv(dfmin_1, "dfmin_1.csv")
    write.csv(dfmax_1, "dfmax_1.csv")
    write.csv(sumsCOUNT_df, "sumsCOUNT_df.csv")
    df_min_cor <- read.csv("dfmin_1.csv")
    df_max_cor <- read.csv("dfmin_1.csv")
    sumsCOUNT_df <- read.csv("sumsCOUNT_df.csv")
    colnames(df_min_cor)[2] <- "min_cor"
    colnames(df_max_cor)[2] <- "max_cor"
    colnames(sumsCOUNT_df)[2] <- "sum_cor"
    
    sumsCOUNT_df_2 <- sumsCOUNT_df[apply(sumsCOUNT_df!=0, 1, all),]
    sumsCOUNT_df_2$X <- gsub("X","",as.character(sumsCOUNT_df_2$X))
    sumsCOUNT_df_2$sum_cor <- sumsCOUNT_df_2$sum_cor / 2
    colnames(sumsCOUNT_df_2)[1] <- "Prot.name(Uniprot)"
    colnames(sumsCOUNT_df_2)[1] <- "EC.number"
    write.csv(sumsCOUNT_df_2, "sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- read.csv("sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- sumsCOUNT_df_2[ , -1]
    library(brendaDb)
    df <- ReadBrenda("brenda_2023_1.txt")
    ID <- df[df$field == "RECOMMENDED_NAME",]
    ID <- ID[ , -c(2)]
    colnames(ID)[1] <- "EC.number"
    sum_coun_w_ID <- merge(sumsCOUNT_df_2, ID)
    sum(sum_coun_w_ID$sum_cor == 1)
    write.csv(sum_coun_w_ID, "prot_name_w_count_corr_EV.csv")
    df222 <- subset(sum_coun_w_ID, sum_coun_w_ID$sum_cor == 1)
    write.csv(df222, "Table_of_only_par_for_EV.csv")
    
    if (Pairs_or_not == "Pairs") {
      test_221 <- sumsCOUNT[sapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x) == 2)]
      test_222 <- lapply(1:length(test_221), function(x) names(test_221[[x]][test_221[[x]] == 1]))
      test_223 <- unlist(test_222)
      test_224 <- as.data.frame(test_223)
      gen_EC <- Evpersent[ , -c(3,4)]
      gen_EC_2 <- gen_EC[!duplicated(gen_EC), ]
      colnames(test_224)[1] <- "gen"
      gen_EC_2$gen <- sub("\\|$", "", gen_EC_2$gen)
      gen_EC_2$gen <- sub("\\gb|", "", gen_EC_2$gen)
      gen_EC_2$gen <- sub("\\|", "", gen_EC_2$gen)
      gen_EC_and_par <- merge(test_224, gen_EC_2)
      gen_EC_and_par$gen <- gsub('gb', '', gen_EC_and_par$gen)
      write.csv(gen_EC_and_par, "gen_EC_and_par_EV_exp_median.csv")
      gen_EC_and_par <- read.csv("gen_EC_and_par_EV_exp_median.csv")
      gen_EC_and_par <- gen_EC_and_par[ , -1]
      library(Biostrings)
      Gl_amino_fastas_after_tranc <- readAAStringSet(trancdecoder_sp)
      Gl_amino_fastas_after_tranc <- as.data.frame(Gl_amino_fastas_after_tranc)
      write.csv(Gl_amino_fastas_after_tranc, "EV_amino_fastas_after_tranc.csv")
      Gl_amino_fastas_after_tranc <- read.csv("EV_amino_fastas_after_tranc.csv")
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p1 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p2 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p3 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p4 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p5 .*","",Gl_amino_fastas_after_tranc$X)
      colnames(Gl_amino_fastas_after_tranc)[1] <- "gen"
      colnames(Gl_amino_fastas_after_tranc)[2] <- "seq_amino"
      EC_gen_seq_GL <- merge(gen_EC_and_par, Gl_amino_fastas_after_tranc)
      df222_222 <- df222[ , -2]
      EC_gen_seq_GL_prot_name <- merge(EC_gen_seq_GL, df222_222)
      aaa <- EC_gen_seq_GL_prot_name[duplicated(EC_gen_seq_GL_prot_name[,1]) | duplicated(EC_gen_seq_GL_prot_name[,1], fromLast=TRUE),]
    } else if (Pairs_or_not == "Unpairs") {
      test_221 <- sumsCOUNT[sapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x) != 0)]
      test_222 <- lapply(1:length(test_221), function(x) names(test_221[[x]][test_221[[x]] != 0]))
          }
    
      } else if (analis == "P") {
    Evpersent <- read.csv(table_sp)
    Evpersent <- Evpersent[ , -c(1, 5) ]
    
    
    Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
    library(tidyverse)
    Evpersent2 <- Evpersent %>%
      group_by(Temperature, gen) %>%
      mutate(medianPercent = median(percent, na.rm = TRUE))
    
    
    Evpersent2 <- Evpersent2[ , -4]
    Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]
    
    Evpersent4 <- Evpersent3 %>% 
      group_by(Temperature, gen) %>%
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = medianPercent)
    Evpersent4 <- Evpersent4[ , -3]
    Evpersent4 <- Evpersent4[ , -c(6, 4, 5)]
    
    Evpersent4 <- Evpersent4[ , c(pattern_Sp)]
  
    Evpersent5 <- Evpersent4[duplicated(Evpersent4[,2]) | duplicated(Evpersent4[,2], fromLast=TRUE),]
    Evpersent5$gen <- sub("\\|$", "", Evpersent5$gen)
    Evpersent5$gen <- sub("\\gb|", "", Evpersent5$gen)
    Evpersent5$gen <- sub("\\|", "", Evpersent5$gen)
    Evpersent5 <- Evpersent5
    
    Evpersent6 <- split(Evpersent5, Evpersent5$EC.number)
    
    
    a <- lapply(Evpersent6[1:length(Evpersent6)], function(x) data.frame(t(x[-c(1, 2)])))
    
    test <- lapply(Evpersent6[c(TRUE)], `[[`, 1)
    colnames(a[[1]]) <- test[[1]]
    
    for (i in 1:length(a)) {
      colnames(a[[i]]) <- test[[i]]
    }
    
    
    
    corr <- lapply(a[1:length(a)], function(x) cor(x))
    diagNA <- lapply(corr[1:length(a)], function(x) `diag<-`(x, NA))
    miiiin <- lapply(diagNA[1:length(a)], function(x) min(x, na.rm = TRUE))
    maaaax <- lapply(diagNA[1:length(a)], function(x) max(x, na.rm = TRUE))
    sumsCOUNT <- lapply(diagNA[1:length(a)], function(x) rowSums(x <= cor_lvl, na.rm = TRUE))
    sumsCOUNT_2 <- lapply(sumsCOUNT[1:length(a)], function(x) sum(x))
    
    
    
    dfmin <- as.data.frame(miiiin)
    dfmax <- as.data.frame(maaaax)
    sumsCOUNT_df <- as.data.frame(sumsCOUNT_2)
    dfmin_1 <- data.frame(t(dfmin[-1]))
    dfmax_1 <- data.frame(t(dfmax[-1]))
    sumsCOUNT_df <- data.frame(t(sumsCOUNT_df[-1]))
    write.csv(dfmin_1, "dfmin_1.csv")
    write.csv(dfmax_1, "dfmax_1.csv")
    write.csv(sumsCOUNT_df, "sumsCOUNT_df.csv")
    df_min_cor <- read.csv("dfmin_1.csv")
    df_max_cor <- read.csv("dfmin_1.csv")
    sumsCOUNT_df <- read.csv("sumsCOUNT_df.csv")
    colnames(df_min_cor)[2] <- "min_cor"
    colnames(df_max_cor)[2] <- "max_cor"
    colnames(sumsCOUNT_df)[2] <- "sum_cor"
    
    sumsCOUNT_df_2 <- sumsCOUNT_df[apply(sumsCOUNT_df!=0, 1, all),]
    sumsCOUNT_df_2$X <- gsub("X","",as.character(sumsCOUNT_df_2$X))
    sumsCOUNT_df_2$sum_cor <- sumsCOUNT_df_2$sum_cor / 2
    colnames(sumsCOUNT_df_2)[1] <- "Prot.name(Uniprot)"
    colnames(sumsCOUNT_df_2)[1] <- "EC.number"
    write.csv(sumsCOUNT_df_2, "sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- read.csv("sumsCOUNT_df_2.csv")
    sumsCOUNT_df_2 <- sumsCOUNT_df_2[ , -1]
    library(brendaDb)
    df <- ReadBrenda("brenda_2023_1.txt")
    ID <- df[df$field == "RECOMMENDED_NAME",]
    ID <- ID[ , -c(2)]
    colnames(ID)[1] <- "EC.number"
    sum_coun_w_ID <- merge(sumsCOUNT_df_2, ID)
    sum(sum_coun_w_ID$sum_cor == 1)
    write.csv(sum_coun_w_ID, "prot_name_w_count_corr_EV.csv")
    df222 <- subset(sum_coun_w_ID, sum_coun_w_ID$sum_cor == 1)
    write.csv(df222, "Table_of_only_par_for_EV.csv")
    
    if (Pairs_or_not == "Pairs") {
      
      test_221 <- sumsCOUNT[sapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x) == 2)]
      test_222 <- lapply(1:length(test_221), function(x) names(test_221[[x]][test_221[[x]] == 1]))
      test_223 <- unlist(test_222)
      test_224 <- as.data.frame(test_223)
      gen_EC <- Evpersent[ , -c(3,4)]
      gen_EC_2 <- gen_EC[!duplicated(gen_EC), ]
      colnames(test_224)[1] <- "gen"
      gen_EC_2$gen <- sub("\\|$", "", gen_EC_2$gen)
      gen_EC_2$gen <- sub("\\gb|", "", gen_EC_2$gen)
      gen_EC_2$gen <- sub("\\|", "", gen_EC_2$gen)
      gen_EC_and_par <- merge(test_224, gen_EC_2)
      gen_EC_and_par <- merge(test_224, gen_EC_2)
      gen_EC_and_par$gen <- gsub('gb', '', gen_EC_and_par$gen)
      write.csv(gen_EC_and_par, "gen_EC_and_par_EV_per.csv")
      gen_EC_and_par <- read.csv("gen_EC_and_par_EV_per.csv")
      gen_EC_and_par <- gen_EC_and_par[ , -1]
      library(Biostrings)
      Gl_amino_fastas_after_tranc <- readAAStringSet(trancdecoder_sp)
      Gl_amino_fastas_after_tranc <- as.data.frame(Gl_amino_fastas_after_tranc)
      write.csv(Gl_amino_fastas_after_tranc, "EV_amino_fastas_after_tranc.csv")
      Gl_amino_fastas_after_tranc <- read.csv("EV_amino_fastas_after_tranc.csv")
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p1 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p2 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p3 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p4 .*","",Gl_amino_fastas_after_tranc$X)
      Gl_amino_fastas_after_tranc$X <- gsub("\\.p5 .*","",Gl_amino_fastas_after_tranc$X)
      colnames(Gl_amino_fastas_after_tranc)[1] <- "gen"
      colnames(Gl_amino_fastas_after_tranc)[2] <- "seq_amino"
      EC_gen_seq_GL <- merge(gen_EC_and_par, Gl_amino_fastas_after_tranc)
      df222_222 <- df222[ , -2]
      EC_gen_seq_GL_prot_name <- merge(EC_gen_seq_GL, df222_222)
      aaa <- EC_gen_seq_GL_prot_name[duplicated(EC_gen_seq_GL_prot_name[,1]) | duplicated(EC_gen_seq_GL_prot_name[,1], fromLast=TRUE),]
      
    } else if (Pairs_or_not == "Unpairs") {
      test_221 <- sumsCOUNT[sapply(sumsCOUNT[1:length(sumsCOUNT)], function(x) sum(x) != 0)]
      test_222 <- lapply(1:length(test_221), function(x) names(test_221[[x]][test_221[[x]] != 0]))
     }
      } else (print("Uncorrected analisis"))
}

pTT_EV <- c("gen", "EC.number", "EveT1.5", "EveB12", "EveB12A", "EveT12", "EveT18", "EveT24")
pTT_EC <- c("gen", "EC.number", "EcyT1.5", "EcyB12", "EcyB12A", "EcyT12", "EcyT18", "EcyT24")
pTT_GL <- c("gen", "EC.number", "GlaT1.5", "GlaB12", "GlaB12A", "GlaT12", "GlaT18", "GlaT24")

EV_exp_cor <- corrFUN("percent_tableEV.csv", pTT_EV, trancdecoder_sp = "EV_exp_median.fa.transdecoder.pep", analis = "E", Pairs_or_not = "Pairs", cor_lvl = -0.9)
EV_per_cor <- corrFUN("percent_tableEV.csv", pTT_EV, trancdecoder_sp = "EV_per_median.fa.transdecoder.pep", analis = "P", Pairs_or_not = "Pairs", cor_lvl = -0.9)
EC_exp_cor <- corrFUN("percent_tableEC.csv", pTT_EC, trancdecoder_sp = "EC_exp_median.fa.transdecoder.pep", analis = "E", Pairs_or_not = "Pairs", cor_lvl = -0.9)
EC_per_cor <- corrFUN("percent_tableEC.csv", pTT_EC, trancdecoder_sp = "EC_per_median.fa.transdecoder.pep", analis = "P", Pairs_or_not = "Pairs", cor_lvl = -0.9)
GL_exp_cor <- corrFUN("percent_tableGL.csv", pTT_GL, trancdecoder_sp = "GL_exp_median.fa.transdecoder.pep", analis = "E", Pairs_or_not = "Pairs", cor_lvl = -0.9)
GL_per_cor <- corrFUN("percent_tableGL.csv", pTT_GL, trancdecoder_sp = "GL_per_17_01_2025.fa.transdecoder.pep", analis = "P", Pairs_or_not = "Pairs", cor_lvl = -0.9)

write.csv(EV_exp_cor, "EV_exp_cor.csv")
write.csv(EV_per_cor, "EV_per_cor.csv")
write.csv(EC_exp_cor, "EC_exp_cor.csv")
write.csv(EC_per_cor, "EC_per_cor.csv")
write.csv(GL_exp_cor, "GL_exp_cor.csv")
write.csv(GL_per_cor, "GL_per_cor.csv")

#Количество коррелирующих изоформ всего
EV_exp_cor_unpairs <- corrFUN("percent_tableEV.csv", pTT_EV, trancdecoder_sp = "EV_exp_median.fa.transdecoder.pep", analis = "E", Pairs_or_not = "Unpairs", cor_lvl = -0.9)
EV_per_cor_unpairs <- corrFUN("percent_tableEV.csv", pTT_EV, trancdecoder_sp = "EV_per_median.fa.transdecoder.pep", analis = "P", Pairs_or_not = "Unpairs", cor_lvl = -0.9)
EC_exp_cor_unpairs <- corrFUN("percent_tableEC.csv", pTT_EC, trancdecoder_sp = "EC_exp_median.fa.transdecoder.pep", analis = "E", Pairs_or_not = "Unpairs", cor_lvl = -0.9)
EC_per_cor_unpairs <- corrFUN("percent_tableEC.csv", pTT_EC, trancdecoder_sp = "EC_per_median.fa.transdecoder.pep", analis = "P", Pairs_or_not = "Unpairs", cor_lvl = -0.9)
GL_exp_cor_unpairs <- corrFUN("percent_tableGL.csv", pTT_GL, trancdecoder_sp = "GL_exp_median.fa.transdecoder.pep", analis = "E", Pairs_or_not = "Unpairs", cor_lvl = -0.9)
GL_per_cor_unpairs <- corrFUN("percent_tableGL.csv", pTT_GL, trancdecoder_sp = "GL_per_17_01_2025.fa.transdecoder.pep", analis = "P", Pairs_or_not = "Unpairs", cor_lvl = -0.9)




###функция для разделения файлов фаста по EC номерам###
sepaFastaby_EC <- function(data_table, name_dir_to, back_dir) {
  library(seqRFLP)
  lenEC <- read.csv(data_table)
  df <- lenEC[ , -c(1, 2, 5)]
  df.fasta = dataframe2fas(df, file="len.fasta")
  IDandEC <- lenEC[ , c(2, 3)]
  IDandEC <- IDandEC[ , c("gen", "EC.number")]
  colnames(IDandEC)[1] <- "ID"
  library(Biostrings)
  library(seqinr)
  library(phylotools)
  myfasta <- read.fasta("len.fasta")
  unique(myfasta$seq.name)
  dir.create(name_dir_to)
  setwd(paste0(getwd(), "/", name_dir_to))
  split_dat(myfasta, IDandEC)
  file.remove("ungrouped.fasta")
  setwd(back_dir)
}

sepaFastaby_EC(data_table = "EV_exp_cor.csv", name_dir_to = "EV_Exp", 
               back_dir = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data")
sepaFastaby_EC(data_table = "EV_per_cor.csv", name_dir_to = "EV_Per", 
               back_dir = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data")
sepaFastaby_EC(data_table = "EC_exp_cor.csv", name_dir_to = "EC_Exp", 
               back_dir = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data")
sepaFastaby_EC(data_table = "EC_per_cor.csv", name_dir_to = "EC_Per", 
               back_dir = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data")
sepaFastaby_EC(data_table = "GL_exp_cor.csv", name_dir_to = "GL_Exp", 
               back_dir = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data")
sepaFastaby_EC(data_table = "GL_per_cor.csv", name_dir_to = "GL_Per", 
               back_dir = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data")

##Цикл для выравнивания mafft#BASH #выравнивание fasta для подсчета процентной идентичности пар последовательностей
#for x in *.fasta; do mafft $x > $x.aln; done

#Функция для подсчета процентной идентичности###

perIdent_FUN <- function(Path_where_fasta, Path_where_df_with_AAseq, df_with_AAseq) {
  setwd(Path_where_fasta)
  library(Biostrings)
  file.list <- list.files(pattern='.aln')
  fastasALNEC <- lapply(file.list, function(x) readAAStringSet(x))
  dfALNEC <- lapply(fastasALNEC, function(x) as.data.frame(x))
  PA <- lapply(dfALNEC[1:length(dfALNEC)], function(x) pairwiseAlignment(x[1 , ], x[2 , ]))
  ID_EC <- lapply(PA[1:length(PA)], function(x) pid(x, type = "PID1"))
  df_ID_EC <- as.data.frame(unlist(ID_EC))
  dffilelist <- as.data.frame(file.list)
  all <- cbind(dffilelist, df_ID_EC)
  all$file.list <- gsub("\\.fa.*","",all$file.list)
  colnames(all)[1] <- "EC.number"
  colnames(all)[2] <- "%ID"
  setwd(Path_where_df_with_AAseq)
  lenEC <- read.csv(df_with_AAseq) 
  lenEC$len <- stringr::str_length(lenEC$seq_amino)
  lenEC <- lenEC[ , -1]
  len_and_ID_EC <- merge(lenEC, all)
  nrow(len_and_ID_EC[len_and_ID_EC$`%ID` == 100 , ])
  WO100 <- len_and_ID_EC[!len_and_ID_EC$`%ID` == 100, ]
}

table_WO100_exp_EC <- perIdent_FUN(Path_where_fasta = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data/test_test",
                               Path_where_df_with_AAseq = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data",
                               df_with_AAseq = "EC_exp_cor.csv")


table_WO100_exp_EV <- perIdent_FUN(Path_where_fasta = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data/EV_Exp",
                              Path_where_df_with_AAseq = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data",
                              df_with_AAseq = "EV_exp_cor.csv")
table_WO100_per_EV <- perIdent_FUN(Path_where_fasta = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data/EV_Per",
                               Path_where_df_with_AAseq = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data",
                               df_with_AAseq = "EV_per_cor.csv")
table_WO100_exp_EC <- perIdent_FUN(Path_where_fasta = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data/EC_Exp",
                               Path_where_df_with_AAseq = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data",
                               df_with_AAseq = "EC_exp_cor.csv")
table_WO100_per_EC <- perIdent_FUN(Path_where_fasta = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data/EC_Per",
                               Path_where_df_with_AAseq = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data",
                               df_with_AAseq = "EC_per_cor.csv")
table_WO100_exp_GL <- perIdent_FUN(Path_where_fasta = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data/GL_Exp",
                               Path_where_df_with_AAseq = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data",
                               df_with_AAseq = "GL_exp_cor.csv")
table_WO100_per_GL <- perIdent_FUN(Path_where_fasta = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data/GL_Per",
                               Path_where_df_with_AAseq = "/media/Data/Work/disser/transcriptomes/EV/clean_scripts_and_data",
                               df_with_AAseq = "GL_per_cor.csv")


write.csv(table_WO100_exp_EV, "table_WO100_exp_EV.csv")
write.csv(table_WO100_per_EV, "table_WO100_per_EV.csv")
write.csv(table_WO100_exp_EC, "table_WO100_exp_EC.csv")
write.csv(table_WO100_per_EC, "table_WO100_per_EC.csv")
write.csv(table_WO100_exp_GL, "table_WO100_exp_GL.csv")
write.csv(table_WO100_per_GL, "table_WO100_per_GL.csv")
###Функция для создания конечной таблицы по анализу корреляции###

finDFcor <- function(table_sp, pattern_Sp, patter_Sp_1, WO100, analis) {
  if (analis == "E") {
    Evpersent <- read.csv(table_sp)
    Evpersent <- Evpersent[ , -c(1, 6) ]
    
    
    Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
    library(tidyverse)
    Evpersent2 <- Evpersent %>%
      group_by(Temperature, gen) %>%
      mutate(meanExpression = median(Expression, na.rm = TRUE))
    
    
    Evpersent2 <- Evpersent2[ , -4]
    Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]
    
    Evpersent4 <- Evpersent3 %>% 
      group_by(Temperature, gen) %>%
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = meanExpression)
    Evpersent4 <- Evpersent4[ , -3]
    Evpersent4 <- Evpersent4[ , -c(6, 4, 5)]
    
    Evpersent4 <- Evpersent4[ , pattern_Sp]
    
    Evpersent4 <- Evpersent4[duplicated(Evpersent4[,2]) | duplicated(Evpersent4[,2], fromLast=TRUE),]
    
    Evpersent2313 <- Evpersent4[ , -c(1, 2)]
    v <- Evpersent4[ , c(1, 2)]
    Evpersent2313$maxmedian <- apply(Evpersent2313, 1, max)
    Evpersent2313 <- cbind(v, Evpersent2313)
    Evpersent2313 <- Evpersent2313[ , c(1, 9)]
    Evpersent2313$gen <- sub("\\|$", "", Evpersent2313$gen)
    Evpersent2313$gen <- sub("\\gb|", "", Evpersent2313$gen)
    Evpersent2313$gen <- sub("\\|", "", Evpersent2313$gen)
    GL <- read.csv(WO100)
    GL_W_maxmedian <- merge(GL, Evpersent2313, by = "gen")
    GL_W_23_max_par_both <- GL_W_maxmedian
    Evpersent4 <- Evpersent4[ , -2 ] 
    Evpersent4$gen <- sub("\\|$", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\gb|", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\|", "", Evpersent4$gen)
    alltemperature <- merge(Evpersent4, GL_W_23_max_par_both, by = "gen")
    
    alltemperature <- alltemperature[ , patter_Sp_1]
  } else if (analis == "P") {
    Evpersent <- read.csv(table_sp)
    Evpersent <- Evpersent[ , -c(1, 5) ]
    
    
    Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
    library(tidyverse)
    Evpersent2 <- Evpersent %>%
      group_by(Temperature, gen) %>%
      mutate(medianPercent = median(percent, na.rm = TRUE))
    
    
    Evpersent2 <- Evpersent2[ , -4]
    Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]
    
    Evpersent4 <- Evpersent3 %>% 
      group_by(Temperature, gen) %>%
      mutate(id = row_number()) %>% 
      pivot_wider(names_from = Temperature, values_from = medianPercent)
    Evpersent4 <- Evpersent4[ , -3]
    Evpersent4 <- Evpersent4[ , -c(6, 4, 5)]
    
    Evpersent4 <- Evpersent4[ , pattern_Sp]
    
    Evpersent4 <- Evpersent4[duplicated(Evpersent4[,2]) | duplicated(Evpersent4[,2], fromLast=TRUE),]
    
    Evpersent2313 <- Evpersent4[ , -c(1, 2)]
    v <- Evpersent4[ , c(1, 2)]
    Evpersent2313$maxmedian <- apply(Evpersent2313, 1, max)
    Evpersent2313 <- cbind(v, Evpersent2313)
    Evpersent2313 <- Evpersent2313[ , c(1, 9)]
    Evpersent2313$gen <- sub("\\|$", "", Evpersent2313$gen)
    Evpersent2313$gen <- sub("\\gb|", "", Evpersent2313$gen)
    Evpersent2313$gen <- sub("\\|", "", Evpersent2313$gen)
    EC <- read.csv(WO100)
    EC_W_maxmedian <- merge(EC, Evpersent2313, by = "gen")
    nrow(EC_W_maxmedian[EC_W_maxmedian$maxmedian >= 50 & EC_W_maxmedian$maxmedian != 100, ])
    EC_W_23 <- EC_W_maxmedian[EC_W_maxmedian$maxmedian > 50 , ]
    EC_W_23_max_par_both <- EC_W_23[duplicated(EC_W_23[,3]) | duplicated(EC_W_23[,3], fromLast=TRUE),]
    hist(EC_W_23_max_par_both$maxmedian, breaks = 30)
    EC_W_23_max_par_both <- EC_W_23_max_par_both[ , -2]
    Evpersent4 <- Evpersent4[ , -2 ] 
    Evpersent4$gen <- sub("\\|$", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\gb|", "", Evpersent4$gen)
    Evpersent4$gen <- sub("\\|", "", Evpersent4$gen)
    alltemperature <- merge(Evpersent4, EC_W_23_max_par_both, by = "gen")
    alltemperature <- alltemperature[ , patter_Sp_1]
    
  } else (print("Uncorrected analisis"))
    
}


pTT_EV <- c("gen", "EC.number", "EveT1.5", "EveB12", "EveB12A", "EveT12", "EveT18", "EveT24")
pTT_EC <- c("gen", "EC.number", "EcyT1.5", "EcyB12", "EcyB12A", "EcyT12", "EcyT18", "EcyT24")
pTT_GL <- c("gen", "EC.number", "GlaT1.5", "GlaB12", "GlaB12A", "GlaT12", "GlaT18", "GlaT24")

pt_1_EV <- c("EC.number", "gen", "description", "len", "X.ID", "maxmedian", "EveT1.5", "EveB12", "EveB12A", "EveT12", "EveT18", "EveT24", "seq_amino")
pt_1_EC <- c("EC.number", "gen", "description", "len", "X.ID", "maxmedian", "EcyT1.5", "EcyB12", "EcyB12A", "EcyT12", "EcyT18", "EcyT24", "seq_amino")
pt_1_GL <- c("EC.number", "gen", "description", "len", "X.ID", "maxmedian", "GlaT1.5", "GlaB12", "GlaB12A", "GlaT12", "GlaT18", "GlaT24", "seq_amino")

Final_table_w_par_exp_EV <- finDFcor(table_sp = "percent_tableEV.csv", pattern_Sp = pTT_EV, patter_Sp_1 = pt_1_EV, WO100 = "table_WO100_exp_EV.csv",
                                     analis = "E")
Final_table_w_par_per_EV <- finDFcor(table_sp = "percent_tableEV.csv", pattern_Sp = pTT_EV, patter_Sp_1 = pt_1_EV, WO100 = "table_WO100_per_EV.csv",
                                     analis = "P")
Final_table_w_par_exp_EC <- finDFcor(table_sp = "percent_tableEC.csv", pattern_Sp = pTT_EC, patter_Sp_1 = pt_1_EC, WO100 = "table_WO100_exp_EC.csv",
                                     analis = "E")
Final_table_w_par_per_EC <- finDFcor(table_sp = "percent_tableEC.csv", pattern_Sp = pTT_EC, patter_Sp_1 = pt_1_EC, WO100 = "table_WO100_per_EC.csv",
                                     analis = "P")
Final_table_w_par_exp_GL <- finDFcor(table_sp = "percent_tableGL.csv", pattern_Sp = pTT_GL, patter_Sp_1 = pt_1_GL, WO100 = "table_WO100_exp_GL.csv",
                                     analis = "E")
Final_table_w_par_per_GL <- finDFcor(table_sp = "percent_tableGL.csv", pattern_Sp = pTT_GL, patter_Sp_1 = pt_1_GL, WO100 = "table_WO100_per_GL.csv",
                                     analis = "P")

library(xlsx)
Final_table_w_par_exp_EV <- Final_table_w_par_exp_EV[order(Final_table_w_par_exp_EV$EC.number),]
write.xlsx(Final_table_w_par_exp_EV, "Final_table_w_par_exp_EV.xlsx")

Final_table_w_par_per_EV <- Final_table_w_par_per_EV[order(Final_table_w_par_per_EV$EC.number),]
write.xlsx(Final_table_w_par_per_EV, "Final_table_w_par_per_EV.xlsx")

Final_table_w_par_exp_EC <- Final_table_w_par_exp_EC[order(Final_table_w_par_exp_EC$EC.number),]
write.xlsx(Final_table_w_par_exp_EC, "Final_table_w_par_exp_EC.xlsx")

Final_table_w_par_per_EC <- Final_table_w_par_per_EC[order(Final_table_w_par_per_EC$EC.number),]
write.xlsx(Final_table_w_par_per_EC, "Final_table_w_par_per_EC.xlsx")

Final_table_w_par_exp_GL <- Final_table_w_par_exp_GL[order(Final_table_w_par_exp_GL$EC.number),]
write.xlsx(Final_table_w_par_exp_GL, "Final_table_w_par_exp_GL.xlsx")

Final_table_w_par_per_GL <- Final_table_w_par_per_GL[order(Final_table_w_par_per_GL$EC.number),]
write.xlsx(Final_table_w_par_per_GL, "Final_table_w_par_per_GL.xlsx")



#Создаем 2 столбца в xlsx файлах (mLen и cor with temperature)



#Function для поиска пар ферментов коррелирующих с температурой
FunCorr_with_temp <- function(df, cor_temp, mLen) {
  library(xlsx)
  SP <- read.xlsx(df, 1)
  SP_1 <- SP[SP$cor_temp <= -cor_temp | SP$cor_temp >= cor_temp, ]
  SP_2 <- SP_1[duplicated(SP_1[,"EC.number"]) | duplicated(SP_1[,"EC.number"], fromLast=TRUE),]
  SP_3 <- SP_2[SP_2$mLen >= -mLen & SP_2$mLen <= mLen, ]
}

EV_cor_temp_per <- FunCorr_with_temp("Final_table_w_par_per_EV.xlsx", 0.4, 0.1)
EC_cor_temp_per <- FunCorr_with_temp("Final_table_w_par_per_EC.xlsx", 0.4, 0.1)
GL_cor_temp_per <- FunCorr_with_temp("Final_table_w_par_per_GL.xlsx", 0.4, 0.1)


dG <- c(1.364, -0.227, 1.091, -0.012, 0.207, 0.052, 0.485, 2.322, -0.296, 1.047, 0.674, 0.638, -0.195)
df_dG <- as.data.frame(dG)
wilcox.test(dG)
t.test(dG, alternative = "greater") 


#выбераем только те пары транскриптов у которых точечные замены после выравнивания



###

#DECRE
All_EV_DECRE <- rbind(DECRE_all_down_EV, DECRE_all_Up_EV)
Pairs_EV_DECRE <- All_EV_DECRE[duplicated(All_EV_DECRE[,3]) | duplicated(All_EV_DECRE[,3], fromLast=TRUE),]
Pairs_EV_DECRE_u_ECnumbers <- unique(Pairs_EV_DECRE$EC.number)

All_EC_DECRE <- rbind(DECRE_all_down_EC, DECRE_all_Up_EC)
Pairs_EC_DECRE <- All_EC_DECRE[duplicated(All_EC_DECRE[,3]) | duplicated(All_EC_DECRE[,3], fromLast=TRUE),]
Pairs_EC_DECRE_u_ECnumbers <- unique(Pairs_EC_DECRE$EC.number)

All_GL_DECRE <- rbind(DECRE_all_down_GL, DECRE_all_Up_GL)
Pairs_GL_DECRE <- All_GL_DECRE[duplicated(All_GL_DECRE[,3]) | duplicated(All_GL_DECRE[,3], fromLast=TRUE),]
Pairs_GL_DECRE_u_ECnumbers <- unique(Pairs_GL_DECRE$EC.number)
#DECPE
All_EV_DECPE <- rbind(DECPE_all_down_EV, DECPE_all_Up_EV)
Pairs_EV_DECPE <- All_EV_DECPE[duplicated(All_EV_DECPE[,3]) | duplicated(All_EV_DECPE[,3], fromLast=TRUE),]
Pairs_EV_DECPE_u_ECnumbers <- unique(Pairs_EV_DECPE$EC.number)

All_EC_DECPE <- rbind(DECPE_all_down_EC, DECPE_all_Up_EC)
Pairs_EC_DECPE <- All_EC_DECPE[duplicated(All_EC_DECPE[,3]) | duplicated(All_EC_DECPE[,3], fromLast=TRUE),]
Pairs_EC_DECPE_u_ECnumbers <- unique(Pairs_EC_DECPE$EC.number)

All_GL_DECPE <- rbind(DECPE_all_down_GL, DECPE_all_Up_GL)
Pairs_GL_DECPE <- All_GL_DECPE[duplicated(All_GL_DECPE[,3]) | duplicated(All_GL_DECPE[,3], fromLast=TRUE),]
Pairs_GL_DECPE_u_ECnumbers <- unique(Pairs_GL_DECPE$EC.number)


#DECRE_reciproc_pairs
All_EV_DECRE <- rbind(DECRE_all_down_EV, DECRE_all_Up_EV)
reciproc_pairs_DECRE_EV <- unique(DECRE_all_Up_EV$EC.number)[unique(DECRE_all_Up_EV$EC.number) %in% unique(DECRE_all_down_EV$EC.number)]

All_EC_DECRE <- rbind(DECRE_all_down_EC, DECRE_all_Up_EC)
reciproc_pairs_DECRE_EC <- unique(DECRE_all_Up_EC$EC.number)[unique(DECRE_all_Up_EC$EC.number) %in% unique(DECRE_all_down_EC$EC.number)]

All_GL_DECRE <- rbind(DECRE_all_down_GL, DECRE_all_Up_GL)
reciproc_pairs_DECRE_GL <- unique(DECRE_all_Up_GL$EC.number)[unique(DECRE_all_Up_GL$EC.number) %in% unique(DECRE_all_down_GL$EC.number)]

#DECPE_reciproc_pairs
All_EV_DECPE <- rbind(DECPE_all_down_EV, DECPE_all_Up_EV)
reciproc_pairs_DECPE_EV <- unique(DECPE_all_Up_EV$EC.number)[unique(DECPE_all_Up_EV$EC.number) %in% unique(DECPE_all_down_EV$EC.number)]

All_EC_DECPE <- rbind(DECPE_all_down_EC, DECPE_all_Up_EC)
reciproc_pairs_DECPE_EC <- unique(DECPE_all_Up_EC$EC.number)[unique(DECPE_all_Up_EC$EC.number) %in% unique(DECPE_all_down_EC$EC.number)]

All_GL_DECPE <- rbind(DECPE_all_down_GL, DECPE_all_Up_GL)
reciproc_pairs_DECPE_GL <- unique(DECPE_all_Up_GL$EC.number)[unique(DECPE_all_Up_GL$EC.number) %in% unique(DECPE_all_down_GL$EC.number)]

type_analisis_v2 <- c("DECRE", "CECRA", "DECPE", "CECPA", "DECRE", "CECRA", "DECPE", "CECPA", "DECRE", "CECRA", "DECPE", "CECPA")
Species_v2 <- c("E. verrucosus", "E. verrucosus", "E. verrucosus", "E. verrucosus", "E. cyaneus", "E. cyaneus", "E. cyaneus", "E. cyaneus", "G. lacustris", "G. lacustris", "G. lacustris", "G. lacustris")
Numbers_v2 <- c(length(reciproc_pairs_DECRE_EV), length(EV_exp_cor_unpairs), length(reciproc_pairs_DECPE_EV), length(EV_per_cor_unpairs),
                length(reciproc_pairs_DECRE_EC), length(EC_exp_cor_unpairs), length(reciproc_pairs_DECPE_EC), length(EC_per_cor_unpairs),
                length(reciproc_pairs_DECRE_GL), length(GL_exp_cor_unpairs), length(reciproc_pairs_DECPE_GL), length(GL_per_cor_unpairs))

df_figure_5_v2 <- data.frame(type_analisis_v2, Species_v2, Numbers_v2) 
write.csv(df_figure_5_v2, "df_figure_5_v2.csv")
df_figure_5_v2 <- read.csv("df_figure_5_v2.csv")
df_figure_5_v2$type_analisis_v2 <- factor(df_figure_5_v2$type_analisis_v2, levels = c("DECRE", "CECRA", "DECPE", "CECPA"))
df_figure_5_v2$Species_v2 <- factor(df_figure_5_v2$Species_v2, levels = c("E. verrucosus", "E. cyaneus", "G. lacustris"))
library(ggplot2)
ggplot(df_figure_5_v2, aes(x = Species_v2, y = Numbers_v2, fill = Species_v2, col = Species_v2, alpha = type_analisis_v2)) +
  geom_bar(position = "dodge", stat="identity") + 
  theme_classic() +
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_alpha_manual(values = c(.1, .35, .65, .95)) +
  scale_x_discrete(name = NULL, 
                   labels = c(substitute(paste(italic("E. verrucosus"))),
                              substitute(paste(italic("E. cyaneus"))), 
                              substitute(paste(italic("G. lacustris"))))) +
  #theme(legend.position="none", legend.title = element_blank()) +
  theme(text = element_text(family="Arial", size = 15, color = "black"),
        axis.text = element_text(family="Arial", size = 13, color = "black")) + 
  theme(legend.position = "none") + 
  ylab("Number of enzymes")
ggsave("Figure_5_V2.png", width = 16, height = 13, units = "cm", dpi = 400)





type_analisis <- c("DECRE", "CECRA", "DECPE", "CECPA", "DECRE", "CECRA", "DECPE", "CECPA", "DECRE", "CECRA", "DECPE", "CECPA")
Species <- c("E. verrucosus", "E. verrucosus", "E. verrucosus", "E. verrucosus", "E. cyaneus", "E. cyaneus", "E. cyaneus", "E. cyaneus", "G. lacustris", "G. lacustris", "G. lacustris", "G. lacustris")
Numbers <- c(length(Pairs_EV_DECRE_u_ECnumbers), length(EV_exp_cor_unpairs), length(Pairs_EV_DECPE_u_ECnumbers), length(EV_per_cor_unpairs),
             length(Pairs_EC_DECRE_u_ECnumbers), length(EC_exp_cor_unpairs), length(Pairs_EC_DECPE_u_ECnumbers), length(EC_per_cor_unpairs),
             length(Pairs_GL_DECRE_u_ECnumbers), length(GL_exp_cor_unpairs), length(Pairs_GL_DECPE_u_ECnumbers), length(GL_per_cor_unpairs))

df_figure_5 <- data.frame(type_analisis, Species, Numbers)
write.csv(df_figure_5, "df_figure_5.csv")
library(ggplot2)
library(ggpattern)
df_figure_5$type_analisis <- factor(df_figure_5$type_analisis, levels = c("DECRE", "CECRA", "DECPE", "CECPA"))
df_figure_5$Species <- factor(df_figure_5$Species, levels = c("E. verrucosus", "E. cyaneus", "G. lacustris"))
cvet <- c("#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9",
          "#009E73", "#009E73", "#009E73", "#009E73", 
          "#E69F00", "#E69F00", "#E69F00")


ggplot(df_figure_5, aes(x = Species, y = Numbers, fill = Species, col = Species, alpha = Species, pattern_fill = Species,
                        pattern_color = Species)) + 
  theme_classic() +
  geom_col_pattern(aes(pattern = type_analisis, pattern_density = Species, pattern_spacing = Species, 
                       pattern_alpha = Species), position = "dodge") + 
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_alpha_manual(values = c(.2, .2, .2)) + 
  scale_pattern_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_pattern_color_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_pattern_manual(values = c("crosshatch", "none", "circle", "stripe")) +
  theme(legend.position = "none") +
  scale_pattern_density_manual(values = c(0.1, 0.1, 0.1)) +
  scale_pattern_spacing_manual(values = c(0.03, 0.03, 0.03)) +
  scale_pattern_alpha_manual(values = c(.5, .5, .5))




ggplot(df_figure_5, aes(x = Species, y = Numbers, fill = Species, col = Species, alpha = type_analisis)) +
  geom_bar(position = "dodge", stat="identity") + 
  theme_classic() +
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#009E73", "#56B4E9", "#E69F00")) +
  scale_alpha_manual(values = c(.1, .35, .65, .95)) +
  scale_x_discrete(name = NULL, 
                   labels = c(substitute(paste(italic("E. verrucosus"))),
                              substitute(paste(italic("E. cyaneus"))), 
                              substitute(paste(italic("G. lacustris"))))) +
  #theme(legend.position="none", legend.title = element_blank()) +
  theme(text = element_text(family="Arial", size = 15, color = "black"),
        axis.text = element_text(family="Arial", size = 13, color = "black")) + 
  theme(legend.position = "none")

ggsave("Figure_5.png", width = 16, height = 13, units = "cm", dpi = 400)








df_for_figure5_DECRE <- df_figure_5[df_figure_5$type_analisis == "DECRE",]
df_for_figure5_DECRE <- df_for_figure5_DECRE[ , -1]
df_for_figure5_DECPE <- df_figure_5[df_figure_5$type_analisis == "DECPE",]
df_for_figure5_DECPE <- df_for_figure5_DECPE[ , -1]
df_for_figure5_CECRA <- df_figure_5[df_figure_5$type_analisis == "CECRA",]
df_for_figure5_CECRA <- df_for_figure5_CECRA[ , -1]
df_for_figure5_CECPA <- df_figure_5[df_figure_5$type_analisis == "CECPA",]
df_for_figure5_CECPA <- df_for_figure5_CECPA[ , -1]

cvet <- c( "#009E73", "#56B4E9", "#E69F00")
df_for_figure5_DECRE$Species <- factor(df_for_figure5_DECRE$Species, levels = c("E. verrucosus", "E. cyaneus", "G. lacustris"))
#p1DECRE
ggplot(df_for_figure5_DECRE, aes(x = fct_inorder(Species), fill=Species)) + 
  geom_bar(aes(y=Numbers), position = "dodge", stat="identity", width = 0.5, col = cvet1, alpha = .5)+ 
  geom_hline(aes(yintercept = 0), col = "white") +
  ylab(NULL)+
  xlab("")+ 
  theme_classic(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=35, family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = "none",
        legend.direction = "horizontal") +
  labs(fill = NULL) +
  scale_fill_manual(values = cvet, labels = c(substitute(paste(italic('E. cyaneus'))),
                                              substitute(paste(italic('E. verrucosus'))), 
                                              substitute(paste(italic('G. lacustris'))))) + 
  ggtitle("DECRE") +
  theme(plot.margin = margin(t=5,1,1,1, "lines"))  +
  theme(plot.title = element_text(vjust = 3.5))






dG <- c(1.364, -0.227, 1.091, -0.012, 0.207, 0.052, 0.485, 2.322, -0.296, 1.047, 0.674, 0.638, -0.195)
df_dG <- as.data.frame(dG)
wilcox.test(dG)
t.test(dG, alternative = "greater")
shapiro.test(dG)
ks.test(dG)
  
library(ggplot2)
ggplot(df_dG, aes(sample = dG)) + 
  stat_qq(col = "black", size = 4, shape = 21, fill = "grey35") + 
  stat_qq_line(linetype = "dashed", col = "black", size = 1.2) + 
  theme_classic() +
  ylab("ΔΔG, kcal/mol") + 
  xlab(NULL) +
  ggtitle("Q-Q plot for ΔΔG of psychrophilic isoforms") +
  theme(text = element_text(family="Arial", size = 17, color = "black"),
        axis.text = element_text(family="Arial", size = 20, color = "black")) +
  theme(axis.title.y = element_text(family="Arial", size = 22, color = "black")) +
  annotate("text", x=-0.5, y=1.7, label= "Shapiro-Wilk normality test:", size = 6, color = "black", family="Arial") +
  annotate("text", x=-0.5, y=1.45, label= substitute(paste(italic("p"), "= 0.2011")), size = 6, color = "black", family="Arial") 
  
ggsave("Figure_5_q-qplot.png", width = 15, height = 13, units = "cm", dpi = 400)
df_dG$name <- "dG"

ggplot(df_dG, aes(y = dG, x = name)) +
  theme_classic() +
  geom_boxplot(fill = "grey80", size =1, col = "black") +
  geom_point(size = 4, shape = 21, col = "black", fill = "grey35") +
  ylab("ΔΔG, kcal/mol") +
  xlab(NULL) +
  theme(text = element_text(family="Arial", size = 17, color = "black"),
        axis.text = element_text(family="Arial", size = 20, color = "black")) +
  ggtitle("ΔΔG mutations of psychophilic isoforms") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(family="Arial", size = 22, color = "black")) +
  theme(axis.line.x = element_line(color = "black")) +
  theme(axis.line.y = element_line(color = "black")) +
  annotate("text", x=1.35, y=1.73, label= "One Sample t-test:", size = 6, color = "black", family="Arial") +
  annotate("text", x=1.35, y=1.57,label= substitute(paste(italic("p"), "= 0.01139")), size = 6, color = "black", family="Arial") 

ggsave("Figure_5_boxplot.png", width = 15, height = 13, units = "cm", dpi = 400)

pTT_EV <- c("gen", "EC.number", "EveT1.5", "EveB12", "EveB12A", "EveT12", "EveT18", "EveT24")
pTT_EC <- c("gen", "EC.number", "EcyT1.5", "EcyB12", "EcyB12A", "EcyT12", "EcyT18", "EcyT24")
pTT_GL <- c("gen", "EC.number", "GlaT1.5", "GlaB12", "GlaB12A", "GlaT12", "GlaT18", "GlaT24")

EV_exp_cor <- corrFUN("percent_tableEV.csv", pTT_EV, trancdecoder_sp = "EV_exp_median.fa.transdecoder.pep", analis = "E")
EV_per_cor <- corrFUN("percent_tableEV.csv", pTT_EV, trancdecoder_sp = "EV_per_median.fa.transdecoder.pep", analis = "P")
EC_exp_cor <- corrFUN("percent_tableEC.csv", pTT_EC, trancdecoder_sp = "EC_exp_median.fa.transdecoder.pep", analis = "E")
write.csv(EC_exp_cor, "EC_exp_cor.csv")
EC_per_cor <- corrFUN("percent_tableEC.csv", pTT_EC, trancdecoder_sp = "EC_per_median.fa.transdecoder.pep", analis = "P")
GL_exp_cor <- corrFUN("percent_tableGL.csv", pTT_GL, trancdecoder_sp = "GL_exp_median.fa.transdecoder.pep", analis = "E")
GL_per_cor <- corrFUN("percent_tableGL.csv", pTT_GL, trancdecoder_sp = "GL_per_median.fa.transdecoder.pep", analis = "P")



###DECRE_EV###

DECRE_resultsEV_1.5 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_1_5_EV.csv", 
                                      nameConrol = "EveB12A_", nameExperiment = "EveT1.5_")

DECRE_resultsEV_6 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_6_EV.csv", 
                                      nameConrol = "EveT12_", nameExperiment = "EveB12_")

DECRE_resultsEV_18 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_18_EV.csv", 
                                      nameConrol = "EveT12_", nameExperiment = "EveT18_")

DECRE_resultsEV_24 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_24_EV.csv", 
                                      nameConrol = "EveT12_", nameExperiment = "EveT24_")

###DECPE_EV###

DECPE_resultsEV_1.5 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_1_5_EV.csv", 
                                            nameConrol = "EveB12A_", nameExperiment = "EveT1.5_")

DECPE_resultsEV_6 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_6_EV.csv", 
                                          nameConrol = "EveT12_", nameExperiment = "EveB12_")

DECPE_resultsEV_18 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_18_EV.csv", 
                                           nameConrol = "EveT12_", nameExperiment = "EveT18_")

DECPE_resultsEV_24 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_24_EV.csv", 
                                           nameConrol = "EveT12_", nameExperiment = "EveT24_")

###DECSE_EV###

DECSE_resultsEV_1.5 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_1_5_EV.csv", 
                                            nameConrol = "EveB12A_", nameExperiment = "EveT1.5_")

DECSE_resultsEV_6 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_6_EV.csv", 
                                          nameConrol = "EveT12_", nameExperiment = "EveB12_")

DECSE_resultsEV_18 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_18_EV.csv", 
                                           nameConrol = "EveT12_", nameExperiment = "EveT18_")

DECSE_resultsEV_24 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EV1.csv", Eve, "Sample_data_24_EV.csv", 
                                           nameConrol = "EveT12_", nameExperiment = "EveT24_")

###DOWN and UP by DECRE EV####

DECRE_Down_EV_1.5 <- filterlfpv(ddsdata = DECRE_resultsEV_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

DECRE_Up_EV_1.5 <- filterlfpv(ddsdata = DECRE_resultsEV_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

DECRE_Down_EV_6 <- filterlfpv(ddsdata = DECRE_resultsEV_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

DECRE_Up_EV_6 <- filterlfpv(ddsdata = DECRE_resultsEV_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

DECRE_Down_EV_18 <- filterlfpv(ddsdata = DECRE_resultsEV_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

DECRE_Up_EV_18 <- filterlfpv(ddsdata = DECRE_resultsEV_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

DECRE_Down_EV_24 <- filterlfpv(ddsdata = DECRE_resultsEV_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

DECRE_Up_EV_24 <- filterlfpv(ddsdata = DECRE_resultsEV_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECRE")

###DOWN and UP by DECPE EV####

DECPE_Down_EV_1.5 <- filterlfpv(ddsdata = DECPE_resultsEV_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

DECPE_Up_EV_1.5 <- filterlfpv(ddsdata = DECPE_resultsEV_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

DECPE_Down_EV_6 <- filterlfpv(ddsdata = DECPE_resultsEV_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

DECPE_Up_EV_6 <- filterlfpv(ddsdata = DECPE_resultsEV_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

DECPE_Down_EV_18 <- filterlfpv(ddsdata = DECPE_resultsEV_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

DECPE_Up_EV_18 <- filterlfpv(ddsdata = DECPE_resultsEV_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

DECPE_Down_EV_24 <- filterlfpv(ddsdata = DECPE_resultsEV_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

DECPE_Up_EV_24 <- filterlfpv(ddsdata = DECPE_resultsEV_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECPE")

###DOWN and UP by DECSE EV####

DECSE_Down_EV_1.5 <- filterlfpv(ddsdata = DECSE_resultsEV_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

DECSE_Up_EV_1.5 <- filterlfpv(ddsdata = DECSE_resultsEV_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

DECSE_Down_EV_6 <- filterlfpv(ddsdata = DECSE_resultsEV_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

DECSE_Up_EV_6 <- filterlfpv(ddsdata = DECSE_resultsEV_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

DECSE_Down_EV_18 <- filterlfpv(ddsdata = DECSE_resultsEV_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

DECSE_Up_EV_18 <- filterlfpv(ddsdata = DECSE_resultsEV_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

DECSE_Down_EV_24 <- filterlfpv(ddsdata = DECSE_resultsEV_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

DECSE_Up_EV_24 <- filterlfpv(ddsdata = DECSE_resultsEV_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EV1.csv", type_analisis = "DECSE")

#Уникальные, значимо изменяющиеся гены DECRE EV

DECRE_all_down_EV <- rbind(DECRE_Down_EV_1.5, DECRE_Down_EV_6, DECRE_Down_EV_18, DECRE_Down_EV_24)

DECRE_all_Up_EV <- rbind(DECRE_Up_EV_1.5, DECRE_Up_EV_6, DECRE_Up_EV_18, DECRE_Up_EV_24)

uniq_DECRE_down_EV <- UniqECnumber(DECRE_all_down_EV, "Down", type_analisis = "DECRE")

uniq_DECRE_Up_EV <- UniqECnumber(DECRE_all_Up_EV, "Up", type_analisis = "DECRE")

DECPE_all_down_EV <- rbind(DECPE_Down_EV_1.5, DECPE_Down_EV_6, DECPE_Down_EV_18, DECPE_Down_EV_24)

DECPE_all_Up_EV <- rbind(DECPE_Up_EV_1.5, DECPE_Up_EV_6, DECPE_Up_EV_18, DECPE_Up_EV_24)

uniq_DECPE_down_EV <- UniqECnumber(DECPE_all_down_EV, "Down", type_analisis = "DECPE")

uniq_DECPE_Up_EV <- UniqECnumber(DECPE_all_Up_EV, "Up", type_analisis = "DECPE")

DECSE_all_down_EV <- rbind(DECSE_Down_EV_1.5, DECSE_Down_EV_6, DECSE_Down_EV_18, DECSE_Down_EV_24)

DECSE_all_Up_EV <- rbind(DECSE_Up_EV_1.5, DECSE_Up_EV_6, DECSE_Up_EV_18, DECSE_Up_EV_24)

uniq_DECSE_down_EV <- UniqECnumber(DECSE_all_down_EV, "Down", type_analisis = "DECSE")

uniq_DECSE_Up_EV <- UniqECnumber(DECSE_all_Up_EV, "Up", type_analisis = "DECSE")


###DECRE_EC###

DECRE_resultsEC_1.5 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_1_5_EC.csv", 
                                            nameConrol = "EcyB12A_", nameExperiment = "EcyT1.5_")

DECRE_resultsEC_6 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_6_EC.csv", 
                                          nameConrol = "EcyT12_", nameExperiment = "EcyB12_")

DECRE_resultsEC_18 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_18_EC.csv", 
                                           nameConrol = "EcyT12_", nameExperiment = "EcyT18_")

DECRE_resultsEC_24 <- dif_express_function(DEC = "DECRE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_24_EC.csv", 
                                           nameConrol = "EcyT12_", nameExperiment = "EcyT24_")

###DECPE_EC###

DECPE_resultsEC_1.5 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_1_5_EC.csv", 
                                            nameConrol = "EcyB12A_", nameExperiment = "EcyT1.5_")

DECPE_resultsEC_6 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_6_EC.csv", 
                                          nameConrol = "EcyT12_", nameExperiment = "EcyB12_")

DECPE_resultsEC_18 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_18_EC.csv", 
                                           nameConrol = "EcyT12_", nameExperiment = "EcyT18_")

DECPE_resultsEC_24 <- dif_express_function(DEC = "DECPE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_24_EC.csv", 
                                           nameConrol = "EcyT12_", nameExperiment = "EcyT24_")

###DECSE_EC###

DECSE_resultsEC_1.5 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_1_5_EC.csv", 
                                            nameConrol = "EcyB12A_", nameExperiment = "EcyT1.5_")

DECSE_resultsEC_6 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_6_EC.csv", 
                                          nameConrol = "EcyT12_", nameExperiment = "EcyB12_")

DECSE_resultsEC_18 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_18_EC.csv", 
                                           nameConrol = "EcyT12_", nameExperiment = "EcyT18_")

DECSE_resultsEC_24 <- dif_express_function(DEC = "DECSE", "30-01-2024_EConlynMATRIX4numbers_EC1.csv", Ecy, "Sample_data_24_EC.csv", 
                                           nameConrol = "EcyT12_", nameExperiment = "EcyT24_")

###DOWN and UP by DECRE EC####

DECRE_Down_EC_1.5 <- filterlfpv(ddsdata = DECRE_resultsEC_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv" , type_analisis = "DECRE")

DECRE_Up_EC_1.5 <- filterlfpv(ddsdata = DECRE_resultsEC_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECRE")

DECRE_Down_EC_6 <- filterlfpv(ddsdata = DECRE_resultsEC_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECRE")

DECRE_Up_EC_6 <- filterlfpv(ddsdata = DECRE_resultsEC_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECRE")

DECRE_Down_EC_18 <- filterlfpv(ddsdata = DECRE_resultsEC_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECRE")

DECRE_Up_EC_18 <- filterlfpv(ddsdata = DECRE_resultsEC_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECRE")

DECRE_Down_EC_24 <- filterlfpv(ddsdata = DECRE_resultsEC_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECRE")

DECRE_Up_EC_24 <- filterlfpv(ddsdata = DECRE_resultsEC_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECRE")


###DOWN and UP by DECPE EC####

DECPE_Down_EC_1.5 <- filterlfpv(ddsdata = DECPE_resultsEC_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")

DECPE_Up_EC_1.5 <- filterlfpv(ddsdata = DECPE_resultsEC_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")

DECPE_Down_EC_6 <- filterlfpv(ddsdata = DECPE_resultsEC_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")

DECPE_Up_EC_6 <- filterlfpv(ddsdata = DECPE_resultsEC_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")

DECPE_Down_EC_18 <- filterlfpv(ddsdata = DECPE_resultsEC_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")

DECPE_Up_EC_18 <- filterlfpv(ddsdata = DECPE_resultsEC_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")

DECPE_Down_EC_24 <- filterlfpv(ddsdata = DECPE_resultsEC_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")

DECPE_Up_EC_24 <- filterlfpv(ddsdata = DECPE_resultsEC_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECPE")


###DOWN and UP by DECSE EC####

DECSE_Down_EC_1.5 <- filterlfpv(ddsdata = DECSE_resultsEC_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")

DECSE_Up_EC_1.5 <- filterlfpv(ddsdata = DECSE_resultsEC_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")

DECSE_Down_EC_6 <- filterlfpv(ddsdata = DECSE_resultsEC_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")

DECSE_Up_EC_6 <- filterlfpv(ddsdata = DECSE_resultsEC_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")

DECSE_Down_EC_18 <- filterlfpv(ddsdata = DECSE_resultsEC_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")

DECSE_Up_EC_18 <- filterlfpv(ddsdata = DECSE_resultsEC_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")

DECSE_Down_EC_24 <- filterlfpv(ddsdata = DECSE_resultsEC_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")

DECSE_Up_EC_24 <- filterlfpv(ddsdata = DECSE_resultsEC_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "30-01-2024_EConlynMATRIX4numbers_EC1.csv", type_analisis = "DECSE")


#Уникальные, значимо изменяющиеся гены DECRE DECPE DECSE EC

DECRE_all_down_EC <- rbind(DECRE_Down_EC_1.5, DECRE_Down_EC_6, DECRE_Down_EC_18, DECRE_Down_EC_24)

DECRE_all_Up_EC <- rbind(DECRE_Up_EC_1.5, DECRE_Up_EC_6, DECRE_Up_EC_18, DECRE_Up_EC_24)

uniq_DECRE_down_EC <- UniqECnumber(DECRE_all_down_EC, "Down", type_analisis = "DECRE")

uniq_DECRE_Up_EC <- UniqECnumber(DECRE_all_Up_EC, "Up", type_analisis = "DECRE")

DECPE_all_down_EC <- rbind(DECPE_Down_EC_1.5, DECPE_Down_EC_6, DECPE_Down_EC_18, DECPE_Down_EC_24)

DECPE_all_Up_EC <- rbind(DECPE_Up_EC_1.5, DECPE_Up_EC_6, DECPE_Up_EC_18, DECPE_Up_EC_24)

uniq_DECPE_down_EC <- UniqECnumber(DECPE_all_down_EC, "Down", type_analisis = "DECPE")

uniq_DECPE_Up_EC <- UniqECnumber(DECPE_all_Up_EC, "Up", type_analisis = "DECPE")

DECSE_all_down_EC <- rbind(DECSE_Down_EC_1.5, DECSE_Down_EC_6, DECSE_Down_EC_18, DECSE_Down_EC_24)

DECSE_all_Up_EC <- rbind(DECSE_Up_EC_1.5, DECSE_Up_EC_6, DECSE_Up_EC_18, DECSE_Up_EC_24)

uniq_DECSE_down_EC <- UniqECnumber(DECSE_all_down_EC, "Down", type_analisis = "DECSE")

uniq_DECSE_Up_EC <- UniqECnumber(DECSE_all_Up_EC, "Up", type_analisis = "DECSE")


###DECRE_GL###

DECRE_resultsGL_1.5 <- dif_express_function(DEC = "DECRE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_1_5_GL.csv", 
                                            nameConrol = "GlaB12A_", nameExperiment = "GlaT1.5_")

DECRE_resultsGL_6 <- dif_express_function(DEC = "DECRE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_6_GL.csv", 
                                          nameConrol = "GlaT12_", nameExperiment = "GlaB12_")

DECRE_resultsGL_18 <- dif_express_function(DEC = "DECRE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_18_GL.csv", 
                                           nameConrol = "GlaT12_", nameExperiment = "GlaT18_")

DECRE_resultsGL_24 <- dif_express_function(DEC = "DECRE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_24_GL.csv", 
                                           nameConrol = "GlaT12_", nameExperiment = "GlaT24_")

###DECPE_GL###

DECPE_resultsGL_1.5 <- dif_express_function(DEC = "DECPE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_1_5_GL.csv", 
                                            nameConrol = "GlaB12A_", nameExperiment = "GlaT1.5_")

DECPE_resultsGL_6 <- dif_express_function(DEC = "DECPE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_6_GL.csv", 
                                          nameConrol = "GlaT12_", nameExperiment = "GlaB12_")

DECPE_resultsGL_18 <- dif_express_function(DEC = "DECPE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_18_GL.csv", 
                                           nameConrol = "GlaT12_", nameExperiment = "GlaT18_")

DECPE_resultsGL_24 <- dif_express_function(DEC = "DECPE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_24_GL.csv", 
                                           nameConrol = "GlaT12_", nameExperiment = "GlaT24_")

###DECSE_GL###

DECSE_resultsGL_1.5 <- dif_express_function(DEC = "DECSE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_1_5_GL.csv", 
                                            nameConrol = "GlaB12A_", nameExperiment = "GlaT1.5_")

DECSE_resultsGL_6 <- dif_express_function(DEC = "DECSE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_6_GL.csv", 
                                          nameConrol = "GlaT12_", nameExperiment = "GlaB12_")

DECSE_resultsGL_18 <- dif_express_function(DEC = "DECSE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_18_GL.csv", 
                                           nameConrol = "GlaT12_", nameExperiment = "GlaT18_")

DECSE_resultsGL_24 <- dif_express_function(DEC = "DECSE", "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", Gla, "Sample_data_24_GL.csv", 
                                           nameConrol = "GlaT12_", nameExperiment = "GlaT24_")

###DOWN and UP by DECRE GL####

DECRE_Down_GL_1.5 <- filterlfpv(ddsdata = DECRE_resultsGL_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

DECRE_Up_GL_1.5 <- filterlfpv(ddsdata = DECRE_resultsGL_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

DECRE_Down_GL_6 <- filterlfpv(ddsdata = DECRE_resultsGL_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

DECRE_Up_GL_6 <- filterlfpv(ddsdata = DECRE_resultsGL_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

DECRE_Down_GL_18 <- filterlfpv(ddsdata = DECRE_resultsGL_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

DECRE_Up_GL_18 <- filterlfpv(ddsdata = DECRE_resultsGL_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

DECRE_Down_GL_24 <- filterlfpv(ddsdata = DECRE_resultsGL_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

DECRE_Up_GL_24 <- filterlfpv(ddsdata = DECRE_resultsGL_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECRE")

###DOWN and UP by DECPE GL####

DECPE_Down_GL_1.5 <- filterlfpv(ddsdata = DECPE_resultsGL_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

DECPE_Up_GL_1.5 <- filterlfpv(ddsdata = DECPE_resultsGL_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

DECPE_Down_GL_6 <- filterlfpv(ddsdata = DECPE_resultsGL_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

DECPE_Up_GL_6 <- filterlfpv(ddsdata = DECPE_resultsGL_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

DECPE_Down_GL_18 <- filterlfpv(ddsdata = DECPE_resultsGL_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

DECPE_Up_GL_18 <- filterlfpv(ddsdata = DECPE_resultsGL_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

DECPE_Down_GL_24 <- filterlfpv(ddsdata = DECPE_resultsGL_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

DECPE_Up_GL_24 <- filterlfpv(ddsdata = DECPE_resultsGL_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECPE")

###DOWN and UP by DECSE GL####

DECSE_Down_GL_1.5 <- filterlfpv(ddsdata = DECSE_resultsGL_1.5, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                                table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")

DECSE_Up_GL_1.5 <- filterlfpv(ddsdata = DECSE_resultsGL_1.5, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                              table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")

DECSE_Down_GL_6 <- filterlfpv(ddsdata = DECSE_resultsGL_6, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                              table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")

DECSE_Up_GL_6 <- filterlfpv(ddsdata = DECSE_resultsGL_6, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                            table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")

DECSE_Down_GL_18 <- filterlfpv(ddsdata = DECSE_resultsGL_18, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")

DECSE_Up_GL_18 <- filterlfpv(ddsdata = DECSE_resultsGL_18, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")

DECSE_Down_GL_24 <- filterlfpv(ddsdata = DECSE_resultsGL_24, logFoldchange = 1, gen_type = "Down", pvalue = 0.001, 
                               table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")

DECSE_Up_GL_24 <- filterlfpv(ddsdata = DECSE_resultsGL_24, logFoldchange = 1, gen_type = "Up", pvalue = 0.001, 
                             table_for_transcript = "09-01-2025_GLonlynMATRIX4numbers_GL1.csv", type_analisis = "DECSE")


#Уникальные, значимо изменяющиеся гены DECRE DECPE DECSE EC

DECRE_all_down_GL <- rbind(DECRE_Down_GL_1.5, DECRE_Down_GL_6, DECRE_Down_GL_18, DECRE_Down_GL_24)

DECRE_all_Up_GL <- rbind(DECRE_Up_GL_1.5, DECRE_Up_GL_6, DECRE_Up_GL_18, DECRE_Up_GL_24)

uniq_DECRE_down_GL <- UniqECnumber(DECRE_all_down_GL, "Down", type_analisis = "DECRE")

uniq_DECRE_Up_GL <- UniqECnumber(DECRE_all_Up_GL, "Up", type_analisis = "DECRE")

DECPE_all_down_GL <- rbind(DECPE_Down_GL_1.5, DECPE_Down_GL_6, DECPE_Down_GL_18, DECPE_Down_GL_24)

DECPE_all_Up_GL <- rbind(DECPE_Up_GL_1.5, DECPE_Up_GL_6, DECPE_Up_GL_18, DECPE_Up_GL_24)

uniq_DECPE_down_GL <- UniqECnumber(DECPE_all_down_GL, "Down", type_analisis = "DECPE")

uniq_DECPE_Up_GL <- UniqECnumber(DECPE_all_Up_GL, "Up", type_analisis = "DECPE")

DECSE_all_down_GL <- rbind(DECSE_Down_GL_1.5, DECSE_Down_GL_6, DECSE_Down_GL_18, DECSE_Down_GL_24)

DECSE_all_Up_GL <- rbind(DECSE_Up_GL_1.5, DECSE_Up_GL_6, DECSE_Up_GL_18, DECSE_Up_GL_24)

uniq_DECSE_down_GL <- UniqECnumber(DECSE_all_down_GL, "Down", type_analisis = "DECSE")

uniq_DECSE_Up_GL <- UniqECnumber(DECSE_all_Up_GL, "Up", type_analisis = "DECSE")


###function for search DECPE and DECSE EC numbers !IN DECRE###
DECnotDECRE <- function(DECRE_sp, DEC_PE_or_SE_sp, length) {
  library(operators)
  if (length == FALSE | length == F) {
    a <- unique(DEC_PE_or_SE_sp)[unique(DEC_PE_or_SE_sp) %!in% unique(DECRE_sp)]
  } else if (length == TRUE | length == T) {
    b <- as.numeric(length(unique(DEC_PE_or_SE_sp)[unique(DEC_PE_or_SE_sp) %!in% unique(DECRE_sp)]))
  } else {print("Ancorrected parametrs")}
}
#DECPE and DECSE !in DECRE for EV#
DECREinDECPE_down_EV <- DECnotDECRE(uniq_DECRE_down_EV, uniq_DECPE_down_EV, length = T)
DECREinDECPE_Up_EV <- DECnotDECRE(uniq_DECRE_Up_EV, uniq_DECPE_Up_EV, length = T)

DECREinDECSE_down_EV <- DECnotDECRE(uniq_DECRE_down_EV, uniq_DECSE_down_EV, length = T)
DECREinDECSE_Up_EV <- DECnotDECRE(uniq_DECRE_Up_EV, uniq_DECSE_Up_EV, length = T)
#DECPE and DECSE !in DECRE for EC#
DECREinDECPE_down_EC <- DECnotDECRE(uniq_DECRE_down_EC, uniq_DECPE_down_EC, length = T)
DECREinDECPE_Up_EC <- DECnotDECRE(uniq_DECRE_Up_EC, uniq_DECPE_Up_EC, length = T)

DECREinDECSE_down_EC <- DECnotDECRE(uniq_DECRE_down_EC, uniq_DECSE_down_EC, length = T)
DECREinDECSE_Up_EC <- DECnotDECRE(uniq_DECRE_Up_EC, uniq_DECSE_Up_EC, length = T)
#DECPE and DECSE !in DECRE for GL#
DECREinDECPE_down_GL <- DECnotDECRE(uniq_DECRE_down_GL, uniq_DECPE_down_GL, length = T)
DECREinDECPE_Up_GL <- DECnotDECRE(uniq_DECRE_down_GL, uniq_DECPE_down_GL, length = T)

DECREinDECSE_down_GL <- DECnotDECRE(uniq_DECRE_down_GL, uniq_DECSE_down_GL, length = T)
DECREinDECSE_Up_GL <- DECnotDECRE(uniq_DECRE_Up_GL, uniq_DECSE_Up_GL, length = T)

#df for figure 3#
Species <- c("E.verrucosus", "E.cyaneus", "G.lacustris")

Up_DECRE <- c(length(uniq_DECRE_Up_EV), length(uniq_DECRE_Up_EC), length(uniq_DECRE_Up_GL))
Down_DECRE <- c(length(uniq_DECRE_down_EV), length(uniq_DECRE_down_EC), length(uniq_DECRE_down_GL))

Up_DECPE <- c(length(uniq_DECPE_Up_EV), length(uniq_DECPE_Up_EC), length(uniq_DECPE_Up_GL))
Down_DECPE <- c(length(uniq_DECPE_down_EV), length(uniq_DECPE_down_EC), length(uniq_DECPE_down_GL))

Up_DECSE <- c(length(uniq_DECSE_Up_EV), length(uniq_DECSE_Up_EC), length(uniq_DECSE_Up_GL))
Down_DECSE <- c(length(uniq_DECSE_down_EV), length(uniq_DECSE_down_EC), length(uniq_DECSE_down_GL))

Up_uniq_DECPE <- c(DECREinDECPE_Up_EV, DECREinDECPE_Up_EC, DECREinDECPE_Up_GL)
Down_uniq_DECPE <- c(DECREinDECPE_down_EV, DECREinDECPE_down_EC, DECREinDECPE_down_GL)

Up_uniq_DECSE <- c(DECREinDECSE_Up_EV, DECREinDECSE_Up_EC, DECREinDECSE_Up_GL)
Down_uniq_DECSE <- c(DECREinDECSE_down_EV, DECREinDECSE_down_EC, DECREinDECSE_down_GL)

df_DECRE <- data.frame(Species, Up_DECRE, Down_DECRE)
df_DECPE <- data.frame(Species, Up_DECPE, Down_DECPE, Up_uniq_DECPE, Down_uniq_DECPE)
df_DECSE <- data.frame(Species, Up_DECSE, Down_DECSE, Up_uniq_DECSE, Down_uniq_DECSE)

##figure 3##

library(ggplot2)
cvet <- c("#56B4E9", "#009E73", "#E69F00")
cvet1 <- c( "#009E73", "#56B4E9", "#E69F00")

E <- ggplot(df_DECRE, aes(x = fct_inorder(Species), fill=Species)) + 
  geom_bar(aes(y=Up_DECRE), position = "dodge", stat="identity", width = 0.5, col = cvet1, alpha = .5)+ 
  geom_bar(aes(y=-Down_DECRE), position = "dodge", stat="identity", width = 0.5,  col = cvet1, alpha = .5)+ 
  geom_hline(aes(yintercept = 0), col = "white") +
  ylab(NULL)+
  xlab("")+ 
  scale_y_continuous(limits = c(-700, 300), breaks = c(-700, -500, -300, -100, 0, 100, 300), labels = abs) + 
  theme_classic(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=35, family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = "none",
        legend.direction = "horizontal") +
  labs(fill = NULL) +
  scale_fill_manual(values = cvet, labels = c(substitute(paste(italic('E. cyaneus'))),
                                              substitute(paste(italic('E. verrucosus'))), 
                                              substitute(paste(italic('G. lacustris'))))) + 
  ggtitle("DECRE") +
  theme(plot.margin = margin(t=5,1,1,1, "lines"))  +
  theme(plot.title = element_text(vjust = 3.5))



P <- ggplot(df_DECPE, aes(x = fct_inorder(Species), fill=Species)) +   
  geom_bar(aes(y=Up_DECPE), position = "dodge", stat="identity", width = 0.5, col = cvet1, alpha = .5)+ 
  geom_bar(aes(y=-Down_DECPE), position = "dodge", stat="identity", width = 0.5,  col = cvet1, alpha = .5)+
  geom_bar(aes(y=Up_uniq_DECPE), position = "dodge", stat="identity", width = 0.5, col = cvet1, alpha = .75)+ 
  geom_bar(aes(y=-Down_uniq_DECPE), position = "dodge", stat="identity", width = 0.5,  col = cvet1, alpha = .75)+
  geom_hline(aes(yintercept = 0), col = "white") +
  ylab(NULL)+
  xlab("")+ 
  scale_y_continuous(limits = c(-700, 300), breaks = c(-700, -500, -300, -100, 0, 100, 300), labels = abs) + 
  theme_classic(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=35, family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = "none",
        legend.direction = "horizontal") +
  labs(fill = NULL) +
  scale_fill_manual(values = cvet, labels = c(substitute(paste(italic('E. cyaneus'))),
                                              substitute(paste(italic('E. verrucosus'))), 
                                              substitute(paste(italic('G. lacustris'))))) + 
  ggtitle("DECPE") +
  theme(plot.margin = margin(t=5,1,1,1, "lines"))  +
  theme(plot.title = element_text(vjust = 3.5)) 

cvet <- c( "#009E73", "#56B4E9", "#E69F00")

cvet1 <- c( "#009E73", "#56B4E9", "#E69F00")

S <- ggplot(df_DECSE, aes(x = fct_inorder(Species), fill=fct_inorder(Species))) +   
  geom_bar(aes(y=Up_DECSE), position = "dodge", stat="identity", width = 0.5, col = cvet1, alpha = .5)+ 
  geom_bar(aes(y=-Down_DECSE), position = "dodge", stat="identity", width = 0.5,  col = cvet1, alpha = .5)+
  geom_bar(aes(y=Up_uniq_DECSE), position = "dodge", stat="identity", width = 0.5, col = cvet1, alpha = .75)+ 
  geom_bar(aes(y=-Down_uniq_DECSE), position = "dodge", stat="identity", width = 0.5,  col = cvet1, alpha = .75)+
  geom_hline(aes(yintercept = 0), col = "white") +
  ylab(NULL)+
  xlab("")+ 
  scale_y_continuous(limits = c(-700, 300), breaks = c(-700, -500, -300, -100, 0, 100, 300), labels = abs) + 
  theme_classic(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=35, family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = c(-0.7, 1.3),
        legend.direction = "horizontal") +
  labs(fill = NULL) +
  scale_fill_manual(values = cvet, breaks = c("E.verrucosus", "E.cyaneus", "G.lacustris"), 
                    labels = c(substitute(paste(italic('E. verrucosus'))),
                               substitute(paste(italic('E. cyaneus'))), 
                               substitute(paste(italic('G. lacustris'))))) + 
  ggtitle("DECSE") +
  theme(plot.margin = margin(t=5,1,1,1, "lines"))  +
  theme(plot.title = element_text(vjust = 3.5))



library(cowplot)
S
library(patchwork)

E + P + S

plot_grid(E, P, S, nrow = 1)

library(ggplot2)

ggsave("Figure_3_13-01-2025.png", width = 43, height = 19, units = "cm", dpi = 400)

#venn diagramm DECRE#
EVall_DECRE <- rbind(DECRE_all_down_EV, DECRE_all_Up_EV)
ECall_DECRE <- rbind(DECRE_all_down_EC, DECRE_all_Up_EC)
GLall_DECRE <- rbind(DECRE_all_down_GL, DECRE_all_Up_GL)

EVall_DECRE <- as.data.frame(unique(EVall_DECRE$EC.number))
ECall_DECRE <- as.data.frame(unique(ECall_DECRE$EC.number))
GLall_DECRE <- as.data.frame(unique(GLall_DECRE$EC.number))

colnames(EVall_DECRE)[1] <- "E.verrucosus"
colnames(ECall_DECRE)[1] <- "E.cyaneus"
colnames(GLall_DECRE)[1] <- "G.lacustis"

library(BioVenn)
svg("DECREvennall.png", width = 10, height = 10)
draw.venn(EVall_DECRE$E.verrucosus, ECall_DECRE$E.cyaneus, GLall_DECRE$G.lacustis,
          title = paste0(""),
          subtitle = "", 
          xtitle = "", ytitle = "", ztitle = "", 
          x_c = "#009E73", y_c="#56B4E9", z_c="#E69F00", 
          xt_f = "Arial", yt_f = "Arial", zt_f = "Arial", 
          xt_fb = 2, yt_fb = 1, zt_fb = 1, 
          xt_s = 1.2, yt_s = 1.2, zt_s = 1.2, nr_s = 3,
          nr_f = "Arial", t_f = "Arial", t_fb = 2, bg_c = "white",
          st_fb = 1, st_f = "Arial")

#svg("Vennd.svg")
dev.off()

#venn diagramm DECPE#

EVall_DECPE <- rbind(DECPE_all_down_EV, DECPE_all_Up_EV)
ECall_DECPE <- rbind(DECPE_all_down_EC, DECPE_all_Up_EC)
GLall_DECPE <- rbind(DECPE_all_down_GL, DECPE_all_Up_GL)

EVall_DECPE <- as.data.frame(unique(EVall_DECPE$EC.number))
ECall_DECPE <- as.data.frame(unique(ECall_DECPE$EC.number))
GLall_DECPE <- as.data.frame(unique(GLall_DECPE$EC.number))

colnames(EVall_DECPE)[1] <- "E.verrucosus"
colnames(ECall_DECPE)[1] <- "E.cyaneus"
colnames(GLall_DECPE)[1] <- "G.lacustis"

library(BioVenn)
svg("DECPEvennall.png", width = 10, height = 10)
draw.venn(EVall_DECPE$E.verrucosus, ECall_DECPE$E.cyaneus, GLall_DECPE$G.lacustis,
          title = paste0(""),
          subtitle = "", 
          xtitle = "", ytitle = "", ztitle = "", 
          x_c = "#009E73", y_c="#56B4E9", z_c="#E69F00", 
          xt_f = "Arial", yt_f = "Arial", zt_f = "Arial", 
          xt_fb = 2, yt_fb = 1, zt_fb = 1, 
          xt_s = 1.2, yt_s = 1.2, zt_s = 1.2, nr_s = 3,
          nr_f = "Arial", t_f = "Arial", t_fb = 2, bg_c = "white",
          st_fb = 1, st_f = "Arial")

#svg("Vennd.svg")
dev.off()


#venn diagramm DECSE#

EVall_DECSE <- rbind(DECSE_all_down_EV, DECSE_all_Up_EV)
ECall_DECSE <- rbind(DECSE_all_down_EC, DECSE_all_Up_EC)
GLall_DECSE <- rbind(DECSE_all_down_GL, DECSE_all_Up_GL)

EVall_DECSE <- as.data.frame(unique(EVall_DECSE$row))
ECall_DECSE <- as.data.frame(unique(ECall_DECSE$row))
GLall_DECSE <- as.data.frame(unique(GLall_DECSE$row))

colnames(EVall_DECSE)[1] <- "E.verrucosus"
colnames(ECall_DECSE)[1] <- "E.cyaneus"
colnames(GLall_DECSE)[1] <- "G.lacustis"

library(BioVenn)
svg("DECSEvennall.png", width = 10, height = 10)
draw.venn(EVall_DECSE$E.verrucosus, ECall_DECSE$E.cyaneus, GLall_DECSE$G.lacustis,
          title = paste0(""),
          subtitle = "", 
          xtitle = "", ytitle = "", ztitle = "", 
          x_c = "#009E73", y_c="#56B4E9", z_c="#E69F00", 
          xt_f = "Arial", yt_f = "Arial", zt_f = "Arial", 
          xt_fb = 2, yt_fb = 1, zt_fb = 1, 
          xt_s = 1.2, yt_s = 1.2, zt_s = 1.2, nr_s = 3,
          nr_f = "Arial", t_f = "Arial", t_fb = 2, bg_c = "white",
          st_fb = 1, st_f = "Arial")

#svg("Vennd.svg")
dev.off()


##Figure 1###

library(xlsx)
library(ggplot2)
library(ggh4x)
library(patchwork)

exp2 <- read.xlsx("fig1.xlsx", 1)
exp1 <- read.xlsx("fig1.xlsx", 2)

quantile(rnorm(100), probs = seq(0, 1, .1))


df12exp2 <- exp2[  exp2$pointy == 12, ]
p1 <-
  ggplot(exp2, aes(x = time, y = temper, col = g1)) + geom_line(linewidth = .8) + 
  geom_point(aes(x = pointx, y = pointy), col = "#2F4F4F", shape = 13, size = 4) +
  geom_point(aes(x = pointxred, y = pointyred), col = "#2F4F4F", shape = 19, size = 2) +
  geom_point(aes(x = pointxred, y = pointyred), col = "#2F4F4F", shape = 13, size = 4) +
  scale_y_continuous(limits = c(0,25), breaks = c(0, 3, 6, 9, 12)) + 
  scale_color_manual(values = c("#50abfb", "#f86b03e5", "#50abfb", "#f86b03e5","#708090")) + 
  theme_classic() + 
  ylab("Temperature, °C") +
  xlab("Time, days") + 
  theme(legend.position = "none") +
  theme(text=element_text(color="black"),axis.text=element_text(color="black")) + 
  theme(text=element_text(size=14, family="Arial")) +
  annotate("text", x = c(7.7), y = c(22), label = substitute(paste(bold("Experiment #2"))), cex=4, col = "black") +
  annotate("text", x = c(7.5), y = c(8), label = substitute(paste("Pre-acclimation")), cex=3, col = "black") + 
  ylim(0, 25)



p2 <-
  ggplot(exp1, aes(x = time, y = temper, col = g1)) + geom_line(linewidth = .8) +
  geom_point(aes(x = pointx, y = pointy), col = "#2F4F4F", shape = 13, size = 4) +
  geom_point(aes(x = pointxred, y = pointyred), col = "#2F4F4F", shape = 19, size = 2) +
  geom_point(aes(x = pointxred, y = pointyred), col = "#2F4F4F", shape = 13, size = 4) +
  scale_y_continuous(limits = c(0,25), breaks = c(0, 6, 12, 18, 24)) +  
  scale_x_continuous(limits = c(0,25), breaks = c(0, 6, 12, 18, 24))  + 
  scale_color_manual(labels = c( "Decreased", "Increased"), values = c("#50abfb", "#f86b03e5")) + 
  theme_classic() + 
  ylab("Temperature, °C") +
  xlab("Time, days") +
  theme(legend.position = "none",
        legend.direction = "horizontal") + 
  labs(col = NULL) +
  theme(text=element_text(color="black"),axis.text=element_text(color="black")) + 
  theme(text=element_text(size=14, family="Arial")) +
  annotate("text", x = c(4), y = c(22), label = substitute(paste(bold("Experiment #1"))), cex=4, col = "black") +
  theme(plot.margin = margin(t=1.5,1,1,1, "lines")) + ylim(0, 25)


p2 / p1

ggsave("fig1_part3.png", width = 13, height = 14, units = "cm", dpi = 600)



#figure 2#
#enzome EV#
library(plyr)
library(ggplot2)
#download start data
startdataEV <- read.csv("Filtered DataEV.csv")
startdataEV <- startdataEV[, c(-1, -3, -5, -6)]
startdataEV <- unique(startdataEV)
#delete all value after ";" 
startdataEV$EC.number <- gsub(";.*","\\1",startdataEV$EC.number)
#remove colon startdate
removeallcolomnGL <- startdataEV[, c(-2, -3, -5, -6)]
#create NA numbers
startdataEV["EC.number"][startdataEV["EC.number"] == ''] <- NA
#delete NA values
deleteNAvalue <- startdataEV[rowSums(is.na(startdataEV)) == 0,]
#ploting with four (1.1.1.1) nomber classification
plotingwithfour <- deleteNAvalue[!grepl("-", deleteNAvalue$EC.number),]
#count duplics
plotingwithfour_countduplics <- ddply(plotingwithfour,.(EC.number),nrow)
#created plot
p1 <- ggplot(plotingwithfour_countduplics, aes(x=V1)) + 
  geom_histogram(bins = 200, fill="#009E73", alpha = 0.7, color="#063303", boundary=1, size = 0.3) + 
  xlim(1, 100) +
  geom_vline(aes(xintercept=median(V1)),
             color="black", linetype="dashed", size=1) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL)  +
  scale_x_continuous(limits = c(1, 50), breaks = c(0, 9, 25, 50, 75, 100), labels = abs) +
  theme(text=element_text(size=28, family="Arial"), 
        axis.title.y = element_text(size = 28)) + 
  theme(text=element_text(color="black"),axis.text=element_text(color="black"))

#created plot
p2 <- ggplot(plotingwithfour_countduplics, aes(x=V1)) + 
  geom_histogram(bins = 185, fill="#009E73", alpha = 0.7, color="#063303", boundary=1, size = 0.3) +
  geom_vline(aes(xintercept=median(V1)),
             color="black", linetype="dashed", size=1) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle(substitute(paste(italic('Eulimnogammarus verrucosus')))) +
  scale_x_continuous(limits = c(50, 1850), breaks = c(50, 350, 600, 850, 1100, 1350, 1600, 1850), 
                     labels = abs) + 
  theme(plot.title = element_text(hjust = -6.2)) +
  theme(text=element_text(size=28, family="Arial"),
        axis.title.y = element_text(vjust = +1.9)) +
  geom_segment(aes(x = 1725, y = 14, xend = 1825, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(1725), y = c(15), 
           label = substitute(paste(bold("EC 2.7.11.1"))), 
           cex=6, col = "black") +
  geom_segment(aes(x = 1206, y = 15, xend = 1206, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(1206), y = c(17), 
           label = substitute(paste(bold("EC 2.7.7.49"))), 
           cex=6, col = "black") +
  geom_segment(aes(x = 750, y = 14, xend = 1025, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(750), y = c(15), 
           label = substitute(paste(bold("EC 2.3.2.27"))), 
           cex=6, col = "black") +
  geom_segment(aes(x = 600, y = 8, xend = 746, yend = 1),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(550), y = c(10), 
           label = substitute(paste(bold("EC 3.6.4.13"))), 
           cex=6, col = "black") + ylim(0, 20) + 
  theme(text=element_text(color="black"),axis.text=element_text(color="black"))
p2
#enzome EC#
#download start data
startdataEC <- read.csv("Filtered DataEC.csv")
startdataEC <- startdataEC[, c(-1, -3, -5, -6)]
startdataEC <- unique(startdataEC)
#delete all value after ";" 
startdataEC$EC.number <- gsub(";.*","\\1",startdataEC$EC.number)
#remove colon startdate
removeallcolomnGL <- startdataEC[, c(-2, -3, -5, -6)]
#create NA numbers
startdataEC["EC.number"][startdataEC["EC.number"] == ''] <- NA
#delete NA values
deleteNAvalue <- startdataEC[rowSums(is.na(startdataEC)) == 0,]
#ploting with four (1.1.1.1) nomber classification
plotingwithfour <- deleteNAvalue[!grepl("-", deleteNAvalue$EC.number),]
#count duplics
plotingwithfour_countduplics <- ddply(plotingwithfour,.(EC.number),nrow)
#created plot
p3 <- ggplot(plotingwithfour_countduplics, aes(x=V1)) + geom_histogram(bins = 200, 
                                                                       fill="#56B4E9", 
                                                                       alpha = 0.7, 
                                                                       color="#0036d1", 
                                                                       boundary=1, 
                                                                       size = 0.3) +
  geom_vline(aes(xintercept=median(V1)),
             color="black", linetype="dashed", size=1) +
  theme_classic() +
  ylab("Number of enzymes") +
  xlab(NULL) +
  scale_x_continuous(limits = c(1, 50), breaks = c(0, 7, 25, 50, 75, 100), labels = abs) +
  theme(text=element_text(size=28, family="Arial")) +
  theme(axis.title.y = element_text(size = 28, vjust = 2)) + 
  theme(text=element_text(color="black"),axis.text=element_text(color="black"))

#created plot
p4 <- ggplot(plotingwithfour_countduplics, aes(x=V1)) + geom_histogram(bins = 185, fill="#56B4E9", 
                                                                       alpha = 0.7, color="#0036d1", 
                                                                       boundary=1, size = 0.3) + 
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  scale_x_continuous(limits = c(50, 1850), breaks = c(50, 350, 600, 850, 1100, 1350, 1600, 1850), 
                     labels = abs) +
  theme(text=element_text(size=28, family="Arial"))  +
  ggtitle(substitute(paste(italic('Eulimnogammarus cyaneus')))) + 
  theme(plot.title = element_text(hjust = -3.74)) +
  annotate("text", x = c(1600), y = c(15), 
           label = substitute(paste(bold("EC 2.7.11.1"))), 
           cex=6, col = "black") + 
  geom_segment(aes(x = 1600, y = 13, xend = 1236, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 350, y = 8, xend = 579, yend = 1.1),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(350), y = c(10), 
           label = substitute(paste(bold("EC 3.6.4.13"))), 
           cex=6, col = "black") +
  geom_segment(aes(x = 600, y = 13, xend = 752, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(600), y = c(15), 
           label = substitute(paste(bold("EC 2.3.2.27"))), 
           cex=6, col = "black")  +
  geom_segment(aes(x = 1100, y = 13, xend = 844, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(1100), y = c(15), 
           label = substitute(paste(bold("EC 2.7.7.49"))), 
           cex=6, col = "black") + ylim(0, 20) + 
  theme(text=element_text(color="black"),axis.text=element_text(color="black"))

#enzome GL#

#download start data
startdataGL <- read.csv("Filtered DataGL.csv")
startdataGL1 <- startdataGL[, c(-1, -3, -5, -6)]
startdataGL2 <- unique(startdataGL1)
#delete all value after ";" 
startdataGL2$EC.number <- gsub(";.*","\\1",startdataGL2$EC.number)
#remove colon startdate
removeallcolomnGL <- startdataGL[, c(-2, -3, -5, -6)]
#create NA numbers
startdataGL2["EC.number"][startdataGL2["EC.number"] == ''] <- NA
#delete NA values
deleteNAvalue <- startdataGL2[rowSums(is.na(startdataGL2)) == 0,]
#ploting with four (1.1.1.1) nomber classification
plotingwithfour <- deleteNAvalue[!grepl("-", deleteNAvalue$EC.number),]
#count duplics
plotingwithfour_countduplics <- ddply(plotingwithfour,.(EC.number),nrow)
#created plot
p5 <- 
  ggplot(plotingwithfour_countduplics, aes(x=V1)) + 
  geom_histogram(bins = 200, fill="#E69F00", alpha = 0.7, color="#82520a", boundary=1, size = 0.3) + 
  xlim(1, 100) + 
  geom_vline(aes(xintercept=median(V1)),
             color="black", linetype="dashed", size=1) + 
  theme_classic() +
  ylab(NULL) +
  xlab("Number of isoforms") +
  scale_x_continuous(limits = c(1, 50), breaks = c(0, 8, 25, 50, 75, 100), labels = abs) +
  theme(text=element_text(size=20, family="Arial"),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_text(size = 28)) + 
  theme(text=element_text(color="black", size = 28),axis.text=element_text(color="black"))

#created plot
p6 <- 
  ggplot(plotingwithfour_countduplics, aes(x=V1)) + 
  geom_histogram(bins = 185, fill="#E69F00", alpha = 0.7, color="#82520a", boundary=1, size = 0.3) +
  scale_x_continuous(limits = c(50, 1850), breaks = c(50, 350, 600, 850, 1100, 1350, 1600, 1850), labels = abs) +
  geom_vline(aes(xintercept=median(V1)),
             color="black", linetype="dashed", size=1) +
  theme_classic() +
  ylab(NULL) +
  xlab("Number of isoforms") +
  ggtitle(substitute(paste(italic("Gammarus lacustris")))) + 
  theme(plot.title = element_text(hjust = -1.74),
        axis.title.x = element_text(size = 28)) +
  theme(text=element_text(size=20, family="Arial")) +
  geom_segment(aes(x = 1660, y = 14, xend = 1660, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(1670), y = c(15), 
           label = substitute(paste(bold("EC 2.7.11.1"))), 
           cex=6, col = "black")  +
  geom_segment(aes(x = 1190, y = 11, xend = 1190, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(1190), y = c(12), 
           label = substitute(paste(bold("EC 2.7.7.49"))), 
           cex=6, col = "black") +
  geom_segment(aes(x = 725, y = 13, xend = 1050, yend = .8),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(725), y = c(15), 
           label = substitute(paste(bold("EC 2.3.2.27"))), 
           cex=6, col = "black") + ylim(0, 20) +
  geom_segment(aes(x = 600, y = 8, xend = 600, yend = 1.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(600), y = c(9), 
           label = substitute(paste(bold("EC 3.6.4.13"))), 
           cex=6, col = "black") + 
  theme(text=element_text(color="black", size = 28),axis.text=element_text(color="black"))


#merge figure 1##
library(patchwork)
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(nrow = 3, widths = c(1, 1.5))
ggsave("Figure 2_V2.png", width = 35, height = 33, units = "cm", dpi = 300)
ggsave("Figure 2_V2.svg", width = 35, height = 33, units = "cm", dpi = 300)


#figure 3###
#EV# Getting table percent_tableEV.csv

#EConlynMATRIX4numbers_EV1forvoolcano <- read.csv("30-01-2024_EConlynMATRIX4numbers_EV1.csv")
#EConlynMATRIX4numbers_EV1forvoolcano <- EConlynMATRIX4numbers_EV1forvoolcano[ , c(-1, -3)]
#colnames(EConlynMATRIX4numbers_EV1forvoolcano)[1] <- "gen"
#detach("package:plyr", unload = TRUE)
#library(dplyr)
#a1 <- EConlynMATRIX4numbers_EV1forvoolcano
##Процентный состав для всей выборки (каждая повторность)
#a2 <- a1 %>%
#  group_by(Temperature, EC.number) %>%
#  mutate(percent = Expression/sum(Expression) *100)
##запись файла с процентным вкладом
#write.csv(a2, "percent_tableEV.csv")


Evpersent <- read.csv("percent_tableEV.csv")
Evpersent <- Evpersent[ , -c(1, 5) ]

Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
library(tidyverse)
Evpersent2 <- Evpersent %>%
  group_by(Temperature, gen) %>%
  mutate(meanPercent = mean(percent, na.rm = TRUE))

Evpersent2 <- Evpersent2[ , -4]
Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]

sum(is.na(Evpersent3))

Evpersent4 <- na.omit(Evpersent3)

Evpersent5 <- Evpersent4[which(Evpersent4$Temperature %in% "EveT12"),]

write.csv(Evpersent5, "ECpersent5.csv")
Evpersent6 <- read.csv("ECpersent5.csv")
library(ggplot2)
ggplot(Evpersent6, aes(x = X, y = meanPercent)) + geom_point()

Evpersent7 <- Evpersent6[which(Evpersent6$meanPercent >= 50),]
Evpersent8 <- Evpersent6[which(Evpersent6$meanPercent < 50),]

uEV7 <- unique(Evpersent7$EC.number)
uEV8 <- unique(Evpersent8$EC.number)


str(Evpersent6)
Evpersent6$EC.number
Evpersent6$X <- as.numeric(Evpersent6$X)

library(tidyverse)
Evpersent6_1 <- Evpersent6 %>%
  arrange(desc(EC.number)) 

Evpersent6_1$number <- 1:23013
ggplot(Evpersent6_1, aes(x = number, y = meanPercent)) + geom_point()
Evpersent7$number <- 1:750
ggplot(Evpersent7, aes(x = number, y = meanPercent)) + geom_point()
Evpersent8$number <- 1:22263
ggplot(Evpersent8, aes(x = number, y = meanPercent)) + geom_point()

uEV7 <- unique(Evpersent7$EC.number)
uEV8 <- unique(Evpersent8$EC.number)


library(plyr)
library(cowplot)
Evpersent9 <- ddply(Evpersent6,.(Evpersent6$meanPercent >= 50),nrow)
Evpersent10 <- ddply(Evpersent7,.(EC.number) , nrow)
Evpersent11 <- ddply(Evpersent6,.(EC.number) , nrow)
detach("package:plyr", unload = TRUE)
Evpersent11$number <- 1:1207
EV12 <- merge(Evpersent11, Evpersent6, by = "EC.number")
ggplot(EV12, aes(x = V1, y = meanPercent)) + geom_point()
library(dplyr, lib.loc = "/usr/lib/R/site-library")
library(operators)
Evpersent14 <- EV12[which(EV12$V1 %!in% "1"),]
detach("package:operators", unload = TRUE)
Evpersent14 <- Evpersent14 %>%
  group_by(EC.number) %>%
  mutate(max = max(meanPercent))
Evpersent14 <- Evpersent14[ , -c(7, 5, 4)]
Evpersent14 <- Evpersent14[!duplicated(Evpersent14), ]
Evpersent14$median <- ave(Evpersent14$max, Evpersent14$V1,FUN=function(y) median(y,na.rm=TRUE))
Evpersent14$ideal <- 100/Evpersent14$V1
p1 <- ggplot(Evpersent14, aes(x = V1, y = max)) + 
  geom_point(col = "#31495a", fill = "#009E73", pch = 21, size = 2.5, alpha = .8) + 
  geom_point(aes(y=median),size=7,colour="#C51D34", pch = 4) +  
  geom_line(aes(y = ideal), col = "#31495a", linetype = "dashed", size = 1.2) +
  scale_x_continuous(trans='log', breaks = c(2, 3, 4, 5, 10, 15, 20, 40, 60, 
                                             80, 120, 175, 250, 400, 700, 1100)) +
  geom_smooth(method = "loess", col = "#C51D34", size = 1.2)  +
  theme_classic() +
  xlab("Number of isoforms per enzyme") + 
  ylab(NULL) + 
  ggtitle(substitute(paste(italic("Eulimnogammarus verrucosus")))) + 
  theme(text=element_text(size=17, family="Arial", colour = "black"),
        axis.text = element_text(size=13, family="Arial", colour = "black"))

p1
#EC# Getting table percent_tableEC.csv

#EConlynMATRIX4numbers_EC1forvoolcano <- read.csv("30-01-2024_EConlynMATRIX4numbers_EC1.csv")
#EConlynMATRIX4numbers_EC1forvoolcano <- EConlynMATRIX4numbers_EC1forvoolcano[ , c(-1, -3)]
#colnames(EConlynMATRIX4numbers_EC1forvoolcano)[1] <- "gen"
#detach("package:plyr", unload = TRUE)
#library(dplyr)
#a1 <- EConlynMATRIX4numbers_EC1forvoolcano
##Процентный состав для всей выборки (каждая повторность)
#a2 <- a1 %>%
#  group_by(Temperature, EC.number) %>%
#  mutate(percent = Expression/sum(Expression) *100)
##запись файла с процентным вкладом
#write.csv(a2, "percent_tableEC.csv")


Evpersent <- read.csv("percent_tableEC.csv")
Evpersent <- Evpersent[ , -c(1, 5) ]

Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
library(tidyverse)
Evpersent2 <- Evpersent %>%
  group_by(Temperature, gen) %>%
  mutate(meanPercent = mean(percent, na.rm = TRUE))

Evpersent2 <- Evpersent2[ , -4]
Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]

sum(is.na(Evpersent3))

Evpersent4 <- na.omit(Evpersent3)

Evpersent5 <- Evpersent4[which(Evpersent4$Temperature %in% "EcyT12"),]

write.csv(Evpersent5, "ECpersent5.csv")
Evpersent6 <- read.csv("ECpersent5.csv")
library(ggplot2)
ggplot(Evpersent6, aes(x = X, y = meanPercent)) + geom_point()

Evpersent7 <- Evpersent6[which(Evpersent6$meanPercent >= 50),]
Evpersent8 <- Evpersent6[which(Evpersent6$meanPercent < 50),]

uEV7 <- unique(Evpersent7$EC.number)
uEV8 <- unique(Evpersent8$EC.number)


str(Evpersent6)
Evpersent6$EC.number
Evpersent6$X <- as.numeric(Evpersent6$X)

library(tidyverse)
Evpersent6_1 <- Evpersent6 %>%
  arrange(desc(EC.number)) 

Evpersent6_1$number <- 1:16888
ggplot(Evpersent6_1, aes(x = number, y = meanPercent)) + geom_point()
Evpersent7$number <- 1:859
ggplot(Evpersent7, aes(x = number, y = meanPercent)) + geom_point()
Evpersent8$number <- 1:16029
ggplot(Evpersent8, aes(x = number, y = meanPercent)) + geom_point()


library(plyr)
Evpersent9 <- ddply(Evpersent6,.(Evpersent6$meanPercent >= 50),nrow)
Evpersent10 <- ddply(Evpersent7,.(EC.number) , nrow)
Evpersent11 <- ddply(Evpersent6,.(EC.number) , nrow)
detach("package:plyr", unload = TRUE)
Evpersent11$number <- 1:1149
EV12 <- merge(Evpersent11, Evpersent6, by = "EC.number")
ggplot(EV12, aes(x = V1, y = meanPercent)) + geom_point()
library(dplyr, lib.loc = "/usr/lib/R/site-library")
library(operators)
Evpersent14 <- EV12[which(EV12$V1 %!in% "1"),]
detach("package:operators", unload = TRUE)
Evpersent14 <- Evpersent14 %>%
  group_by(EC.number) %>%
  mutate(max = max(meanPercent))
Evpersent14 <- Evpersent14[ , -c(7, 5, 4)]
Evpersent14 <- Evpersent14[!duplicated(Evpersent14), ]
Evpersent14$median <- ave(Evpersent14$max, Evpersent14$V1,FUN=function(y) median(y,na.rm=TRUE))
Evpersent14$ideal <- 100/Evpersent14$V1
p2 <- ggplot(Evpersent14, aes(x = V1, y = max)) + 
  geom_point(col = "#31495a", fill = "#56B4E9", pch = 21, size = 2.5, alpha = .8) + 
  geom_point(aes(y=median),size=7,colour="#C51D34", pch = 4) +   
  geom_line(aes(y = ideal), col = "#31495a", linetype = "dashed", size = 1.2) +
  scale_x_continuous(trans='log', breaks = c(2, 3, 4, 5, 10, 15, 20, 40, 60, 
                                             80, 120, 175, 250, 400, 700, 1100)) +
  theme_classic() +
  geom_smooth(method = "loess", col = "#C51D34", size = 1.2) +
  xlab("Number of isoforms per enzyme") + 
  ylab(NULL) + 
  ggtitle(substitute(paste(italic("Eulimnogammarus cyaneus")))) + 
  theme(text=element_text(size=17, family="Arial"),
        axis.text = element_text(size=13, family="Arial", colour = "black"))
p2


#GL# Getting table percent_tableGL.csv

#EConlynMATRIX4numbers_EC1forvoolcano <- read.csv("09-01-2025_GLonlynMATRIX4numbers_GL1.csv")
#EConlynMATRIX4numbers_EC1forvoolcano <- EConlynMATRIX4numbers_EC1forvoolcano[ , c(-1, -3)]
#colnames(EConlynMATRIX4numbers_EC1forvoolcano)[1] <- "gen"
#detach("package:plyr", unload = TRUE)
#library(dplyr)
#a1 <- EConlynMATRIX4numbers_EC1forvoolcano
##Процентный состав для всей выборки (каждая повторность)
#a2 <- a1 %>%
#  group_by(Temperature, EC.number) %>%
#  mutate(percent = Expression/sum(Expression) *100)
##запись файла с процентным вкладом
#write.csv(a2, "percent_tableGL.csv")


Evpersent <- read.csv("percent_tableGL.csv")
Evpersent <- Evpersent[ , -c(1, 5) ]

Evpersent$Temperature <- gsub("\\_.*","",Evpersent$Temperature)
library(tidyverse)
Evpersent2 <- Evpersent %>%
  group_by(Temperature, gen) %>%
  mutate(meanPercent = mean(percent))

Evpersent2 <- Evpersent2[ , -4]
Evpersent3 <- Evpersent2[!duplicated(Evpersent2), ]

sum(is.na(Evpersent3))

Evpersent4 <- na.omit(Evpersent3)

Evpersent5 <- Evpersent4[which(Evpersent4$Temperature %in% "GlaT12"),]

write.csv(Evpersent5, "glpersent5.csv")
Evpersent6 <- read.csv("glpersent5.csv")
library(ggplot2)
ggplot(Evpersent6, aes(x = X, y = meanPercent)) + geom_point()

Evpersent7 <- Evpersent6[which(Evpersent6$meanPercent >= 50),]
Evpersent8 <- Evpersent6[which(Evpersent6$meanPercent < 50),]


uEV7 <- unique(Evpersent7$EC.number)
uEV8 <- unique(Evpersent8$EC.number)


str(Evpersent6)
Evpersent6$EC.number
Evpersent6$X <- as.numeric(Evpersent6$X)

library(tidyverse)
Evpersent6_1 <- Evpersent6 %>%
  arrange(desc(EC.number)) 

Evpersent6_1$number <- 1:18354
ggplot(Evpersent6_1, aes(x = number, y = meanPercent)) + geom_point()
Evpersent7$number <- 1:681
ggplot(Evpersent7, aes(x = number, y = meanPercent)) + geom_point()
Evpersent8$number <- 1:17673
ggplot(Evpersent8, aes(x = number, y = meanPercent)) + geom_point()


library(plyr)
Evpersent9 <- ddply(Evpersent6,.(Evpersent6$meanPercent >= 50),nrow)
Evpersent10 <- ddply(Evpersent7,.(EC.number) , nrow)
Evpersent11 <- ddply(Evpersent6,.(EC.number) , nrow)
detach("package:plyr", unload = TRUE)
Evpersent11$number <- 1:1099
EV12 <- merge(Evpersent11, Evpersent6, by = "EC.number")
ggplot(EV12, aes(x = V1, y = meanPercent)) + geom_point()
library(dplyr, lib.loc = "/usr/lib/R/site-library")
library(operators)
Evpersent14 <- EV12[which(EV12$V1 %!in% "1"),]
detach("package:operators", unload = TRUE)
Evpersent14 <- Evpersent14 %>%
  group_by(EC.number) %>%
  mutate(max = max(meanPercent))
Evpersent14 <- Evpersent14[ , -c(7, 5, 4)]
Evpersent14 <- Evpersent14[!duplicated(Evpersent14), ]
Evpersent14$median <- ave(Evpersent14$max, Evpersent14$V1,FUN=function(y) median(y,na.rm=TRUE))
Evpersent14$ideal <- 100/Evpersent14$V1
p3 <- ggplot(Evpersent14, aes(x = V1, y = max)) + 
  geom_point(col = "#31495a", fill = "#E69F00", pch = 21, size = 2.5, alpha = .8) + 
  geom_point(aes(y=median),size=7,colour="#C51D34", pch = 4) +    
  geom_line(aes(y = ideal), col = "#31495a", linetype = "dashed", size = 1.2) +
  scale_x_continuous(trans='log', breaks = c(2, 3, 4, 5, 10, 15, 20, 40, 60, 
                                             80, 120, 175, 250, 400, 700, 1100)) +
  geom_smooth(method = "loess", col = "#C51D34", size = 1.2) +
  theme_classic() +
  xlab("Number of isoforms per enzyme") + 
  ylab("Contribution of leading isoform, % of total enzyme expression") + 
  ggtitle(substitute(paste(italic("Gammarus lacustris")))) + 
  theme(text=element_text(size=17, family="Arial"))  + 
  theme(axis.title.y = element_text(hjust = -.7),
        axis.text = element_text(size=13, family="Arial", colour = "black"))
p3


library(patchwork)

p1 / p2/ p3
ggsave("figure 4.png", width = 22, height = 32, units = "cm", dpi = 400)

library(dplyr)
