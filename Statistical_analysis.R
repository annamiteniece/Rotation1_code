#Prerequisite: the csv files have to be stored in the manner outlined in Python script for this to work


#PACKAGES (might need to install the packages first with install.package("package_name"))
library(lme4)
library(lmerTest)
library(lsmeans) 
library(pbkrtest)

#FUNCTIONS

p_value_code_generator <- function(df, p_value_column){
  codes <- c()
  
  for (i in seq(1, length(p_value_column))) {
    
    if (is.na(p_value_column[i])) {
      codes[i] <- "NA"
    } else if (p_value_column[i] <= 1 & p_value_column[i] > 0.1) {
      codes[i] <- "ns"
    } else if (p_value_column[i] <= 0.1 & p_value_column[i] > 0.05) {
      codes[i] <- "."
    } else if (p_value_column[i] <= 0.05 & p_value_column[i] > 0.01) {
      codes[i] <- "*"
    } else if (p_value_column[i] <= 0.01 & p_value_column[i] > 0.001) {
      codes[i] <- "**"
    } else if (p_value_column[i] <= 0.001) {
      codes[i] <- "***"
    }
  }
  
  
  df$codes <- codes
  return(df)
  
}


#TOTAL_DF ANALYSIS FROM EACH OF THE SESSIONS (e.g., nov2023 OF_reg or aug2022 habituation)

# Output is total_comparison_df_list containing individual sublists for each session
## Each sublist contains:
### normalcy_df, which shows p-value from Shapiro test for CFA/Saline mice for each column (variable, e.g., average_speed_total)
### total_comparison_df, which shows statistic from t-test or wilcoxon-test (depending on the normality of the distribution), degrees of freedom (only available for t-test); p_value from t_test or wilcoxon-test, type of test performed for each column/variable


base_path = "/Users/annam/Desktop/Rotation1/STATS/total_df/"

file_names <- list.files(base_path, full.names = FALSE)

total_comparison_df_list <- list()

for (file_name in file_names) {
  print(file_name)
  total_comparison_df_list[[file_name]] <- list(total_comparison_df = data.frame(), normalcy_df = data.frame())
}


for (file_name in file_names) {
  
  path = paste(base_path, file_name, sep = "")
  df <- read.csv(path, header = TRUE, stringsAsFactors = TRUE)
  groups = c("CFA", "Saline")
  
  column_names <- colnames(df)
  column_remove <- c("mouse_names", "group_names", "session")
  column_names <- column_names[!column_names %in% column_remove]
  
  column_list = list()
  
  
  for (column_name in column_names) {
    column_list[[column_name]] <- integer(0)
  }
  
  temp <- data.frame()
  group_names <- c()
  
  for (group in groups) {
    
    group_names <- append(group_names, group) 
    temp <- df[df$group_names == group, ]
    
    for (column_name in column_names) {
      
      shaptest <- shapiro.test(temp[[column_name]])
      column_list[[column_name]] <- append(column_list[[column_name]], shaptest$p.value)
    
    }
  }
  
  normalcy_df <- data.frame(groups)
  
  for (column_name in column_names) {
    normalcy_df[[column_name]] <- column_list[[column_name]]
  }
  
  total_comparison_df_list[[file_name]]$normalcy_df <- normalcy_df
  
  column_names_df <- c()
  stat <- c()
  deg_free <- c()
  p_values <- c()
  test_type <- c()
  
  for (column_name in column_names) {
    
    column_names_df <- append(column_names_df, column_name)
    
    if (normalcy_df[normalcy_df$groups == "CFA", ][[column_name]] >= 0.05 & 
        normalcy_df[normalcy_df$groups == "Saline", ][[column_name]] >= 0.05) {
      
      test_type <- append(test_type, "t-test")
      t_test <- t.test(df[df$group_names == "CFA", ][[column_name]], df[df$group_names == "Saline", ][[column_name]], paired = TRUE)
      stat <- append(stat, t_test$statistic)
      deg_free <- append(deg_free, t_test$parameter)
      p_values <- append(p_values, t_test$p.value)
    }
    
    if (normalcy_df[normalcy_df$groups == "CFA", ][[column_name]] < 0.05 | normalcy_df[normalcy_df$groups == "Saline", ][[column_name]] < 0.05) {
      
      test_type <- append(test_type, "wilcox-test")
      w_test <- wilcox.test(df[df$group_names == "CFA", ][[column_name]], df[df$group_names == "Saline", ][[column_name]], paired = TRUE)
      stat <- append(stat, w_test$statistic)
      deg_free <- append(deg_free, "NA")
      p_values <- append(p_values, w_test$p.value)
    }
  }
  
  total_comparison_df <- data.frame(column_names_df, stat, deg_free, p_values, test_type) 
  total_comparison_df <- p_value_code_generator(total_comparison_df, total_comparison_df$p_values)
  
  total_comparison_df_list[[file_name]]$total_comparison_df <- total_comparison_df
  
  print(file_name)
}

#The total_comparison_df for OF_reg can be accessed like this: total_comparison_df_list[["nov2023_OF_reg_total_df"]][["total_comparison_df"]] or this total_comparison_df_list$nov2023_OF_reg_total_df$total_comparison_df





#INTERVAL_DF ANALYSIS FOR EACH OF THE SESSIONS 

# Output is interval_comparison_df_list containing 2 individual sublists for each session: 
## normalcy_df, which shows p-value from Shapiro test for CFA/Saline mice for each interval for each column (variable, e.g., interval_speed_total) --> can be accessed: interval_comparison_df_list$nov2023_OF_reg_interval_df$normalcy_df
## column_list_comparisons, which contains 17 sublists, each named after the column on which analysis is being performed (e.g., interval_speed_total) containing: 

### anova_summary sublist, which contains summary of the ANOVA run for the column (factors: group_names * interval_names) --> can be accessed: interval_comparison_df_list$nov2023_OF_reg_interval_df$column_list_comparisons$interval_speed_total$anova_summary
### lmm_summary sublist, which contains summary of the LMM run for the column (factors: group_names * interval_names + 1/mouse_names) --> can be accessed: interval_comparison_df_list$nov2023_OF_reg_interval_df$column_list_comparisons$interval_speed_total$lmm_summary

### treatment sublist, which containts two dataframes:
#### anova, which contains p values from TukeyHSD pairwise comparisons from the ANOVA fit between averages of the column for CFA and saline in specific intervals (e.g., interval_speed_total in interval 0:5 CFA vs 0:5 Saline) --> can be accessed: interval_comparison_df_list$nov2023_OF_reg_interval_df$column_list_comparisons$interval_speed_total$treatment$anova
### lmm, which contains p values from TukeyHSD pairwise comparisons from the LMM fit between averages of the column for CFA and saline in specific intervals (e.g., interval_speed_total in interval 0:5 CFA vs 0:5 Saline) --> can be accessed: interval_comparison_df_list$nov2023_OF_reg_interval_df$column_list_comparisons$interval_speed_total$treatment$lmm

### time sublist, which contains 2 sublists: 
### CFA, which contains anova and lmm df that are similar as those described for treatment$anova/lmm, but the pairwise comparisons are within group between time intervals (e.g., interval_speed_total in interval 0:5 CFA vs interval 5:10 CFA) --> interval_comparison_df_list$nov2023_OF_reg_interval_df$column_list_comparisons$interval_speed_total$time$lmm
### Saline, which contains anova and lmm as for CFA 

base_path = "/Users/annam/Desktop/Rotation1/STATS/interval_df/"

file_names <- list.files(base_path, full.names = FALSE)

interval_comparison_df_list <- list()

for (file_name in file_names) {
  print(file_name)
  interval_comparison_df_list[[file_name]] <- list(normalcy_df = data.frame() , column_list_comparisons = list())
}


for (file_name in file_names) {
  path = paste(base_path, file_name, sep = "")
  df <- read.csv(path, header = TRUE, stringsAsFactors = TRUE)
  column_names <- colnames(df)
  column_remove <- c("mouse_names", "group_names", "interval_names")
  column_names <- column_names[!column_names %in% column_remove]
  
  interval_names <- levels(df$interval_names)
  groups = c("CFA", "Saline")
  
  #Testing normality 
  column_list = list()
  for (column_name in column_names) {
    for (interval_name in interval_names) {
      for (group in groups) {
        column_list[[column_name]][[interval_name]][[group]] <- integer(0)
      }
    }
  }
  
  for (interval_name in interval_names) {
    
    for (group in groups)  {
      
      temp <- df[df$group_names == group & df$interval_names == interval_name, ]
      
      for (column_name in column_names) {
        
        cleaned <- as.vector(na.omit(temp[[column_name]]))
        
        if (length(cleaned) >=3 & length(unique(temp[[column_name]])) > 1) {
          shaptest <- shapiro.test(temp[[column_name]])
          column_list[[column_name]][[interval_name]][[group]] <- append(column_list[[column_name]][[interval_name]][[group]], shaptest$p.value)
        }
        
        if (length(cleaned) < 3 | length(unique(temp[[column_name]])) == 1)  {
          column_list[[column_name]][[interval_name]][[group]] <- append(column_list[[column_name]][[interval_name]][[group]], NA)
        }
        
      }
    }
  }
  
  group_names_mini <- c()
  interval_names_mini <- c()
  
  for (interval_name in interval_names) {
    for (group in groups) {
      group_names_mini <- append(group_names_mini, group)
      interval_names_mini <- append(interval_names_mini, interval_name)
    }
  }
  
  normalcy_df <- data.frame(interval_names_mini, group_names_mini)

  shap_values_mini <- c()
  for (column_name in column_names){
    for (interval_name in interval_names) {
      for (group in groups) {
        shap_values_mini <- append(column_list[[column_name]][[interval_name]][[group]], shap_values_mini)
      }
    }
    normalcy_df[[column_name]] <- shap_values_mini
    shap_values_mini <- c()
  }
  
  interval_comparison_df_list[[file_name]]$normalcy_df <- normalcy_df
  
  
  #FITTING ANOVA 
  
  for (column_name in column_names) {
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$anova_summary <- list()
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$lmm_summary <- list()
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$treatment <-list(anova = data.frame(), lmm = data.frame())
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$time <- list(CFA = list(anova = data.frame(), lmm = data.frame()), Saline = list(anova = data.frame(), lmm = data.frame()))
  }
  
  for (column_name in column_names) {
    column_df <- df[, c("interval_names", "group_names", "mouse_names", column_name)]
    column_df <- na.omit(column_df)
    
    anova_result <- aov(column_df[[column_name]] ~ column_df$interval_names * column_df$group_names)
    anova_summary = summary(anova_result)
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$anova_summary <- anova_summary
    
    
    tukey_results <- TukeyHSD(anova_result)
    
    tukey_df <- data.frame(tukey_results$`column_df$interval_names:column_df$group_names`)
    
    #FINDING CFA VS SALINE ANOVA COMPARISONS
    row_names <- c()
    for (interval in interval_names) {
      row_name <- paste(interval, ":Saline-", interval, ":CFA", sep = "")
      row_names <- append(row_names, row_name)
    }
    
    interval_comparison_df_anova <- tukey_df[row_names, ]
    interval_comparison_df_anova <- p_value_code_generator(interval_comparison_df_anova, interval_comparison_df_anova$p.adj)
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$treatment$anova <- interval_comparison_df_anova
    
    #FINDING INTERVAL COMPARISONS
    
    CFA_rownames <- c()
    Saline_rownames <- c()
    for (i in 1:(length(interval_names) - 1)) {
      for (j in (i + 1):length(interval_names)) {
        CFA_rowname <- paste(interval_names[j], ":CFA-", interval_names[i], ":CFA", sep = "")
        Saline_rowname <- paste(interval_names[j], ":Saline-", interval_names[i], ":Saline", sep = "")
        
        CFA_rownames <- c(CFA_rownames, CFA_rowname)
        Saline_rownames <- c(Saline_rownames, Saline_rowname)
      }
    }
    
    
    CFA_tukey_df <- tukey_df[CFA_rownames, ] 
    CFA_tukey_df <- p_value_code_generator(CFA_tukey_df, CFA_tukey_df$p.adj)
    CFA_tukey_df <- CFA_tukey_df[, (ncol(CFA_tukey_df) - 1):ncol(CFA_tukey_df)]
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$time$CFA$anova <- CFA_tukey_df
    
    Saline_tukey_df <- tukey_df[Saline_rownames, ]
    Saline_tukey_df <- p_value_code_generator(Saline_tukey_df, Saline_tukey_df$p.adj)
    Saline_tukey_df <- Saline_tukey_df[, (ncol(Saline_tukey_df) - 1):ncol(Saline_tukey_df)]
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$time$Saline$anova <- Saline_tukey_df
    
    
    #FITTING Linear mixed model (LMM) 

    
    mixed_model <- lmer(column_df[[column_name]] ~ column_df$group_names * column_df$interval_names + (1| column_df$mouse_names))
    lmm_summary <- summary(mixed_model)
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$lmm_summary <- lmm_summary
    
    
    emmeans_results <- emmeans(mixed_model, ~ group_names * interval_names)
    pairwise_comp <- pairs(emmeans_results)
    pairwise_comp <- data.frame(pairwise_comp)
    rownames(pairwise_comp) <- pairwise_comp$contrast
    pairwise_comp$contrast <- NULL
    
    #FINDING CFA VS SALINE LMM COMPARISONS
    row_names <- c()
    for (interval in interval_names) {
      row_name <- paste("CFA ", interval, " - Saline ", interval, sep = "")
      
      row_names <- append(row_names, row_name)
    }
  
    interval_comparison_df_lmm <- pairwise_comp[row_names, ]
    interval_comparison_df_lmm <- p_value_code_generator(interval_comparison_df_lmm, interval_comparison_df_lmm$p.value)
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$treatment$lmm <- interval_comparison_df_lmm
    
    #FINDING INTERVAL LMM COMPARISONS
    
    CFA_rownames <- c()
    Saline_rownames <- c()
    for (i in 1:(length(interval_names) - 1)) {
      for (j in (i + 1):length(interval_names)) {
        CFA_rowname <- paste("CFA", interval_names[i], "-","CFA", interval_names[j], sep = " ")
        Saline_rowname <- paste("Saline", interval_names[i], "-","Saline", interval_names[j], sep = " ")
        
        CFA_rownames <- c(CFA_rownames, CFA_rowname)
        Saline_rownames <- c(Saline_rownames, Saline_rowname)
      }
    }
    
    CFA_lmm_df <- pairwise_comp[CFA_rownames, ] 
    CFA_lmm_df <- p_value_code_generator(CFA_lmm_df, CFA_lmm_df$p.value)
    CFA_lmm_df <- CFA_lmm_df[, (ncol(CFA_lmm_df) - 1):ncol(CFA_lmm_df)]
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$time$CFA$lmm <- CFA_lmm_df
    
    Saline_lmm_df <- pairwise_comp[Saline_rownames, ]
    Saline_lmm_df <- p_value_code_generator(Saline_lmm_df, Saline_lmm_df$p.value)
    Saline_lmm_df <- Saline_lmm_df[, (ncol(Saline_lmm_df) - 1):ncol(Saline_lmm_df)]
    interval_comparison_df_list[[file_name]]$column_list_comparisons[[column_name]]$time$Saline$lmm <- Saline_lmm_df
    
  }
}






#WITHIN OF_THREAT SUBSESSION ANALYSIS 

# Output is column_list containing 4 individual sublists for each session: 
### anova_summary sublist, which contains summary of the ANOVA run for the column (factors: group_names * session_names) --> can be accessed: column_list$average_speed_total$anova_summary
### lmm_summary sublist, which contains summary of the LMM run for the column (factors: group_names * session_names + 1/mouse_names) --> can be accessed: column_list$average_speed_total$lmm_summary

### treatment sublist, which containts two dataframes:
#### anova, which contains p values from TukeyHSD pairwise comparisons from the ANOVA fit between averages of the column for CFA and saline in specific subsessions of OF_threat (e.g., average_speed_total in open_arena CFA vs open_arena Saline) --> can be accessed: column_list$average_speed_total$treatment$anova
### lmm, which contains p values from TukeyHSD pairwise comparisons from the LMM fit between averages of the column for CFA and saline in specific subsessions of OF_threat (e.g., interval_speed_total in interval 0:5 CFA vs 0:5 Saline) --> can be accessed: column_list$average_speed_total$treatment$lmm

### time sublist, which contains 2 sublists: 
### CFA, which contains anova and lmm df that are similar as those described for treatment$anova/lmm, but the pairwise comparisons are within group between subessions of OF_threat (e.g., average_speed_total in open_arena CFA vs novel_object  CFA) --> column_list$average_speed_total$time$anova
### Saline, which contains anova and lmm as for CFA 



open_arena_df <- read.csv("/Users/annam/Desktop/Rotation1/STATS/total_df/nov2023_OF_threat_total_df", header = TRUE, stringsAsFactors = TRUE)
session <- rep("open_arena", nrow(open_arena_df))
open_arena_df$session <- session

novel_object_df <- read.csv("/Users/annam/Desktop/Rotation1/STATS/total_df/nov2023_OF_threat_bodypart_df_novel_total_df", header = TRUE, stringsAsFactors = TRUE)
session <- rep("novel_object", nrow(novel_object_df))
novel_object_df$session <- session

post_object_df <- read.csv("/Users/annam/Desktop/Rotation1/STATS/total_df/nov2023_OF_threat_bodypart_df_post_object_total_df", header = TRUE, stringsAsFactors = TRUE)
session <- rep("post_object", nrow(post_object_df))
post_object_df$session <- session

threat_df <- read.csv("/Users/annam/Desktop/Rotation1/STATS/total_df/nov2023_OF_threat_bodypart_df_threat_total_df", header = TRUE, stringsAsFactors = TRUE)
session <- rep("threat", nrow(threat_df))
threat_df$session <- session

threat_session_df <- rbind(open_arena_df, novel_object_df, post_object_df, threat_df)
threat_session_df$session <- as.factor(threat_session_df$session)


column_names <- colnames(threat_session_df)
column_remove <- c("mouse_names", "group_names", "session")
column_names <- column_names[!column_names %in% column_remove]

session_names <- unique(threat_session_df$session)


column_list <- list()

for (column_name in column_names) {
  
  column_list[[column_name]] <- list(anova_summary = data.frame(), lmm_summary = data.frame(), treatment = list(anova = data.frame(), lmm = data.frame), time = list(CFA = list(anova = data.frame(), lmm = data.frame), Saline = list(anova = data.frame(), lmm = data.frame)))
  
  column_df <- threat_session_df[, c("session", "group_names", "mouse_names", column_name)]
  column_df <- na.omit(column_df)
  
  #FITTING ANOVA 
  
  anova_result <- aov(column_df[[column_name]] ~ column_df$session * column_df$group_names)
  anova_summary <- summary(anova_result)
  column_list[[column_name]]$anova_summary <- anova_summary
  
  tukey_results <- TukeyHSD(anova_result)
  tukey_results
  tukey_df <- data.frame(tukey_results$`column_df$session:column_df$group_names`)
  tukey_df <- p_value_code_generator(tukey_df, tukey_df$p.adj)
  
  #FINDING CFA VS SALINE ANOVA COMPARISONS 
  
  row_names <- c()
  for (session in session_names) {
    row_name <- paste(session, ":Saline-", session, ":CFA", sep = "")
    row_names <- append(row_names, row_name)
  }
  
  interval_comparison_df_anova <- tukey_df[row_names, ]
  column_list[[column_name]]$treatment$anova <- interval_comparison_df_anova
  
  #FINDING INTERVAL ANOVA COMPARISONS 
  
  CFA_rownames <- c()
  Saline_rownames <- c()
  for (i in 1:(length(session_names) - 1)) {
    for (j in (i + 1):length(session_names)) {
      # Construct the string with the word inserted between the letters
      CFA_rowname <- paste(session_names[j], ":CFA-", session_names[i], ":CFA", sep = "")
      Saline_rowname <- paste(session_names[j], ":Saline-", session_names[i], ":Saline", sep = "")
      
      # Add the new string to the result vector
      CFA_rownames <- c(CFA_rownames, CFA_rowname)
      Saline_rownames <- c(Saline_rownames, Saline_rowname)
    }
  }
  CFA_tukey_df <- tukey_df[CFA_rownames, ] 
  CFA_tukey_df <- CFA_tukey_df[, (ncol(CFA_tukey_df) - 1):ncol(CFA_tukey_df)]
  column_list[[column_name]]$time$CFA$anova <- CFA_tukey_df
  
  Saline_tukey_df <- tukey_df[Saline_rownames, ]
  Saline_tukey_df <- Saline_tukey_df[, (ncol(Saline_tukey_df) - 1):ncol(Saline_tukey_df)]
  column_list[[column_name]]$time$Saline$anova <- Saline_tukey_df
  
  
  # FITTING Linear mixed model (LMM) 
  
  
  mixed_model <- lmer(column_df[[column_name]] ~ column_df$group_names * column_df$session + (1| column_df$mouse_names))
  lmm_summary <- summary(mixed_model)
  column_list[[column_name]]$lmm_summary <- lmm_summary
  
  emmeans_results <- emmeans(mixed_model, ~ group_names * session)
  pairwise_comp <- pairs(emmeans_results)
  pairwise_comp <- data.frame(pairwise_comp)
  rownames(pairwise_comp) <- pairwise_comp$contrast
  pairwise_comp$contrast <- NULL
  pairwise_comp <- p_value_code_generator(pairwise_comp, pairwise_comp$p.value)
  
  # FINDING CFA VS SALINE LMM COMPARISONS
  row_names <- c()
  for (session in session_names) {
    row_name <- paste("CFA ", session, " - Saline ", session, sep = "")
    row_names <- append(row_names, row_name)
  }
  
  interval_comparison_df_lmm <- pairwise_comp[row_names, ]
  column_list[[column_name]]$treatment$lmm <- interval_comparison_df_lmm
  
  # FINDING INTERVAL LMM COMPARISONS
  
  CFA_rownames <- c()
  Saline_rownames <- c()
  for (i in 1:(length(session_names) - 1)) {
    for (j in (i + 1):length(session_names)) {
      # Construct the string with the word inserted between the letters
      CFA_rowname <- paste("CFA", session_names[i], "-","CFA", session_names[j], sep = " ")
      Saline_rowname <- paste("Saline", session_names[i], "-","Saline", session_names[j], sep = " ")
      
      # Add the new string to the result vector
      CFA_rownames <- c(CFA_rownames, CFA_rowname)
      Saline_rownames <- c(Saline_rownames, Saline_rowname)
    }
  }
  
  CFA_rownames[1] <- "CFA novel_object - CFA open_arena" 
  Saline_rownames[1] <- "Saline novel_object - Saline open_arena" 
  
  CFA_lmm_df <- pairwise_comp[CFA_rownames, ] 
  CFA_lmm_df <- CFA_lmm_df[, (ncol(CFA_lmm_df) - 1):ncol(CFA_lmm_df)]
  column_list[[column_name]]$time$CFA$lmm <- CFA_lmm_df
  
  Saline_lmm_df <- pairwise_comp[Saline_rownames, ]
  Saline_lmm_df <- Saline_lmm_df[, (ncol(Saline_lmm_df) - 1):ncol(Saline_lmm_df)]
  column_list[[column_name]]$time$Saline$lmm <- Saline_lmm_df
  
}




#OBJECT APPROACH ANALYSIS 

#Running this code will print:

#For total approaches: 
## Shapiro test p values for CFA vs Saline
#t-test results for CFA vs Saline 

#For interval approaches: 
## Shapiro test p values for CFA intervals and Saline intervals 
## summary of ANOVA (factors: group_names*interval_names) 
## pairwise ANOVA comparisons for intervals CFA vs Saline (0:1 CFA vs 0:1 Saline)
## pairwise ANOVA comparisons for within group between intervals (0:1 CFA vs 1:2 CFA)
## summary of LMM (factors: group_names*interval_names + 1/mouse_names) 
## pairwise LMM comparisons for intervals CFA vs Saline (0:1 CFA vs 0:1 Saline)
## pairwise LMM comparisons for within group between intervals (0:1 CFA vs 1:2 CFA)



interval_df <- read.csv("/Users/annam/Desktop/Rotation1/object_approach_csv.csv", header = TRUE, stringsAsFactors = TRUE)
mouse_names <- levels(interval_df$Mouse_ID)

CFA <- c("90L", "90RR", "89R", "89LL", "88LL", "88R", "87RR", "87L") 
Saline <- c("90R", "90LL", "89RR", "89L", "88RR", "88L", "87R", "87LL")

group_names_interval <- c()
for (i in seq(nrow(interval_df))) {
  if (interval_df$Mouse_ID[i] %in% CFA) {
    group_names_interval <- append(group_names_interval, "CFA")
  }
  else {
    group_names_interval <- append(group_names_interval, "Saline")
  }
}

interval_df$group_names <- group_names_interval
  
total_approach_nr <- c()
group_names_total <- c()

for (i in seq(length(mouse_names))) {
  mouse_df <- interval_df[interval_df$Mouse_ID == mouse_names[i], ]
  total_approach_nr <- append(total_approach_nr, sum(mouse_df$Object_approaches))
  if (mouse_names[i] %in% CFA) {
    group_names_total <- append(group_names_total, "CFA")
  } 
  
  else {
    group_names_total <- append(group_names_total, "Saline")
    }
}

total_df <- data.frame(mouse_names, group_names_total, total_approach_nr)

total_df_g1 <- total_df[total_df$group_names_total == "CFA", ]
shaptest_g1 <- shapiro.test(total_df_g1$total_approach_nr)
print(paste("CFA shapiro test p-value is", shaptest_g1$p.value))


total_df_g2 <- total_df[total_df$group_names_total == "Saline", ]
shaptest_g2 <- shapiro.test(total_df_g2$total_approach_nr)
print(paste("CFA shapiro test p-value is", shaptest_g2$p.value))


t_test <- t.test(total_df[total_df$group_names_total == "CFA", ]$total_approach_nr, total_df[total_df$group_names_total == "Saline", ]$total_approach_nr, paired = TRUE)
print(t_test)


interval_names <- unique(interval_df$Until_minute)

#Testing normality 

interval_df_g1 = interval_df[interval_df$group_names == "CFA", ]
temp <- data.frame()
shap_p_value_g1 <- c()
for (interval_name in interval_names) {
  temp <- interval_df_g1[interval_df_g1$Until_minute == interval_name, ]
  shaptest <- shapiro.test(temp$Object_approaches)
  shap_p_value_g1 <- append(shap_p_value_g1, shaptest$p.value)
}

interval_shap_g1 <- data.frame(interval_names, shap_p_value_g1)
print("CFA interval normality")
print(interval_shap_g1)

interval_df_g2 = interval_df[interval_df$group_names == "Saline", ]
temp <- data.frame()
shap_p_value_g2 <- c()

for (interval_name in interval_names) {
  temp <- interval_df_g2[interval_df_g2$Until_minute == interval_name, ]
  shaptest <- shapiro.test(temp$Object_approaches)
  shap_p_value_g2 <- append(shap_p_value_g2, shaptest$p.value)
}
interval_shap_g2 <- data.frame(interval_names, shap_p_value_g2)
print("Saline interval normality")
print(interval_shap_g2)

interval_df$Until_minute <- factor(interval_df$Until_minute)
interval_df$group_names <- factor(interval_df$group_names)

#ANOVA

anova_result <- aov(Object_approaches ~ Until_minute * group_names, data = interval_df)
print("Interval ANOVA summary")
summary(anova_result)

tukey_results <- TukeyHSD(anova_result)

tukey_df <- data.frame(tukey_results$`Until_minute:group_names`)
  
row_names <- c()
for (interval in interval_names) {
  row_name <- paste(interval, ":Saline-", interval, ":CFA", sep = "")
  row_names <- append(row_names, row_name)
}

print("CFA vs Saline ANOVA comparisons")
interval_comparison_df_anova_OBJECT <- tukey_df[row_names, ]
interval_comparison_df_anova_OBJECT <- p_value_code_generator(interval_comparison_df_anova_OBJECT, interval_comparison_df_anova_OBJECT$p.adj)

CFA_rownames <- c()
Saline_rownames <- c()
for (i in 1:(length(interval_names) - 1)) {
  for (j in (i + 1):length(interval_names)) {
    # Construct the string with the word inserted between the letters
    CFA_rowname <- paste(interval_names[j], ":CFA-", interval_names[i], ":CFA", sep = "")
    Saline_rowname <- paste(interval_names[j], ":Saline-", interval_names[i], ":Saline", sep = "")
    
    # Add the new string to the result vector
    CFA_rownames <- c(CFA_rownames, CFA_rowname)
    Saline_rownames <- c(Saline_rownames, Saline_rowname)
  }
}

CFA_tukey_df <- tukey_df[CFA_rownames, ] 
CFA_tukey_df <- p_value_code_generator(CFA_tukey_df, CFA_tukey_df$p.adj)
CFA_tukey_df <- CFA_tukey_df[, (ncol(CFA_tukey_df) - 1):ncol(CFA_tukey_df)]
print("CFA interval ANOVA comparisons")
print(CFA_tukey_df)

Saline_tukey_df <- tukey_df[Saline_rownames, ]
Saline_tukey_df <- p_value_code_generator(Saline_tukey_df, Saline_tukey_df$p.adj)
Saline_tukey_df <- Saline_tukey_df[, (ncol(Saline_tukey_df) - 1):ncol(Saline_tukey_df)]
print("Saline interval ANOVA comparisons")
print(Saline_tukey_df)


# Linear mixed model (LMM) 

mixed_model <- lmer(Object_approaches ~ group_names * Until_minute + (1| Mouse_ID), data = interval_df)
print("Interval LMM summary")
summary(mixed_model)

emmeans_results <- emmeans(mixed_model, ~ group_names * Until_minute)
pairwise_comp <- pairs(emmeans_results)
pairwise_comp <- data.frame(pairwise_comp)
rownames(pairwise_comp) <- pairwise_comp$contrast
pairwise_comp$contrast <- NULL

row_names <- c()
for (interval in interval_names) {
  row_name <- paste("CFA Until_minute", interval, " - Saline Until_minute", interval, sep = "")
  
  row_names <- append(row_names, row_name)
}

interval_comparison_df_lmm_OBJECT <- pairwise_comp[row_names, ]
interval_comparison_df_lmm_OBJECT <- p_value_code_generator(interval_comparison_df_lmm_OBJECT, interval_comparison_df_lmm_OBJECT$p.value)

print("CFA vs Saline LMM comparisons")
print(interval_comparison_df_lmm_OBJECT)

CFA_rownames <- c()
Saline_rownames <- c()
for (i in 1:(length(interval_names) - 1)) {
  for (j in (i + 1):length(interval_names)) {
    # Construct the string with the word inserted between the letters
    CFA_rowname <- paste("CFA Until_minute", interval_names[i], " - CFA Until_minute", interval_names[j], sep = "")
    Saline_rowname <- paste("Saline Until_minute", interval_names[i], " - Saline Until_minute", interval_names[j], sep = "")
    
    # Add the new string to the result vector
    CFA_rownames <- c(CFA_rownames, CFA_rowname)
    Saline_rownames <- c(Saline_rownames, Saline_rowname)
  }
}

CFA_lmm_df <- pairwise_comp[CFA_rownames, ] 
CFA_lmm_df <- p_value_code_generator(CFA_lmm_df, CFA_lmm_df$p.value)
CFA_lmm_df <- CFA_lmm_df[, (ncol(CFA_lmm_df) - 1):ncol(CFA_lmm_df)]
print("CFA interval LMM comparisons")


Saline_lmm_df <- pairwise_comp[Saline_rownames, ]
Saline_lmm_df <- p_value_code_generator(Saline_lmm_df, Saline_lmm_df$p.value)
Saline_lmm_df <- Saline_lmm_df[, (ncol(Saline_lmm_df) - 1):ncol(Saline_lmm_df)]
print("Saline interval LMM comparisons")






