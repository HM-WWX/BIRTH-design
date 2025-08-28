library(readxl)
library(openxlsx)
library(ggplot2)
library(tidyr)
library(gridExtra)

data_readxl <- read_excel("C://Users//dell//Desktop//结果汇总4.xlsx", sheet = "Sheet1")
data <- as.matrix(data_readxl)
delete <- seq(1,136,by=8)

data <- data[-delete, ]
data <- data[, colSums(is.na(data)) == 0]

con_bias_row <- seq(1,112,by=7)
trt_theta_row <- seq(2,112,by=7) 
MSE_row <- seq(3,112,by=7) 
trt_95CI_row <- seq(4,112,by=7) 
con_95CI_row <- seq(5,112,by=7) 
Power_row <- seq(6,112,by=7) 
ESS_row <- seq(7,112,by=7) 
label <- rep(c("NB","MAP","CPP","DAW-CPP-MAP","1:1 trial"),16)

Bias <- matrix(NA,ncol=5,nrow=1)
for (i in con_bias_row){
  Bias_matrix <- matrix(data[i, ], ncol = 5, byrow = TRUE)
  Bias_matrix <- sweep(Bias_matrix, 2, Bias_matrix[1, ], FUN = "-")
  Bias <- rbind(Bias,Bias_matrix)
}
Bias <- Bias[-1,]
Bias <- apply(Bias,c(1,2), as.numeric)
Bias <- apply(Bias, 2, function(x) sprintf("%.4f", x))
Bias <- cbind(label,Bias)


Bias_wwx <- matrix(NA,ncol=5,nrow=1)
for (i in con_bias_row){
  Bias_matrix_wwx <- matrix(data[i, ], ncol = 5, byrow = TRUE)
  Bias_wwx <- rbind(Bias_wwx,Bias_matrix_wwx)
}
Bias_wwx <- Bias_wwx[-1,]
Bias_wwx <- apply(Bias_wwx,c(1,2), as.numeric)
Bias_wwx <- apply(Bias_wwx, 2, function(x) sprintf("%.4f", x))
Bias_wwx <- cbind(label,Bias_wwx)

MSE <- matrix(NA,ncol=5,nrow=1)
for (i in MSE_row){
  MSE_matrix <- matrix(data[i, ], ncol = 5, byrow = TRUE)
  MSE <- rbind(MSE,MSE_matrix)
}
MSE <- MSE[-1,]
MSE <- apply(MSE,c(1,2), as.numeric)
MSE <- apply(MSE, 2, function(x) sprintf("%.4f", x))
MSE <- cbind(label,MSE)

con_95CI <- matrix(NA,ncol=5,nrow=1)
for (i in con_95CI_row){
  con_95CI_matrix <- matrix(data[i, ], ncol = 5, byrow = TRUE)
  con_95CI <- rbind(con_95CI,con_95CI_matrix)
}
con_95CI <- con_95CI[-1,]
con_95CI <- apply(con_95CI,c(1,2), as.numeric)
con_95CI <- apply(con_95CI, 2, function(x) sprintf("%.4f", x))
con_95CI <- cbind(label,con_95CI)

Power <- matrix(NA,ncol=5,nrow=1)
for (i in Power_row){
  Power_matrix <- matrix(data[i, ], ncol = 5, byrow = TRUE)
  Power <- rbind(Power,Power_matrix)
}
Power <- Power[-1,]
Power <- apply(Power,c(1,2), as.numeric)
Power <- apply(Power, 2, function(x) sprintf("%.3f", x))
Power <- cbind(label,Power)

ESS <- matrix(NA,ncol=5,nrow=1)
for (i in ESS_row){
  ESS_matrix <- matrix(data[i, ], ncol = 5, byrow = TRUE)
  ESS <- rbind(ESS,ESS_matrix)
}
ESS <- ESS[-1,]
ESS <- apply(ESS,c(1,2), as.numeric)
ESS <- apply(ESS, 2, function(x) sprintf("%.1f", x))
ESS <- cbind(label,ESS)

ESS <- as.data.frame(ESS)
write.xlsx(ESS, file = "C:/Users/dell/Desktop/Work/result_table/ESS.xlsx")
Power <- as.data.frame(Power)
write.xlsx(Power, file = "C:/Users/dell/Desktop/Work/result_table/Power.xlsx")
con_95CI <- as.data.frame(con_95CI)
write.xlsx(con_95CI, file = "C:/Users/dell/Desktop/Work/result_table/con_95CI.xlsx")
MSE <- as.data.frame(MSE)
write.xlsx(MSE, file = "C:/Users/dell/Desktop/Work/result_table/MSE.xlsx")
Bias_wwx <- as.data.frame(Bias_wwx)
write.xlsx(Bias_wwx, file = "C:/Users/dell/Desktop/Work/result_table/Bias.xlsx")

# 假设 MSE 数据框已经存在并包含所需的数据
# 初始化列表以存储图形
Bias_list <- list()
for (i in 1:16) {
  start <- 1 + 5 * (i - 1)
  end <- 5 + 5 * (i - 1)
  V2 <- as.numeric(Bias[start:end, 2])
  V3 <- as.numeric(Bias[start:end, 3])
  V4 <- as.numeric(Bias[start:end, 4])
  V5 <- as.numeric(Bias[start:end, 5])
  V6 <- as.numeric(Bias[start:end, 6])
  
  # 创建示例数据框
  data <- data.frame(
    label = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"),
    V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6
  )
  
  # 将数据转换为长格式
  data_long <- pivot_longer(data, cols = V2:V6, names_to = "variable", values_to = "value")
  data_long$label <- factor(data_long$label, levels = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"))
  data_long$variable <- factor(data_long$variable, levels = c("V2", "V3", "V4", "V5", "V6"), labels = c("Sub 1", "Sub 2", "Sub 3", "Sub 4", "Sub 5"))
  
  colors <- c("NB" = "#449945", "MAP" = "#EA7827", "CPP" = "#83639F", "DAW-CPP-MAP" = "#C22F2F", "1:1 trial" = "#1F70A9")
  
  # 绘制折线图
  results_Bias <- ggplot(data_long, aes(x = variable, y = value, color = label, group = label)) +
    geom_line(size = 0.5) +  # 设置线条粗细
    geom_point(size= 1) +
    scale_y_continuous(limits = c(-0.1, 0.1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(face = "bold", size = 8),
      axis.text.y = element_text(face = "bold", size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 8),
      axis.title.y = element_text(face = "bold", size = 8),
      panel.grid = element_blank(),
      legend.position = ifelse(i == 16, "right", "none")
    ) +
    labs(x = "Subgroup", y = "Relative Bias", title = paste("Scenario", i), color = "Method") +
    scale_color_manual(values = colors)
  
  Bias_list[[i]] <- results_Bias
}

# 将图像排列成所需的布局
  Bias_1 <- do.call(grid.arrange, c(Bias_list[1:3], ncol = 3, nrow = 1))
  Bias_2 <- do.call(grid.arrange, c(Bias_list[4:6], ncol = 3, nrow = 1))
  Bias_3 <- do.call(grid.arrange, c(Bias_list[7:9], ncol = 3, nrow = 1))
  Bias_4 <- do.call(grid.arrange, c(Bias_list[10:12], ncol = 3, nrow = 1))
  Bias_5 <- do.call(grid.arrange, c(Bias_list[13:15], ncol = 3, nrow = 1))
  Bias_6 <- do.call(grid.arrange, c(Bias_list[16], ncol = 2, nrow = 1))
  ggsave("C:/Users/dell/Desktop/Work/result_picture/Bias/Bias_1.png", plot = Bias_1, device = "png", width = 24.5, height = 6.4, units = "cm")
  ggsave("C:/Users/dell/Desktop/Work/result_picture/Bias/Bias_2.png", plot = Bias_2, device = "png", width = 24.5, height = 6.4, units = "cm")
  ggsave("C:/Users/dell/Desktop/Work/result_picture/Bias/Bias_3.png", plot = Bias_3, device = "png", width = 24.5, height = 6.4, units = "cm")
  ggsave("C:/Users/dell/Desktop/Work/result_picture/Bias/Bias_4.png", plot = Bias_4, device = "png", width = 24.5, height = 6.4, units = "cm")
  ggsave("C:/Users/dell/Desktop/Work/result_picture/Bias/Bias_5.png", plot = Bias_5, device = "png", width = 24.5, height = 6.4, units = "cm")
  ggsave("C:/Users/dell/Desktop/Work/result_picture/Bias/Bias_6.png", plot = Bias_6, device = "png", width = 24.5, height = 6.4, units = "cm")

MSE_list <- list()
for (i in 1:16){
  start <- 1+5*(i-1)
  end <- 5+5*(i-1)
  V2 <- as.numeric(MSE[start:end,2])
  V3 <- as.numeric(MSE[start:end,3])
  V4 <- as.numeric(MSE[start:end,4])
  V5 <- as.numeric(MSE[start:end,5])
  V6 <- as.numeric(MSE[start:end,6])
  # 创建示例数据框
  data <- data.frame(
    label = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"),
    V2 = V2,V3 = V3,V4 = V4,V5 = V5,V6 = V6
  )
  # 将数据转换为长格式
  data_long <- pivot_longer(data, cols = V2:V6, names_to = "variable", values_to = "value")
  data_long$label <- factor(data_long$label, levels = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"))
  data_long$variable <- factor(data_long$variable, levels = c("V2", "V3", "V4", "V5", "V6"), labels = c("Sub 1", "Sub 2", "Sub 3", "Sub 4", "Sub 5"))
  colors <- c("NB" = "#449945", "MAP" = "#EA7827", "CPP" = "#83639F", "DAW-CPP-MAP" = "#C22F2F", "1:1 trial" = "#1F70A9")
  # 绘制条形图
  results_MSE <- ggplot(data_long, aes(x = variable, y = value, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", color = NA, width = 0.8) +
    scale_y_continuous(limits = c(0, 0.017)) +
    theme_minimal() +
    theme(axis.text.x = element_text(face="bold", size = 8),
          axis.text.y=element_text(face="bold", size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 8),  # 横轴标签加粗，12号字体
          axis.title.y = element_text(face = "bold", size = 8),
          panel.grid = element_blank(),
          legend.position = ifelse(i %in% c(16), "right", "none"),
          ) +
    labs(x = "Subgroup", y = "MSE", title = paste("Scenario", i),fill="Method") +
    scale_fill_manual(values = colors)
    MSE_list[[i]] <- results_MSE
  }
    MSE_1<-do.call(grid.arrange, c(MSE_list[1:3], ncol = 3, nrow = 1))
    MSE_2<-do.call(grid.arrange, c(MSE_list[4:6], ncol = 3, nrow = 1))
    MSE_3<-do.call(grid.arrange, c(MSE_list[7:9], ncol = 3, nrow = 1))
    MSE_4<-do.call(grid.arrange, c(MSE_list[10:12], ncol = 3, nrow = 1))
    MSE_5<-do.call(grid.arrange, c(MSE_list[13:15], ncol = 3, nrow = 1))
    MSE_6<-do.call(grid.arrange, c(MSE_list[16], ncol = 2, nrow = 1))
    ggsave("C:/Users/dell/Desktop/Work/result_picture/MSE/MSE_1.png", plot = MSE_1, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/MSE/MSE_2.png", plot = MSE_2, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/MSE/MSE_3.png", plot = MSE_3, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/MSE/MSE_4.png", plot = MSE_4, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/MSE/MSE_5.png", plot = MSE_5, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/MSE/MSE_6.png", plot = MSE_6, device = "png", width = 24.5, height = 6.4, units = "cm")
    
con_95CI_list <- list()
for (i in 1:16){
  start <- 1+5*(i-1)
  end <- 5+5*(i-1)
  V2 <- as.numeric(con_95CI[start:end,2])
  V3 <- as.numeric(con_95CI[start:end,3])
  V4 <- as.numeric(con_95CI[start:end,4])
  V5 <- as.numeric(con_95CI[start:end,5])
  V6 <- as.numeric(con_95CI[start:end,6])
  # 创建示例数据框
  data <- data.frame(
    label = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"),
    V2 = V2,V3 = V3,V4 = V4,V5 = V5,V6 = V6
  )
  # 将数据转换为长格式
  data_long <- pivot_longer(data, cols = V2:V6, names_to = "variable", values_to = "value")
  data_long$label <- factor(data_long$label, levels = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"))
  data_long$variable <- factor(data_long$variable, levels = c("V2", "V3", "V4", "V5", "V6"), labels = c("Sub 1", "Sub 2", "Sub 3", "Sub 4", "Sub 5"))
  colors <- c("NB" = "#449945", "MAP" = "#EA7827", "CPP" = "#83639F", "DAW-CPP-MAP" = "#C22F2F", "1:1 trial" = "#1F70A9")
  # 绘制条形图
  results_con_95CI <- ggplot(data_long, aes(x = variable, y = value, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", color = NA, width = 0.8) +
    scale_y_continuous(limits = c(0, 0.015)) +
    theme_minimal() +
    theme(axis.text.x = element_text(face="bold", size = 8),
          axis.text.y=element_text(face="bold", size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 8),  # 横轴标签加粗，12号字体
          axis.title.y = element_text(face = "bold", size = 8),
          panel.grid = element_blank(),
          legend.position = ifelse(i %in% c(16), "right", "none"),
          ) +
    labs(x = "Subgroup", y = "Width 95%CI", title = paste("Scenario", i),fill="Method") +
    scale_fill_manual(values = colors)
    con_95CI_list[[i]] <- results_con_95CI
  }
    con_95CI_1<-do.call(grid.arrange, c(con_95CI_list[1:3], ncol = 3, nrow = 1))
    con_95CI_2<-do.call(grid.arrange, c(con_95CI_list[4:6], ncol = 3, nrow = 1))
    con_95CI_3<-do.call(grid.arrange, c(con_95CI_list[7:9], ncol = 3, nrow = 1))
    con_95CI_4<-do.call(grid.arrange, c(con_95CI_list[10:12], ncol = 3, nrow = 1))
    con_95CI_5<-do.call(grid.arrange, c(con_95CI_list[13:15], ncol = 3, nrow = 1))
    con_95CI_6<-do.call(grid.arrange, c(con_95CI_list[16], ncol = 2, nrow = 1))
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Width/con_95CI_1.png", plot = con_95CI_1, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Width/con_95CI_2.png", plot = con_95CI_2, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Width/con_95CI_3.png", plot = con_95CI_3, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Width/con_95CI_4.png", plot = con_95CI_4, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Width/con_95CI_5.png", plot = con_95CI_5, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Width/con_95CI_6.png", plot = con_95CI_6, device = "png", width = 24.5, height = 6.4, units = "cm")
    
Power_list <- list()
for (i in 1:16){
  start <- 1+5*(i-1)
  end <- 5+5*(i-1)
  V2 <- as.numeric(Power[start:end,2])*100
  V3 <- as.numeric(Power[start:end,3])*100
  V4 <- as.numeric(Power[start:end,4])*100
  V5 <- as.numeric(Power[start:end,5])*100
  V6 <- as.numeric(Power[start:end,6])*100
  # 创建示例数据框
  data <- data.frame(
    label = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"),
    V2 = V2,V3 = V3,V4 = V4,V5 = V5,V6 = V6
  )
  # 将数据转换为长格式
  data_long <- pivot_longer(data, cols = V2:V6, names_to = "variable", values_to = "value")
  data_long$label <- factor(data_long$label, levels = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"))
  data_long$variable <- factor(data_long$variable, levels = c("V2", "V3", "V4", "V5", "V6"), labels = c("Sub 1", "Sub 2", "Sub 3", "Sub 4", "Sub 5"))
  colors <- c("NB" = "#449945", "MAP" = "#EA7827", "CPP" = "#83639F", "DAW-CPP-MAP" = "#C22F2F", "1:1 trial" = "#1F70A9")
  # 绘制条形图
  results_Power <- ggplot(data_long, aes(x = variable, y = value, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", color = NA, width = 0.8) +
    geom_text(aes(label = round(value, 1)), position = position_dodge(width = 0.8), vjust = -0.8, size = 1.25) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(face="bold", size = 8),
          axis.text.y=element_text(face="bold", size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 8),  # 横轴标签加粗，12号字体
          axis.title.y = element_text(face = "bold", size = 8),
          panel.grid = element_blank(),
          legend.position = ifelse(i %in% c(16), "right", "none"),
          ) +
    labs(x = "Subgroup", y = "Power", title = paste("Scenario", i),fill="Method") +
    scale_fill_manual(values = colors)
    if (i %in% c(4, 7, 8,9,10,11,12,13,14,15,16)) {
    results_Power <- results_Power + geom_hline(yintercept = 5, linetype = "dashed", color = "red")
    }
    if (i %in% c(4, 7, 8,9,10,11,12,13,14,15,16)) {
    results_Power <- results_Power + geom_hline(yintercept = 10, linetype = "dashed", color = "blue")
    }
    Power_list[[i]] <- results_Power
 }
    Power_1<-do.call(grid.arrange, c(Power_list[1:3], ncol = 3, nrow = 1))
    Power_2<-do.call(grid.arrange, c(Power_list[4:6], ncol = 3, nrow = 1))
    Power_3<-do.call(grid.arrange, c(Power_list[7:9], ncol = 3, nrow = 1))
    Power_4<-do.call(grid.arrange, c(Power_list[10:12], ncol = 3, nrow = 1))
    Power_5<-do.call(grid.arrange, c(Power_list[13:15], ncol = 3, nrow = 1))
    Power_6<-do.call(grid.arrange, c(Power_list[16], ncol = 2, nrow = 1))    
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Power/Power_1.png", plot = Power_1, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Power/Power_2.png", plot = Power_2, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Power/Power_3.png", plot = Power_3, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Power/Power_4.png", plot = Power_4, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Power/Power_5.png", plot = Power_5, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/Power/Power_6.png", plot = Power_6, device = "png", width = 24.5, height = 6.4, units = "cm")
    
ESS_list <- list()
for (i in 1:16){
  start <- 1+5*(i-1)
  end <- 5+5*(i-1)
  V2 <- as.numeric(ESS[start:end,2])
  V3 <- as.numeric(ESS[start:end,3])
  V4 <- as.numeric(ESS[start:end,4])
  V5 <- as.numeric(ESS[start:end,5])
  V6 <- as.numeric(ESS[start:end,6])
  # 创建示例数据框
  data <- data.frame(
    label = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"),
    V2 = V2,V3 = V3,V4 = V4,V5 = V5,V6 = V6
  )
  # 将数据转换为长格式
  data_long <- pivot_longer(data, cols = V2:V6, names_to = "variable", values_to = "value")
  data_long$label <- factor(data_long$label, levels = c("NB", "MAP", "CPP", "DAW-CPP-MAP", "1:1 trial"))
  data_long$variable <- factor(data_long$variable, levels = c("V2", "V3", "V4", "V5", "V6"), labels = c("Sub 1", "Sub 2", "Sub 3", "Sub 4", "Sub 5"))
  colors <- c("NB" = "#449945", "MAP" = "#EA7827", "CPP" = "#83639F", "DAW-CPP-MAP" = "#C22F2F", "1:1 trial" = "#1F70A9")
  # 绘制条形图
  results_ESS <- ggplot(data_long, aes(x = variable, y = value, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", color = NA, width = 0.8) +
    geom_text(aes(label = round(value, 1)), position = position_dodge(width = 0.8), vjust = -0.8, size = 1.25) +
    scale_y_continuous(limits = c(0, 102)) +
    theme_minimal() +
    theme(axis.text.x = element_text(face="bold", size = 8),
          axis.text.y=element_text(face="bold", size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 8),  # 横轴标签加粗，12号字体
          axis.title.y = element_text(face = "bold", size = 8),
          panel.grid = element_blank(),
          legend.position = ifelse(i %in% c(16), "right", "none"),
          ) +
    labs(x = "Subgroup", y = "ESS", title = paste("Scenario", i),fill="Method") +
    scale_fill_manual(values = colors)
    ESS_list[[i]] <- results_ESS
  }
    ESS_1<-do.call(grid.arrange, c(ESS_list[1:3], ncol = 3, nrow = 1))
    ESS_2<-do.call(grid.arrange, c(ESS_list[4:6], ncol = 3, nrow = 1))
    ESS_3<-do.call(grid.arrange, c(ESS_list[7:9], ncol = 3, nrow = 1))
    ESS_4<-do.call(grid.arrange, c(ESS_list[10:12], ncol = 3, nrow = 1))
    ESS_5<-do.call(grid.arrange, c(ESS_list[13:15], ncol = 3, nrow = 1))
    ESS_6<-do.call(grid.arrange, c(ESS_list[16], ncol = 2, nrow = 1))    
    ggsave("C:/Users/dell/Desktop/Work/result_picture/ESS/ESS_1.png", plot = ESS_1, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/ESS/ESS_2.png", plot = ESS_2, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/ESS/ESS_3.png", plot = ESS_3, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/ESS/ESS_4.png", plot = ESS_4, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/ESS/ESS_5.png", plot = ESS_5, device = "png", width = 24.5, height = 6.4, units = "cm")
    ggsave("C:/Users/dell/Desktop/Work/result_picture/ESS/ESS_6.png", plot = ESS_6, device = "png", width = 24.5, height = 6.4, units = "cm")
    
    
    
    
    
