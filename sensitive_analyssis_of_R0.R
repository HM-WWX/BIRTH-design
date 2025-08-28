library(readxl)
library(tidyverse)

# 最终版数据读取函数（支持动态列数）
read_bias_data <- function(path, method_name) {
  # 读取前6列数据（跳过后续列）
  df <- read_excel(path, 
                   col_names = paste0("col", 1:7),  # 命名前6列
                   range = cell_cols(1:7)) %>%      # 强制读取前6列
    select(1:7)  # 确保只保留前6列
  
  df_clean <- df %>%
    mutate(
      row_id = row_number(),
      # 识别场景行（兼容中英文）
      is_scenario = str_detect(col1, "Scenario|场景"),
      # 识别数据行（支持R0和NB）
      is_data_row = str_detect(col2, "R0|^NB$")
    ) %>%
    # 填充场景标签
    mutate(scenario = ifelse(is_scenario, col1, NA_character_)) %>%
    fill(scenario, .direction = "down") %>%
    # 过滤保留数据行
    filter(is_data_row) %>%
    # 类型安全转换
    mutate(
      across(col3:col7, as.numeric),  # 转换数值列
      R0 = parse_number(str_remove(col2, "R0=")),        # 提取R0数值
      Method = method_name
    ) %>%
    # 清理最终数据
    select(Method, scenario, R0, 
           subgroup1 = col3, subgroup2 = col4, 
           subgroup3 = col5, subgroup4 = col6, subgroup5 = col7) %>%
    pivot_longer(
      cols = subgroup1:subgroup5,
      names_to = "subgroup",
      values_to = "bias"
    ) %>%
    filter(!is.na(R0), !is.na(bias))  # 确保有效数据
  
  return(df_clean)
}


# 后续统计和绘图代码保持不变...

# 重新读取数据
daw <- read_bias_data("C:/Users/dell/Desktop/Work/文章整理结果/result_table/Bias_DAW-MAP-CPP.xlsx", "BIRTH design")
cpp <- read_bias_data("C:/Users/dell/Desktop/Work/文章整理结果/result_table/Bias_CPP.xlsx", "1:1 trial")
combined <- bind_rows(daw, cpp)

# 场景顺序控制（关键修改点）----------------------------------------------------
scenario_levels <- str_sort(unique(combined$scenario), numeric = TRUE)
combined <- combined %>%
  mutate(scenario = factor(scenario, levels = scenario_levels))

# 敏感性分析
sensitivity_data <- combined %>%
  filter(R0 > 0) %>%
  mutate(logR0 = log10(R0)) %>%
  group_by(scenario, subgroup, Method) %>%
  summarise(
    sensitivity = lm(bias ~ logR0)$coefficients[2],
    .groups = "drop"
  ) %>%
  mutate(scenario = factor(scenario, levels = scenario_levels))

# 综合统计量
overall_sensitivity <- sensitivity_data %>%
  group_by(Method, scenario) %>%
  summarise(
    mean_sensitivity = mean(sensitivity),
    se_sensitivity = sd(sensitivity) / sqrt(n()),
    .groups = "drop"
  )

# 统计检验（修正版）
sensitivity_test <- sensitivity_data %>%
  pivot_wider(
    names_from = Method,
    values_from = sensitivity
  ) %>%
  group_by(scenario) %>%
  summarise(
    t_test = t.test(`BIRTH design`, `1:1 trial`, paired = TRUE) %>% broom::tidy(),
    .groups = "drop"
  ) %>%
  unnest(t_test) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    scenario = factor(scenario, levels = scenario_levels)
  )

ggplot(overall_sensitivity, aes(x = scenario, y = mean_sensitivity, fill = Method)) +
  geom_col(
    position = position_dodge(width = 0.55),
    width = 0.5,
    color = "black",
    alpha = 0.9,
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(ymin = mean_sensitivity - 1.96*se_sensitivity,
        ymax = mean_sensitivity + 1.96*se_sensitivity),
    position = position_dodge(width = 0.55),
    width = 0.2,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = sensitivity_test,
    aes(x = scenario, y = 0.033, label = significance),
    inherit.aes = FALSE,
    size = 3.2,
    color = "red",
    fontface = "bold",
    vjust = 0.5
  ) +
  scale_fill_manual(values = c("BIRTH design" = "#D53E4F", "1:1 trial" = "#3288BD")) +
  scale_x_discrete(
    name = "Scenario",
    position = "bottom"  # 确保标签在底部
  ) +
  scale_y_continuous(
    name = "Mean Sensitivity Coefficient\n(ΔBias/ΔlogR0)",  # 换行处理
    limits = c(-0.04, 0.04),
    breaks = seq(-0.04, 0.04, 0.01),
    expand = c(0, 0),
    labels = scales::number_format(accuracy = 0.01)  # 优化标签显示
  ) +
  labs(
    title = "Method Sensitivity Comparison",
    subtitle = "BIRTH design demonstrates significantly higher parameter sensitivity",
    caption = "Error bars: 95% CI | Significance: ***p<0.001, **p<0.01, *p<0.05"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(t = -0.2, b = -0.4, unit = "cm"),  # 压缩图例边距
    legend.text = element_text(margin = margin(r = 0.2, unit = "cm")),
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 0,  # 不倾斜
      hjust = 0.5,  # 水平居中对齐
      vjust = 0.5,  # 垂直居中对齐
      face = "bold",
      margin = margin(t = 0.2, unit = "cm")  # 调整与轴的距离
    ),
    axis.title.y = element_text(
      margin = margin(r = 0.2, unit = "cm"),
      vjust = 2  # 调整标题位置
    ),
    axis.title.x = element_text(
      margin = margin(t = 0.2, unit = "cm")
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = rel(1.1),
      margin = margin(b = 0.1, unit = "cm")
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = rel(0.9),
      margin = margin(b = 0.3, unit = "cm")
    ),
    plot.caption = element_text(
      hjust = 0,
      lineheight = 0.9,
      size = rel(0.75),
      margin = margin(t = 0.2, unit = "cm")
    ),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm")  # 极简边距
  )

# 保存设置
ggsave("C:/Users/dell/Desktop/Work/文章整理结果/result_picture/sensitivity_1.png", 
       device = "png",
       width = 24.5, 
       height = 6.4,
       units = "cm",
       dpi = 600)  # 提高DPI保证清晰度
