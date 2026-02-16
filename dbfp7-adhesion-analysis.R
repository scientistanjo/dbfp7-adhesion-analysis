# AFM Adhesion Analysis and Plotting for Obille-2026.

# === DATA LOADING ===

file_list <- list.files(pattern = "\\.csv$")

file_df <- data.frame(filenames = file_list)

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(tidyverse)

file_df <- file_df %>%
  separate(filenames, into = c("Date", "SampleID", "ImageNum", "ExportID", "DataType"), 
           sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
  mutate(DataType = gsub("\\.csv$", "", DataType))

filenames_without_suffix <- sub("_[ahlm]\\.csv", "", file_df$filenames)

unique_filenames <- unique(filenames_without_suffix)

table.names.a <- c("x","y","a")
table.names.h <- c("x","y","h")
table.names.m <- c("x","y","m")

df_list <- list()

# Generate full table including height, adhesion, and modulus channels, joined by (x,y)
for (i in seq_along(unique_filenames)) {
  file_name <- unique_filenames[i]
  
  file_name_a <- read.csv(paste0(file_name, "_a.csv"), header = FALSE, sep = "")
  file_name_h <- read.csv(paste0(file_name, "_h.csv"), header = FALSE, sep = "")
  file_name_m <- read.csv(paste0(file_name, "_m.csv"), header = FALSE, sep = "")
  
  print(paste("file_name_a - Rows:", nrow(file_name_a), "Columns:", ncol(file_name_a)))
  print(paste("file_name_h - Rows:", nrow(file_name_h), "Columns:", ncol(file_name_h)))
  print(paste("file_name_m - Rows:", nrow(file_name_m), "Columns:", ncol(file_name_m)))
  
  colnames(file_name_a) <- table.names.a
  colnames(file_name_h) <- table.names.h
  colnames(file_name_m) <- table.names.m
  
  y <- full_join(file_name_h, file_name_a)
  alpha <- full_join(y, file_name_m)
  
  df_list[[file_name]] <- alpha
}

for (file_name in unique_filenames) {
  current_df <- df_list[[file_name]]
  current_df <- current_df %>% 
    mutate(hnew = h - min(h)) %>% 
    mutate(anew = a - min(a)) %>% 
    mutate(mnew = m - min(m))
  df_list[[file_name]] <- current_df
}

# === SUMMARY TABLES ===
list2env(df_list, envir = .GlobalEnv)

summary_tables <- list()

summary_tables <- lapply(unique_filenames, function(fname) {
  table <- get(fname) 
  summary_df <- as.data.frame(as.matrix(summary(table)))
  summary_df$unique_filename <- fname
  
  return(summary_df)
})

wlist <- lapply(summary_tables, function(df) {
  split_text <- strsplit(df$Freq, ":")
  
  df %>%
    mutate(statistic = sapply(split_text, `[`, 1),
           value = sapply(split_text, `[`, 2))
})

vtable.names <- c("variable", "statistic", "value", "unique_filename")

vlist <- lapply(wlist, function(w) {
  data.frame(variable = w$Var2, statistic = w$statistic, 
             value = as.numeric(w$value), unique_filename = w$unique_filename,
             stringsAsFactors = FALSE) %>%
    setNames(vtable.names)
})

vlist <- lapply(vlist, function(df) {
  df %>% filter(!is.na(statistic))
})

vlist_wd <- lapply(vlist, function(wd){
  wd %>% spread(key = statistic, value = value)
})

names(vlist_wd) <- unique_filenames

new_vlist_wd_colnames <- c("variable", "filename", "q1", "q3", "max", "mean", "median", "min") 

vlist_wd <- lapply(vlist_wd, function(df) {
  setNames(df, new_vlist_wd_colnames)
})

vlist_wd <- lapply(vlist_wd, function(df) {
  df %>% mutate(variable = str_trim(df$variable))
})

# === BOXPLOTS ===
box_a_datalist <- lapply(vlist_wd, function(df) {
  df %>% filter(variable == "a")
})

combined_box_a_data <- bind_rows(box_a_datalist, .id = "filename")

box_a_table <- data.frame(
  variable = combined_box_a_data$variable,
  filename = combined_box_a_data$filename,
  q1 = combined_box_a_data$q1,
  q3 = combined_box_a_data$q3,
  max = combined_box_a_data$max,
  mean = combined_box_a_data$mean,
  median = combined_box_a_data$median,
  min = combined_box_a_data$min
)

combined_image_data <- do.call(rbind, Map(cbind, df_list, index = seq_along(df_list)))
combined_image_data <- rownames_to_column(combined_image_data, var = "imagename")
combined_image_data$imagename <- sub("\\..*$", "", combined_image_data$imagename)
combined_image_data$imagename <- sub("^[^_]*_", "", combined_image_data$imagename)
combined_image_data$imagename <- factor(combined_image_data$imagename)
combined_image_data <- combined_image_data %>% 
  mutate(apn = a*1e12, hnm = h*1e9, apnnew = anew*1e12, hnmnew = hnew*1e9)
combined_image_data <- combined_image_data %>%
  mutate(sampleID = sub("_.*", "", imagename))
combined_image_data <- combined_image_data %>%
  mutate(exportID = sub(".*_", "", imagename))
combined_image_data <- combined_image_data %>%
  mutate(group = case_when(
    sampleID == "Dbfp7*" ~ "Dbfp7",
    sampleID == "BSA*" ~ "BSA",
    sampleID == "Blank" ~ "Mica",
    sampleID == "CT*" ~ "CellTak"
    TRUE ~ sampleID # keep original values if no match
  ))
unique_group <- unique(combined_image_data$group)
unique_exportID <- unique(combined_image_data$exportID)
cid_colnames <- colnames(combined_image_data)

adhf_viol <- ggplot(combined_image_data, aes(x = group, y = apn, fill = group)) + 
  geom_violin(position = position_dodge(width = 0.8)) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, coef = 0.5, fill = NA, color = NA) +
  labs(title = "Adhesion Force Violin Plots", x = "Sample", y = expression("Adhesion Force"~(pN))) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "gray"),
        axis.ticks = element_line(size = 0.5),
        text = element_text(size = 14, family = "Arial"),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
adhf_viol

custom_colors <- c('Dbfp7' = '#298791', 'BSA' = '#913329', 'CellTak' = '#F18F01', 'Mica' = 'grey')

# === ADHESION ENERGY CALC ===
# Need to iterate adhesion energy calculation according to the modulus
# First step is calculating the Tabor parameter
# Then for that pixel, take the corresponding DMT or JKR value and put into a column called adhe_fin

# Calculate the Tabor parameter
# Constants
R <- 20e-9 #m
W <- ((2*0.52e-9)/(3*pi*R)) #N/m
epsilon <- 0.15e-9 #nm
E1 <- 250e9 #Pa
nu1 <- 0.25
nu2 <- 0.5

# Tabor function
calc_tabor <- function(E2) {
  E_star_inv <- (1 - nu1^2) / E1 + (1 - nu2^2) / E2
  E_star <- 1 / E_star_inv
  x <- (R * W^2) / (E_star^2 * epsilon^3)
  tabor <- sqrt(sqrt(x))
  return(tabor)
}

# Apply to dataframe
combined_image_data <- combined_image_data %>%
  mutate(tabor = calc_tabor(mnew))

# calculate energy per contact area
calc_W2 <- function(E2){
  E_star_inv <- (1 - nu1^2) / E1 + (1 - nu2^2) / E2
  E_star <- 1 / E_star_inv
  a_cubed <- (0.75*(R/E_star)) * (0.52e-9 + (3*pi*R*W) + sqrt((6*pi*R*W*0.52e-9)+(3*pi*R*W)^2))
  a1 <- sqrt(sqrt(a_cubed))
  area <- pi*a1^2
  W2 <- 1.67e-17/area
  return(W2)
}

combined_image_data <- combined_image_data %>% 
  mutate(W2 = calc_W2(mnew))

calc_tabor2 <- Vectorize(function(E2, W2) {
  E_star_inv <- (1 - nu1^2) / E1 + (1 - nu2^2) / E2
  E_star <- 1 / E_star_inv
  x <- (R * W2) / (E_star^2 * epsilon^3)
  tabor2 <- sqrt(sqrt(x))
  return(tabor2)
})

combined_image_data <- combined_image_data %>%
  mutate(tabor2 = calc_tabor2(mnew, W2))

combined_long_tabor <- combined_image_data %>%
  pivot_longer(cols = c(tabor, tabor2), names_to = "TaborIter", values_to = "Value")

df_long_clean <- combined_long_tabor %>%
  filter(is.finite(Value))

taborviol <- ggplot(combined_image_data, aes(x = imagename, y = apn, group = imagename)) + 
  geom_violin(position = position_dodge(width = 0.8)) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA, color = NA) +
  labs(title = "Adhesion Box Plots", x = "Image Name", y = "Adhesion (pN)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "gray"),
        axis.ticks = element_line(size = 0.5),
        text = element_text(size = 14, family = "Arial"),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_hline(yintercept = c(-1e3), linetype = "solid", color = "red") +
  geom_hline(yintercept = c(1e3), linetype = "solid", color = "red") +
  geom_hline(yintercept = c(0), linetype = "solid", color = "green") +
  ylim(-2.5e3, 2.5e3)
taborviol

# number of points in JKR regime
sum(combined_image_data$tabor > 1)
sum(combined_image_data$tabor2 > 1, na.rm = TRUE)

# number of points in DMT regime
sum(combined_image_data$tabor < 0.1)
sum(combined_image_data$tabor2 < 0.1, na.rm = TRUE)

# proportion JKR
sum(combined_image_data$tabor > 1) / totalnumpoints
sum(combined_image_data$tabor2 > 1, na.rm = TRUE) / totalnumpoints

# proportion DMT
sum(combined_image_data$tabor < 0.1) / totalnumpoints
sum(combined_image_data$tabor2 < 0.1, na.rm = TRUE) / totalnumpoints

# proportion in transition regime
1 - ((sum(combined_image_data$tabor > 1) / totalnumpoints) + (sum(combined_image_data$tabor < 0.1) / totalnumpoints))

1 - ((sum(combined_image_data$tabor2 > 1, na.rm = TRUE) / totalnumpoints) + (sum(combined_image_data$tabor2 < 0.1, na.rm = TRUE) / totalnumpoints))

combined_image_data <- combined_image_data %>%
  mutate(model = case_when(
    tabor2 < 0.1 ~ "DMT",
    tabor2 > 1 ~ "JKR",
    TRUE ~ "Transition"
  ))

hnmnew_by_model <- ggplot(combined_image_data, aes(x = model, y = hnmnew, group = model)) + 
  geom_violin(position = position_dodge(width = 0.8)) + 
  labs(x = "Model", y = "Height (nm)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "gray"),
        axis.ticks = element_line(size = 0.5),
        text = element_text(size = 14, family = "Arial"),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")) +
  ylim(min(combined_image_data$hnmnew), 25)
hnmnew_by_model

# calculate DMT and JKR adhe

DMT20_function <- function(df) {
  df$Adh_DMT20 <- 1000*(df$apn*1e-12)/(2*pi*20e-9)
  return(df)
}

JKR20_function <- function(df) {
  df$Adh_JKR20 <- 1000*(df$apn*2e-12)/(3*pi*20e-9)
  return(df)
}

Maugis_function <- function(df) {
  df$Adh_Maugis <- 1000*(df$apn*1e-12)/(((7/4)+(-1/4)*(((4.04*(1.1570*df$tabor2)^1.4)-1)/((4.04*(1.1570*df$tabor2)^1.4)+1)) )*pi*20e-9)
  return(df)
}

functions_list <- list(DMT20_function, JKR20_function)
adhe_functions <- function(df){
  for (f in functions_list){
    df <- f(df)
  }
  return(df)
}

combined_image_data <- adhe_functions(combined_image_data)

combined_image_data <- combined_image_data %>%
  mutate(adhe = case_when(
    model == "DMT" ~ Adh_DMT20,
    model == "JKR" ~ Adh_JKR20,
    TRUE ~ NA_real_
  ))

adhe_viol <- ggplot(combined_image_data, aes(x = group, y = adhe)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(x = NULL, y = "Adhesion Energy") +
  ylim(-2.5, 10)
adhe_viol

combined_image_data <- Maugis_function(combined_image_data)

adhe_viol_Maugis <- ggplot(combined_image_data, aes(x = group, y = Adh_Maugis, fill = group)) +
  geom_violin(trim = FALSE) +
  labs(x = NULL, y = expression("Adhesion Energy"~(mJ/m^2))) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(text = element_text(size = 12, family = "Arial"),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  ylim(-2.5, 12.5)
adhe_viol_Maugis

adhf_plot <- ggplot(combined_image_data, aes(x = group, y = apn, fill = group)) +
  geom_violin(trim = FALSE) +
  labs(x = NULL, y = "Adhesion Force (pN)") +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(text = element_text(size = 12, family = "Arial"),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
        axis.line = element_line(colour = "grey"),
        axis.ticks = element_line(size = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  ylim(-200, 1200)
adhf_plot

adhe_box_Maugis <- ggplot(combined_image_data, aes(x = group, y = Adh_Maugis, fill = group)) +
  geom_boxplot(width = 0.8, outlier.shape = NA) +
  labs(x = NULL, y = expression("Adhesion Energy"~(mJ/m^2))) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    limits = c(-1.5, 12.5),
    breaks = c(0, 4, 8, 12),
    minor_breaks = c(2, 6, 10)
  ) +
  theme_classic() +
  guides(
    y = guide_axis(minor.ticks = TRUE) 
  ) +
  theme(text = element_text(size = 12, family = "Arial"),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
        axis.line = element_line(colour = "grey"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks.y.right = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none")
adhe_box_Maugis

df_long <- combined_image_data %>%
  pivot_longer(cols = c(Adh_DMT20, Adh_JKR20, Adh_Maugis),
               names_to = "AdhesionType",
               values_to = "Adhesion")

box_allgroup_alltype <- ggplot(df_long, aes(x = AdhesionType, y = Adhesion, fill = AdhesionType)) +
  geom_boxplot(width = 0.8, outlier.shape = NA) +
  facet_wrap(~ group) +
  theme_classic() +
  theme(text = element_text(size = 12, family = "Arial"),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
        axis.line = element_line(colour = "grey"),
        axis.ticks = element_line(size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Adhesion Model", y = expression("Adhesion Energy"~(mJ/m^2))) +
  ylim(-1.5, 12.5)
box_allgroup_alltype

height_on_adhe <- ggplot(combined_image_data, aes(x = Adh_Maugis, y = hnmnew, color = group)) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  theme(text = element_text(size = 12, family = "Arial"),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
        axis.line = element_line(colour = "grey"),
        axis.ticks = element_line(size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right") +
  xlim(-1.5, 12.5)
height_on_adhe

combined_image_data <- combined_image_data %>% mutate(
  normforce = (1000*(apn*1e-12))/(20e-9)
)

adhf_plot <- ggplot(combined_image_data,
                    aes(x = group, y = apn, fill = group)) +
  geom_violin(trim = FALSE) +
  labs(
    x = NULL,
    y = "Adhesion Force (pN)"
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    limits = c(-200, 1200),
    sec.axis = sec_axis(
      trans = ~ . * (1e-12*1000/20e-9),
      name  = "Normalized Force (mN/m)"
    )
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
    axis.line = element_line(colour = "grey"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none"
  )
adhf_plot