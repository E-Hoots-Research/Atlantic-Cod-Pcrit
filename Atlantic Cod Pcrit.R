rm(list=ls())


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pkgs <- c(
  # ----- Data Visualisation -----
  "ggthemes", "bayesplot", "gt", "gtsummary", "plotly", "qqplotr", "gridExtra", 
  "ggplot2",
  
  # ----- Tidy Data and Wrangling -----
  "tidyverse", "janitor", "readxl", "broom.mixed", "data.table", "hms", "devtools",
  "mclust", "readxl", "dplyr",
  
  # ----- Modelling and Statistical Analysis -----
  "brms", "rstan", "marginaleffects", "performance", "emmeans",
  "tidybayes", "respirometry", "future", "furrr"
)

# ---- Install and load all packages using pacman ----
suppressPackageStartupMessages(
  pacman::p_load(char = pkgs, install = TRUE)
)

#read in data
rawdata <- read_xlsx("C:/Users/bhoot/OneDrive - Deakin University/Shared Lab Space/K-berg 2025/pcrit_tidydata.xlsx", sheet = "master")
metadata <- read_xlsx("C:/Users/bhoot/OneDrive - Deakin University/Shared Lab Space/K-berg 2025/pcrit_tidydata.xlsx", sheet = "meta")

#calcSMR function
calcSMR <- function(Y) {
  u <- sort(Y)
  tenpc <- round(0.1 * length(u))
  SD10pc <- sd(u[1:tenpc])
  low10pc = mean(u[(which((u > (mean(u[1:tenpc])-SD10pc)))):(tenpc+which((u > (mean(u[1:tenpc])-SD10pc -u[1]))))])
  return(list(low10pc = low10pc))
}

#visualize data
ggplot(data = rawdata, aes(x = PO2_mg_perL, y = MO2_mg_perL_persec, color = factor(Stage))) +
  geom_point() +
  geom_line(aes(group = ID_fish)) +
  facet_grid(.~Temp) +
  scale_y_reverse() +
  scale_x_reverse() +
  theme_bw()

#mass range
ggplot(data = metadata, aes(x = 1, y = mass, fill = factor(treatment))) +
  geom_boxplot() +
  geom_point() +
  facet_grid(.~treatment) +
  scale_fill_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  theme_bw()

#merge raw and metadata
data <- left_join(rawdata, metadata %>% dplyr::select(ID_fish, mass, length, chamber_vol, backresp), by = "ID_fish") %>%
  mutate(Stage = as.factor(Stage),
         MO2_mg_permin = (abs(MO2_mg_perL_persec)-abs(backresp))*((chamber_vol - mass)/1000)*60)

data_mmr <- data %>%
  group_by(ID_fish, mass) %>%
  arrange(ID_fish) %>%
  summarise(
    MMR = max(MO2_mg_permin)
  ) 

data_smr <- data[which(data$Stage == "smr"),] %>%
  group_by(ID_fish, mass) %>%
  arrange(ID_fish) %>%
  summarise(
    SMR = calcSMR(MO2_mg_permin)$low10pc %>% unname()
    )

data_rmr <- data[which(data$Stage == "smr"),] %>%
  group_by(ID_fish) %>%
  arrange(DateTime) %>%
  slice(9:(n() - 1)) %>%
  ungroup()  %>%
  group_by(
    ID_fish, mass) %>% 
  arrange(ID_fish) %>%
  summarise(
    RMR = mean(MO2_mg_permin)
  )

data1 <- left_join(data, data_mmr %>% dplyr::select(ID_fish, MMR), by = "ID_fish")
data2 <- left_join(data1, data_smr %>% dplyr::select(ID_fish, SMR), by = "ID_fish")
data3 <- left_join(data2, data_rmr %>% dplyr::select(ID_fish, RMR), by = "ID_fish")

data3_subset <- data3[!data3$ID_fish %in% c("A2_20h_2", "B1_20h_2", "B2_20h_2"), ]

ggplot(data = data3_subset, aes(x = PO2_mg_perL, y = MO2_mg_permin, color = factor(Treatment))) +
  geom_point() +
  facet_wrap(~ID_fish, scales = "free_y") +
  geom_line(aes(y = SMR), color = "black", linetype = "dashed") +
  geom_line(aes(y = RMR), color = "red", linetype = "dashed") +
  theme_bw() +
  scale_x_reverse() +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1"))+
  geom_text(aes(x = 6.5, y = 0.1, label = paste(mass, "g")),
            hjust = 1, vjust = 0, size = 3, color = "black", show.legend = FALSE) +
  labs(x = "Tank PO2 (mg/L)", y = "Metabolic Rate (mgO2/min)", color = "Temp (degC)")

##########################################################################################

# Set output directory
output_figs <- "C:/Users/bhoot/OneDrive - Deakin University/Shared Lab Space/K-berg 2025/Cod pcrit plots"
setwd(output_figs)

fish_dfs <- split(data3_subset, data3_subset$ID_fish)

# Initialize results data frame with correct column names
pcrit_results <- data.frame(ID_fish = character(),
                            Alpha = numeric(),
                            Breakpoint = numeric(),
                            LLO = numeric(),
                            NLR = numeric(),
                            stringsAsFactors = FALSE)

# Loop over each fish dataset
for (i in seq_along(fish_dfs)) {
  test <- fish_dfs[[i]]
  fish_id <- unique(test$ID_fish)
  
  # Run pcrit calculation
  pcrit_vals <- calc_pcrit(po2 = test$PO2_mg_perL, 
                           mo2 = test$MO2_mg_permin, 
                           MR = test$SMR[1], #sets the MR for the estimate (Alpha and LLO only) to SMR
                           method = "All") #shows BSR, LLO, NLR, Alpha, and sub-PI estimation methods
  
  # Store results
  result_row <- data.frame(ID_fish = fish_id,
                           Alpha = pcrit_vals["Alpha"], #Alpha-based PO2
                           Breakpoint = pcrit_vals["Breakpoint"], #Broken stick regression
                           LLO = pcrit_vals["LLO"], #Linear low O2
                           NLR = pcrit_vals["NLR"]) #Non-linear regression
  pcrit_results <- rbind(pcrit_results, result_row)
  
  # Save plot
  png(filename = paste0("plot_", fish_id, ".png"), width = 1600, height = 1600, res = 150)
  plot_pcrit(po2 = test$PO2_mg_perL, mo2 = test$MO2_mg_permin, MR = test$SMR[1], method = "All")
  dev.off()
}

#calc_pcrit does not work perfectly with all of these data. If error received, run from here:

data4 <- left_join(data3, pcrit_results, by = "ID_fish") %>% 
  dplyr::select(-LLO)

data5 <- data4 %>% 
  dplyr::select(-MO2_mg_perL_persec, -MO2_mg_permin, -PO2_mg_perL, -DateTime, -Stage, -chamber_vol) %>% 
  distinct(ID_fish, .keep_all = TRUE)

data4_subset <- data4[!data4$ID_fish %in% c("A2_20h_2", "B1_20h_2", "B2_20h_2"), ]

output_dir <- "C:/Users/bhoot/OneDrive - Deakin University/Shared Lab Space/K-berg 2025"
setwd(output_dir)

write.csv(data5, file = "pcrit_results.csv", col.names = TRUE, row.names = FALSE)

###########################
# Pcrit visual assessment 
###########################

ggplot(data = data4_subset, aes(x = PO2_mg_perL, y = MO2_mg_permin, color = factor(Treatment))) +
  geom_point() +
  facet_wrap(~ID_fish, scales = "free") +
  geom_line(aes(y = SMR), color = "red", linetype = "dashed") +
  geom_line(aes(y = RMR), color = "orangered", linetype = "dashed") +
  geom_line(aes(y = MMR), color = "hotpink", linetype = "dashed") +
  geom_line(aes(x = Alpha), color = "limegreen") +
  theme_bw() +
  scale_x_reverse() +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1"))+
  geom_text(aes(x = 6.5, y = 0.1, label = paste(mass, "g")),
            hjust = 0.5, vjust = 0, size = 3, color = "black", show.legend = FALSE) +
  labs(x = "Tank PO2 (mg/L)", y = "Metabolic Rate (mgO2/min)", color = "Temp (degC)")

#by treatment, with MMR, RMR, SMR, and pcrit indicated:

#15C normoxia
Figure1 <- ggplot(data = subset(data4_subset, Treatment == "15_norm"), aes(x = PO2_mg_perL, y = MO2_mg_permin)) +
  geom_point(color = "skyblue3") +
  facet_wrap(~ID_fish, scales = "free") +
  geom_line(aes(y = SMR), color = "red", linetype = "dashed") +
  geom_line(aes(y = RMR), color = "orangered", linetype = "dashed") +
  geom_line(aes(y = MMR), color = "hotpink", linetype = "dashed") +
  geom_line(aes(x = Alpha), color = "limegreen") +
  theme_bw() +
  scale_x_reverse() +
  geom_text(aes(x = 6.5, y = 0.1, label = paste(mass, "g")),
            hjust = 0.5, vjust = 0, size = 3, color = "black", show.legend = FALSE) +
  labs(x = expression(Tank~PO[2]~(mg~O[2]~l^{-1})), 
       y = expression(Metabolic~rate~(mg~O[2]~min^{-1})))

#15C hyperoxia
Figure2 <- ggplot(data = subset(data4_subset, Treatment == "15_hyper"), aes(x = PO2_mg_perL, y = MO2_mg_permin)) +
  geom_point(color = "cadetblue3") +
  facet_wrap(~ID_fish, scales = "free") +
  geom_line(aes(y = SMR), color = "red", linetype = "dashed") +
  geom_line(aes(y = RMR), color = "orangered", linetype = "dashed") +
  geom_line(aes(y = MMR), color = "hotpink", linetype = "dashed") +
  geom_line(aes(x = Alpha), color = "limegreen") +
  theme_bw() +
  scale_x_reverse() +
  geom_text(aes(x = 6.5, y = 0.1, label = paste(mass, "g")),
            hjust = 0.5, vjust = 0, size = 3, color = "black", show.legend = FALSE) +
  labs(x = expression(Tank~PO[2]~(mg~O[2]~l^{-1})), 
       y = expression(Metabolic~rate~(mg~O[2]~min^{-1})))

#20C normoxia
Figure3 <- ggplot(data = subset(data4_subset, Treatment == "20_norm"), aes(x = PO2_mg_perL, y = MO2_mg_permin)) +
  geom_point(color = "coral1") +
  facet_wrap(~ID_fish, scales = "free") +
  geom_line(aes(y = SMR), color = "red", linetype = "dashed") +
  geom_line(aes(y = RMR), color = "orangered", linetype = "dashed") +
  geom_line(aes(y = MMR), color = "hotpink", linetype = "dashed") +
  geom_line(aes(x = Alpha), color = "limegreen") +
  theme_bw() +
  scale_x_reverse() +
  geom_text(aes(x = 6.5, y = 0.1, label = paste(mass, "g")),
            hjust = 0.5, vjust = 0, size = 3, color = "black", show.legend = FALSE) +
  labs(x = expression(Tank~PO[2]~(mg~O[2]~l^{-1})), 
       y = expression(Metabolic~rate~(mg~O[2]~min^{-1})))

#20C hyperoxia
Figure4 <- ggplot(data = subset(data4_subset, Treatment == "20_hyper"), aes(x = PO2_mg_perL, y = MO2_mg_permin)) +
  geom_point(color = "pink") +
  facet_wrap(~ID_fish, scales = "free") +
  geom_line(aes(y = SMR), color = "red", linetype = "dashed") +
  geom_line(aes(y = RMR), color = "orangered", linetype = "dashed") +
  geom_line(aes(y = MMR), color = "hotpink", linetype = "dashed") +
  geom_line(aes(x = Alpha), color = "limegreen") +
  theme_bw() +
  scale_x_reverse() +
  geom_text(aes(x = 6.5, y = 0.1, label = paste(mass, "g")),
            hjust = 0.5, vjust = 0, size = 3, color = "black", show.legend = FALSE) +
  labs(x = expression(Tank~PO[2]~(mg~O[2]~l^{-1})), 
       y = expression(Metabolic~rate~(mg~O[2]~min^{-1})))

####################
# Allometric plots
####################

#plot of SMR vs. mass
Figure5 <- ggplot(data = data5, aes(x = mass, y = SMR, color = factor(Treatment), fill = factor(Treatment))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(.~Temp) +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  scale_fill_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  theme_bw() +
  guides(fill = FALSE) +
  labs(
    x = expression(paste("Mass (g)")),
    y = expression(SMR~(mg~O[2]~min^{-1})),
    color = expression(paste("Treatment (", degree, "C, ", PO[2], ")"))
  )+
  scale_x_log10() +
  scale_y_log10()

#smr scaling equation:
lm_smr_15n <- lm(log10(SMR) ~ log10(mass), data = subset(data5, Treatment == "15_norm"))
lm_smr_15h <- lm(log10(SMR) ~ log10(mass), data = subset(data5, Treatment == "15_hyper"))
lm_smr_20n <- lm(log10(SMR) ~ log10(mass), data = subset(data5, Treatment == "20_norm"))
lm_smr_20h <- lm(log10(SMR) ~ log10(mass), data = subset(data5, Treatment == "20_hyper"))

#plot of RMR vs. mass
Figure6 <- ggplot(data = data5, aes(x = mass, y = RMR, color = factor(Treatment), fill = factor(Treatment))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(.~Temp) +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  scale_fill_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  theme_bw() +
  guides(fill = FALSE) +
  labs(
    x = expression(paste("Mass (g)")),
    y = expression(RMR~(mg~O[2]~min^{-1})),
    color = expression(paste("Treatment (", degree, "C, ", PO[2], ")"))
  )+
  scale_x_log10() +
  scale_y_log10()

#rmr scaling equation:
lm_rmr_15n <- lm(log10(RMR) ~ log10(mass), data = subset(data5, Treatment == "15_norm"))
lm_rmr_15h <- lm(log10(RMR) ~ log10(mass), data = subset(data5, Treatment == "15_hyper"))
lm_rmr_20n <- lm(log10(RMR) ~ log10(mass), data = subset(data5, Treatment == "20_norm"))
lm_rmr_20h <- lm(log10(RMR) ~ log10(mass), data = subset(data5, Treatment == "20_hyper"))

#plot of MMR vs. mass
Figure7 <- ggplot(data = data5, aes(x = mass, y = MMR, color = factor(Treatment), fill = factor(Treatment))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(.~Temp) +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  scale_fill_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  theme_bw() +
  guides(fill = FALSE) +
 labs(x = expression(paste("Mass (g)")),
      y = expression(MMR~(mg~O[2]~min^{-1})),
      color = expression(paste("Treatment (", degree, "C, ", PO[2], ")"))) +
  scale_x_log10() +
  scale_y_log10()

#mmr scaling equation:
lm_mmr_15n <- lm(log10(MMR) ~ log10(mass), data = subset(data5, Treatment == "15_norm"))
lm_mmr_15h <- lm(log10(MMR) ~ log10(mass), data = subset(data5, Treatment == "15_hyper"))
lm_mmr_20n <- lm(log10(MMR) ~ log10(mass), data = subset(data5, Treatment == "20_norm"))
lm_mmr_20h <- lm(log10(MMR) ~ log10(mass), data = subset(data5, Treatment == "20_hyper"))

#plot of pcrit (breakpoint) vs. mass
ggplot(data = data5, aes(x = mass, y = Breakpoint, color = factor(Treatment), fill = factor(Treatment))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(.~Temp) +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  scale_fill_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  theme_bw() +
  guides(fill = FALSE) +
  labs(
    x = expression(mass~(g)),
    y = expression(P[crit]~"by breakpoint"~(mg~l^{-1})),
    color = expression(paste("Treatment (", degree, "C, ", PO[2], ")"))) +
  scale_x_log10() +
  scale_y_log10()


#plot of pcrit (alpha) vs. mass
Figure8 <- ggplot(data = data5, aes(x = mass, y = Alpha, color = factor(Treatment), fill = factor(Treatment))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(.~Temp) +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  scale_fill_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  theme_bw() +
  guides(fill = FALSE) +
  labs(
    x = expression(mass~(g)),
    y = expression(P[CRIT]~"by alpha"~(mg~l^{-1})),
    color = expression(paste("Treatment (", degree, "C, ", PO[2], ")"))) +  
  scale_x_log10() +
  scale_y_log10()

#add LOE column for each fish, giving the PO2 where the pcrit trial ended
data6 <- data4_subset %>%
  group_by(ID_fish, mass, Treatment, Temp) %>%
  arrange(ID_fish) %>%
  summarise(
    LOE = min(PO2_mg_perL)
  ) 

#plot of loe vs. mass
Figure9 <- ggplot(data = data6, aes(x = mass, y = LOE, color = factor(Treatment), fill = factor(Treatment))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(.~Temp) +
  scale_color_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  scale_fill_manual(values = c("cadetblue3", "skyblue3", "pink", "coral1")) +
  theme_bw() +
  guides(fill = FALSE) +
  labs(
    x = expression(mass~(g)),
    y = expression(LOE~(mg~l^{-1})),
    color = expression(Treatment~"(" * degree * "C, " * PO[2] * ")")) +  
  scale_x_log10() +
  scale_y_log10()
