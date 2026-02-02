
###################################################################################################################
# The following script was written by Jiang Li
# The purpose is to generate a Kaplan-Meier curve in the specific NCI-requested layout style. 

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(gridExtra)
})

km_cuminc_save_survminer <- function(
    data, time_col, event_col, group_col,
    xlim = c(0, 90),
    outfile_prefix = "km_plot",
    legend_title = group_col,
    legend_labels = NULL,
    title = "Time-to-Event",
    palette = c('#E7B900', '#2E9FDF')
) {
  stopifnot(is.data.frame(data))
  df <- data
  # Prepare columns (days -> years)
  df$time_years <- as.numeric(df[[time_col]]) / 365.25
  df$event_ind  <- as.numeric(df[[event_col]])
  if (!is.factor(df[[group_col]])) df[[group_col]] <- factor(df[[group_col]])
  df$group_var  <- df[[group_col]]
  if (is.null(legend_labels)) legend_labels <- levels(df$group_var)
  
  # Fit KM
  fit <- survival::survfit(Surv(time_years, event_ind) ~ group_var, data = df)
  
  # Log-rank p-value
  sd <- survival::survdiff(Surv(time_years, event_ind) ~ group_var, data = df)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  
  # Try modern fun="event"
  fun_arg <- "event"
  test <- try(survminer::ggsurvplot(fit, data = df, fun = fun_arg, pval = FALSE), silent = TRUE)
  if (inherits(test, "try-error")) fun_arg <- function(s) 1 - s
  
  gsp <- survminer::ggsurvplot(
    fit,
    data = df,
    fun = fun_arg,                     # cumulative probability
    xlim = xlim,
    xlab = "Age (Years)",
    ylab = "Cumulative Probability",
    break.x.by = 10,
    break.y.by = 0.2,
    surv.scale = "percent",
    palette = palette,
    censor = FALSE,
    conf.int = TRUE,
    pval = FALSE,                      # we annotate our own p-value
    legend.title = legend_title,
    legend.labs  = legend_labels,
    font.legend  = c(18, "plain", "black"),
    risk.table = TRUE,                 # simple, version-stable
    risk.tables.height = 0.30,
    risk.tables.y.text = TRUE,
    risk.tables.fontsize = 2,
    title = title,
    ggtheme = theme_light(),
    font.x = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black")
  )
  
  # Annotate p-value near top-right
  xy <- c(xlim[2] - 0.05 * diff(xlim), 0.95)
  gsp$plot <- gsp$plot +
    annotate("text", x = xy[1], y = xy[2],
             label = paste("p =", format(pval, digits = 3)),
             size = 5, color = "black")
  
  # Combine plot + table and save
  combined <- gridExtra::grid.arrange(gsp$plot, gsp$table, ncol = 1, heights = c(3, 1))
  tiff_file <- sprintf("%s_%dto%d.tiff", outfile_prefix, xlim[1], xlim[2])
  pdf_file  <- sprintf("%s_%dto%d.pdf",  outfile_prefix, xlim[1], xlim[2])
  ggsave(tiff_file, plot = combined, device = "tiff", dpi = 150, width = 14, height = 8)
  ggsave(pdf_file,  plot = combined, device = "pdf",  dpi = 150, width = 14, height = 8)
  
  invisible(list(fit = fit, pval = pval, gsp = gsp, files = c(tiff_file, pdf_file)))
}


##################################################################################################

# # ---------- Generate the dataset ----------
# TERT_syn <- generate_TERT_dataset(
#   n = 160864,
#   carriers_n = 138,         # matches your earlier cohort
#   seed = 20260317,
#   mean_followup_age = 66,
#   sd_followup_age   = 14,
#   base_hazard = 0.0016,
#   HR_carrier  = 1.3
# )

# Peek
load("/path-to-synthetic-data/Example_TERT_synthetic_dataframe.RData")
head(TERT_syn, 10)

# PT_ID carrier current_age    sex death_yn   bmi tobacco_use_ever alcohol_use_ever Genetic_Race frz175 deceased_yn Male
# 1    PT1       0       79.97   Male       no 30.64                0                0          EUR     NA           0    1
# 2    PT2       0       59.13 Female       no 35.74                1                1          EUR     NA           0    0
# 3    PT3       0       59.95   Male       no 19.44                0                0          EUR    442           0    1
# 4    PT4       0       60.26   Male       no 17.63                0                1          EUR     NA           0    1
# 5    PT5       0       94.86 Female       no 25.37                1                1          EUR   8330           0    0
# 6    PT6       0       79.49 Female       no 29.20                1                1          EUR     NA           0    0
# 7    PT7       0       78.07 Female       no 28.98                1                0          EUR     NA           0    0
# 8    PT8       0       67.69 Female       no 23.00                0                0          AMR     NA           0    0
# 9    PT9       0       51.53   Male       no 20.88                0                1          EUR     NA           0    1
# 10  PT10       0       52.24 Female       no 32.39                0                1          EUR  29110           0    0
# timediff_death_censor event_death_censor timediff_death_censor_80yr event_death_censor_80yr timediff_death_censor_85yr
# 1               29207.59                  0                   29207.59                       0                   29207.59
# 2               21596.71                  0                   21596.71                       0                   21596.71
# 3               21897.44                  0                   21897.44                       0                   21897.44
# 4               22010.83                  0                   22010.83                       0                   22010.83
# 5               34647.42                  0                   29220.00                       0                   31046.25
# 6               29033.44                  0                   29033.44                       0                   29033.44
# 7               28514.21                  0                   28514.21                       0                   28514.21
# 8               24722.20                  0                   24722.20                       0                   24722.20
# 9               18819.82                  0                   18819.82                       0                   18819.82
# 10              19080.50                  0                   19080.50                       0                   19080.50
# event_death_censor_85yr timediff_death_censor_90yr event_death_censor_90yr
# 1                        0                   29207.59                       0
# 2                        0                   21596.71                       0
# 3                        0                   21897.44                       0
# 4                        0                   22010.83                       0
# 5                        0                   32872.50                       0
# 6                        0                   29033.44                       0
# 7                        0                   28514.21                       0
# 8                        0                   24722.20                       0
# 9                        0                   18819.82                       0
# 10                       0                   19080.50                       0

# 0–90 panel
km_cuminc_save_survminer(
  data = TERT_syn,
  time_col = "timediff_death_censor_90yr",
  event_col = "event_death_censor_90yr",
  group_col = "carrier",
  xlim = c(0, 90),
  outfile_prefix = "km_90yrdeath_TERT_alltiersvscontrol_03172025",
  legend_title = "TERT",
  legend_labels = c("noncarrier", "carrier"),
  title = "Time-to-death"
)

# 30–90 panel
km_cuminc_save_survminer(
  data = TERT_syn,
  time_col = "timediff_death_censor_90yr",
  event_col = "event_death_censor_90yr",
  group_col = "carrier",
  xlim = c(30, 90),
  outfile_prefix = "km_90yrdeath_TERT_alltiersvscontrol_03172025",
  legend_title = "TERT",
  legend_labels = c("noncarrier", "carrier"),
  title = "Time-to-death"
)

# combined <- km_cuminc_save_survminer()
# save the file in two formats with the designated figure layout 
ggsave(tiff_file, plot = combined, device = "tiff", dpi = 300, width = 14, height = 8)
ggsave(pdf_file,  plot = combined, device = "pdf",  dpi = 300, width = 14, height = 8)
