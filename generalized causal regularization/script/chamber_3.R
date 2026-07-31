# Alternative Causal Chamber analysis using a strong t_vis_2 intervention and
# an intervention on the direct cause red. Outputs use suffix _3.

Sys.setenv(
  CHAMBER_ENVIRONMENT_1 = "uniform_t_vis_2_strong.csv",
  CHAMBER_ENVIRONMENT_2 = "uniform_red_mid.csv",
  CHAMBER_OUTPUT_SUFFIX = "_3",
  CHAMBER_EXCLUDE_PREDICTORS = "t_vis_2",
  CHAMBER_DIRECT_Y_MIN = "-1800"
)
source("script/chamber.R")
