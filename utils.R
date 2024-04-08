library(compositions)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

# Constants
short_sleep_hours <- 6
hrs_in_day <- 24
mins_in_day <- 1440
mins_in_hour <- 60
sub_steps <- 6
sub_step_mins <- mins_in_hour / sub_steps
ncpus <- as.integer(Sys.getenv("NCPUS"))
bootstrap_iterations <- as.integer(Sys.getenv("BOOT_ITRS"))
maxit <- as.integer(Sys.getenv("MAXIT"))


## Define SBP
sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)

v <- gsi.buildilrBase(t(sbp))


## strip unneccessary model components  

strip_glm <- function(cm) {
  cm$y <- c()
  cm$model <- c()
  
  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()
  
  cm$family$variance <- c()
  cm$family$dev.resids <- c()
  cm$family$aic <- c()
  cm$family$validmu <- c()
  cm$family$simulate <- c()
  
  return(cm)
}

