library(compositions)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

# Constants
short_sleep_hours <- 6
hrs_in_day <- 24
mins_in_day <- 1440
sub_steps <- 12
sub_step_mins <- 5
ncpus <- as.integer(Sys.getenv("NCPUS"))
bootstrap_iterations <- as.integer(Sys.getenv("BOOT_ITRS"))

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
