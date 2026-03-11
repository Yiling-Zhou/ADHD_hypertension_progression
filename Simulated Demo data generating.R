# demo_generate_small_data.R
# Creates a very small simulated demo dataset for Step2 and microsimulation
# Authour: GPT 5.4

library(data.table)

set.seed(20260311)

dir.create("processData", showWarnings = FALSE)
dir.create("processData/adhd", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

n <- 30
pin <- sprintf("P%04d", 1:n)

# -----------------------------
# cohort_admin.csv
# -----------------------------
cohort_admin <- data.table(
  pin = pin,
  sex = sample(c(0L, 1L), n, replace = TRUE),
  birth_yr = sample(1945:2002, n, replace = TRUE),
  birth_mn = sample(1:12, n, replace = TRUE),
  emigration_date = as.IDate(NA),
  death_date = as.IDate(NA),
  SES = sample(c("low", "middle", "high"), n, replace = TRUE),
  otherIndication = 0L,
  
  # psychiatric comorbidity indicators
  asd = rbinom(n, 1, 0.08),
  sud = rbinom(n, 1, 0.10),
  scz = rbinom(n, 1, 0.03),
  dep = rbinom(n, 1, 0.15),
  anx = rbinom(n, 1, 0.18),
  bip = rbinom(n, 1, 0.04),
  ed  = rbinom(n, 1, 0.03),
  pd  = rbinom(n, 1, 0.06),
  intellectual = rbinom(n, 1, 0.03),
  
  # medication covariates
  add_med = rbinom(n, 1, 0.05),
  anxio_med = rbinom(n, 1, 0.08),
  dep_med = rbinom(n, 1, 0.12),
  ep_med = rbinom(n, 1, 0.04),
  hy_se_med = rbinom(n, 1, 0.05),
  psycho_med = rbinom(n, 1, 0.06),
  
  beta = rbinom(n, 1, 0.15),
  t2d = rbinom(n, 1, 0.10),
  hyperlipid = rbinom(n, 1, 0.20)
)

# some emigration/death after baseline period for realism
cohort_admin[sample(.N, 2), emigration_date := as.IDate(sample(seq(as.Date("2019-01-01"), as.Date("2022-06-30"), by = "day"), 2))]
cohort_admin[sample(.N, 2), death_date := as.IDate(sample(seq(as.Date("2020-01-01"), as.Date("2022-12-15"), by = "day"), 2))]

fwrite(cohort_admin, "processData/cohort_admin.csv")

# -----------------------------
# rx_hyper_med.csv
# -----------------------------
eligible_atc <- c("C03A", "C03B", "C08C", "C08G", "C09A", "C09B", "C09C", "C09D")
all_hyper <- rbindlist(lapply(pin, function(id) {
  baseline <- sample(seq(as.Date("2013-01-01"), as.Date("2020-12-31"), by = "day"), 1)
  
  # 1 baseline fill + 1 to 5 follow-up fills
  k <- sample(2:6, 1)
  dates <- sort(c(
    baseline,
    baseline + cumsum(sample(45:180, k - 1, replace = TRUE))
  ))
  
  # occasional pre-baseline antihypertensive history for exclusion
  if (runif(1) < 0.15) {
    dates <- c(baseline - sample(200:1200, 1), dates)
  }
  
  data.table(
    pin = id,
    disp_date = as.IDate(dates),
    atc = sample(eligible_atc, length(dates), replace = TRUE)
  )
}))

fwrite(all_hyper, "processData/rx_hyper_med.csv")

# -----------------------------
# rx_adhd_med.csv
# -----------------------------
adhd_people <- sample(pin, 8)

rx_adhd_med <- rbindlist(lapply(adhd_people, function(id) {
  first_date <- sample(seq(as.Date("2008-01-01"), as.Date("2022-06-30"), by = "day"), 1)
  k <- sample(1:4, 1)
  dates <- sort(first_date + cumsum(c(0, sample(30:200, k - 1, replace = TRUE))))
  data.table(pin = id, adhd_disp_date = as.IDate(dates))
}))

fwrite(rx_adhd_med, "processData/rx_adhd_med.csv")

# -----------------------------
# dx_adhd_diag.csv
# -----------------------------
dx_adhd_diag <- rbindlist(lapply(sample(pin, 10), function(id) {
  data.table(
    pin = id,
    adhd_diag_date = as.IDate(sample(seq(as.Date("2011-01-01"), as.Date("2022-06-30"), by = "day"), 1))
  )
}))
dx_adhd_diag <- unique(dx_adhd_diag, by = c("pin", "adhd_diag_date"))

fwrite(dx_adhd_diag, "processData/dx_adhd_diag.csv")

# -----------------------------
# dx_outcomes.csv
# -----------------------------
make_optional_date <- function(prob, start = "2013-01-01", end = "2022-12-31") {
  if (runif(1) < prob) {
    as.IDate(sample(seq(as.Date(start), as.Date(end), by = "day"), 1))
  } else {
    as.IDate(NA)
  }
}

dx_outcomes <- data.table(
  pin = pin,
  stroke_date = as.IDate(NA),
  AMI_date = as.IDate(NA),
  HHF_date = as.IDate(NA),
  CKD_date = as.IDate(NA),
  ESKD_date = as.IDate(NA),
  unnaturaldeath_date = as.IDate(NA),
  CRdeath_date = as.IDate(NA),
  noCRdeath_date = as.IDate(NA)
)

# generate some post-baseline events
for (i in 1:n) {
  base_i <- min(all_hyper[pin == dx_outcomes$pin[i], disp_date])
  if (!is.finite(base_i)) next
  
  if (runif(1) < 0.15) dx_outcomes[i, CKD_date := as.IDate(base_i + sample(100:1200, 1))]
  if (runif(1) < 0.12) dx_outcomes[i, stroke_date := as.IDate(base_i + sample(120:1200, 1))]
  if (runif(1) < 0.10) dx_outcomes[i, AMI_date := as.IDate(base_i + sample(120:1200, 1))]
  if (runif(1) < 0.10) dx_outcomes[i, HHF_date := as.IDate(base_i + sample(120:1200, 1))]
}

# death can occur after an event
for (i in 1:n) {
  evt_dates <- unlist(dx_outcomes[i, .(stroke_date, AMI_date, HHF_date, CKD_date)])
  evt_dates <- evt_dates[!is.na(evt_dates)]
  if (length(evt_dates) > 0 && runif(1) < 0.15) {
    dx_outcomes[i, CRdeath_date := min(evt_dates) + sample(30:600, 1)]
  } else if (runif(1) < 0.12) {
    dx_outcomes[i, noCRdeath_date := as.IDate(sample(seq(as.Date("2015-01-01"), as.Date("2022-12-15"), by = "day"), 1))]
  }
}

fwrite(dx_outcomes, "processData/dx_outcomes.csv")

message("Demo files written to processData/")