#' Simple Version of GuanRank
#'
#' Simple version of GuanRank (https://www.nature.com/articles/s43588-021-00083-2)
## the code here is similar to the code found in https://github.com/GuanLab/GuanRank-R
## but makes use of survival::survfit for computing survival probability for simplicity
#'
#' @param surv_data input should have two columns, the first column is time,
#' the second is status an indicator with 0 = censored, 1 = event
#'
#' @return
#' Returns `surv_data` data with an addition guan_rank variable of the GuanRank
#'
#' @export
#'
#' @examples
#' dat <- UCSCXenaTools::getTCGAdata(project = 'LUAD', clinical = TRUE, download = TRUE)
#' clin <- data.table::fread(dat$destfiles, data.table = FALSE)
#'
#' survData <- data.frame(
#'   time = ifelse(!is.na(clin$days_to_last_followup),
#'                 as.numeric(clin$days_to_last_followup),
#'                 as.numeric(clin$days_to_death)),
#'   status = dplyr::recode(
#'     clin$vital_status, 'LIVING' = 0,'DECEASED' = 1, .default = NA_real_
#'   )
#' )
#'
#' gr <- calculate_guan_rank(surv_data = survData)
#'
#' # Plotting code
#' Col <- as.vector(factor(gr$status, labels = c('slategrey', 'cyan')))
#' par(mfrow = c(2,2))
#' plot(survival::survfit(survival::Surv(gr$time, gr$status)~1), xlab = "Days",
#'      ylab = "Survival %", conf.int = FALSE, mark.time = TRUE, col = "slategrey",
#'      cex = .5, main = 'LUAD')
#' plot(gr$time,gr$guan_rank, pch = 19, col = Col, cex = .5, xlab = "Days",
#'      ylab = "Guan Rank", main = "Guan Rank vs Time")
#' boxplot(gr$guan_rank ~ gr$status, xlab = "Status", ylab = "Guan Rank",
#'         main = "Guan Rank vs Status")
#' hist(gr$guan_rank, breaks = 30, xlab = "Guan Rank",
#'      main = "Guan Rank Distribution")
#'
calculate_guan_rank <- function(surv_data){

  if (!all(c('time', 'status') %in% colnames(surv_data))) {
    stop('expecting columns "time" and "status"')
  }
  if (!all(surv_data$status %in% 0:1 | is.na((surv_data$status)))) {
    stop('status must be 0, 1, or missing')
  }
  if (!is.numeric(surv_data$status)) {
    stop('status must a numeric vector')
  }

  surv_data           <- data.frame(time = surv_data[,1],
                                    status = surv_data[,2],
                                    # adding dummy id to enable mapping back to orignal order
                                    id = 1:nrow(surv_data))
  # removing NA's since survfit will break with them
  tData               <- stats::na.omit(surv_data);
  # reorder to facilitate calculation below
  tData               <- tData[order(tData[,"time"]),]

  #compute survival probabilities by time
  km_fit              <- survival::survfit(survival::Surv(tData$time, tData$status)~1)
  tKM                 <- data.frame(time = km_fit$time, surv_prob = km_fit$surv)
  km_curve            <- dplyr::left_join(tData, tKM, "time")
  km_rank_mat         <- cbind(km_curve,"guan_rank" = rep(0,nrow(tData)))
  time                <- km_rank_mat[,"time"]
  status              <- km_rank_mat[,"status"]

  # note: the section below does not break down cleanly into cases 1 - 6 from
  #https://www.nature.com/articles/s43588-021-00083-2
  # instead it is a clever computation of rank using the same conditionals but
  #in a different order to ease the rank calculation.
  # it taken from https://github.com/GuanLab/GuanRank-R
  for (i in 1:nrow(tData)) {
    tA <- km_rank_mat[i,"time"]
    rA <- km_rank_mat[i,"surv_prob"]
    sA <- km_rank_mat[i,"status"]

    if (sA == 1) {
      tBgttA <- km_rank_mat[time > tA, "surv_prob"]
      tBletA_sBeq0 <- km_rank_mat[time <= tA & status == 0, "surv_prob"]
      tBeqtA_sBeq1 <- km_rank_mat[time == tA & status == 1, "surv_prob"]
      km_rank_mat[i, "guan_rank"] <- ifelse(length(tBgttA) == 0, 0,
                                           1 * length(tBgttA)) +
        ifelse(length(tBletA_sBeq0) == 0, 0, sum(rA/tBletA_sBeq0)) +
        ifelse(length(tBeqtA_sBeq1) == 0, 0, 0.5 * length(tBeqtA_sBeq1))
    }else if (sA == 0) {
      tBgetA_sBeq0 <- km_rank_mat[time >= tA & status == 0, "surv_prob"]
      tBgetA_sBeq1 <- km_rank_mat[time >= tA & status == 1, "surv_prob"]
      tBlttA_sBeq0 <- km_rank_mat[time < tA & status == 0, "surv_prob"]
      km_rank_mat[i, "guan_rank"] <- ifelse(length(tBgetA_sBeq0) == 0, 0,
                                           sum(1 - 0.5*tBgetA_sBeq0/rA)) +
        ifelse(length(tBgetA_sBeq1) == 0, 0, sum(1 - tBgetA_sBeq1/rA)) +
        ifelse(length(tBlttA_sBeq0) == 0, 0, sum(0.5*rA/tBlttA_sBeq0))
    }
  }

  # 0.5 is the correction for self-comparison
  km_rank_mat[,"guan_rank"] <- km_rank_mat[,"guan_rank"] - 0.5
  # normalization to [0,1]
  km_rank_mat[,"guan_rank"] <- km_rank_mat[,"guan_rank"] /
    max(km_rank_mat[,"guan_rank"])

  # mapping back to the order of the input data
  gRank <- dplyr::left_join(surv_data,
                            km_rank_mat[,c(3,5)],
                            by = "id", keep = FALSE)
  # return without dummy id
  gRank[,-3]
}
