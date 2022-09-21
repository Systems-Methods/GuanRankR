## simple version of GuanRank (https://www.nature.com/articles/s43588-021-00083-2) 
## the code here is similar to the code found in https://github.com/GuanLab/GuanRank-R but
## it makes use of survival::survfit for computing survival probability for simplicity

# surv_data input should have two columns, the first column is time, the second is an indicater with 0 = censored, 1 = event
# this code is a simplified version of code from https://github.com/GuanLab/GuanRank-R
calc_guanRank <- function(surv_data){
  surv_data           <- data.frame(time = surv_data[,1],status = surv_data[,2], id = 1:nrow(surv_data))
  tData               <- na.omit(surv_data); 
  tData               <- tData[order(tData[,"time"]),] # reorder to facilitate calculation below

  #compute survival probabilities by time
  km_fit              <- survival::survfit(survival::Surv(tData$time, tData$status)~1)
  tKM                 <- data.frame(time = km_fit$time, surv_prob = km_fit$surv)
  km_curve            <- dplyr::left_join(tData, tKM, "time")
  km_rank_mat         <- cbind(km_curve,"guan_rank"=rep(0,nrow(tData)))
  time                <- km_rank_mat[,"time"]
  status              <- km_rank_mat[,"status"]
  
  # note: the section below does not break down cleanly into cases 1 - 6 from https://www.nature.com/articles/s43588-021-00083-2
  # instead it is a clever computation of rank using the same conditionals but in a different order to ease the rank calculation.
  # it taken from https://github.com/GuanLab/GuanRank-R
  for(i in 1:nrow(tData)){
    tA <- km_rank_mat[i,"time"]
    rA <- km_rank_mat[i,"surv_prob"]
    sA <- km_rank_mat[i,"status"]
    
    if(sA==1){
      tBgttA <- km_rank_mat[time > tA,"surv_prob"]                   
      tBletA_sBeq0 <- km_rank_mat[time <= tA & status==0,"surv_prob"]  
      tBeqtA_sBeq1 <- km_rank_mat[time == tA & status==1,"surv_prob"]  
      km_rank_mat[i,"guan_rank"] <- ifelse(length(tBgttA) == 0, 0, 1 * length(tBgttA)) +
        ifelse(length(tBletA_sBeq0) == 0, 0, sum(rA/tBletA_sBeq0)) +
        ifelse(length(tBeqtA_sBeq1) == 0, 0, 0.5 * length(tBeqtA_sBeq1))
    }
    
    if(sA==0){
      tBgetA_sBeq0 <- km_rank_mat[time >= tA & status == 0,"surv_prob"]
      tBgetA_sBeq1 <- km_rank_mat[time >= tA & status == 1,"surv_prob"]
      tBlttA_sBeq0 <- km_rank_mat[time < tA & status == 0,"surv_prob"]
      km_rank_mat[i,"guan_rank"] <- ifelse(length(tBgetA_sBeq0) == 0, 0, sum(1 - 0.5*tBgetA_sBeq0/rA)) +
        ifelse(length(tBgetA_sBeq1) == 0, 0, sum(1 - tBgetA_sBeq1/rA)) +
        ifelse(length(tBlttA_sBeq0) == 0, 0, sum(0.5*rA/tBlttA_sBeq0))
    }
  }
  
  km_rank_mat[,"guan_rank"] <- km_rank_mat[,"guan_rank"]-0.5 # 0.5 is the correction for self-comparison
  km_rank_mat[,"guan_rank"] <- km_rank_mat[,"guan_rank"]/max(km_rank_mat[,"guan_rank"]) # normalization to [0,1]
  
  gRank                <- dplyr::left_join(surv_data, km_rank_mat[,c(3,5)], by="id", keep=F)
  return(gRank[,-3]) # return without dummy id
}

