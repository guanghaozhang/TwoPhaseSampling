####Final Simulation Step 2: merge data
#March 10, 2021
#Guanghao Zhang
#ghzhang@umich.edu

index_length=50000
###################close
#big np
idexist <- matrix(ncol = 3, nrow = index_length)
for (j in 1:3) {
  for (i in 1:index_length){#1:index_length){
    idexist[i,j] <- file.exists(paste0("intermediate_simdata11/sim_beta_bignp_", i, "_", j, ".Rdata"))
  }
}

beta_bignp <- matrix(ncol = 7, nrow = 3)
RE_bignp <- matrix(ncol = 7, nrow = 3)

for (j in 1:3) {
  bootbetafull <- matrix(ncol = 7, nrow = index_length)
  for (i in 1:index_length){#1:index_length){
    if (idexist[i,j] == T){
      load(paste0("intermediate_simdata18/sim_beta_bignp_", i, "_", j, ".Rdata"))
      bootbetafull <- rbind(bootbetafull, beta)} else {
        bootbetafull <- rbind(bootbetafull, matrix(data=NA,nrow=1,ncol=7))
      }
  }
  mean_beta=apply(bootbetafull, 2, mean, na.rm=TRUE)
  var_beta=apply(bootbetafull, 2, var, na.rm=TRUE)
  beta_bignp[j,] <- mean_beta
  RE_bignp[j,] <- var_beta/var_beta[1]
  }

beta_bignp
RE_bignp

RE_bignp[,3]/RE_bignp[,2]
RE_bignp[,6]/RE_bignp[,2]
RE_bignp[,7]/RE_bignp[,2]

save(beta_bignp, RE_bignp, file = "sim_bignp_close_final.RData")

#scp ghzhang@armis2-xfer.arc-ts.umich.edu:/home/ghzhang/ghzhang/sim_bignp.RData /Users/zhangguanghao/Downloads

#small np
idexist <- matrix(ncol = 3, nrow = index_length)
for (j in 1:3) {
  for (i in 1:index_length){#1:index_length){
    idexist[i,j] <- file.exists(paste0("intermediate_simdata11/sim_beta_smallnp_", i, "_", j, ".Rdata"))
  }
}

beta_smallnp <- matrix(ncol = 7, nrow = 3)
RE_smallnp <- matrix(ncol = 7, nrow = 3)

for (j in 1:3) {
  bootbetafull2 <- matrix(ncol = 7, nrow = index_length)
  for (i in 1:index_length){#1:index_length){
    if (idexist[i,j] == T){
    load(paste0("intermediate_simdata18/sim_beta_smallnp_", i, "_", j, ".Rdata"))
    bootbetafull2 <- rbind(bootbetafull2, beta)} else {
      bootbetafull2 <- rbind(bootbetafull2, matrix(data=NA,nrow=1,ncol=7))
    }
  }
  mean_beta=apply(bootbetafull2, 2, mean, na.rm=TRUE)
  var_beta=apply(bootbetafull2, 2, var, na.rm=TRUE)
  beta_smallnp[j,] <- mean_beta
  RE_smallnp[j,] <- var_beta/var_beta[1]
}

beta_smallnp
RE_smallnp

RE_smallnp[,3]/RE_smallnp[,2]
RE_smallnp[,6]/RE_smallnp[,2]
RE_smallnp[,7]/RE_smallnp[,2]

save(beta_smallnp, RE_smallnp, file = "sim_smallnp_close_final.RData")

#############################far
#big np
idexist <- matrix(ncol = 3, nrow = index_length)
for (j in 1:3) {
  for (i in 1:index_length){#1:index_length){
    idexist[i,j] <- file.exists(paste0("intermediate_simdata22/sim_beta_bignp_", i, "_", j, ".Rdata"))
  }
}

beta_bignp <- matrix(ncol = 7, nrow = 3)
RE_bignp <- matrix(ncol = 7, nrow = 3)

for (j in 1:3) {
  bootbetafull <- matrix(ncol = 7, nrow = index_length)
  for (i in 1:index_length){#1:index_length){
    if (idexist[i,j] == T){
      load(paste0("intermediate_simdata18/sim_beta_bignp_", i, "_", j, ".Rdata"))
      bootbetafull <- rbind(bootbetafull, beta)} else {
        bootbetafull <- rbind(bootbetafull, matrix(data=NA,nrow=1,ncol=7))
      }
  }
  mean_beta=apply(bootbetafull, 2, mean, na.rm=TRUE)
  var_beta=apply(bootbetafull, 2, var, na.rm=TRUE)
  beta_bignp[j,] <- mean_beta
  RE_bignp[j,] <- var_beta/var_beta[1]
}

beta_bignp
RE_bignp

RE_bignp[,3]/RE_bignp[,2]
RE_bignp[,6]/RE_bignp[,2]
RE_bignp[,7]/RE_bignp[,2]

save(beta_bignp, RE_bignp, file = "sim_bignp_far_final.RData")

#scp ghzhang@armis2-xfer.arc-ts.umich.edu:/home/ghzhang/ghzhang/sim_bignp.RData /Users/zhangguanghao/Downloads

#small np
idexist <- matrix(ncol = 3, nrow = index_length)
for (j in 1:3) {
  for (i in 1:index_length){#1:index_length){
    idexist[i,j] <- file.exists(paste0("intermediate_simdata22/sim_beta_smallnp_", i, "_", j, ".Rdata"))
  }
}

beta_smallnp <- matrix(ncol = 7, nrow = 3)
RE_smallnp <- matrix(ncol = 7, nrow = 3)

for (j in 1:3) {
  bootbetafull2 <- matrix(ncol = 7, nrow = index_length)
  for (i in 1:index_length){#1:index_length){
    if (idexist[i,j] == T){
      load(paste0("intermediate_simdata18/sim_beta_smallnp_", i, "_", j, ".Rdata"))
      bootbetafull2 <- rbind(bootbetafull2, beta)} else {
        bootbetafull2 <- rbind(bootbetafull2, matrix(data=NA,nrow=1,ncol=7))
      }
  }
  mean_beta=apply(bootbetafull2, 2, mean, na.rm=TRUE)
  var_beta=apply(bootbetafull2, 2, var, na.rm=TRUE)
  beta_smallnp[j,] <- mean_beta
  RE_smallnp[j,] <- var_beta/var_beta[1]
}

beta_smallnp
RE_smallnp

RE_smallnp[,3]/RE_smallnp[,2]
RE_smallnp[,6]/RE_smallnp[,2]
RE_smallnp[,7]/RE_smallnp[,2]

save(beta_smallnp, RE_smallnp, file = "sim_smallnp_far_final.RData")

#scp ghzhang@armis2-xfer.arc-ts.umich.edu:/home/ghzhang/ghzhang/sim_smallnp.RData /Users/zhangguanghao/Downloads