library(mvtnorm)
library(corpcor)
library(Hmisc)


inference_pg <- function(mu, sigma){
  inference <- t(rmvnorm(1, mu, sigma))
  return(inference)
}


reward_pg <- function(inference, radius){
  if (sqrt((treasure_location[1]-inference[1])^2 + (treasure_location[2]-inference[2])^2) < radius) return(1)
  else return(0)
}

policy_update <- function(inference, mu, sigma, reward, alpha_m, beta_m, alpha_s, beta_s){
  s <- solve(sigma)
  mu <- alpha_m * (reward - beta_m) * s %*% (inference - mu) + mu
  d_sig <- (-0.5 * alpha_s *(reward - beta_s) * (s - s%*%(inference-mu)%*%t(inference-mu)%*%s))
  sigma <- d_sig + sigma
  return(list("mu" = mu, "sigma" = sigma))
}


policy_gradient <- function(n){
  rho <- 0
  mu1 <- 0.5*nx #32
  s1 <- 10
  mu2 <- 0.5*ny #24
  s2 <- 10
  # Parameters for bivariate normal distribution
  mu <- matrix(c(mu1,mu2), nrow=2) # Mean 
  sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),nrow=2) # Covariance matrix
  
  alpha_m = 15
  beta_m = 0.1#### beta??? ????????? ?????? ??????: reward zone??? shirinking?????? ?????????
  alpha_s = 15  #### reward??? ???????????? ???????????? ?????????????????? ???????????? ?????? ???????????? ??????
  beta_s = 0.05
  
  infer_history <- c(NA, NA)  ########
  mu_history <- t(mu)
  sigma_history <- sigma
  rad = 32
  rew_history = 0
  rew_history2 <- NA
  rad_history <- NA
  
  
    for (t in 1:n){
      
      
      infer <- inference_pg(mu,sigma)
      
      ###truncate
      if (infer[1] > nx) {infer[1] = nx
      } else if (infer[1] < 0) {infer[1] = 0}
      
      if (infer[2] > ny) {infer[2] = ny
      } else if (infer[2] < 0) {infer[2] = 0}
      
      rew <- reward_pg(infer, rad)
      rad_history <- rbind(rad_history, rad)
      
      if (rew == 1) rew_history <- rew_history + 1
      rew_history2 <- rbind(rew_history2, rew)
      rad <- 32 - (rew_history) * (6.24)
      if (rew_history > 4) rad = 7.04 
      
      param_pg <- policy_update(infer, mu, sigma, rew, alpha_m, beta_m, alpha_s, beta_s)
      mu <- param_pg$mu
      sigma <- param_pg$sigma
      
      
      ##############################
      ##############################
      ##############################
      sigma <- make.positive.definite(sigma)
      sigma[1,2] <- sigma[2,1] #####
      ##############################
      ##############################
      ##############################
      
      mu_history <- rbind(mu_history, t(mu))
      sigma_history <- rbind(sigma_history, sigma)
      infer_history <- rbind(infer_history, t(infer))
      #print(sprintf("N of the trial: %d", t))
    }
  treasure_found <- rew_history2[nrow(rew_history2)] 
  # compute distances btw the target and every points
  a<-rep(NA,n)
  for (i in 1:n){
    a[i] <- sqrt((treasure_location[1]-infer_history[i+1,1])^2 + (treasure_location[2]-infer_history[i+1,2])^2)}
  
  # set the outputs    
  dataset <- list("treasure location" = treasure_location, "inference" = infer_history[-1,], "distance" = a, "reward history" = rew_history2[-1], "radius" = rad_history[-1], "result" = treasure_found, "mu" = mu_history, "sigma" = sigma_history)  
  return(dataset)  
  
}


nx <<- 64
ny <<- 48
treasure_location <<- c(sample(1:nx, 1), sample(1:ny, 1))
data <- policy_gradient(100)
plottingthepath(data,64,48,100)







iteration_pg <- function(n){
  distances <- matrix(NA, nrow=n, ncol=100)
  reward_history <- matrix(NA, nrow=n, ncol=100)
  for (i in 1:n){
    treasure_location <<- c(sample(1:nx, 1), sample(1:ny, 1))
    data <- policy_gradient(100)
    distances[i,] <- data$distance
    reward_history[i,] <- data$`reward history`
    print(sprintf("iteration %s", i))
  }
  
  dataset <- list("distance" = distances, "reward history" = reward_history)
  return(dataset)}


iterationdata_1000 <- iteration_pg(1000)




iterationplot_with_errbar <- function(data,n){
  b <- matrix(rep(NA, n*2), nrow=n, ncol=2)
  for (i in 1:100){
    b[i,] <- quantile(data$distance[,i], probs = c(0.25,0.75), na.rm = FALSE,names = TRUE, type = 7)}
  errbar(x = 1:100, y = apply(data$distance, 2, mean), yplus=b[,2], yminus=b[,1],xlim = c(1,100), ylim = c(0,60),xlab = "Number of Trials", ylab = "Average Distance from the Treasure")
}


iterationplot_pg <- function(data,n){
    par(mfrow = c(1,2), mar=c(4,4,3,4), oma=c(0.5, 0.5, 2, 0.5))
    
    iterationplot_with_errbar(iterationdata_1000,100)

    plot(x = 1, y = mean(data$distance[,1]), pch=19, cex=0.7, xlim = c(1,100), ylim = c(0,1000), type ="n", xlab = "Number of Trials", ylab = "Number of the Rewarded Trials")
    for (i in 1:n){
      points(x = i, y = sum(data$`reward history`[,i]), pch=19, cex=0.1, type="p")
      lines(1:100, apply(data$`reward history`, 2, sum))}
    mtext("1000 iteration of Policy Gradient", outer = TRUE, cex = 1.5)}

iterationplot_pg(iterationdata_1000, 100)



