library(mvtnorm)
library(corpcor)
library(Hmisc)

bayesian_update <- function(prior, inference, sigma, reward){
  
  likelihood <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  posterior <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  
  for(i in 1:nx){
    for(j in 1:ny){
      mu <- matrix(c(i,j), nrow=2)
      likelihood[i,j] <- dmvnorm(t(inference), mu, sigma)}}
      maxl <- max(likelihood)
        if (reward == 0){
          likelihood <- maxl-likelihood}
  for(i in 1:nx){
    for(j in 1:ny){
      posterior[i,j] <- likelihood[i,j] * prior[i,j]}}
  return(posterior)
  
}


inference_bs <- function(prior){ ##
  inference <- which(prior == max(prior), arr.ind=TRUE) #
  inference <- inference[sample(1:nrow(inference), 1),]
  return(inference)
}


reward_bs <- function(inference, radius, treasure_location){ ##
  if (sqrt((treasure_location[1]-inference[1])^2 + (treasure_location[2]-inference[2])^2) < radius) return(1)
  else return(0)
}

bayesian_MAP <- function(n){
  rho <- 0 # ?????? x
  rad = 400/16
  s <- rad # 
  lambda <- 3.75
  sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2) #
  prior <- matrix(data = 1/(nx*ny), nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  prior_history <- prior
  infer_history <- c(NA, NA)  
  treasure_location <<- c(sample(1:nx, 1), sample(1:ny, 1))
  
  
  rew_history = NA
  rad_history <- NA
  
  
  for (t in 1:n){
    
    infer <- inference_bs(prior)
    ###truncate
    if (infer[1] > nx) {infer[1] = nx
    } else if (infer[1] < 0) {infer[1] = 0}
    
    if (infer[2] > ny) {infer[2] = ny
    } else if (infer[2] < 0) {infer[2] = 0}
    
    infer_history <- rbind(infer_history, infer)
    
    
    rew <- reward_bs(infer, rad, treasure_location)
    rad_history[t] <- rad
    
    rew_history[t] <- rew
    rad <- (400 - (sum(rew_history)) * (60))/16
    if (sum(rew_history) > 5) rad = 100/16
    
    posterior <- bayesian_update(prior, infer, sigma, rew)
    
    if (rew == 1 & sum(rew_history) <= 5) 
      {s <- s - lambda}
    else {s <- s} # == s2 + decay
    print(s)
    sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2) #
    
    prior <- posterior
    prior_history <- rbind(prior_history, data.matrix(prior))
    
    
    #print(sprintf("N of the trial: %d", t))
  }
  treasure_found <- rew_history[nrow(rew_history)] 
  # compute distances btw the target and every points
  a<-rep(NA,n)
  for (i in 1:n){
    a[i] <- sqrt((treasure_location[1]-infer_history[i+1,1])^2 + (treasure_location[2]-infer_history[i+1,2])^2)}
  
  # set the outputs    
  dataset <- list("prior history" = prior_history, "treasure location" = treasure_location, "inference" = infer_history[-1,], "distance" = a, "reward history" = rew_history, "radius" = rad_history, "result" = treasure_found)  
  return(dataset)  
  
}

nx <<- 64
ny <<- 48
data <- bayesian_MAP(30)
plottingthepath(data,64,48,30)



i <- 6
prior <- data$`prior history`[(nx*(i-1)+1):(nx*i),]



gg_multiheatmap_10 <- function(matrix, col){
  p1 <- gg_heatmap(matrix, 1, col)
  p2 <- gg_heatmap(matrix, 2, col)
  p3 <- gg_heatmap(matrix, 3, col)
  p4 <- gg_heatmap(matrix, 4, col)
  p5 <- gg_heatmap(matrix, 5, col)
  p6 <- gg_heatmap(matrix, 6, col)
  p7 <- gg_heatmap(matrix, 7, col)
  p8 <- gg_heatmap(matrix, 8, col)
  p9 <- gg_heatmap(matrix, 9, col)
  p10 <- gg_heatmap(matrix, 10, col)
  # <- gg_heatmap(matrix, 11, col)
  #p12 <- gg_heatmap(matrix, 12, col)
  multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, 
            #p11, p12, 
            cols=2)}


n <<- 10
nx <<- 64
gg_heatmap(data$`prior history`, 10, "red")
gg_multiheatmap_10(data$`prior history`, "red")







iteration_bs <- function(n){
  distances <- matrix(NA, nrow=n, ncol=100)
  reward_history <- matrix(NA, nrow=n, ncol=100)
  for (i in 1:n){
    treasure_location <<- c(sample(1:nx, 1), sample(1:ny, 1))
    data <- bayesian_MAP(100)
    distances[i,] <- data$distance
    reward_history[i,] <- data$`reward history`
    print(sprintf("iteration %s", i))
  }
  
  dataset <- list("distance" = distances, "reward history" = reward_history)
  return(dataset)}


iterationdata_1000 <- iteration_bs(100)




iterationplot_with_errbar <- function(data,n){
  b <- matrix(rep(NA, n*2), nrow=n, ncol=2)
  for (i in 1:100){
    b[i,] <- quantile(data$distance[,i], probs = c(0.25,0.75), na.rm = FALSE,names = TRUE, type = 7)}
  errbar(x = 1:100, y = apply(data$distance, 2, mean), yplus=b[,2], yminus=b[,1],xlim = c(1,100), ylim = c(0,60),xlab = "Number of Trials", ylab = "Average Distance from the Treasure")
}


iterationplot_bs <- function(data,n){
  par(mfrow = c(1,2), mar=c(4,4,3,4), oma=c(0.5, 0.5, 2, 0.5))
  
  iterationplot_with_errbar(iterationdata_1000,100)
  
  plot(x = 1, y = mean(data$distance[,1]), pch=19, cex=0.7, xlim = c(1,100), ylim = c(0,1), type ="n", xlab = "Number of Trials", ylab = "Percentage of the Rewarded Trials")
  for (i in 1:n){
    points(x = i, y = sum(data$`reward history`[,i])/n, pch=19, cex=0.1, type="p")
    lines(1:100, apply(data$`reward history`/n, 2, sum))}
  mtext("1000 iteration of Bayesian MAP", outer = TRUE, cex = 1.5)}

iterationplot_bs(iterationdata_1000, 100)


