

#########################
#### greedy inference ###
#########################
inference_greedy <- function(val_fn){
  a = which(val_fn == max(val_fn), arr.ind=TRUE)  ### GREEDY
  infer <- a[sample(nrow(a),1),]
  return(infer)
}


##########################
#### softmax inference ###
##########################
inference <- function(beta, val_fn){
  softmax_prob <- matrix(data = 0, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  for(i in 1:nx){
    for(j in 1:ny){
      softmax_prob[i,j] <- (exp(beta * val_fn[i,j])) / (sum(exp(beta * val_fn)))
    }}
  
  softmax_matrix <- matrix(NA, nrow=nx*ny, ncol=3)
  softmax_matrix[,1] <- rep(1:nx, times = 1, length.out = NA, each = ny)
  softmax_matrix[,2] <- rep(1:ny, times = nx, length.out = NA, each = 1)
  softmax_matrix[,3] <- as.vector(t(softmax_prob))
  a <- sample(x = c(1:(nx*ny)), 1, replace = T, prob = as.vector(t(softmax_prob))) 
  b <- c(softmax_matrix[a,1],softmax_matrix[a,2])
  return(b)
}


####################
#### eligibility ###
####################
eligibility <- function(mu, infer_history){
  current_point <- infer_history[nrow(infer_history)-1,]
  eligibility <- matrix(data = 0, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  for(i in 1:nx){
    for(j in 1:ny){
      eligibility[i,j] <- (sqrt((i-current_point[1])^2 + (j-current_point[2])^2))
    }}
  eligibility <- exp(eligibility/(-mu)) 
  return(eligibility)}

eligibility2 <- function(mu, infer_history){
  next_point <- infer_history[nrow(infer_history),]
  eligibility <- matrix(data = 0, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  for(i in 1:nx){
    for(j in 1:ny){
      eligibility[i,j] <- (sqrt((i-next_point[1])^2 + (j-next_point[2])^2))
    }}
  eligibility <- exp(eligibility/(-mu)) 
  return(eligibility)}


######################################################################
#### reward ### whether the inferred point is included in reward zone
######################################################################
reward_td <- function(infer_history, radius){
  current_point <- infer_history[nrow(infer_history)-1,]
  if (sqrt((treasure_location[1]-current_point[1])^2 + (treasure_location[2]-current_point[2])^2) < radius) return(1)
  else return(0)
}


########################
#### value update TD ###
########################
value_update <- function(alpha, gamma, val_fn, elig, reward, infer_history){ ##
  for(i in 1:nx){
    for(j in 1:ny){
      next_point <- infer_history[nrow(infer_history),]
      current_point <- infer_history[nrow(infer_history)-1,]
      #print(dim(val_fn))
      #print(next_point)
      val_fn[i,j] <- alpha * (reward + gamma*val_fn[next_point[1],next_point[2]] - val_fn[current_point[1],current_point[2]]) * elig[i,j] + val_fn[i,j]
      
      #if (val_fn[i,j] < 0) val_fn[i,j] <- 0
      #else if(val_fn[i,j] > 1) val_fn[i,j] <- 1 ##truncate?
    }}
  return(val_fn)}

###########################
#### eligibility update ###
###########################
eligibility_update <- function(gamma, mu, lambda, infer_history, elig){
  for(i in 1:nx){
    for(j in 1:ny){
      elig[i,j] <- (gamma*lambda*elig[i,j] + eligibility2(mu, infer_history))[i,j] 
    }}
  return(elig)}


############################
#### temporal difference ###
############################
temp_diff <- function(n){
  #### hyperparameters
  mu <- 32 
  alpha <- 0.7
  gamma <- 0.9
  beta <- 3
  lambda <- 0.5
  
  #### settings and prerequisites
  rad = 32
  rew_history = 0
  rew_history2 <- NA
  rad_history <- NA
  treasure_found <- 0
  val_fn <- matrix(data = 0, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  initial_point <- c(sample(1:nx, 1), sample(1:ny, 1)) # centroid first
  infer_history <- initial_point ###############
  elig_matrix <- rep(NA, ny)
  val_matrix <- rep(NA, ny)
  
  #### run the learning
  for (t in 1:n){
    rad_history <- rbind(rad_history, rad)
    infer <- inference(beta, val_fn) # infer new points (next point)
    infer_history <- rbind(infer_history, infer)
    
    if (t==1) elig <- eligibility(mu, infer_history) ## ????????? ??????????????? ????????? ??? # compute eligibility based on current point
    
    rew <- reward_td(infer_history, rad) # compute reward of current point
    if (rew == 1) rew_history <- rew_history + 1
    rew_history2 <- rbind(rew_history2, rew)
              rad <- 32 - (rew_history) * (6.24)
              if (rew_history > 4) rad = 7.04 
              mu <- rad# radius shirinks as reward increases
    
    val_fn <- value_update(alpha, gamma, val_fn, elig, rew, infer_history)
    elig_matrix <- rbind(elig_matrix,data.matrix(elig))
    val_matrix <- rbind(val_matrix, data.matrix(val_fn))
    
    elig <- eligibility_update(gamma, mu, lambda, infer_history, elig)
    print(sprintf("N of the trial: %d", t))} 
  # value fn is updated based on current point, current value of current and next point, eligibility, reward.
  
  #### after the learning
  
  elig_matrix <- elig_matrix[-1,]
  val_matrix <- val_matrix[-1,]
  
  
  # whether he FINALLY has found the treasure
  treasure_found <- rew_history2[nrow(rew_history2)] 
  
  # let the dataset seem better
  #rownames(infer_history) = c("initial point","infer_1","infer_2","infer_3","infer_4","infer_5","infer_6","infer_7","infer_8","infer_9","infer_10")
  
  # compute distances btw the target and every points
  a<-rep(NA,n)
  for (i in 1:n){
    a[i] <- sqrt((treasure_location[1]-infer_history[i,1])^2 + (treasure_location[2]-infer_history[i,2])^2)}
  
  # set the outputs    
  dataset <- list("treasure location" = treasure_location, "inference" = infer_history[-1,], "distance" = a, "reward history" = rew_history2[-1], "radius" = rad_history[-1], "result" = treasure_found, "value function" = val_matrix, "eligibility" = elig_matrix)  
  return(dataset)  
}

#############
#### run ####
#############
nx <<- 64
ny <<- 48
treasure_location <<- c(sample(1:nx, 1), sample(1:ny, 1))
data <- temp_diff(100)







####### ITERATION #########

iteration <- function(n){
  prob_of_win <- 0
  distances <- matrix(NA, nrow=n, ncol=100)
  reward_history <- matrix(NA, nrow=n, ncol=100)
  for (i in 1:n){
    treasure_location <<- c(sample(1:nx, 1), sample(1:ny, 1))
    data <- temp_diff(100)
    prob_of_win <- prob_of_win + data$result
    
    distances[i,] <- data$distance
    reward_history[i,] <- data$`reward history`
    print(sprintf("iteration %s", i))
  }
  prob_of_win <- prob_of_win / n
  
  dataset <- list("probability to win" = prob_of_win, "distance" = distances, "reward history" = reward_history)
  return(dataset)}


iterationdata_1000 <- iteration(10)


iterationplot <- function(n){
  par(mfrow = c(1,1), mar=c(4,4,3,4), oma=c(0.5, 0.5, 2, 0.5))
  plot(x = 1, y = mean(iterationdata$distance[,1]), pch=19, cex=0.7, xlim = c(1,100), ylim = c(0,80), type ="n", xlab = "Number of Trials", ylab = "AVERAGE Distance from the treasure")
  for (i in 1:n){
    points(x = i, y = mean(iterationdata$distance[,i]), pch=19, cex=0.3, type="p")
    
    #b <- quantile(iterationdata$distance[,i], probs = c(0.25,0.75), na.rm = FALSE,
    #              names = TRUE, type = 7)
    #errbar(x = i, y = iterationdata$distance[,i], yplus=b[2], yminusb[1])
    }}

iterationplot(100)
#
#data_compile <- function(data){
#  val_fn2 <- t(data$`value function`)
#  val_fn2 <- data.matrix(val_fn2)
#  #data_heatmap <- heatmap(val_fn, Rowv=NA, Colv=NA, col=cm.colors(256), scale="column")
#
#  data_matrix <- matrix(NA, nrow=nx*ny, ncol=3)
#  data_matrix[,1] <- rep(1:nx, times = 1, length.out = NA, each = ny)
#  data_matrix[,2] <- rep(1:ny, times = nx, length.out = NA, each = 1)
#  data_matrix[,3] <- as.vector(t(val_fn2))
#  
#}
#data_compile(data)
