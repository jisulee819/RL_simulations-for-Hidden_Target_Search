library(ggplot2)

#plot( x = rep(target_x, each=6), y= rep(target_y, 6), xlim=c(0,1024), ylim=c(0,768))
#rect( 0, 0, 1024, 768)

Random_Simulation <- function(n){
  rho <- 0
  ratio <- 1024/nx
  
  pr <- read.csv("/Users/jisu/Documents/Hidden Target/Behavioral Data/raw_data/priorreduced.csv")
  pr <- pr[,-1]
  prior <- matrix(NA, nrow=nx, ncol=ny)
  for (a in 1:nx){
    for (b in 1:ny){
      prior[a,b] <- pr[a,b]
    }
  }

  prior_history <- prior
  target_x <- c(-352.0000, -211.2000, -70.4000, 70.4000, 211.2000, 352.0000)
  target_y <- c(-224.0000, -134.4000, -44.8000, 44.8000, 134.4000, 224.0000)
  target_x <- (target_x + 1024/2)/ratio
  target_y <- (target_y + 768/2)/ratio
  rad <- (ny/2)
  
  target_location <- cbind(rep(target_x, each=6), rep(target_y, 6))
  rew <- NA
  rew_history = NA
  rad_history = NA
  infer_history = matrix(NA,ncol=2,nrow=10)

  for (t in 1:n){

    infer <- c(sample(1:nx, 1), sample(1:ny, 1))
    infer_history[t,] <- infer
    
    dist.from.target <- sqrt((target_location[t,1]-infer[1])^2 + ((target_location[t,1]-infer[2])^2))
    if (dist.from.target < rad) {rew <- 1
    } else {rew <- 0
    }
    
    rew_history[t] <- rew
    rad <- (ny/2) * (0.78^sum(rew_history))
    rad_history[t] <- rad
  }
  
  dataset <- data.frame(infer_history, rew_history, rad_history)
  return(dataset)
}



nx <- 64
ny <- 48

for (i in 1:100000){
  data <- Random_Simulation(10)
  data <- cbind(data, i)
  if (i == 1){
    random.simul.data <- data
  } else {
  random.simul.data <- rbind(random.simul.data, data)
  }
}
random.simul.data$trial <- rep(1:10)
random.simul.data$block <- rep(1:36, each=10)
# Average Reward
sum(random.simul.data$rew_history)/nrow(random.simul.data)
# 0.253614

# Reward rate per trial
rew.rate.trial <- NA
for (i in 1:10){
  rew.rate.trial[i] <- mean(random.simul.data$rew_history[random.simul.data$trial==i])
}
data <- data.frame("Trial" = c(1:10), "rew.rate" = rew.rate.trial)  
plot(c(1:10), rew.rate.trial)

pdf("/Users/jisu/Documents/Hidden Target/Visualizations/1906/Reward Rate per Trial of Random Searcher(n=100,000)", width=10, height = 7)

ggplot(data = data, aes(x = Trial, y = rew.rate.trial)) + 
  scale_y_continuous(name="Reward Rate", limits=c(0, 1)) + 
  scale_x_continuous("Trial",limits=c(1,10),breaks=c(1:10)) +
  geom_point() +
  ggtitle("Reward Rate per Trial of Random Searcher (n=100,000)") + 
  theme(plot.title=element_text(face="bold", size=20, hjust=0.5))

dev.off()

#for i in 
#random.simul.data$rad_history <- i

# Reward zone per trial
rad.trial <- NA
for (i in 1:10){
  rad.trial[i] <- mean(random.simul.data$rad_history[random.simul.data$trial==i])
}
data <- data.frame("Trial" = c(1:10), "rad" = rad.trial)  
plot(c(1:10), rad.trial)

pdf("/Users/jisu/Documents/Hidden Target/Visualizations/1906/Radius per Trial of Random Searcher(n=100,000)", width=10, height = 7)

ggplot(data = data, aes(x = Trial, y = rad.trial)) + 
  scale_y_continuous(name="Reward Rate", limits=c(0, 25)) + 
  scale_x_continuous("Trial",limits=c(1,10),breaks=c(1:10)) +
  geom_point() +
  ggtitle("Radius per Trial of Random Searcher (n=100,000)") + 
  theme(plot.title=element_text(face="bold", size=20, hjust=0.5))

dev.off()


### feedback -> 얻은 정보의 양
### target별로, specific한 target의 위치들에 대해 (averaged by 10 trial)
### ttest: 정보량 > 0
### Correlation: beta series (averaged by target) and information gain estimated by this
### KL-divergence distance








