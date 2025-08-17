######################################################################
# Propagating uncertainties from species identification to occupancy #
#                           with MONTE CARLO                         #
######################################################################

# CASE STUDY 2 : Occupancy modelling with pictures from camtraps 
       

library(purrr) ; library(dplyr)
library(ggplot2) ; library(ggExtra) ; library(RColorBrewer)
library(unmarked)
# "C:/Users/monchy/Desktop/Projets/Uncertainties_Occupancy/MonteCarlo"

# Import data  ----
pics <- read.csv(file="./cleaned_data2.csv",
                row.names=1)
rownames(pics) <- NULL


# Prepare data for unmarked package  ---- 
n_point <- length(unique(pics$point))

# Aggregate pictures by month to create visit samplings
occas2017 <- paste("2017", sort(unique(pics$month[pics$year==2017]),decreasing=F), sep=".")
occas2018 <- paste("2018", sort(unique(pics$month[pics$year==2018]),decreasing=F), sep=".")
occas2019 <- paste("2019", sort(unique(pics$month[pics$year==2019]),decreasing=F), sep=".")

n_occas <- sum(n_distinct(pics$month[pics$year==2017]), 
               n_distinct(pics$month[pics$year==2018]),
               n_distinct(pics$month[pics$year==2019]))

# get the sites
point <- unique(pics$point)

# Build the detection matrix from thresholding the classification scores
detect_simulate <- function(point, W, thresh){
  point %>%
    map_dfr(function(s){
      mat <- rep(NA, n_occas)
      for(o in 1:n_occas){
        if( nrow(W[W$point==s & W$occas==o,4])!=0 ){
          if(any(W[W$point==s & W$occas==o,4] > thresh)){
            mat[o]=1
          } else{ mat[o]=0 }
        }
      }
      return(c(mat=mat))
    })
}


# LYNX Occupancy modelling  ----
## Preparing identification matrix... ---- 
### ... from ground truth data ----
W1 <- pics[,c(1,3,4,7)] %>%
  group_by(point,year,month)

W1 <- W1 %>% 
  mutate(lynx=ifelse(espece=="lynx",1,0), .before = espece)

# Replace year and month by the sampling occasion
W1$occas <- map2_dbl(W1$year, W1$month, function(x,y){
  if(x==2017){
    occas = y
  } 
  if(x==2018){
    occas = y+12
  }
  if(x==2019){
    occas = y+24
  }
  return(occas=occas)
})

thresh <- 0.5
Y_manual <- as.data.frame(detect_simulate(point, W1, thresh))
colnames(Y_manual) <- c(occas2017, occas2018, occas2019)
rownames(Y_manual) <- unique(pics$point)
Y_manual


### ... from DL predictions ----
W0 <- pics[,c(1,3,4,19)] %>%
  group_by(point,year,month)

# Replace year and month by the sampling occasion
W0$occas <- map2_dbl(W0$year, W0$month, function(x,y){
  if(x==2017){
    occas = y
  } 
  if(x==2018){
    occas = y+12
  }
  if(x==2019){
    occas = y+24
  }
  return(occas=occas)
})


# score threshold <- 0.5
Y_DL_50 <- as.data.frame(detect_simulate(point, W0, 0.5))
colnames(Y_DL_50) <- c(occas2017, occas2018, occas2019)
rownames(Y_DL_50) <- unique(pics$point)
Y_DL_50

# score threshold <- 0.95
Y_DL_95 <- as.data.frame(detect_simulate(point, W0, 0.95))
colnames(Y_DL_95) <- c(occas2017, occas2018, occas2019)
rownames(Y_DL_95) <- unique(pics$point)
Y_DL_95

## Occupancy estimates with unmarked package ----

### Estimates from manually labelled data ---
# Get data for 2017 into unmarked format
frame_lynx <- unmarked::unmarkedFrameOccu(y=Y_manual[,1:12]) # FOR 2017
summary(frame_lynx)
# fit
mod_lynx <- unmarked::occu(~1~1, frame_lynx, linkPsi = "logit")
# estimates
summary(mod_lynx)
unmarked::predict(mod_lynx, 'state')[1,]
unmarked::predict(mod_lynx, 'det')[1,] ; backTransform(mod_lynx, "det")
# save estimates and confidence intervals (90%)
unmrkd_lynx_manu <- mod_lynx@opt$par
unmrkd_lynx_manu <- plogis(c(unmrkd_lynx_manu, 
                             confint(mod_lynx,type="state", level=0.95, method="profile"),
                             confint(mod_lynx,type="det", level=0.95, method="profile")))

### Estimates from DL output with threshold at 0.50
frame_lynx <- unmarked::unmarkedFrameOccu(y=Y_DL_50[,1:12]) # FOR 2017
summary(frame_lynx)
# fit
mod_lynx <- unmarked::occu(~1~1, frame_lynx, linkPsi = "logit")
# estimates
summary(mod_lynx)
unmarked::predict(mod_lynx, 'state')[1,] ; unmarked::predict(mod_lynx, 'det')[1,]
# save estimates and confidence intervals (90%)
unmrkd_lynx_50 <- mod_lynx@opt$par
unmrkd_lynx_50 <- plogis(c(unmrkd_lynx_50, 
                           confint(mod_lynx,type="state", level=0.95, method="profile"),
                           confint(mod_lynx,type="det", level=0.95, method="profile")))

### Estimates from DL output with threshold at 0.95
frame_lynx <- unmarked::unmarkedFrameOccu(y=Y_DL_95[,1:12]) # FOR 2017
summary(frame_lynx)
# fit
mod_lynx <- unmarked::occu(~1~1, frame_lynx, linkPsi = "logit")
# estimates
summary(mod_lynx)
unmarked::predict(mod_lynx, 'state')[1,] ; unmarked::predict(mod_lynx, 'det')[1,]
# save estimates and confidence intervals (90%)
unmrkd_lynx_95 <- mod_lynx@opt$par
unmrkd_lynx_95 <- plogis(c(unmrkd_lynx_95, 
                           confint(mod_lynx,type="state", level=0.95, method="profile"),
                           confint(mod_lynx,type="det", level=0.95, method="profile")))


## MONTE CARLO simulations (x100)  ----
k <- 100 # number of simulations

# using scores as probability
recap <- data.frame()
Y_map <- matrix(0,nrow=n_point, ncol=12)
ptm <- proc.time()[3]

set.seed(123)
for(i in 1:k){
  # Simulate the identification matrix from scores
  Wbis <- subset(W0, year==2017) # W0
  Wbis$lynx <- rbinom(nrow(Wbis), 1, W0$lynx[W0$year==2017])
  # Simulate the detection matrix
  Ybis <- as.data.frame(detect_simulate(point, Wbis, 0.5))
  
  # save the details of detection matrix
  # Y_map <- Ybis
  Y_map <- as.matrix(Y_map) + as.matrix(Ybis[,1:12])
  
  # put data in unmarked format
  frame <- unmarked::unmarkedFrameOccu(y=Ybis[,1:12])
  # fit
  mod_lynx <- unmarked::occu(~1~1, frame, linkPsi = "logit")
  
  # Save estimates in CI
  recap[i,1] <- plogis(mod_lynx@opt$par[1])
  recap[i,2:3] <- plogis(confint(mod_lynx, type="state", level=0.95, method="profile"))
  recap[i,4] <- plogis(mod_lynx@opt$par[2])
  recap[i,5:6] <- plogis(confint(mod_lynx, type="det", level=0.95, method="profile"))
  
  # Save the number of occupied sites (at least one identification)
  recap[i,7] <-  sum(unlist(lapply(apply(Ybis, 1,sum, na.rm=T),
                                   function(x){ifelse(x==0, return(0), return(1))})))
}
print(proc.time()[3]-ptm)
colnames(recap) <- c("psi","psi_inf", "psi_sup","p", "p_inf", "p_sup", "N_occ")
beepr::beep(1)

write.csv(recap, file="./recap_lynx_100.csv")

recap_lynx <- recap
# recap_lynx <- read.csv(file="./recap_lynx_100_095.csv", row.names=1)


## Fig4 - estimates LYNX ----
### occupancy ----
recap_lynx$psi_inf[recap_lynx$psi_inf==0] <- NA
psi_lynx_df <- data.frame(method=c(c("1_95","1_50"), rep("2",100)), sample=c(rep(1,2),1:100),
                          estimates=c(unmrkd_lynx_95[1], unmrkd_lynx_50[1], recap_lynx$psi),
                          CI_inf=c(unmrkd_lynx_95[3], unmrkd_lynx_50[3], recap_lynx$psi_inf),
                          CI_sup=c(unmrkd_lynx_95[4], unmrkd_lynx_50[4], recap_lynx$psi_sup))
psi_lynx_df <- psi_lynx_df %>%
  mutate(method=factor(method, levels=c("2","1_50", "1_95")))

fig4_psi_lynx <- ggplot(data=psi_lynx_df, aes(x=method,y=estimates, group=method, col=method)) +
  geom_linerange(aes(ymin=CI_inf, ymax=CI_sup),position=position_dodge2(width=1), linewidth=0.8, color="gray60") +
  geom_point(position=position_dodge2(width=1), size=2.5) + 
  scale_color_manual(values=c("#F26419","#86BBD8","#2F4858"), 
                     labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) +
  scale_x_discrete(labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) + 
  theme_light() + ylim(0.5,1) +
  labs(x=NULL,y=NULL) + # , title="Lynx"
  theme(panel.grid.major.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(color = "gray20",size = 10))

### detection ----
recap_lynx$p_inf[recap_lynx$p_inf==0] <- NA
p_lynx_df <- data.frame(method=c(c("1_95","1_50"), rep("2",100)), sample=c(rep(1,2),1:100),
                        estimates=c(unmrkd_lynx_95[2], unmrkd_lynx_50[2], recap_lynx$p),
                        CI_inf=c(unmrkd_lynx_95[5], unmrkd_lynx_50[5], recap_lynx$p_inf),
                        CI_sup=c(unmrkd_lynx_95[6], unmrkd_lynx_50[6], recap_lynx$p_sup))
p_lynx_df <- p_lynx_df %>%
  mutate(method=factor(method, levels=c("2","1_50", "1_95")))

fig4_p_lynx <- ggplot(data=p_lynx_df, aes(x=method,y=estimates, group=method, col=method)) +
  geom_linerange(aes(ymin=CI_inf, ymax=CI_sup),position=position_dodge2(width=1), linewidth=0.8, color="gray60") +
  geom_point(position=position_dodge2(width=1), size=2.5) + # shape=21,  stroke=2,  fill="gray60"
  scale_color_manual(values=c("#F26419","#86BBD8","#2F4858"), labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) +
  scale_x_discrete(labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) + 
  theme_light() + ylim(0,1) +
  labs(x=NULL,y=NULL) +
  theme(panel.grid.major.x = element_blank(),legend.position = "none",
        axis.text.x = element_text(color = "gray20",size = 10))



# LAGOMORPHE Occupancy  ----
## Preparing identification matrix ... ----
### ... from ground truth data ----
W1 <- pics[,c(1,3,4,7)] %>%
  group_by(point,year,month)

W1 <- W1 %>% 
  mutate(lagomorph=ifelse(espece=="lagomorph",1,0), .before = espece)

# Replace year and month by the sampling occasion
W1$occas <- map2_dbl(W1$year, W1$month, function(x,y){
  if(x==2017){
    occas = y
  } 
  if(x==2018){
    occas = y+12
  }
  if(x==2019){
    occas = y+24
  }
  return(occas=occas)
})

thresh <- 0.5
Y_manual <- as.data.frame(detect_simulate(point, W1, thresh))
colnames(Y_manual) <- c(occas2017, occas2018, occas2019)
rownames(Y_manual) <- unique(pics$point)
Y_manual


### ...from DL predictions ----
W0 <- pics[,c(1,3,4,17)] %>%
  group_by(point,year,month)

# Replace year and month by the sampling occasion
W0$occas <- map2_dbl(W0$year, W0$month, function(x,y){
  if(x==2017){
    occas = y
  } 
  if(x==2018){
    occas = y+12
  }
  if(x==2019){
    occas = y+24
  }
  return(occas=occas)
})


# score threshold <- 0.5
Y_DL_50 <- as.data.frame(detect_simulate(point, W0, 0.5))
colnames(Y_DL_50) <- c(occas2017, occas2018, occas2019)
rownames(Y_DL_50) <- unique(pics$point)
Y_DL_50

# score threshold <- 0.95
Y_DL_95 <- as.data.frame(detect_simulate(point, W0, 0.95))
colnames(Y_DL_95) <- c(occas2017, occas2018, occas2019)
rownames(Y_DL_95) <- unique(pics$point)
Y_DL_95

## Occupancy estimates with unmarked package ----
### Estimates from manually labelled data
# Get data for 2017 into unmarked format
frame_lagomorph <- unmarked::unmarkedFrameOccu(y=Y_manual[,1:12]) # FOR 2017
summary(frame_lagomorph)
# fit
mod_lagomorph <- unmarked::occu(~1~1, frame_lagomorph, linkPsi = "logit")
# estimates
summary(mod_lagomorph)
unmarked::predict(mod_lagomorph, 'state')[1,]
unmarked::predict(mod_lagomorph, 'det')[1,] ; backTransform(mod_lagomorph, "det")
# save estimates and confidence intervals (90%)
unmrkd_lagomorph_manu <- mod_lagomorph@opt$par
unmrkd_lagomorph_manu <- plogis(c(unmrkd_lagomorph_manu, 
                                  confint(mod_lagomorph,type="state", level=0.95, method="profile"),
                                  confint(mod_lagomorph,type="det", level=0.95, method="profile")))

### Estimates from DL output with threshold at 0.50
frame_lagomorph <- unmarked::unmarkedFrameOccu(y=Y_DL_50[,1:12]) # FOR 2017
summary(frame_lagomorph)
# fit
mod_lagomorph <- unmarked::occu(~1~1, frame_lagomorph, linkPsi = "logit")
# estimates
summary(mod_lagomorph)
unmarked::predict(mod_lagomorph, 'state')[1,] ; unmarked::predict(mod_lagomorph, 'det')[1,]
# save estimates and confidence intervals (90%)
unmrkd_lagomorph_50 <- mod_lagomorph@opt$par
unmrkd_lagomorph_50 <- plogis(c(unmrkd_lagomorph_50, 
                                confint(mod_lagomorph,type="state", level=0.95, method="profile"),
                                confint(mod_lagomorph,type="det", level=0.95, method="profile")))

### Estimates from DL output with threshold at 0.95
frame_lagomorph <- unmarked::unmarkedFrameOccu(y=Y_DL_95[,1:12]) # FOR 2017
summary(frame_lagomorph)
# fit
mod_lagomorph <- unmarked::occu(~1~1, frame_lagomorph, linkPsi = "logit")
# estimates
summary(mod_lagomorph)
unmarked::predict(mod_lagomorph, 'state')[1,] ; unmarked::predict(mod_lagomorph, 'det')[1,]
# save estimates and confidence intervals (90%)
unmrkd_lagomorph_95 <- mod_lagomorph@opt$par
unmrkd_lagomorph_95 <- plogis(c(unmrkd_lagomorph_95, 
                                confint(mod_lagomorph,type="state", level=0.95, method="profile"),
                                confint(mod_lagomorph,type="det", level=0.95, method="profile")))


## MONTE CARLO simulations (x100)  ----
k <- 100 # number of simulations

# the probability is the score
recap <- data.frame()
Y_map <- matrix(0,nrow=n_point, ncol=12)
ptm <- proc.time()[3]

set.seed(123)
for(i in 1:k){
  # Simulate the identification matrix from scores
  Wbis <- subset(W0, year==2017) # W0
  Wbis$lagomorphe <- rbinom(nrow(Wbis), 1, W0$lagomorphe[W0$year==2017])
  # Simulate the detection matrix
  Ybis <- as.data.frame(detect_simulate(point, Wbis, 0.5))
  
  # save the details of detection matrix
  # Y_map <- Ybis
  Y_map <- as.matrix(Y_map) + as.matrix(Ybis[,1:12])
  
  # put data in unmarked format
  frame <- unmarked::unmarkedFrameOccu(y=Ybis[,1:12])
  # fit
  mod_lagomorph <- unmarked::occu(~1~1, frame, linkPsi = "logit")
  
  # Save estimates in CI
  recap[i,1] <- plogis(mod_lagomorph@opt$par[1])
  recap[i,2:3] <- plogis(confint(mod_lagomorph, type="state", level=0.95, method="profile"))
  recap[i,4] <- plogis(mod_lagomorph@opt$par[2])
  recap[i,5:6] <- plogis(confint(mod_lagomorph, type="det", level=0.95, method="profile"))
  
  # Save the number of occupied sites (at least one identification)
  recap[i,7] <-  sum(unlist(lapply(apply(Ybis, 1,sum, na.rm=T),
                                   function(x){ifelse(x==0, return(0), return(1))})))
}
print(proc.time()[3]-ptm)
colnames(recap) <- c("psi","psi_inf", "psi_sup","p", "p_inf", "p_sup", "N_occ")
beepr::beep(1)

write.csv(recap, file="./recap_lagomorphe_100.csv")

recap_lago <- recap
# recap_lago <- read.csv(file="./recap_lagomorphe_100_095.csv", row.names=1)


## Fig4 - estimates LAGOMORPHE ----
### occupancy ----
# recap$psi_inf[recap$psi_inf==0] <- NA
psi_lago_df <- data.frame(method=c(c("1_95","1_50"), rep("2",100)), sample=c(rep(1,2),1:100),
                          estimates=c(unmrkd_lagomorph_95[1], unmrkd_lagomorph_50[1], recap_lago$psi),
                          CI_inf=c(unmrkd_lagomorph_95[3], unmrkd_lagomorph_50[3], recap_lago$psi_inf),
                          CI_sup=c(unmrkd_lagomorph_95[4], unmrkd_lagomorph_50[4], recap_lago$psi_sup))
psi_lago_df <- psi_lago_df %>%
  mutate(method=factor(method, levels=c("2","1_50", "1_95")))
# fct_relevel

fig4_psi_lago <- ggplot(data=psi_lago_df, aes(x=method,y=estimates, group=method, col=method)) +
  geom_linerange(aes(ymin=CI_inf, ymax=CI_sup),position=position_dodge2(width=1), linewidth=0.8, color="gray60") +
  geom_point(position=position_dodge2(width=1), size=2.5) + # shape=21,  stroke=2,  fill="gray60"
  scale_color_manual(values=c("#F26419","#86BBD8","#2F4858"), labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) +
  # scale_y_continuous(breaks=seq(1.0150,1.0250, 0.0010)) +
  scale_x_discrete(labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) + # expand=c(0.1,0.1),
  theme_light() + ylim(0.5,1) +
  labs(x=NULL,y=NULL) + #, title="Lagomorph"
  theme(panel.grid.major.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(color = "gray20",size = 10))

### detection ----
# recap$p_inf[recap$p_inf==0] <- NA
p_lago_df <- data.frame(method=c(c("1_95","1_50"), rep("2",100)), sample=c(rep(1,2),1:100),
                        # data.frame(method=c(c("1_50","1_95"), rep(2,100)), sample=c(rep(1,2),1:100), 
                        estimates=c(unmrkd_lagomorph_95[2], unmrkd_lagomorph_50[2], recap_lago$p),
                        CI_inf=c(unmrkd_lagomorph_95[5], unmrkd_lagomorph_50[5], recap_lago$p_inf),
                        CI_sup=c(unmrkd_lagomorph_95[6], unmrkd_lagomorph_50[6], recap_lago$p_sup))
p_lago_df <- p_lago_df %>%
  mutate(method=factor(method, levels=c("2","1_50", "1_95")))

fig4_p_lago <- ggplot(data=p_lago_df, aes(x=method,y=estimates, group=method, col=method)) +
  geom_linerange(aes(ymin=CI_inf, ymax=CI_sup),position=position_dodge2(width=1), linewidth=0.8, color="gray60") +
  geom_point(position=position_dodge2(width=1), size=2.5) + # shape=21,  stroke=2,  fill="gray60"
  scale_color_manual(values=c("#F26419","#86BBD8","#2F4858"), labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) +
  # scale_y_continuous(breaks=seq(1.0150,1.0250, 0.0010)) +
  scale_x_discrete(labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) + # expand=c(0.1,0.1),
  theme_light() + ylim(0,1) +
  labs(x=NULL,y=NULL) +
  theme(panel.grid.major.x = element_blank(), legend.position="none",
        axis.text.x = element_text(color = "gray20",size = 10))



# CHAMOIS Occupancy  ----
## Preparing identification matrix ... ----
### ... from ground truth data ----
W1 <- pics[,c(1,3,4,7)] %>%
  group_by(point,year,month)

W1 <- W1 %>% 
  mutate(chamois=ifelse(espece=="chamois",1,0), .before = espece)

# Replace year and month by the sampling occasion
W1$occas <- map2_dbl(W1$year, W1$month, function(x,y){
  if(x==2017){
    occas = y
  } 
  if(x==2018){
    occas = y+12
  }
  if(x==2019){
    occas = y+24
  }
  return(occas=occas)
})

thresh <- 0.5
Y_manual <- as.data.frame(detect_simulate(point, W1, thresh))
colnames(Y_manual) <- c(occas2017, occas2018, occas2019)
rownames(Y_manual) <- unique(pics$point)
Y_manual


### ...from DL predictions ----
W0 <- ref[,c(1,3,4,11)] %>%
  group_by(point,year,month)

# Replace year and month by the sampling occasion
W0$occas <- map2_dbl(W0$year, W0$month, function(x,y){
  if(x==2017){
    occas = y
  } 
  if(x==2018){
    occas = y+12
  }
  if(x==2019){
    occas = y+24
  }
  return(occas=occas)
})


# score threshold <- 0.5
Y_DL_50 <- as.data.frame(detect_simulate(point, W0, 0.5))
colnames(Y_DL_50) <- c(occas2017, occas2018, occas2019)
rownames(Y_DL_50) <- unique(pics$point)
Y_DL_50

# score threshold <- 0.95
Y_DL_95 <- as.data.frame(detect_simulate(point, W0, 0.95))
colnames(Y_DL_95) <- c(occas2017, occas2018, occas2019)
rownames(Y_DL_95) <- unique(pics$point)
Y_DL_95

## Occupancy estimates with unmarked package ----
### Estimates from manually labelled data
# Get data for 2017 into unmarked format
frame_chamois <- unmarked::unmarkedFrameOccu(y=Y_manual[,1:12]) # FOR 2017
summary(frame_chamois)
# fit
mod_chamois <- unmarked::occu(~1~1, frame_chamois, linkPsi = "logit")
# estimates
summary(mod_chamois)
unmarked::predict(mod_chamois, 'state')[1,]
unmarked::predict(mod_chamois, 'det')[1,] ; backTransform(mod_chamois, "det")
# save estimates and confidence intervals (90%)
unmrkd_chamois_manu <- mod_chamois@opt$par
unmrkd_chamois_manu <- plogis(c(unmrkd_chamois_manu, 
                                confint(mod_chamois,type="state", level=0.95, method="profile"),
                                confint(mod_chamois,type="det", level=0.95, method="profile")))

### Estimates from DL output with threshold at 0.50
frame_chamois <- unmarked::unmarkedFrameOccu(y=Y_DL_50[,1:12]) # FOR 2017
summary(frame_chamois)
# fit
mod_chamois <- unmarked::occu(~1~1, frame_chamois, linkPsi = "logit")
# estimates
summary(mod_chamois)
unmarked::predict(mod_chamois, 'state')[1,] ; unmarked::predict(mod_chamois, 'det')[1,]
# save estimates and confidence intervals (90%)
unmrkd_chamois_50 <- mod_chamois@opt$par
unmrkd_chamois_50 <- plogis(c(unmrkd_chamois_50, 
                              confint(mod_chamois,type="state", level=0.95, method="profile"),
                              confint(mod_chamois,type="det", level=0.95, method="profile")))

### Estimates from DL output with threshold at 0.95
frame_chamois <- unmarked::unmarkedFrameOccu(y=Y_DL_95[,1:12]) # FOR 2017
summary(frame_chamois)
# fit
mod_chamois <- unmarked::occu(~1~1, frame_chamois, linkPsi = "logit")
# estimates
summary(mod_chamois)
unmarked::predict(mod_chamois, 'state')[1,] ; unmarked::predict(mod_chamois, 'det')[1,]
# save estimates and confidence intervals (90%)
unmrkd_chamois_95 <- mod_chamois@opt$par
unmrkd_chamois_95 <- plogis(c(unmrkd_chamois_95, 
                              confint(mod_chamois,type="state", level=0.95, method="profile"),
                              confint(mod_chamois,type="det", level=0.95, method="profile")))


## MONTE CARLO simulations (x1000)  ----
k <- 100 # number of simulations

# the probability is the score
recap <- data.frame()
Y_map <- matrix(0,nrow=n_point, ncol=12)
ptm <- proc.time()[3]

set.seed(123)
for(i in 1:k){
  # Simulate the identification matrix from scores
  Wbis <- subset(W0, year==2017) # W0
  Wbis$chamois <- rbinom(nrow(Wbis), 1, W0$chamois[W0$year==2017])
  # Simulate the detection matrix
  Ybis <- as.data.frame(detect_simulate(point, Wbis, 0.5))
  
  # save the details of detection matrix
  # Y_map <- Ybis
  Y_map <- as.matrix(Y_map) + as.matrix(Ybis[,1:12])
  
  # put data in unmarked format
  frame <- unmarked::unmarkedFrameOccu(y=Ybis[,1:12])
  # fit
  mod_chamois <- unmarked::occu(~1~1, frame, linkPsi = "logit")
  
  # Save estimates in CI
  recap[i,1] <- plogis(mod_chamois@opt$par[1])
  recap[i,2:3] <- plogis(confint(mod_chamois, type="state", level=0.95, method="profile"))
  recap[i,4] <- plogis(mod_chamois@opt$par[2])
  recap[i,5:6] <- plogis(confint(mod_chamois, type="det", level=0.95, method="profile"))
  
  # Save the number of occupied sites (at least one identification)
  recap[i,7] <-  sum(unlist(lapply(apply(Ybis, 1,sum, na.rm=T),
                                   function(x){ifelse(x==0, return(0), return(1))})))
}
print(proc.time()[3]-ptm)
colnames(recap) <- c("psi","psi_inf", "psi_sup","p", "p_inf", "p_sup", "N_occ")
beepr::beep(1)


write.csv(recap, file="./recap_chamois_100.csv")

# recap_chamois <- read.csv(file="./recap_chamois_100.csv",
                          row.names=1)
recap_chamois <- recap

## Fig4 - estimates CHAMOIS ----
### occupancy ----
recap_chamois$psi_inf[recap_chamois$psi_inf==0] <- NA
psi_chamois_df <- data.frame(method=c(c("1_95","1_50"), rep("2",100)), sample=c(rep(1,2),1:100),
                             estimates=c(unmrkd_chamois_95[1], unmrkd_chamois_50[1], recap_chamois$psi),
                             CI_inf=c(unmrkd_chamois_95[3], unmrkd_chamois_50[3], recap_chamois$psi_inf),
                             CI_sup=c(unmrkd_chamois_95[4], unmrkd_chamois_50[4], recap_chamois$psi_sup))
psi_chamois_df <- psi_chamois_df %>%
  mutate(method=factor(method, levels=c("2","1_50", "1_95")))

fig4_psi_chamois <- ggplot(data=psi_chamois_df, aes(x=method,y=estimates, group=method, col=method)) +
  geom_linerange(aes(ymin=CI_inf, ymax=CI_sup),position=position_dodge2(width=1), linewidth=0.8, color="gray60") +
  geom_point(position=position_dodge2(width=1), size=2.5) + # shape=21,  stroke=2,  fill="gray60"
  scale_color_manual(values=c("#F26419","#86BBD8","#2F4858"), labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) +
  # scale_y_continuous(breaks=seq(1.0150,1.0250, 0.0010)) +
  scale_x_discrete(labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) + # expand=c(0.1,0.1),
  theme_light() + ylim(0.5,1) +
  labs(x=NULL,y=NULL) + # , title="Chamois", y="Occupancy proabability (CI 95%)"
  theme(panel.grid.major.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(color = "gray20",size = 10))

### detection ----
# recap$p_inf[recap$p_inf==0] <- NA
p_chamois_df <- data.frame(method=c(c("1_95","1_50"), rep("2",100)), sample=c(rep(1,2),1:100),
                           # data.frame(method=c(c("1_50","1_95"), rep(2,100)), sample=c(rep(1,2),1:100), 
                           estimates=c(unmrkd_chamois_95[2], unmrkd_chamois_50[2], recap$p),
                           CI_inf=c(unmrkd_chamois_95[5], unmrkd_chamois_50[5], recap$p_inf),
                           CI_sup=c(unmrkd_chamois_95[6], unmrkd_chamois_50[6], recap$p_sup))
p_chamois_df <- p_chamois_df %>%
  mutate(method=factor(method, levels=c("2","1_50", "1_95")))

fig4_p_chamois <- ggplot(data=p_chamois_df, aes(x=method,y=estimates, group=method, col=method)) +
  geom_linerange(aes(ymin=CI_inf, ymax=CI_sup),position=position_dodge2(width=1), linewidth=0.8, color="gray60") +
  geom_point(position=position_dodge2(width=1), size=2.5) + # shape=21,  stroke=2,  fill="gray60"
  scale_color_manual(values=c("#F26419","#86BBD8","#2F4858"), labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) +
  # scale_y_continuous(breaks=seq(1.0150,1.0250, 0.0010)) +
  scale_x_discrete(labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) + # expand=c(0.1,0.1),
  theme_light() + ylim(0,1) +
  labs(x=NULL,y=NULL) + # y="Detection proabability (CI 95%)"
  theme(panel.grid.major.x = element_blank(), legend.position="none",
        axis.text.x = element_text(color = "gray20",size = 10))


# Save ESTIMATES ----
estimates <- as.data.frame(rbind(unmrkd_lynx_manu, unmrkd_lynx_50, unmrkd_lynx_95,
                                 unmrkd_lagomorph_manu, unmrkd_lagomorph_50, unmrkd_lagomorph_95,
                                 unmrkd_chamois_manu, unmrkd_chamois_50, unmrkd_chamois_95))
estimates <- estimates %>% 
  mutate(type=rep(c("manu","DL50","DL95"),3), .before=V1) %>%
  mutate(esp=rep(c("lynx","lagomorphe","chamois"),each=3), .after=type) %>%
  tibble::remove_rownames()
colnames(estimates)[3:8] <- c("psi","p","psi_inf", "psi_sup","p_inf", "p_sup")
estimates

# write.csv(estimates, file="./estimates.csv")
estimates <- read.csv(file="./estimates.csv", row.names=1)


# PLOT :Figure3 ----
col_fig3 <- data.frame(n_sites=18:29, 
                       col=c("#0C1E2F","#22403B", "#4A976C", 
                             "#B9C968","#ECE08B","#ECBE6B", "#F9A53C","#D67D21",#DE8E72
                             "#f06a45","#BD3432", "#762936", "#401744"))


## LYNX ----
p0_lynx <- ggplot(recap_lynx, aes(x=p, y=psi, col=as.factor(N_occ))) + 
  geom_point(alpha=1/3, size=4) +
  scale_color_manual(values=col_fig3$col[col_fig3$n_sites %in% 25:27]) +
  scale_x_continuous(breaks=seq(0.54,0.62,0.02)) +
  scale_y_continuous(breaks=seq(0.75,1,0.05), minor_breaks=seq(0.73,1,0.01)) +
  theme_bw() + theme(legend.direction = "horizontal")+ 
  theme(legend.position="none") + labs(x=NULL, y=NULL) +
  theme(axis.text=element_text(size=10, color="black"))
p0_lynx
# with marginal histogram 
p4_lynx <- ggMarginal(p0_lynx, type='density', fill="gray70",col=NA, bw=0.005)
p4_lynx

## LAGO ----
p0_lago <- ggplot(recap_lago, aes(x=p, y=psi, col=as.factor(N_occ))) +
  geom_point(alpha=1/3, size=4) + 
  scale_color_manual(values=col_fig3$col[col_fig3$n_sites %in% 18:23]) +
  scale_y_continuous(breaks=seq(0.50,1,0.05), minor_breaks=seq(0.5,1,0.01)) +
  scale_x_continuous(breaks=seq(0.26,0.38,0.02)) +
  theme_bw() + theme(legend.position="none") + labs(x=NULL, y=NULL) +
  theme(axis.text=element_text(size=10, color="black"))
p0_lago
# with marginal histogram 
p4_lago <- ggMarginal(p0_lago, type='density', fill="gray70", col=NA,bw=0.005)
p4_lago

## Chamois ----
p0_chamois <- ggplot(recap_chamois, aes(x=p, y=psi, col=as.factor(N_occ))) + 
  geom_point(alpha=1/3, size=4) + 
  scale_color_manual(values=col_fig3$col[col_fig3$n_sites %in% 24:29]) +
  scale_x_continuous(breaks=seq(0.5,0.6,0.02)) +
  scale_y_continuous(breaks=seq(0.75,1,0.05), minor_breaks=seq(0.73,1,0.01)) +
  theme_bw() + theme(legend.direction = "horizontal")+ # labs(col="N sites occupés") 
  theme(legend.position="none") + labs(x=NULL, y=NULL) +
  theme(axis.text=element_text(size=10, color="black"))
p0_chamois
# with marginal histogram 
p4_chamois <- ggMarginal(p0_chamois, type='density', fill="gray70",col=NA, bw=0.005)
p4_chamois

cowplot::plot_grid(p4_lynx, p4_lago, p4_chamois, nrow=1)


# PLOT :Figure4 ----
cowplot::plot_grid(fig4_psi_lynx,fig4_psi_lago,fig4_psi_chamois,
                   fig4_p_lynx,fig4_p_lago,fig4_p_chamois,
                   nrow=2, byrow=T)



# RELATIVE uncertainties COVERAGE statistics ----
recouvrement <- function(simu,pics){
  inter <- max(0, (min(simu[2], pics[2], na.rm=T)-max(simu[1], pics[1], na.rm=T)))
  long <- simu[2] - simu[1]# pics[2]-pics[1]
  return((inter/long)*100)
}

# /!\ to execute one after another /!\
recap <- recap_lynx ; est <- unmrkd_lynx_50
recap <- recap_lago; est <- unmrkd_lagomorph_50
recap <- recap_chamois ; est <- unmrkd_chamois_50

psi_pics <- c(est[3], est[4])
total_uncert_psi <- c(min(recap$psi_inf, na.rm=T), max(recap$psi_sup, na.rm=T))
recouvrement(total_uncert_psi,psi_pics)

p_pics <- c(est[5], est[6])
total_uncert_p <- c(min(recap$p_inf, na.rm=T), max(recap$p_sup, na.rm=T))
recouvrement(total_uncert_p, p_pics)
total_uncert_p

# DL uncertainties
nrow(recap[recap$psi>=unmrkd_lynx_manu[3] & recap$psi<=unmrkd_lynx_manu[4],])
nrow(recap[recap$p>=unmrkd_lynx_manu[5] & recap$p<=unmrkd_lynx_manu[6],])

nrow(recap_chamois[recap_chamois$psi<=unmrkd_chamois_50[4],])
