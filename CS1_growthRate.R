######################################################################
# Propagating uncertainties from species identification to occupancy #
#                           with MONTE CARLO                         #
######################################################################

# CASE STUDY 1 : Stochastic growth rate prediction from RFID data 
  
library(tidyverse)
library(ggridges)
library(popbio)

# Import data  ----
# Breeding data given by Deep Learning (RFIDeep)
breeding_data <- read.csv("pred_breeding_output_sex_SCORES.csv")

breeding_data['breeding_outcome'] = breeding_data$pred_S_F
breeding_data[which(breeding_data['pred_NB_B'] == 'NB'),'breeding_outcome'] = 'NB'
head(breeding_data)


# Data on birds (selection from the database, only birds tagged as chicks from Crozet)
cmr_data =  read.csv("./Gael/birds_data_processed_27_09_2023.csv")

cmr_data[which(cmr_data$rfid == 'A 00000 0 964 001013464988'),'birth_year'] = 2016 # correction of a mistake in the database
cmr_data_chicks <- subset(cmr_data, rfid_stage=='Chick')
head(cmr_data_chicks)

# Import detection data as time series
# Detections corrected (same algorithm than the one used in data augmentation step of RFIDeep)
detections = read.csv('Gael/corrected_detections_25_04_2023.csv')

detections$dtime = as.Date(detections$dtime, "%Y-%m-%d %H:%M:%S")
head(detections)


# Building of the 'CMR' table ----
presenceTable <- data.frame(rfid=cmr_data_chicks$rfid, birth_year=cmr_data_chicks$birth_year)  
list_year = 1999:2022

for (year in list_year){
  detect_year = unique(detections[which((detections$dtime > paste(year-1, "-09-01", sep='')) &
                                          (detections$dtime < paste(year+2, "-09-01", sep=''))),'rfid'])
  presence_vector = c()
  for (i in 1:nrow(cmr_data_chicks)){
    if (cmr_data_chicks[i, 'rfid_year'] > year){
      presence_vector = c(presence_vector, -1)
    } else if ((cmr_data_chicks[i, 'rfid_year'] == year) | (cmr_data_chicks[i, 'rfid'] %in% detect_year)){
      presence_vector = c(presence_vector, 1)
    }else{
      presence_vector = c(presence_vector, 0)
    }
  }
  df_year = data.frame(presence_vector)
  colnames(df_year)= year
  presenceTable = cbind(presenceTable, df_year)
}
head(presenceTable)

# indiv in the presence table for which we have no data on reproduction: 284
not_in_breeding_data <- setdiff(presenceTable$rfid, breeding_data$rfid)
subset(presenceTable, rfid %in% not_in_breeding_data)
# they never came back after their first departure, in other words they died between 1 and 2 years ago 

# indiv for which we have data on reproduction but not in the presence table: 4033   
not_in_presenceTable <- setdiff(breeding_data$rfid, presenceTable$rfid)
unique(subset(cmr_data, rfid %in% not_in_presenceTable, birth_year))


## Population matrix by age group and year ----
# Population size calculation
calcul_effectif <- function(pop, thresh){
  dta <- data.frame()
  c <- 1
  for (yr in 1999:2022){
    cohort <- presenceTable %>% 
      filter(birth_year < yr & .data[[as.character(yr)]]==1)
    if(pop=="B"){
      cohort <- cohort %>%
        left_join(breeding_data[,c(1,2,6,8)], by = "rfid") %>%
        filter(year==yr & score_B>=thresh)
    }
    if(pop=="S"){
      cohort <- cohort %>%
        left_join(breeding_data[,c(1,2,6,8)], by = "rfid") %>%
        filter(year== yr & score_B>=thresh & score_S>=thresh)
    }
    # else{}
    cohort <- cohort %>%
      group_by(birth_year) %>% 
      summarise(n=n()) %>%
      mutate(age=yr - birth_year)
    # }
    dta[cohort$age,c] <- cohort$n
    c <- c+1
  }
  colnames(dta) <- 1999:2022
  # Appliquer l'hypothèse des 10 ans max
  if(pop %in% c("B", "S")){
    # if(pop == "S"){
    dta[1:4,] <- dta[1:4,] %>% mutate_at(-(1:3), ~replace(., is.na(.), 0))
    dta10 <- dta[1:10,]
    dta10[10,] <- colSums(dta[10:nrow(dta),], na.rm=T)
  }
  else{
    # dta <- dta %>% mutate_all(~replace(., is.nan(.), 0))
    dta10 <- dta[1:11,]
    dta10[11,] <- colSums(dta[11:nrow(dta),], na.rm=T)
  }
  return(dta10)
}


# Survival no longer changes from the age of 10 (11th year)
Nt10 <- calcul_effectif(pop="non",thresh=NA)

## SURVIVAL ----

surv <- matrix(nrow=10, ncol=23)
for(i in 1:10){
  for(j in 1:23){
    if(i==10){
      # if(j=>10){
      surv[i,j] <- (Nt10[10,j+1]+Nt10[11,j+1]-(Nt10[9,j]*surv[9,j])) / (Nt10[10,j]+Nt10[11,j])
    } else{
      surv[i,j] <- Nt10[i+1,j+1] / Nt10[i,j]
    }
  }
}
colnames(surv) <- 1999:2021
surv

taux_surv <- data.frame(Year=2008:2021, 
                        mean=colMeans(surv[,10:23], na.rm=T), 
                        sd=apply(surv[,10:23], 2, sd))
# PLOT
# ggplot(taux_surv) +
#   geom_line(aes(Year, mean), size=1.5, color="gray80") +
#   geom_point(aes(Year, mean), size=3, color="gray20") +
#   geom_errorbar(aes(x=Year, ymin = (mean - 0.5*sd), ymax = (mean + 0.5*sd)), color="gray20") + 
#   scale_x_continuous(breaks=taux_surv$Year) +
#   labs(x=NULL, y='Mean Survival Rate', 
#        title="Average survival per year (all age class combined)") + theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# taux_surv_age <- data.frame(Age=1:10, 
#                             mean=rowMeans(surv[,10:23], na.rm=T), 
#                             sd=apply(surv[,10:23], 1, sd))
# ggplot(taux_surv_age) +
#   geom_line(aes(Age, mean), size=1.5, color="gray80") +
#   geom_point(aes(Age, mean), size=3, color="gray20") +
#   geom_errorbar(aes(x=Age, ymin = (mean - 0.5*sd), ymax = (mean + 0.5*sd)), color="gray20") + 
#   scale_x_continuous(breaks=taux_surv_age$Age) +
#   labs(x=NULL, y='Mean Survival Rate',
#        title="Average survival by age class (all years combined)") + theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


## FECUNDITY ----

# Breeding individuals matrix
Ft05 <- calcul_effectif(pop="S", thresh=0.5)
Ft05

for(i in 2008:2020){
  print(nrow(breeding_data[breeding_data$year==i & breeding_data$breeding_outcome=="S",]))
  print(length(unique(breeding_data$rfid[breeding_data$year==i & breeding_data$breeding_outcome=="S"])))
}

calcul_fecond <- function(num, denom){
  fecond <- num[1:9,]/ denom[1:9,] # Ft05_10[1:9,]/ Nt10[1:9,]
  if(nrow(denom)>10){
    fecond <- rbind(fecond, num[10,]/(denom[10,]+denom[11,]))
  }
  else{
    fecond <- rbind(fecond, num[10,]/(denom[10,]))
    fecond <- fecond %>% mutate_all(~replace(., is.nan(.), 0))
  }
  colnames(fecond) <- 1999:2022
  return(fecond)
}

# Calculation of fertility (threshold=0.5) taking only breeding individuals
NBreeder_05 <- calcul_effectif(pop="B", 0.5)

fecond05 <- calcul_fecond(Ft05, NBreeder_05)
round(fecond05,3)

# Applying a threshold of 0.95
Ft095 <- calcul_effectif(pop="S", 0.95)
NBreeder_095 <- calcul_effectif(pop="B", 0.95)

fecond095 <- calcul_fecond(Ft095, NBreeder_095)
round(fecond095,3)

taux_repro <- data.frame(Year=2008:2021, 
                         mean05=colMeans(fecond05[,10:23], na.rm=T), 
                         sd05=apply(fecond05[,10:23],2,sd),
                         mean095=colMeans(fecond095[,10:23], na.rm=T), 
                         sd095=apply(fecond095[,10:23],2,sd))
# # PLOT 
# ggplot(taux_repro) +                                      
#   geom_line(aes(Year, mean05), size=1.2, color="#86BBD8", alpha=0.3) +
#   geom_point(aes(Year, mean05), size=3, color="#86BBD8") +
#   geom_errorbar(aes(x=Year, ymin = (mean05 - 0.5*sd05), ymax = (mean05 + 0.5*sd05)), color="#86BBD8", width=0.5) + 
#   geom_line(aes(Year, mean095), size=1.2, color="#2F4858", alpha=0.3) +
#   geom_point(aes(Year, mean095), size=3, color="#2F4858") +
#   geom_errorbar(aes(x=Year, ymin = (mean095 - 0.5*sd095), ymax = (mean095 + 0.5*sd095)), color="#2F4858", width=0.5) + 
#   scale_x_continuous(breaks=2007:2022, minor_breaks =NULL) + # ylim(c(0,1)) +
#   # scale_color_manual(values=c("data1"="dark green", "data2"="purple")) +
#   labs(x=NULL, y='Mean fecundity', title="Average annual fecundity") + theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# taux_repro_age <- data.frame(Age=1:10, 
#                              mean05=rowMeans(fecond05[,10:23], na.rm=T), 
#                              sd05=apply(fecond05[,10:23],1,sd),
#                              mean095=rowMeans(fecond095[,10:23], na.rm=T), 
#                              sd095=apply(fecond095[,10:23],1,sd))
# ggplot(taux_repro_age) +                                      
#   geom_line(aes(Age, mean05), size=1.2, color="#86BBD8", alpha=0.3) + # color="gray80"
#   geom_point(aes(Age, mean05), size=3, color="#86BBD8") + # color="gray20"
#   geom_errorbar(aes(x=Age, ymin = (mean05 - 0.5*sd05), ymax = (mean05 + 0.5*sd05)), color="#86BBD8", width=0.5) + 
#   geom_line(aes(Age, mean095), size=1.2, color="#2F4858", alpha=0.3) + # color="pink"
#   geom_point(aes(Age, mean095), size=3, color="#2F4858") + # color="darkred"
#   geom_errorbar(aes(x=Age, ymin = (mean095 - 0.5*sd095), ymax = (mean095 + 0.5*sd095)), color="#2F4858", width=0.5) + 
#   scale_x_continuous(breaks=1:10, minor_breaks=NULL) +
#   # scale_color_manual(name="bla", breaks=c("0.5","0.95"), values=c("#86BBD8", "#2F4858")) +
#   # scale_color_manual(values=c("data1"="dark green", "data2"="purple")) +
#   labs(x=NULL, y='Average fecundity', title="Average fecundity by age class") + theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



## GROWTH Rate estimation ----
# Building transtion matrices between 2008 and 2021
years <- 2008:2021

# threshold 50
list_matrices05 <- map(years, function(year) {
  pop_matrix <- matrix(0, nrow = 10, ncol = 10)
  # pop_matrix[1,] <- fecond[, as.character(year)]
  pop_matrix[1,] <- fecond05[, as.character(year)]
  pop_matrix[2:10, 1:9] <- diag(surv[1:9, as.character(year)])
  pop_matrix[10,10] <- surv[10, as.character(year)]
  pop_matrix[is.na(pop_matrix)] <- 0
  return(pop_matrix)
})
list_matrices05[1:2]

# threshold 95
list_matrices095 <- map(years, function(year) {
  pop_matrix <- matrix(0, nrow = 10, ncol = 10)
  pop_matrix[1,] <- fecond095[, as.character(year)]
  pop_matrix[2:10, 1:9] <- diag(surv[1:9, as.character(year)])
  pop_matrix[10,10] <- surv[10, as.character(year)]
  pop_matrix[is.na(pop_matrix)] <- 0
  return(pop_matrix)
})
list_matrices095[1:2]

stoch_lambda05 <- popbio::stoch.growth.rate(list_matrices05, verbose=T)
stoch_lambda05
stoch_lambda05$approx ; stoch_lambda05$sim.CI

stoch_lambda095 <- popbio::stoch.growth.rate(list_matrices095, verbose=T)
stoch_lambda095
stoch_lambda095$approx ; stoch_lambda095$sim.CI

sens <- popbio::stoch.sens(list_matrices05)

# MONTE CARLO simulations ----
# Prepare inputs
fecond05_Nt <- calcul_fecond(Ft05, Nt10)
fecond095_Nt <- calcul_fecond(Ft095, Nt10)
taux_repro_Nt <- data.frame(Year=2008:2021, 
                            mean05=colMeans(fecond05_Nt[,10:23], na.rm=T), 
                            sd05=apply(fecond05_Nt[,10:23],2,sd),
                            mean095=colMeans(fecond095_Nt[,10:23], na.rm=T), 
                            sd095=apply(fecond095_Nt[,10:23],2,sd))

# initialize tables
k = 100
success <- data.frame() ; breeders <- data.frame()
taux_repro_MC <- data.frame(Year=2008:2021)
growth_rates_MC <- data.frame()

set.seed(123)
for(i in 1:k){
  # Matrice du nombre de reproducteurs par cohorte par an
  c <- 1
  for (yr in 1999:2022){
    cohort <- presenceTable %>% 
      filter(birth_year < yr & .data[[as.character(yr)]]==1) %>%
      left_join(breeding_data[,c(1,2,6,8,9)], by = "rfid") #%>%
    # Breeding success is conditionned by status
    cohort$B = rbinom(nrow(cohort),1, cohort$score_B)
    cohort$B_S_drawing = ifelse(cohort$B==1, rbinom(nrow(cohort),1, cohort$score_S) ,0)
    # N successful breeders
    suc <- cohort %>%
      filter(year==yr & B_S_drawing==1) %>%
      mutate(age=yr - birth_year) %>%
      group_by(age) %>%
      summarize(n=n())
    success[suc$age,c] <- suc$n
    # N breeders
    bre <- cohort %>%
      filter(year==yr & B==1) %>%
      mutate(age=yr - birth_year) %>%
      group_by(age) %>%
      summarize(n=n())# %>%
    # pull(n)
    breeders[bre$age,c] <- bre$n
    c <- c+1
  }
  success[1:4,] <- success[1:4,] %>% 
    mutate_at(-(1:3), ~ replace(., is.na(.), 0))
  breeders[1:4,] <- breeders[1:4,] %>%
    mutate_at(-(1:3), ~ replace(., is.na(.), 0))
  
  # constant after 10 years old
  success10 <- success[1:10,]
  success10[10,] <- colSums(success[10:nrow(success),], na.rm=T)
  breeders10 <- breeders[1:10,]
  breeders10[10,] <- colSums(breeders[10:nrow(breeders),], na.rm=T)
  
  fecond_MC <- calcul_fecond(success10, breeders10)
  
  taux_repro_MC <- cbind(taux_repro_MC,
                         colMeans(fecond_MC[,10:23], na.rm=T))
  
  list_matrices_MC <- map(2008:2021, function(year) {
    pop_matrix <- matrix(0, nrow = 10, ncol = 10)
    pop_matrix[1,] <- fecond_MC[, as.character(year)]
    pop_matrix[2:10, 1:9] <- diag(surv[1:9, as.character(year)])
    pop_matrix[10,10] <- surv[10, as.character(year)]
    pop_matrix[is.na(pop_matrix)] <- 0
    return(pop_matrix)
  })
  
  stoch_lambda_MC <- popbio::stoch.growth.rate(list_matrices_MC)

  growth_rates_MC = rbind(growth_rates_MC,
                          c(log(Re(eigen(mean(list_matrices_MC))$values[1])),
                            stoch_lambda_MC$approx,stoch_lambda_MC$sim.CI))
}
colnames(growth_rates_MC) <- c("average","approx", "IC_inf", "IC_sup")
head(growth_rates_MC)

fecond_simus <- taux_repro_MC %>% pivot_longer(cols = starts_with('col'))

# PLOT Fig 1 ----
ggplot() +                           
  geom_violin(data=fecond_simus, aes(x=Year, y=value, group=Year), 
              width=0.8, col="#F26419") + # col="gray30""#F6AE2D"
  stat_summary(data=fecond_simus, aes(x=Year, y=value, group=Year), 
               fun.y=mean, geom="point", size=2, shape=20, col="#F26419") + 
  geom_point(data=taux_repro, aes(x=Year, y=mean05), col="#86BBD8", size=3)+
  geom_line(data=taux_repro, aes(x=Year, y=mean05), col="#86BBD8") + 
  geom_point(data=taux_repro, aes(x=Year, y=mean095), col="#2F4858", size=3)+
  geom_line(data=taux_repro, aes(x=Year, y=mean095), col="#2F4858") + 
  scale_x_continuous(breaks=2007:2022, minor_breaks = NULL) +
  labs(x=NULL, y='Average fecundity') + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Range of variation between simulations in 2017 and 2008 :
max(fecond_simus$value[fecond_simus$Year==2010]) - min(fecond_simus$value[fecond_simus$Year==2010]) ;
max(fecond_simus$value[fecond_simus$Year==2018]) - min(fecond_simus$value[fecond_simus$Year==2018]) 
max(fecond_simus$value[fecond_simus$Year==2019]) - min(fecond_simus$value[fecond_simus$Year==2019])

avg_fecond <- fecond_simus %>%
  group_by(Year) %>%
  summarize(avg=mean(value))
mean(avg_fecond$avg - taux_repro$mean05)
mean(avg_fecond$avg - taux_repro$mean095)

# PLOT FIG 2 ----
growth_rates_MC2 <- growth_rates_MC[sample(1:100, size=100, replace=F),]

plot_df <- data.frame(method=c(c("1_95","1_50"), rep("2",100)), sample=c(rep(1,2),1:100),
                      # data.frame(method=c(c("1_50","1_95"), rep(2,100)), sample=c(rep(1,2),1:100), 
                      estimates=c(stoch_lambda095$approx, stoch_lambda05$approx, growth_rates_MC2$approx),
                      CI_inf=c(stoch_lambda095$sim.CI[1], stoch_lambda05$sim.CI[1], growth_rates_MC2$IC_inf),
                      CI_sup=c(stoch_lambda095$sim.CI[2], stoch_lambda05$sim.CI[2], growth_rates_MC2$IC_sup))
plot_df <- plot_df %>%
  mutate(method=factor(method, levels=c("2","1_50", "1_95")))
# fct_relevel

ggplot(data=plot_df, aes(x=method,y=estimates, group=method, col=method)) +
  geom_point(position=position_dodge2(width=1), size=4) + # shape=21,  stroke=2,  fill="gray60"
  scale_color_manual(values=c("#F26419","#86BBD8","#2F4858"), labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) +
  scale_y_continuous(breaks=seq(0.05,0.3, 0.05), minor_breaks=seq(0,0.3, 0.01)) +
  scale_x_discrete(labels=c("MC simulations","Single-class 0.5", "Single-class 0.95")) + # expand=c(0.1,0.1),
  theme_light() + labs(x=NULL,y="Stochastic growth rate (CI 95%)") +
  theme(panel.grid.major.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(color = "gray20",size = 10))
