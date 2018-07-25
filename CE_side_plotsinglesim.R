# MR, 16.07.2018 
########### Side: plot the output of sovling ODE with 1 type, 1 scenario 


# copied from R_ce folder: CE_plot_.R

# Plot the outputs from core script

sim <- as.data.frame(ode(init, time, HPV_1, parms=params))

colnames(sim) <- col_nam
sim[1:10,1:10]

#Use R script: CE_core_script.R to make the simulations which can be plotted here


#define de denominator depending on which subgroup we are interested in
w_l_a1 <- (sr*agr[1]* rgr[1]) # women, low risk group, age 1
w_h_a1 <- (sr*agr[1]* rgr[2]) #women, high risk group, age 1
w_lh_a <- c( sr * agr)

#color sets
col1 <- colorRampPalette(c("grey20", "lightgrey"))(length(agegroups))
col2 <- colorRampPalette(c("red", "orange"))(length(agegroups))
col3 <- colorRampPalette(c("darkblue", "lightblue"))(length(agegroups))
col4 <- colorRampPalette(c("darkgreen", "lightgreen"))(length(agegroups))

colS <- colorRampPalette(c("orange", "yellow"))(length(agegroups))
colI <- colorRampPalette(c("darkred", "darkorange"))(length(agegroups))
colR <- colorRampPalette(c("darkblue", "lightblue"))(length(agegroups))
colV <- colorRampPalette(c("darkgreen", "lightgreen"))(length(agegroups))

timepl <- c(0,200)
timepl <- c(0,100)

# sim <- read_csv("O:/Test_R/R/2017_juin_HPV_cost-effectiveness/sim_vacc_05_younggirls.csv")
# sim <- as.data.frame(sim[,2:(length(sim[1,] ) )] )
summary(sim)

sim[,1:2]

############################### General plots to test how the model behaves #################################

####### PLOT 1 (test plot)
#Women, first age group, high and low risk groups (straight and dashed lines, respectively)
#######
plot(NA, xlim=timepl, ylim=c(0,1), 
     xlab= "years",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)

lines(1:max(time+1),  (sim$`S_f_l_15-19` ) / w_l_a1 ,  lty=1, col="yellow", lwd=3)
lines(1:max(time+1),  (sim$`I_f_l_15-19` ) / w_l_a1 ,  lty=1, col="red", lwd=3)
lines(1:max(time+1),  (sim$`R_f_l_15-19` ) / w_l_a1 ,  lty=1, col="blue", lwd=3)
lines(1:max(time+1),  (sim$`V4v_f_l_15-19` ) / w_l_a1 ,  lty=1, col="darkgreen", lwd=3)

lines(1:max(time+1),  (sim$`S_f_h_15-19` ) / w_h_a1 ,  lty=2, col="yellow", lwd=2)
lines(1:max(time+1),  (sim$`I_f_h_15-19` ) / w_h_a1  ,  lty=2, col="red", lwd=2)
lines(1:max(time+1),  (sim$`R_f_h_15-19` ) / w_h_a1  ,  lty=2, col="blue", lwd=2)
lines(1:max(time+1),  (sim$`V4v_f_h_15-19` ) / w_h_a1  ,  lty=2, col="darkgreen", lwd=2)

####### PLOT 2 (test plot)
#women, low and high risk groups together, different age groups (darker lines= younger age groups, lighter = older)
#######
plot(NA, xlim=timepl, ylim=c(0,1), 
     xlab= "years",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)

for(i in 1:4){
  lines(1:max(time+1),  (sim[,1+ ( (i-1)*(length(comp) * 4)  + 1)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 1)+1 ) ) ]  ) / w_lh_a[i] ,  lty=1, col=colS[i], lwd=2)
  lines(1:max(time+1),  (sim[,1+ ( (i-1)*(length(comp) * 4)  + 2)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 1)+2 ) ) ] ) / w_lh_a[i] ,  lty=1, col=colI[i], lwd=2)
  lines(1:max(time+1),  (sim[,1+ ( (i-1)*(length(comp) * 4)  + 12)] + sim[,1+ (  (i-1)*(length(comp) * 4)  +((length(comp) * 1)+12 ) ) ] ) / w_lh_a[i] ,  lty=1, col=colR[i], lwd=2)
  lines(1:max(time+1),  (sim[,1+ ( (i-1)*(length(comp) * 4)  + 13)]+ sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 1)) ) ]  ) / w_lh_a[i] ,  lty=1, col=colV[i], lwd=2)
}

####### PLOT 3 (test plot)
#women, low and high risk groups together, different age groups (thicker lines= younger age groups, thiner = older), for CIN 
####### 
plot(NA, xlim=timepl, ylim=c(0,0.003), 
     xlab= "years",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)

#CIN2
for(i in 1:4){
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 3 )] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+3) ) ] ) / w_lh_a[i] ,  lty=1, col=col1[i], lwd=1)
}
#CIN3
for(i in 1:4){
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 4)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+4) ) ] ) / w_lh_a[i] ,  lty=1, col=col2[i], lwd=1)
}

####### PLOT 4 (test plot)
# And now for cancer stages
####### 

1e-4 * 100000

plot(NA, xlim=timepl, ylim=c(0,3e-4), 
     xlab= "years",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)
#LCC (undiagn and diagn)
for(i in 1:4){
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 5)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+5) ) ] ) / w_lh_a[i] ,  lty=2, col=col3[i], lwd=1)
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 6)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+6) ) ] ) / w_lh_a[i] ,  lty=1, col=col3[i], lwd=2)
}

#DCC
for(i in 1:4){
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 7)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+7) ) ] ) / w_lh_a[i] ,  lty=2, col=col4[i], lwd=1)
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 8)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+8) ) ] ) / w_lh_a[i] ,  lty=1, col=col4[i], lwd=2)
}

#DCC
for(i in 1:4){
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 9)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+9) ) ] ) / w_lh_a[i] ,  lty=2, col=col1[i], lwd=1)
  lines(1:max(time+1),   (sim[,1+ ( (i-1)*(length(comp) * 4)  + 10)] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2)+10) ) ] ) / w_lh_a[i] ,  lty=1, col=col1[i], lwd=2)
}

####### PLOT 5 (test plot)
#### Total diagnosed cancer per age groups:
####### 
plot(NA, xlim=timepl, ylim=c(0,0.0005), 
     xlab= "years",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)

# lines(1:max(time+1),  sim$`ULCC_f_h_42-64`  / w_lh_a[i] ,  lty=1, col=col3[i], lwd=2 )

# for(i in 1:4){
#   lines(1:max(time+1),   (sim[,1+ ( (i-1)*88  + 5)] + sim[,1+ (  (i-1)*88  + (44+5) ) ] +
#                             sim[,1+ ( (i-1)*88  + 7)] + sim[,1+ (  (i-1)*88  + (44+7) ) ] +
#                             sim[,1+ ( (i-1)*88  + 9)] + sim[,1+ (  (i-1)*88  + (44+9) ) ]) 
#         / w_lh_a[i] ,  lty=2, col=col3[i], lwd=1) }

for(i in 1:4){
  lines(1:max(time+1),   (  sim[,1+ ( (i-1)*(length(comp) * 4)  + 6 ) ] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2) + 6) )  ] +
                              sim[,1+ ( (i-1)*(length(comp) * 4)  + 8 ) ] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2) + 8) )  ] +
                              sim[,1+ ( (i-1)*(length(comp) * 4)  + 10) ] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2) + 10) ) ]) 
        / w_lh_a[i] ,  lty=1, col=col3[i], lwd=1)
}


###Plot number of incident cancers   
plot(NA, xlim=timepl, ylim=c(0,0.0005), 
     xlab= "years",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)

lines(1:max(time+1),   ( ( sim[,1+ ( (i-1)*(length(comp) * 4)  + 6 ) ] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2) + 6) )  ] / sr )  +
                           ( sim[,1+ ( (i-1)*(length(comp) * 4)  + 8 ) ] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2) + 8) )  ] ) /sr +
                           (sim[,1+ ( (i-1)*(length(comp) * 4)  + 10) ] + sim[,1+ (  (i-1)*(length(comp) * 4)  + ((length(comp) * 2) + 10) ) ] ) /sr  ) 
      / w_lh_a[i] ,  lty=1, col=col3[i], lwd=1)

( z[1,y]*ssoa[5] + z[2, y]*ssoa[7] +  z[3, y]*ssoa[9] ) * 4.0*10^6 



############################### Plots to be compared over different simulations (with/withoug vaccine etc. #############################


sim1 <- read_csv("O:/Test_R/R/2017_juin_HPV_cost-effectiveness/sim_no_vacc.csv")
sim1 <- as.data.frame(sim1[,2:(length(sim1[1,] ) )] )

sim2 <- read_csv("O:/Test_R/R/2017_juin_HPV_cost-effectiveness/sim_vacc_05_younggirls.csv")
sim2 <- as.data.frame(sim2[,2:(length(sim2[1,] ) )] )

sim3 <- read_csv("O:/Test_R/R/2017_juin_HPV_cost-effectiveness/sim_vacc_075_younggirls.csv")
sim3 <- as.data.frame(sim3[,2:(length(sim3[1,] ) )] )

dataList <- list(sim1 = sim1, sim2=sim2, sim3=sim3)

# dataList[[1]][,2]

####### PLOT 6 (test plot)
### total cancer cases, not per age
####### 
plot(NA, xlim=timepl, ylim=c(0,0.0001), 
     xlab= "years",
     ylab="Overall HPV-16 cancers",
     cex.lab=2, frame.plot=FALSE)

for(i in 1:3){
  lines(1:max(time+1),   (dataList[[i]][,1+ ( (1-1)*(length(comp) * 4)  + 5)]+ dataList[[i]][,1+ (  (1-1)*(length(comp) * 4)  + ((length(comp) * 2)+5) ) ] + 
                            dataList[[i]][,1+ ( (1-1)*(length(comp) * 4)  + 7)] + dataList[[i]][,1+ (  (1-1)*(length(comp) * 4)  + ((length(comp) * 2)+7) ) ] +
                            dataList[[i]][,1+ ( (1-1)*(length(comp) * 4)  + 9)] + dataList[[i]][,1+ (  (1-1)*(length(comp) * 4)  + ((length(comp) * 2)+9) ) ] +
                            
                            dataList[[i]][,1+ ( (2-1)*(length(comp) * 4)  + 5)] + dataList[[i]][,1+ (  (2-1)*(length(comp) * 4)  + ((length(comp) * 2)+5) ) ] +
                            dataList[[i]][,1+ ( (2-1)*(length(comp) * 4)  + 7) ] + dataList[[i]][,1+ (  (2-1)*(length(comp) * 4)  + ((length(comp) * 2)+7) ) ] +
                            dataList[[i]][,1+ ( (2-1)*(length(comp) * 4)  + 9)] + dataList[[i]][,1+ (  (2-1)*(length(comp) * 4)  + ((length(comp) * 2)+9) ) ] +
                            
                            dataList[[i]][,1+ ( (3-1)*(length(comp) * 4)  + 5)] + dataList[[i]][,1+ (  (3-1)*(length(comp) * 4)  + ((length(comp) * 2)+5) ) ] +
                            dataList[[i]][,1+ ( (3-1)*(length(comp) * 4)  + 7)] + dataList[[i]][,1+ (  (3-1)*(length(comp) * 4)  + ((length(comp) * 2)+7) ) ] +
                            dataList[[i]][,1+ ( (3-1)*(length(comp) * 4)  + 9)] + dataList[[i]][,1+ (  (3-1)*(length(comp) * 4)  + ((length(comp) * 2)+9) ) ] +
                            
                            dataList[[i]][,1+ ( (4-1)*(length(comp) * 4)  + 5)] + dataList[[i]][,1+ (  (4-1)*(length(comp) * 4)  + ((length(comp) * 2)+5) ) ] +
                            dataList[[i]][,1+ ( (4-1)*(length(comp) * 4)  + 7) ] + dataList[[i]][,1+ (  (4-1)*(length(comp) * 4)  + ((length(comp) * 2)+7) ) ] +
                            dataList[[i]][,1+ ( (4-1)*(length(comp) * 4)  + 9) ] + dataList[[i]][,1+ (  (4-1)*(length(comp) * 4)  + ((length(comp) * 2)+9) ) ] )  / sr,  lty=1, col=col1[i], lwd=2)
}




####### PLOT 7 ICER
### total cancer cases, not per age
####### 
dim(totcostqualy) #time, vacc strat and cost/qaly
colce <- colorRampPalette(c("grey20", "lightgrey"))(3)

plot(NA, xlim=timepl, ylim=c(min(totcostqualy[,,1]),max(totcostqualy[,,1])), 
     xlab= "years",
     ylab="Costs for specific type",
     cex.lab=2, frame.plot=FALSE)

for(i in 1: length(vaccstrat)){
  lines(time, totcostqualy[,i,1] , lty=1, col=colce[i], lwd=c(1,1,1,2,2,2) )
}

for(i in 1: 3){
  lines(time, totcostqualy[,i,1] , lty=1, col=colce[i], lwd=c(1,1,1,2,2,2) )
}


plot(NA, xlim=timepl, ylim=c(min(totcostqualy[,,2]),max(totcostqualy[,,2])), 
     xlab= "years",
     ylab="Costs for specific type",
     cex.lab=2, frame.plot=FALSE)

for(i in 1: length(vaccstrat)){
  lines(time, totcostqualy[,i,2] , lty=1, col=colce[i], lwd=c(1,1,1,2,2,2) )
}

icerut <- (totcostqualy[,1:6,1] - totcostqualy[,1,1]) / (totcostqualy[,1:6,2] - totcostqualy[,1,2])

timepl1 <- c(0,30)

plot(NA, xlim=timepl1, ylim=c (min(icerut[5:201,c(2,3,5,6)]) ,max(icerut[5:201,c(2,3,5,6)]) ), 
     xlab= "years",
     ylab="Costs for specific type",
     cex.lab=2, frame.plot=FALSE)

for(i in 1: ( length(vaccstrat))  ){
  lines(time, icerut[,i] , lty=1, col=colce[i], lwd=c(1,1,1,2,2,2) )
}


min(totcostqualy[,,1])
max(totcostqualy[,,1])

dim(totcostqualy)
plot( time, totcostqualy[,4,1] )


lines(time, totcostqualy[,5,1])


lines( time, totcostqualy[,2,1] )
plot( time, totcostqualy[,3,1] )
lines( time, totcostqualy[,6,1] )



