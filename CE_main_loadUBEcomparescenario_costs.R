# MR, 13.07.2018 
########### For loading results of UBELIX cluster scripts: load and show COSTS from CE_UBE_comparescenario.R


# copied from R_ce folder: CE_loadandshowcosts.R

# Load Final data from models for specific outcomes and vaccination strategies but for all different hpv strains : COSTS aspects

#directories
dirl <- "tab_params/" #where to load parameter data from
plotdir <- "plots/" #where to save plots
loaddir <- "outp/m11/p" #where the UBE data was stored

### Load optimbeta from ubelix
source("CE_main_model_parameters.R")

time <- seq(0, 200, 1)
vaccstrat <- c("V4v_05","V4v_06" , "V4v_07", "V9v_05","V9v_06" , "V9v_07" )
diffv <- c(0.5,0.6,0.7)
y <-  which(dis_outcome=="cervix")





############## LOAD dataset
dfinal <- array(0, dim=c(length(time),length(col_nam),length(vaccstrat), length(hpvtype)), dimnames = list(time,col_nam, vaccstrat, hpvtype) )  
for(i in 1:13){
  load( paste(loaddir,i,".RData", sep=""))
  dfinal[,,,i] <- dperHPVt
}
dim(dfinal)
##############

#### add costs aspect to model calculated without the costs
source("CE_main_model_costsparameters.R")

#### function which groups previous compartements in new compartements (SIR together, disease stages are added for all types, and vaccination coverages)
########## Function that loads the costs for the different vaccination strategies or parameters, and the ICER

#give yp which is the weight to be applied on the incoming population in order to scale it to the European or Swiss population structure
popmod <- 8.4 # population size end of 2016

# in 2016, from the total population residing in Switzerland (all nationality, both sexes, permanently residing population)
# in the age group 11-14 years (which is the vaccination target age) the population sizes are 11y:81149, 12y:81083, 13y: 79781, 14y:81032, mean: 80761.25
yp <- 80761.25 / (g[length(agegroups)]* sum( as.vector(N[,,length(agegroups) ]) ) * 1  * popmod*1e6 )
g[length(agegroups)]* sum( as.vector(N[,,length(agegroups) ]) ) * yp  * popmod*1e6
#check with ESP
(agw * popmod*1e6)[1] / 15 

#load function script:
source("CE_main_loadUBEcomparescenario_costs_functions.R") 

#the function is called: fcostsvaccstrat_varyp, needs value for C9v and gives, as a list: list(cqoutcome,icer, icer2, icer100 )
modpop <- popmod/2 * 1e6
#calculates also costs and qaly, for overall population (weighted according to ESP)
diffcosts <- c("pap_costs", "CIN2CIN3_costs", "cancer_costs", "vacc_costs", "AE_costs", "AE_QALY", "TOT_QALY", "TOT_costs")

########### new compartements: Groups SIR together, and disease stages remain and will be added for each types, (vaccination stays)
#a pre-function gives dfinal2, which is a data frame that groups the compartements in new compartements (SIR together, disease stages are added for all types, and vaccination coverages)
newcomp
col_nam2
dim(dfinal2)

############################################# Check if cancer incidence is correct
# in total, in women:
ncomb <- length(newcomp)

totmodinc <- array(0, dim=c(length(time), length(vaccstrat), 1),  dimnames = list(time, vaccstrat , "incidence")  )
tempin <- c()
for(j in 1:length(time) ){
  for( i in 1:length(vaccstrat)){
    for(k in 0:(length(agegroups)-1) ){
      tempin[k+1] <- sum(dfinal2[j, 1+c(k*ncomb*4+4,                k*ncomb*4+6 ,         k*ncomb*4+8, 
                                        k*ncomb*4+ncomb*2+4, k*ncomb*4+ncomb*2+6, k*ncomb*4+ncomb*2+8) ,i]   * rep( z[,y], 2) / (sr*agr[1+k])  * agw[2+k])
    }
    totmodinc[j,i,] <- sum(tempin  )
  }}

totmodinc[1,1,]
totmodinc[1,1,] * popmod*1e6*sr
totmodinc[1,1,] * 1e5

######## seems to be ok because it gives me the same incidence as before, however, we can't be sure for the SIR yet.


#look, after 100y vaccination start, the ICER of different vaccination strategies (compared to baseline/current  situation)
# integrate over defined time period (area under curve)
#for the vaccination cohort life time (15-84y.o) i.e 69 years
endt <- 102
startt <-2
diffcosts9v <- seq(from=200, to=500, by=10)
icerpricevar <- array(0, dim=c(length(diffcosts9v), length(vaccstrat), 2),  dimnames = list(diffcosts9v, vaccstrat , c("compared_to_baseline", "compared_to_novacc") )  )


for(k in 1:length(diffcosts9v)){
  
  tempd <- fcostsvaccstrat_varyp(C4v, diffcosts9v[k])
  # tempd <- fcostsvaccstrat_varyp(diffcosts9v[k], diffcosts9v[k])
  
  oCosts_b <- sum(tempd[startt:endt,8,1]) # compared with current vaccination strategy as baseline (50% V4v) 
  oQALY_b <- sum(tempd[startt:endt,7,1])
  
  oCosts_nv <- sum(tempd[startt,1:3,1])*(endt-startt) #compared to no vaccination (starting values)
  oQALY_nv <- (tempd[startt,7,1]+ tempd[startt,6,1])*(endt-startt) #total but add the QALY losts ba vacc AE
  
  iCosts <- apply(tempd[,8,], 2, function(x) sum(x[startt:endt]) )
  iQALY <- apply(tempd[,7,], 2, function(x) sum(x[startt:endt] ) )
  
  icerpricevar[k,,1] <- (iCosts - oCosts_b) / (iQALY - oQALY_b)
  icerpricevar[k,,2] <- (iCosts - oCosts_nv) / (iQALY - oQALY_nv)
}

tempd[10,,] 

icerpricevar[,,1]-0


col3 <- colorRampPalette(c("darkblue", "lightblue"))(3)
col4 <- colorRampPalette(c("darkgreen", "lightgreen"))(3)
col1 <- c( col3, col4)
maintles <- c("ICER compared with current situation", "ICER compared with no vaccination")

# pdf(paste(plotdir, "CC_ICER_100y_diffvaccprices_100y", Sys.Date(), ".pdf", sep=""), width=13, height=6)
pdf(paste(plotdir, "CC_ICER_100y_diffvaccprices_fivedC4v", endt-1, Sys.Date(), ".pdf", sep=""), width=13, height=6)
par(mfrow=c(1,2))
par(oma = c(4, 3, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(2, 3, 3, 1)) # make the plots be closer together

for(j in 1:2){
  # par(mar = c(4, 6, 3, 3))
  plot(NA, xlim=c(min(diffcosts9v),max(diffcosts9v) ), ylim=c(min(icerpricevar[,,j], na.rm=TRUE),max(icerpricevar[,,j], na.rm=TRUE) ), 
       xlab= "",
       ylab="",
       cex.lab=2, frame.plot=FALSE, main=maintles[j])
  
  for(i in 1:length(vaccstrat)){
    lines(diffcosts9v, icerpricevar[,i,j], col=col1[i], lwd=2)
  }
  abline(h=0, lty=2, lwd=0.5)
  if(j==1){
    legend(200, max(icerpricevar, na.rm=TRUE) , legend=c( paste("V4v:", (diffv*100), "%", c("(baseline)","", ""), sep="") ,
                                                          paste("V9v:", (diffv*100), "%", sep="") ), 
           pch=16, col=col1, lwd=2, bty="n", title="Vaccine type and coverage" )
  }}
################## Overall x and y labels
mtext('vaccination cost (CHF)', side = 1, outer = TRUE, line = 1, cex=1.5)
mtext('ICER (CHF/QALY)', side = 2, outer = TRUE, line = 1, cex=1.5)

dev.off()



cqoutcome <-  fcostsvaccstrat_varyp(C4v, 310 )



########################
########### -----------------------------------------    plot different costs over time
modpop <- popmod/2  #population modelled here is Swiss women pop in age of being screened or having cancer and less thant 85

plott <- 60

pdf(paste(plotdir, "cdiffcosts", Sys.Date(), ".pdf", sep=""), width=15, height=6)
par(mfrow=c(1,2))
par(oma = c(3, 3, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together

par(mar=c(5,5,1,2))
########################### Cervical cancer
plot(NA, xlim=c(0,plott ), ylim=c(min(cqoutcome[,3,])*modpop ,max(cqoutcome[,3,]*modpop)), 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, main="cervical cancer")

for(i in 1:length(vaccstrat)){
  lines(time, cqoutcome[,3,i]*modpop, col=col1[i], lwd=2)
}

legend(35, max(cqoutcome[,3,]*modpop, na.rm=TRUE) , legend=c( paste("V4v:", (diffv*100), "%", c("(baseline)","", ""), sep="") ,
                                                              paste("V9v:", (diffv*100), "%", sep="") ), 
       pch=16, col=col1, lwd=2, bty="n", title="Vaccine type and coverage" )

################################ CIN
# par(mar=c(5,7,7,2))
plot(NA, xlim=c(0,plott ), ylim=c(min(cqoutcome[,2,])*modpop ,max(cqoutcome[,2,]*modpop)), 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, main="CIN2 and CIN3")

for(i in 1:length(vaccstrat)){
  lines(time, cqoutcome[,2,i]*modpop, col=col1[i], lwd=2)
}

mtext('years after vaccination onset', side = 1, outer = TRUE, line = -1, cex=1.5)
mtext('costs in mio. CHF', side = 2, outer = TRUE, line = -1, cex=1.5)

dev.off()

# legend(60, max(cqoutcome[,2,]*modpop, na.rm=TRUE) , legend=c( paste("V4v:", (diffv*100), "%", c("(baseline)","", ""), sep="") ,
#                                                               paste("V9v:", (diffv*100), "%", sep="") ), 
#        pch=16, col=col1, lwd=2, bty="n", title="Vaccine type and coverage" )
################################ HPV screening
# par(mar=c(5,7,7,2))
plot(NA, xlim=c(0,plott ), ylim=c(min(cqoutcome[,1,])*modpop ,max(cqoutcome[,1,]*modpop)), 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, main="HPV screening")

for(i in 1:length(vaccstrat)){
  lines(time, cqoutcome[,1,i]*modpop, col=col1[i], lwd=2)
}
################## Overall x and y labels
# mtext('years after vaccination onset', side = 1, outer = TRUE, line = -1, cex=1.5)
# mtext('costs in mio. CHF', side = 2, outer = TRUE, line = -1, cex=1.5)
# 
# dev.off()

#vacc costs
plot(NA, xlim=c(0,200 ), ylim=c(min(cqoutcome[,4,])*modpop ,max(cqoutcome[,4,]*modpop)), 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, main="HPV vaccination")

for(i in 1:length(vaccstrat)){
  lines(time, cqoutcome[,4,i]*modpop, col=col1[i], lwd=2)
}

####################### ---------------------------------costs and qualy over time
########### ----  plot different costs over time
# modpop <- popmod/2  #population modelled here is Swiss women pop
plott <- 100

pdf(paste(plotdir, "totcostsandqaly", Sys.Date(), ".pdf", sep=""), width=12, height=6)
par(mfrow=c(1,2))
par(oma = c(3, 3, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together

par(mar=c(5,5,3,2))
########################### Total costs
plot(NA, xlim=c(0,plott ), ylim=c(min(cqoutcome[,8,])*modpop ,max(cqoutcome[,8,]*modpop)), 
     xlab= "",
     ylab="CHF (mio.)",
     cex.lab=2, frame.plot=FALSE, main="total costs")

for(i in 1:length(vaccstrat)){
  lines(time, cqoutcome[,8,i]*modpop, col=col1[i], lwd=2)
}

legend(30, max(cqoutcome[,8,]*modpop, na.rm=TRUE) , legend=c( paste("V4v:", (diffv*100), "%", c("(baseline)","", ""), sep="") ,
                                                              paste("V9v:", (diffv*100), "%", sep="") ), 
       pch=16, col=col1, lwd=2, bty="n", title="Vaccine type and coverage" )

################################ Total QALYS
# QALY gained with vaccination
qalydiff <- t(apply(cqoutcome[,7,],1, function(x) x*modpop*1e6 -cqoutcome[1,7,1]*modpop*1e6  ) )



par(mar=c(5,5,3,2))
plot(NA, xlim=c(0,plott ), ylim=c(min(qalydiff) ,max(qalydiff)), 
     xlab= "",
     ylab="QALY",
     cex.lab=2, frame.plot=FALSE, main="quality of life adjusted years gained \n compared with no vaccination")

for(i in 1:length(vaccstrat)){
  lines(time, qalydiff[,i], col=col1[i], lwd=2)
}

################## Overall x and y labels
mtext('years after vaccination onset', side = 1, outer = TRUE, line = -1, cex=1.5)
# mtext('costs in mio. CHF', side = 2, outer = TRUE, line = -1, cex=1.5)

dev.off()




