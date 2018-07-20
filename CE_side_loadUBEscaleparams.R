# MR, 13.07.2018 
########### For loading results of UBELIX cluster scripts: scale beta and pr3c parameter to cervical cancer and HPV prevalence data, CE_UBE_scale_betaORpr3c.R


# copied from R_ce folder: CE_load_varybeta.R

#load varybeta from ubelix results, for varying beta in a firts step and then vary pr3c

#date
todaysdate <- gsub("-", "_", Sys.Date() )

#CE_UBE_scale_betaORpr3c
sex <- c("f", "m")
riskgroups <- c("l", "h")
hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "other_HR")
time <- seq(0, 200, 1)

y <-  which(dis_outcome=="cervix")

vbeta <- seq(0.55,1,0.05) #has to be the same as the one used in the UBE script

#directory where the results from UBE file are saved ############# CHOOSE DIRECTORY TO LOAD DATA FROM; CHANGES FOR BETA OR PR3C
loaddir <- "param_scaling/beta/m3/p"
loaddir <- "param_scaling/pr3c/m6/p"

#directory where to save optimised params
svdir <- "tab_params/"

#directory where to save table used for params
svdirtab <- "O:/Cost_effectiveness project/Preparation_drafts/tab/"

#where to save plots
plotdir <- "plots/"

######################################  LOAD results from CE_UBE_scale_betaORpr3c: FOR BETA

#first take only for HPV-16, and look at the output
i <- 1
load( paste(loaddir,i,".RData", sep=""))

tosave


plot(tosave[[1]]~ vbeta)
abline(tosave[[2]])

abline( l)
# (- tosave[[2]]$coefficients[1] + prev_pt_2024[10] ) / tosave[[2]]$coefficients[2] 

####### Load for all HPV types
betap_ub <- array(0, dim=c(length(sex),length(hpvtype) ), dimnames = list(sex,hpvtype) )

for(i in (1:length(hpvtype) ) ) {
  load( paste(loaddir,i,".RData", sep=""))
  betap_ub[,i] <- tosave[[5]]
}

# write optimised beta csv.
write.csv(betap_ub, paste(svdir, "optimised_beta", todaysdate,  ".csv", sep="" ) )


#### more information about the optimisation 
# betap_proof <- array(0, dim=c(length(hpvtype), length(vbeta)+3), dimnames = list(hpvtype, c(vbeta, c("min", "opt", "line") )) ) 
# 
# for(i in 1:(length(hpvtype) -3) ){
#   load( paste(loaddir,i,".RData", sep=""))
#   betap_proof[c(1:7,9:10,13)[i],] <- unlist(tosave[c(1,4,5,3)])
# }
# betap_proof
lincorr <- round(tosave[[1]][3],2) ==  round(tosave[[2]]$coefficients[2] * vbeta[3] + tosave[[2]]$coefficients[1],2)




############################################# Vary pr3c based on previous beta calibration
ccpt1  <- read.csv(paste(svdir, "cc_inc_pertype4e6.csv" , sep=""))
ccpt <- as.data.frame(ccpt1[,2])
row.names(ccpt) <- hpvtype

loaddir

# i <- 1
# load( paste(loaddir,i,".RData", sep=""))
# # if(i==2){
# #   vpr3c <- seq(0.4,1,0.1)
# # } else {
# #   vpr3c <- seq(0,0.15,0.05) } #has to be the same as the one used in UBE
# tosave

pdf(paste(plotdir, "pr3cvary", Sys.Date(), ".pdf", sep=""), width=9, height=9)
par(mfrow=c(5,2))
par(oma = c(4, 4, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together

for(i in (1:10)){
  wt <- c(1:7,9:10,13)[i]
  load( paste(loaddir,wt,".RData", sep=""))
  vpr3c <- vptprf(wt)
  plot(NA, xlim=c(min(vpr3c),max(vpr3c) ), ylim=c(0,max(tosave[[1]])) ,
       # plot(NA, xlim=c(min(vpr3c),max(vpr3c) ), ylim=c(0,4) , 
       xlab="",
       ylab="",
       cex.lab=2, frame.plot=FALSE,  main=hpvtype[wt])
  points(tosave[[1]]~ vpr3c, main=hpvtype[wt])
  lines(tosave[[1]]~ vpr3c)
  abline(h=ccpt$`ccpt1`[wt], lty=2, lwd=0.5)
  mtext('pr3c', side = 1, outer = TRUE, line = 2)
  mtext('cevical cancer incidence', side = 2, outer = TRUE, line = 2)
  
}

dev.off()

plot(tosave[[1]]~ vpr3c)
lines(tosave[[1]]~ vpr3c)

pr3c_ub <- array(NA, dim=c(length(hpvtype), 30, 2 ), dimnames = list(hpvtype, 1:30, c( "cc_inc", "vpr3c") ) )
opt_pr3c <- rep(NA, length(hpvtype))

# pr3c_ub <- rep(0,13)
# names(pr3c_ub) <- hpvtype

for(i in 1:(length(hpvtype)-2 ) ) {
  wt <- c(1:7,9:10,13)[i]
  wt <- c(1:10,13)[i]
  load( paste(loaddir,wt,".RData", sep=""))
  vpr3c <- vptprf(wt)
  pr3c_ub[wt, 1:length(vpr3c) ,1] <- tosave[[1]]
  pr3c_ub[wt, 1:length(vpr3c) ,2] <- vpr3c
  clos_smal <- which.min((ccpt[wt,] - pr3c_ub[wt, , 1])[(ccpt[wt,] - pr3c_ub[wt, , 1])>0])  #closest point to the given CC_inc which is smaller than CC_inc
  lr <- lm(pr3c_ub[wt, c( clos_smal, pr3c_ub+1) ,1]~ pr3c_ub[wt, c( clos_smal, pr3c_ub+1) ,2])  #linear regression between the two closest points estimates
  opt_pr3c[wt] <- (-lr$coefficients[1]+ccpt[wt,] ) /lr$coefficients[2] 
}
names(opt_pr3c) <- hpvtype

opt_pr3c[which(is.na(opt_pr3c))] <- 0

todaysdate <- gsub("-", "_", Sys.Date() )
write.csv(opt_pr3c, paste(svdir,"optimised_pr3c", todaysdate,  ".csv", sep="" ) )

i <- 1
load( paste(loaddir,i,".RData", sep=""))
tosave

pr3c_ub
###################################### plot beta variation
tbeta <- rep(NA, length(hpvtype))

pdf(paste(plotdir, "betavariationpertype", Sys.Date(), ".pdf", sep=""), width=9, height=9)
par(mfrow=c(5,2))
par(oma = c(4, 4, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together

for(i in 1:10){
  wt <- c(1:7,9:10,13)[i]
  load( paste(loaddir,wt,".RData", sep=""))
  
  vbeta <- seq(0.55,1,0.05)
  
  plot(NA, xlim=c(0.5,max(vbeta) ), ylim=c(0,max(tosave[[1]])) ,
       # plot(NA, xlim=c(min(vpr3c),max(vpr3c) ), ylim=c(0,4) , 
       xlab="",
       ylab="",
       cex.lab=2, frame.plot=FALSE,  main=hpvtype[wt])
  points(tosave[[1]]~ vbeta, main=hpvtype[wt])
  lines(tosave[[1]]~ vbeta)
  abline(h=prev_pt_2024[wt], lty=2, lwd=0.5)
  mtext('beta', side = 1, outer = TRUE, line = 2)
  mtext('HPV prevalence', side = 2, outer = TRUE, line = 2)
  clos_smal <- which.min((prev_pt_2024[wt] - tosave[[1]] )[(prev_pt_2024[wt] - tosave[[1]])>0])  #closest point to the given prev which is smaller than prev
  lr <- lm(tosave[[1]][c( clos_smal, clos_smal+1) ]~ vbeta[ c( clos_smal, clos_smal+1) ])  #linear regression between the two closest points estimates
  abline(lr, col="red")
  tbeta[wt] <- (-lr$coefficients[1]+prev_pt_2024[wt] ) /lr$coefficients[2]
}

dev.off()

#Write the parameter set used for the beta optimisation:
#re-write the new (optimised) beta in table P_per_hpvtype

load( paste(loaddir,1,".RData", sep=""))
pr3c <- as.vector( read.csv(paste(svdir, "optimised_pr3c2018_07_17.csv" , sep=""))[,2])
betap_ub <- read.csv(paste(svdir, "optimised_beta2018_07_16.csv" , sep=""))
pr2c <- 0.2*pr3c

write.table(tosave[[2+4]], paste(svdirtab,"P_general_",Sys.Date(),".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
# tosave[[3+4]][1,2:14] <- round(betap_ub[1,], 3)
write.table(tosave[[3+4]][,2:14], paste(svdirtab, "P_perhtype_",Sys.Date(),".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tosave[[4+4]], paste(svdirtab,"P_perdisease_",Sys.Date(),".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tosave[[5+4]], paste(svdirtab,"P_economic_",Sys.Date(),".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tosave[[6+4]], paste(svdirtab, "P_c_inc_",Sys.Date(),".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)

