# MR, 13.07.2018 
########### For UBELIX cluster scripts: scale beta and pr3c parameter to cervical cancer and HPV prevalence data


# copied from R_ce folder: CE_varybeta_ube2.R

############################## For the UBELIX .job job
# Args from array in bash
args=(commandArgs(TRUE))
i=as.numeric(unlist(args))
##############################

#Decide wether run over all 10 HR hpv types which have a CC incidence >0, or on all (for beta for all, for pr3c only for the 10)
# wt <- c(1:7,9:10,13)[i]
wt <- c(1:13)[i]

hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

print(hpvtype[wt])

#define directory to take the external data from (in ubelix it is different)
# dirl <- "tab_params/"
dirl <- ""

#define date of betaparams to be used 
dobeta <- "2018_07_23"
dopr3c <- "2018_07_23"

################################ Choose which parameter to scale first (beta or pr3c). First the beta has to be scaled on the HPV prev per type, then 
# the pr3c on the cervical cancer incidence per type. Thus the scaling or optimisation of the two parameters has to be done in two steps.

ptv <- c("vbeta", "vpr3c")[2] # 1 for vaying beta, 2 for varying pr3c

vbeta <- seq(0.55,1,0.05) #first step: says the range in which beta can be

#second step (with different value for hpv types) #range of pr3c (differs for types, to be more precise)

# vpr3c <- seq(0.0, 0.02,0.001)

vptprf <- function(wt){
  if(wt==which(hpvtype=="hpv16") ){
  vpr3c <- seq(0.033,0.038,0.001)
} else if(wt==which(hpvtype=="hpv18")  ){
  vpr3c <- seq(0.04,0.044,0.001)
} else if(wt==which(hpvtype=="otherHR")  ){
  vpr3c <- seq(0.03,0.035,0.001)
} else if(wt==which(hpvtype=="hpv45")  ){
  vpr3c <- seq(0.06,0.065,0.001)
} else if(wt==which(hpvtype=="hpv35")  ){
  vpr3c <- seq(0.055,0.060,0.001)
}else {vpr3c <- seq(0.007,0.023,0.002)}
  return(vpr3c)
}


# vptprf <- function(wt){
#   if(wt==which(hpvtype=="hpv16") ){
#     vpr3c <- seq(0.01,0.05,0.005)
#   } else if(wt==which(hpvtype=="hpv18")  ){
#     vpr3c <- seq(0.02,0.12,0.01)
#   }else {vpr3c <- seq(0.0,0.05,0.005)}
#   return(vpr3c)
# }

# vptprf <- function(wt){
#   if(wt==which(hpvtype=="hpv35") |wt==which(hpvtype=="hpv45") |wt==which(hpvtype=="hpv16")  ){
#     vpr3c <- seq(0.02,0.1,0.005)
#    }else {vpr3c <- seq(0.0,0.05,0.005) }
#   return(vpr3c)
# }

vpr3c <- vptprf(wt)
# vpr3c <- seq(0.0,0.05,0.005)


## requires paramters 
source("CE_main_model_parameters.R")

#and the model script
source("CE_main_model.R")


#use optimised data for beta if varying pr3c
betap <- array(0, dim=c(length(sex),length(hpvtype) ), dimnames = list(sex,hpvtype) )
txtti <-  paste(dirl, "optimised_beta", dobeta,".csv" , sep="")

betap[1,] <- as.numeric(read.csv(txtti)[1,2:14] )
betap[2,] <- as.numeric(read.csv(txtti)[2,2:14] )

#load cervical cancer incidence per type 
ccpt1  <- read.csv( paste(dirl,  "cc_inc_pertype4e6.csv", sep="") )
ccpt <- as.data.frame(ccpt1[,2])

#load type ratio for normal cytology
nGUAN_typerationormal <- read.csv(paste(dirl, "nGUAN_typerationormal.csv", sep=""))[,2]
#calculate type specific prevalence based on Natsal'3 estimates of HR HPV in 15-19 years olds (from CE_compl_hpvprevnatsal.R)
#attribute this HR prevalence to each HPV type according to the type ratio from Guan et al. and assumes that this corresponds to the 
#HPV HR prevalence in 20-24 year old Swiss women
prev_pt_2024 <- nGUAN_typerationormal * 0.2435849
hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "other_HR")
names(prev_pt_2024) <- hpvtype


#has to be on if vbeta
pr3c <- read.csv(paste(dirl, "optimised_pr3c", dopr3c ,".csv",sep=""))[,2]
# pr3c[c(8,11,12)] <- pr3c[3]
pr2c <- pr3c*0.2


# w <- which(hpvtype=="hpv16")
w <- wt
y <-  which(dis_outcome=="cervix")
prev_ipag2<- c()
inc_cc <- c()
print(wt)

############################################ Calculates, for each individual HPV type at time, model outputs when varying either beta or pr3c
# pr2c is considered to be 0.2*pr3c. Gives, the HPV prevalence for beta and cervical cancer incidence for pr3c, the linear regression between the two closest points estimates,
#says if the rounded resulting estimates is the same as the data point estimate, indicates the lowest value possible (= 0), and the correct estimate for the data point to be fitted on.


########################### for varying beta
if(ptv == "vbeta"){
  for(l in 1:length(vbeta)){
    params <- c(w= wt, y = y, beta_fm=vbeta[l], beta_mf=vbeta[l] ) 
    RS <- runsteady(y = I0, fun = HPV_1,  parms = params, times = c(0, 1e5)  ) 
    init <- RS$y
    spw_esp <- array(0, dim=c(  length(agegroups), length(comp)  ),
                     dimnames =list(  agegroups, comp) )
    
    
    for(i in 1: length(agegroups)) {
      for(j in 1: length(comp) ){
        spw_esp[i,j] <- ( (init[(i-1)*(length(comp)*4) + j] + init[(i-1)*(length(comp)*4) + (length(comp)*2)+j] ) / (sr*agr[i]) )* agw[1+i]
      } }
    
    prev_ipag2[l] <- spw_esp[2,2] / agw[3]
  } 
  clos_smal <- which.min((prev_pt_2024[wt] - prev_ipag2 )[(prev_pt_2024[wt] - prev_ipag2)>0])  #closest point to the given prev which is smaller than prev
  lr <- lm(prev_ipag2[c( clos_smal, clos_smal+1) ]~ vbeta[ c( clos_smal, clos_smal+1) ])  #linear regression between the two closest points estimates
  lincorr <- round(prev_ipag2[3],2) ==  round(lr$coefficients[2] * vbeta[3] + lr$coefficients[1],2)
  minbeta <- -lr$coefficients[1] /lr$coefficients[2]
  betadp <- (-lr$coefficients[1]+prev_pt_2024[wt] ) /lr$coefficients[2]
  inc <- prev_ipag2
  
} else {
  ########################### for varying pr3c and pr2c  
  for(l in 1:length(vpr3c)){
    pr3c[wt] =vpr3c[l]
    
    pr2c[wt] <- pr3c[wt]*0.2
    params <- c(w= wt, y = y, beta_fm = betap[1,wt], beta_mf= betap[2,wt]  ) 
    RS <- runsteady(y = I0, fun = HPV_1,  parms = params, times = c(0, 1e5)  ) 
    init <- RS$y
    spw_esp <- array(0, dim=c(  length(agegroups), length(comp)  ),
                     dimnames =list(  agegroups, comp) )
    
    
    for(i in 1: length(agegroups)) {
      for(j in 1: length(comp) ){
        spw_esp[i,j] <- ( (init[(i-1)*(length(comp)*4) + j] + init[(i-1)*(length(comp)*4) + (length(comp)*2)+j] ) / (sr*agr[i]) )* agw[1+i]
      } }
    inc_cc[l] <- ( ( z[1,y]*sum(spw_esp[,5]) + z[2,y]*sum(spw_esp[,7]) + z[3,y]*sum(spw_esp[,9] ) ) * 1e5 )
  }
  # lr <- lm(inc_cc[4:5]~vpr3c[4:5]) 
  dp <- ccpt$`ccpt1[, 2]`[wt]
  clos_smal <- which.min((dp - inc_cc )[(dp - inc_cc)>0])  #closest point to the given prev which is smaller than prev
  lr <- lm(inc_cc[c( clos_smal, clos_smal+1) ]~ vpr3c[ c( clos_smal, clos_smal+1) ]) 
  lincorr <- round(inc_cc[3],2) ==  round(lr$coefficients[2] * vpr3c[3] + lr$coefficients[1],2)
  minbeta <- -lr$coefficients[1] /lr$coefficients[2]
  betadp <- (-lr$coefficients[1]+dp ) /lr$coefficients[2] 
  inc <- inc_cc
  
}


print(prev_ipag2)
print(vbeta)
print(ptv)
print(lr)


# tosave <- list(inc, lr,lincorr, minbeta, betadp )
#for HPV-16, it also saves the parameter tables which were used to retrieve optimised parameters

# savingf <- function(wt){
#   if(wt==1) {
#     source("CE_main_model_costsparameters.R")
#     source("CE_side_loadLateXtable.R")
#     tosave <- list(inc, lr,lincorr, minbeta, betadp,tab_generalparams, tab_hpvtparams,tab_hpvdparams, tab_economicev, tab_cinc )
#     return(save(tosave, file="p.RData"))
#   }
#   else
#     tosave <- list(inc, lr,lincorr, minbeta, betadp )
#   return(save(tosave, file="p.RData"))
# }
# savingf(wt)

tosave <- list(inc, lr,lincorr, minbeta, betadp )
save(tosave, file="p.RData")
