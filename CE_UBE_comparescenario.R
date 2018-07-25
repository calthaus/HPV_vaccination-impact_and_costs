# MR, 13.07.2018 
########### For UBELIX cluster scripts: Model CE for different HPV types and vaccination scenarios with ubelix

# copied from R_ce folder: CE_model_scenariwithcost_ube_v2.R

############################## For the UBELIX .job job
# Args from array in bash
args=(commandArgs(TRUE))
i=as.numeric(unlist(args))
##############################
# w <- c(1:7, 9:10, 13)[i]
w <- c(1:13)[i]

#define directory to take the external data from (in ubelix it is different)
# dirl <- "tab_params/"
dirl <- ""

#define date of betaparams to be used 
# dobeta <- "2018_07_16"
# dopr3c <- "2018_07_17"
dobeta <- "2018_07_23"
dopr3c <- "2018_07_23"


######## make initial params for specific type
source("CE_main_model_parameters.R")
#and the model script
source("CE_main_model.R")
#make the sexual mixing matrix
rho <- rho_sra(eps_r, eps_a)

#call some of the parameters
#before running the  model, define hpv type and disease outcome:
# w <- which(hpvtype=="hpv16")
print(w)
y <-  which(dis_outcome=="cervix")

betap <- array(0, dim=c(length(sex),length(hpvtype) ), dimnames = list(sex,hpvtype) )
txttibeta <- paste( dirl, "optimised_beta", dobeta , ".csv", sep="")
txttibpr3c <- paste(dirl, "optimised_pr3c", dopr3c, ".csv", sep="")

betap[1,] <- as.numeric(read.csv(txttibeta)[1,1+(1:length(hpvtype))] )
betap[2,] <- as.numeric(read.csv(txttibeta)[2,1+(1:length(hpvtype))] )

pr3c <- as.numeric(read.csv(txttibpr3c)[,2] )
pr2c[w] <- pr3c[w]*0.2


print(betap[1,w])
params <- c(w= w, y = y, beta_fm = betap[1,w], beta_mf= betap[2,w] )

######################################### STEADY state for all HPV types, before starting preventions intervention (vaccination)

#set vaccination to 0 (here for quadrivalent vaccination, all other already are at 0)
p4v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,
                                                                                         agegroups)  )
p9v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,
                                                                                         agegroups)  )
#this runs it until steady state and gives the initial vector at steady state takes about 2 min
RS <- runsteady(y = I0, fun = HPV_1,  parms = params, times = c(0, 1e5)  )   


#set the new initial values based on steady-state
init <- RS$y


#load initial values for specific hpv type and outcome 
# load("init_t1_o1.RData")
# init_t1_o1

#Set time after vaccine introduction
time <- seq(0, 200, 1)

#define different vaccination scenarios:
diffv <- c(0.5,0.6,0.7)
vaccstrat <- c("V4v_05","V4v_06" , "V4v_07", "V9v_05","V9v_06" , "V9v_07" )

#start with setting vaccination to 0 
p4v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups, agegroups)  )
p9v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,agegroups)  )

dperHPVt <- array(0, dim=c(length(time),length(col_nam),length(vaccstrat)), dimnames = list(time,col_nam,  vaccstrat)  )

# pbr <- txtProgressBar(1,length(diffv), style=3)
print(w)
print(betap[1,w])

# First three vaccination scenario with V4v vaccination
for (j in 1:length(diffv) ){
  p4v[1,,1] <- diffv[j]
  sim <- as.data.frame(ode(init, time, HPV_1, parms=params))
  dperHPVt[,,j] <- as.matrix(sim)
  # setTxtProgressBar(pbr, j)
}

p4v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,
                                                                                         agegroups)  )
p9v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,
                                                                                         agegroups)  )

for (j in 1:length(diffv) ){
  p9v[1,,1] <- diffv[j]
  sim <- as.data.frame(ode(init, time, HPV_1, parms=params))
  dperHPVt[,,3+j] <- as.matrix(sim)
  # setTxtProgressBar(pbr, j)
}


save(dperHPVt, file="p.RData")



