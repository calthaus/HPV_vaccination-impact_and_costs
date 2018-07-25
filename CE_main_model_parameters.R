# MR, 12.07.2018

########### Parameters for the main HPV model (CE_main_model.R) transmission dynamics


# copied from R_ce folder: CE_model_params_nocost_v10.R

# each HPV strain has be run individually


#--------------------------- Parameters ------------------------------


## dimensions **********

sex <- c("f", "m")

riskgroups <- c("l", "h")

agegroups <- c( "15-19", "20-24","25-29", "30-39", "40-64", "65-84") #max 72 for natsal3, and 99 for shs

comp <- c("S", "I", "CIN2", "CIN3", 
          "ULCC", "DLCC", "URCC", "DRCC",  "UDCC", "DDCC", "DCC", "R", "V2v", "V4v", "V9v")

hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

dis_outcome <- c("cervix", "anus", "penile", "oropharynx", "warts")

#directory to take data from
#dirl <- "tab_params/"


# FIXED PARAMS
## population **********

##### Vaccination
#vaccination parameters are set to 0, this can be changed in the core script before running the model, according to prevention strategy choosen

#vaccination coverage rates p[sex,riskgrous,agegroups] for bivalent vaccine
p2v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,  agegroups)  )
# proportion young women vaccinated (same between risk groups) is : p2v[1,,2] <- 0.0

#vaccination coverage rates p[sex,riskgrous,agegroups] for quadrivalent vaccine
p4v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups, agegroups)  )

#vaccination coverage rates p[sex,riskgrous,agegroups] for nonavalent vaccine
p9v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,  agegroups)  )

# array of the same size as vaccination arrays, but with only 0. This then used for non vaccine targeted types
pv0 <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups, agegroups)  )
p2v_p <- pv0
p4v_p <- pv0
p9v_p <- pv0


#population size N[sex,riskgroup,agegroup], the sum of all the population (including deaths) is 1
N <- array(1, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,
                                                                                       agegroups)  )
#sex ratio
sr <- 0.5

N[1,,] <- N[1,,]*sr 
N[2,,] <- N[2,,]*sr 

#risk groups % in high or low sexual activity groups. This is calculated in the CE_compl_partnerchangerate_1_retrieveSwissdata.R and CE_compl_partnerchangerate_2_MLEactivitygroups.R
# script and saved as sbehav_ch3.csv
sexbehavdata <- paste(dirl, "sbehav_ch3.csv", sep="")
sbehav <- as.data.frame(read.csv(sexbehavdata) )

rgr <- c(sbehav[1,2],1-sbehav[1,2]) # % in low or high sexual activity group

N[,1,] <- N[,1,]*rgr[1]
N[,2,] <- N[,2,]*rgr[2]


#makes age groups and adapt size according to number of year spent in each group 
agr1 <- c((19-15) , 24-20, 29-25, 39-30, 64-40,84-65)
agr1 <- agr1 +1
agr <- agr1/ mean(agr1)

agr <- agr* (1/ length(agegroups) )


for( i in 1:length(agegroups) ){
  N[,,i] <- N[,,i]*agr[i]
}

#rate of entering/exiting system, according to the min and max age included
mu <- 1/(84-15)

#rates of changes between agegroups
g <- mu/agr

g <- matrix(g, nrow=1, ncol=length(agegroups))
colnames(g) <- agegroups

#number of contacts C[sex, riskgroups, age]
C <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups, agegroups)  )

# sexbehavdata <- "sbehav_ch3.csv"
sbehav <- as.data.frame(read.csv(sexbehavdata) )
rownames(sbehav) <- sbehav[,1]
sbehav1 <- sbehav[-1,2:3]

for( i in 1: (length(agegroups)) ){
  C[,1,i] <- sbehav1[c( (i*2) -1 , (i*2) ),1]
  C[,2,i] <- sbehav1[c( ( ( (i*2) -1) + (length(agegroups))*2), ( (i*2) + (length(agegroups))*2) ), 1]
}



#sexual mixing between risk groups (eps_r) and age groups (eps_a), assortativity (0= random, 1=assortative), 
# eps_r is based on the litterature, and eps_a was retrieved from information hold in Natsal'3 dataset (CE_ageassortativity.R)
eps_r <- 0.5
eps_a <- 0.5242733

m <- 1

################################################### PER HPV TYPE VARYING PARAMETERS

hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

## ********** infection dynamics **********

#transmissibility probablity
betap <- array(0, dim=c(length(sex),length(hpvtype) ), dimnames = list(sex,hpvtype) )

#scaled to fit to prevalence data. Here is only an example. Has to be called before running the model
# can be adapted to women to men and men to women, however, we assume the same for both sexes
betap[1,1:13] <-c(0.95, 0.93, 0.62, 0.89, 0.66, 0.61, 0.68, NA, 0.63, 0.65,NA, NA,  0.89)
betap[2,1:13] <- betap[1,1:13] 

#Proportion of infections that develop HPV-type specific antibodies #from Durham and for 6 11 and other it is the mean value
alpha <- array(0, dim=c(length(hpvtype) ), dimnames = list(hpvtype) )
alpha[1:13] <- c(56.7, 41.0, 56.6, 40.7, 88.9, 47.5, 37.9, 51.33, 15.4, 57.9, 63.6, 51.33, 51.33 ) #for hpv51, 11 and OHR, the average was used
alpha[1:13] <- alpha[1:13]/100


###################################################  PER DISEASE OUTCOME VARYING

dis_outcome <- c("cervix", "anus", "penile", "oropharynx", "warts")

#rate of progression to different stage: pri1 = from infectious status to CIN2, pri2 = from infectious to CIN3
pri1 <- array(0, dim=c(length(sex),length(hpvtype), length(dis_outcome) ), dimnames = list(sex,hpvtype, dis_outcome) )
pri2 <- array(0, dim=c(length(sex),length(hpvtype), length(dis_outcome) ), dimnames = list(sex,hpvtype, dis_outcome) )

#for cervical cancer, per type (from Durham's params ) #-------------------------------------------------------CHECK OK
# pri1[1,,1] <- c(0.04853086, 0.01812457, 0.02612056, 0.04792516, 0.01332602, 0.01332602, 0.0185395, 0.01332602, 
#                 0.01332602, 0.01332602, 0.01332602, 0.01332602, 0.01332602) #from Jaisamrarn et al 2013 ref. recalculated (similar to Durham but not exactly same)

#proportion values from paper, with CI, (multivariate analysis Hazard Ratio ): hpv16:9.25 (95%CI:6.84-12.51), hpv18: 3.53(2.4-5.27),
# hpv31: 5.09(3.56-7.29), hpv33: 9.14 (6.34-13.18), hpv45: 3.64 (2.08-6.40) other oncogenic types: 2.63 (1.94-3.57)

#from infected directly to CIN3 (from Durham's params ) #-------------------------------------------------------CHECK OK
# pri2[1,,1] <- c(0.020411, 0.005025168, 0.01010135, 0.02564665, 0.0178654, 0.0178654, 0.02276026, 0.0178654, 0.0178654, 0.0178654, 0.0178654, 0.0178654, 0.0178654)

#proportion values from paper, with CI, (probability of develop. CIN3+ at month 24 and 95%CI): 
# hpv16: 0.04 (0.03-0.06), hpv18: 0.01 (0.01-0.02),
# hpv31: 0.02 (0.01-0.03), hpv33: 0.05 (0.03-0.07), 
# hpv45: 0.0445(0.0117-0.1697) other oncogenic types: 0.0351(0.0162-0.0759) # for these two, multivariate analysis hazard ratios were used

#from Durham only
pri1[1,,1] <- c(0.042, 0.018, 0.026, 0.045, 0.011, 0.016, 0.011, 0.011, 0.011, 0.011, 0, 0, 0.011)
pri2[1,,1] <- c(0.015, 0.003, 0.006, 0.016, 0.0022, 0.003, 0.0022, 0.0022,
               0.0022, 0.0022, 0, 0, 0.0022)

#for cervix, 0 for men
pri1[2,,1] <- 0
#for penis, 0 for women
pri1[1,,3] <- 0

#for anus , per type (random number) - to be defined
pri1[,,2] <- c(0.001)
#for penis , per type (random number)  - to be defined
pri1[2,,3] <- c(0.001)
#for oropharynx , per type (random number)  - to be defined
pri1[,,4] <- c(0.001)
#for anogenital warts , per type (random number) in women  - to be defined
pri1[1,,5] <- c(rep(0,10), 0.05, 0.05, 0)
#for anogenital warts , per type (random number) in men     - to be defined
pri1[2,,5] <- c(rep(0,10), 0.05, 0.05, 0)

#for hpv types 6 and 11 (low risk, the progression to cancer is 0)
pri1[,c(11,12),1:4] <- 0

#no progression to warts
pri2[,,5] <- 0

# progression from CIN2 to CIN3, From Durham reference: Moscicki A-B et al. 2010
pr23 <- array(0.05417298, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) )
pr23[c(11,12),1:4] <- 0 # for low risk types for other outcomes than genital warts
pr23[-c(11,12),5] <- 0

# progression from CIN2 to cancer is dependent on progression CIN3 to cancer, which is scaled to cancer incidence (pr2c = 0.2*pr3c), From Campos et al. 2014
# This is just an example, data has to be loaded before running the model
pr2c <- c(0.040, 0.080 ,0.095, 0.005, 0.020 ,0.175,0.105, 0.0, 0.045, 0.040, 0,0, 0.010)*0.2
names(pr2c) <- hpvtype
# progression from CIN3 to cancer is scaled in incidence data. This is just an example, data has to be loaded before running the model
pr3c <- c(0.040, 0.080 ,0.095, 0.005, 0.020 ,0.175,0.105, 0.0, 0.045, 0.040, 0,0, 0.010)
names(pr3c) <- hpvtype

#rate of regression from CIN3 and CIN2 to susceptibles # from Moscicki 2010
pr3s <- array(0, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) )
pr2s <- array(0, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) ) 

pr3s <- array(0.1, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) )
pr3s[c(11,12),2:5] <- 0

pr2s <- array(0.5047092, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) ) #Moscicki 2010
pr2s[1:2,1] <- 0.2669108

#This was assumed to be an unnecessary parameters since we scale beta and pr2c/pr3c. Also, the estimates were not very relieable in my opinon (McCredie et al 2008)
pr32 <- 0

#progression localised to regional and regional to distant, from Campos et al (2014), slightly different from Durham's
prlcrc <- 0.2424 #
prrcdc <- 0.3038

#relative risk of re-infection following clearance, from Castellsagué et al. 2014
tau <- array(0, dim=c(length(hpvtype), length(dis_outcome) ), 
             dimnames = list(  hpvtype, dis_outcome) )

tau[,1:5] <- c(0.64, rep(0.96,12) ) #only known for hpv16 and 18. All other types are considered to be like hpv18
# CI: 0.64 (0.53-0.78), 0.94 (0.75-1.19)

#Vaccine efficacy #hpv16:0.91-0.96, 18:0.91-0.97): from Lu et al 2011. 
#for the other hpv types: 0.967 CI 0.809-0.998) from Joura et al 2015
ve <- rep(0, length(hpvtype))
ve <- c(0.94, 0.95, 0.967, 0.967, 0, 0,0.967,0,0.967, 0.967, 0.967, 0.967, 0)
names(ve) <- hpvtype

#screening rates per age group and disease (only women)
#take % of women who had pap smear in the last year from SHS 2012, retrieved from SHS 2012 through: CE_screeningCH.R, saved under: sceeningrates_2012SHS_V2.csv
s1  <- read.csv( paste(dirl, "sceeningrates_2012SHS_v2.csv", sep="") )
s1 <- as.data.frame(s1[,2:13])

s <- array(0, dim=c(length(sex), length(agegroups), length(dis_outcome) ), 
           dimnames = list( sex, agegroups, dis_outcome) )
#screening only occurs for cervix outcomes and only in women # for the first age group there was no indication, I took half of the second age group
s[1,,1] <-  c( s1[1,6] /2 , s1[1:5,6] )

#transform proportions in rate.
s[1,,1] <- (-log(1-s[1,,1]) )/1


#previously used screening rates:
# s[1,,1] <- c(0.28,0.56,0.66,0.68,0.61,0.38)
# s[1,,1] <- c(0.2780318, 0.5560636, 0.6612978, 0.6818641, 0.6084159, 0.3820578)

#screening sensitivity
# ssens <- 0.96 #0.95-0.98
ssens  <- 0.587 #from Bigras 2005, CI48.6-68.2

#screening specificity
# sspec <- 0.91 #0.90-0.93
sspec <- 0.969 #from Bigras 2005, CI 96.6-97.2

################################################### PER DISEASE STAGE SPECIFIC (for cervical cancer)

#stage specific probablity of diagn.(?), from Campos (2014) Web appendix, table 1
dstages <- c("cin2", "cin3", "local", "regional", "distant")
dstages[3:5]

z <- array(0, dim=c(length(dstages[3:5]), length(dis_outcome) ), 
           dimnames = list( dstages[3:5],  dis_outcome) )
z[,1:5] <- c(0.2106379, 0.9160948,2.302646)


#stage specific probability of cure with treatement, from Durham who cites: Elbasha EH et al. 2007 #----------------------------------CHECK value from modelling paper..
pi <- array(0, dim=c(length(dstages[3:5]), length(dis_outcome) ), 
            dimnames = list( dstages[3:5],  dis_outcome) )
pi[,1:5] <-  c(0.92, 0.55, 0.17)

#mortality of undiagnosed and diagnosed cancers (stage specific), from Campos et al. (2014), slightly different from Durham but very close
um <- c(0.01921,0.11455,0.35685)
dm <- c(0.01080486,0.04327795,0.09154833)


#infection clearance rate (1/gamma = duration of the infectiousness), or disease clearance rate
#gamma[infection clearance, CIN2 clearance, CIN3 clearane]
gamma <- array(0, dim=c(length(sex),length(hpvtype)), dimnames = list(sex,hpvtype) ) 

# per hpv type from Jaisamararn et al (2013). Different from Durham #----------------------------------CHECK OK, median duration, no multivariate analysis
# gamma[1,] <- c(0.49, 0.66, 0.57, 0.66, 0.69,0.69, 0.66, 0.69, 0.69 ,0.69, 0.69, 0.69, 0.69 ) #these are Durhams
gamma[1,] <- c(17.11, 11.84, 13.8, 12.00, 11.77,11.77, 11.48, 11.77, 11.77 ,11.77, 8.26, 8.26, 11.77 )
gamma[1,] <- 1/ ( gamma[1,] /12 )
#for men               
# gamma[2,] <- c(1.08, 1.33, rep(1.41,11 )  )
gamma[2,] <- c(7.7, 6.2, rep(6.1,11 )  )
gamma[2,] <- 1/ ( gamma[2,] /12 )

#immunity duration (1/omega= duration for becoming susceptible again)
######### New from Johnson (2012) paper because I could not retrieve the values in the Durham reference (Syrjänen 2009)#----------------------------------CHECK OK
omega <- array(0, dim=c(length(sex),length(hpvtype)), dimnames = list(sex,hpvtype) ) 

tempo <- c(0.024, 0.017, 0.018, 0.018 , 0.016, 0.021, 0.016, 0.019, 0.027, 0.017,0.020,0.017  ) # these are all available HPV HR values
mean(tempo)
omega[1,] <- c(0.024, 0.017, 0.018, 0.018 , 0.016, 0.021, 0.016, 0.019, 0.027, 0.020,mean(tempo), mean(tempo), mean(tempo)  ) 
omega[2,] <- omega[1,]

##################################### END of parameters for infection dynamics. 

#initial params I0[compartement, sex, riskgroups, agegroups]
I0_1<- array(0, dim=c( length(comp), length(sex), length(riskgroups), length(agegroups)),
             dimnames =list(comp,sex,riskgroups,
                            agegroups))

#set initial infection proportion in the population for the different risk groups and age group
#here I assume it is the same for both age groups and for male and female

# women & men
# Suseptibles is N-I0
I0_1[1,,1,] <-  N[,1,]- 0.01*N[,1,] 
# highrisk group
I0_1[1,,2,] <- N[,2,]- 0.1*N[,2,]


# lowrisk group
I0_1[2,,1,] <- 0.01*N[,1,] 
# highrisk group
I0_1[2,,2,] <- 0.1*N[,2,]


I0 <- I0_1
I0[,1,,1] / (sr* agr[1])


##### Output column names vector
# gives column names to output matrix
compnamef <- function(compartements){
  c_name <- rep(compartements, length(sex)*length(riskgroups)*length(agegroups))
  s_name <- rep( rep(sex, each=length(compartements)) , length(riskgroups)*length(agegroups))
  r_name <- rep( rep(riskgroups, each=length(compartements)*length(sex) ), length(agegroups))
  a_name <- rep(agegroups, each=length(compartements)*length(sex)*length(riskgroups) )
  col_nam <- c("time", paste(c_name, s_name, r_name, a_name, sep="_") )
  return(col_nam)
}

col_nam <- compnamef(comp) 


#age group weights according to the European Standard population (ESP 2013)
#agw <- c( sum(0.01,0.04,0.055,0.055),0.055,0.06,0.06, sum(0.065,0.07), sum(0.07,0.07,0.07,0.065,0.06), sum(0.055,0.05,0.04,0.025), 0.025 )
#names(agw) <- c("0:14",agegroups, "+85")

# from http://ec.europa.eu/eurostat/documents/3859598/5926869/KS-RA-13-028-EN.PDF/e713fa79-1add-44e8-b23d-5e8fa09b3f8f p.31, projections 2011-2020
nav <- c(1078.641, 4373.749, 5410.346, 5252.859,  5410.049, 6066.914,  6711.973,  7023.97, 7135.495, 
         7126.248, 7087.804, 6938.434, 6595.514, 6095.677, 5307.002, 4328.78, 3419.627, 2492.941,1452.548,555.307,136.12 )
names(nav) <- c("0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
                 "55-59", "60-64","65-69", "70-74","75-79", "80-84" , "85-89","90-94", "95+")

nav /sum(nav)

#make a vector with age group weights with our age groups and the younger and older pop
agw <- c( sum(nav[1:4]) ,nav[5:7], sum(nav[8:9]), sum(nav[10:14]), sum(nav[15:18]), sum(nav[19:21]) )
names(agw) <- c("0:14",agegroups, "+85")
agw <- agw/sum(agw)


