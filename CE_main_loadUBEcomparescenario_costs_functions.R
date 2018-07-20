# MR, 13.07.2018 
########### For loading and showing costs: function adding HPV types and calculating costs


# copied from R_ce folder: CE_loadandshowcosts_functions.R


################## Pre functions and codes to load the function which gives the costs and qalys -----------------------------------------------
##################### 1st step is to combine the different compartements from dfinal across the different hpv types
##################### The disease stages from the different types will be added up. Vaccination uptake has to be accounted once, and the rest of the
##################### compartements, except the deaths, will be groupes (SIR) and calculated as 1 - (disease stages + vaccination uptake)
#function which gives the emplacement of specific compartements over age groupes sex risk groups. 
fagmf <- function(comptla, compartements){
  nbcompt <- rep(NA, length(agegroups)*length(riskgroups))
  for(i in 0: (length(agegroups)-1) ){
    nbcompt[ ( i*(length(comptla)*4) +1 ) : ( ( i*(length(comptla)*4) +1 )+ 3 ) ] <- 1+ c( (i*length(compartements)*4 + comptla), #female low
                                                                                           (i*length(compartements)*4 + length(compartements)*1 + comptla), #male low
                                                                                           (i*length(compartements)*4 + length(compartements)*2 + comptla), #female high
                                                                                           (i*length(compartements)*4 + length(compartements)*3 + comptla))  #male high
  }
  return(nbcompt)
}
# fagmf(1, length(newcomp) )

#same only for women
fag <- function(comptla, compartements){
  nbcompt <- rep(NA, length(agegroups)*length(riskgroups))
  for(i in 0: (length(agegroups)-1) ){
    nbcompt[ ( i*(length(comptla)*2) +1 ) : ( ( i*(length(comptla)*2) +1 )+ 1 ) ] <- 1+ c( (i*length(compartements)*4 + comptla),                                                                                        (i*length(compartements)*4 + length(compartements)*2 + comptla)) 
  }
  return(nbcompt)
}

# fag(3, newcomp)
# col_nam2[ c( fag(3, newcomp)) ]

dim(dfinal)
dim(dfinal[1,,1,])
x <- dfinal[100,,1,]
temp <- dfinal[,,1,]
dim(temp)

########### new compartements: Groups SIR together, and disease stages remain and will be added for each types, (vaccination stays)
newcomp <- c("SIR", comp[3:11], comp[13:15])
col_nam2 <- compnamef(newcomp) 

#vector of element which are not SIR, over sex, riskgroups and age groups
vnsir <- array(0, dim=c(length(3:14),24), dimnames = list(col_nam2[3:14], col_nam2[fagmf(2,newcomp)] ) ) 
for(i in 1:12){
  vnsir[i,] <- fagmf(i, newcomp)
}
vnsir

## These for loops, take together, for each sex, riskgroups, agegroups, the disease stages for all different types, the % vaccinated and the number of deaths and makes 1- the sum
#of all this to have the SIR grouped, for the costs calculations. It groups it according to the newcomp vector

dfinal2 <- array(0, dim=c(length(time),length(col_nam2),length(vaccstrat)), dimnames = list(time,col_nam2, vaccstrat) )  

for(t in 0:length(time+1) ){
  for(i in 1:length(vaccstrat) ) {
    for( j in c(3:11 )){
      dfinal2[t,fagmf(j-1, newcomp),i] <-  rowSums(dfinal[t,fagmf(j, comp),i,] )  
    }# "CIN2" "CIN3" "ULCC" "DLCC" "URCC" "DRCC" "UDCC" "DDCC" "DCC" (comp[3:11])
    for(k in c(13:15)){
      dfinal2[t,fagmf(k-2, newcomp),i] <- (dfinal[t,fagmf(k, comp),i,1] )  
    }  #vaccinated    
    for(l in 1:length(vnsir[1,]) ){
      dfinal2[t,fagmf(1, newcomp)[l],i] <-   as.vector(N)[l] - sum( dfinal2[t, as.vector(vnsir[,l])+1 ,i  ] ) #for sex, risk group and age groups, SIR is
    }  #for sex, risk group and age groups, SIR is N[sex, riskgroup, agegroup] - sum(CIN2, CIN3, ULCC, etc.)
    } }

dfinal2[1,,1]

################### Function for load and show costs

#Big function giving, for women, cervical cancer and women vaccination the costs and icer

diffcosts <- c("pap_costs", "CIN2CIN3_costs", "cancer_costs", "vacc_costs", "AE_costs", "AE_QALY", "TOT_QALY", "TOT_costs")

cqoutcome <- array(0, dim=c(length(time),length(diffcosts), length(vaccstrat) ), dimnames = list(time,diffcosts, vaccstrat) )
icer  <- matrix(0, nrow= length(time), ncol= length(vaccstrat) )
icer2 <- matrix(0, nrow= length(time), ncol= length(vaccstrat) )
icer100 <- matrix(NA, ncol=6, nrow=2); colnames(icer100) <- vaccstrat;rownames(icer100) <- c("50%V4v_baseline", "novacc_baseline")

#function that gives the costs for the different vaccination strategies or parameters, and the ICER. Only women vaccination. FOr men vaccination, another function is needed
# yp <- 1

fcostsvaccstrat_varyp <- function(C4v_v, C9v_v){
  for(i in 1:length(vaccstrat) ){
    p4v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups, agegroups)  )
    p9v <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups, agegroups)  )
    if(i==1 ){
      p4v[1,,1] <- diffv[i] }else if (i==2){
        p4v[1,,1] <- diffv[i]} else if (i==3){ 
          p4v[1,,1] <- diffv[i]}else if (i==4 ){
            p9v[1,,1] <- diffv[i-3]} else if (i==5){
              p9v[1,,1] <- diffv[i-3]} else if (i==6){
                p9v[1,,1] <- diffv[i-3]}
    
    dataf <- dfinal2[,,i]
    # pap_costs
    cqoutcome[,1,i] <- apply(dataf, 1, function(x) (  sum( rep( rep(s[1,,1]*agw[2:7], each=length(riskgroups)),2) *  Cpt*(1-ssens)* x[ c( fag(2, newcomp), fag(3, newcomp) ) ]  )   #costs of false negative screens  * ()
                                                      +  sum( rep( rep(s[1,,1]*agw[2:7], each=length(riskgroups)),4) * (Cpt + (1-sspec) *Cfp ) * (x[c( fag(1, newcomp), fag(11, newcomp), #for SIR and vaccinated
                                                                                                                                              fag(12, newcomp), fag(13, newcomp) )] ) ) ) )#cost of follow up of false positives + negative screens
    # CIN2CIN3_costs 
    cqoutcome[,2,i] <- apply(dataf, 1, function(x) (  sum( rep( rep(s[1,,1]*agw[2:7], each=2),2)*ssens* 
                                                             rep( Ct[1:2], each=length(riskgroups)*length(agegroups) )* x[c( fag(2, newcomp), fag(3, newcomp))]))) #costs true positive treatement costs
    # cancer_costs 
    cqoutcome[,3,i] <- apply(dataf, 1, function(x) (  sum( rep( Ct[3:5,y], each=length(riskgroups)*length(agegroups) )*
                                                             rep( z[1:3,y], each=length(riskgroups)*length(agegroups) )*
                                                             x[c( fag(4, newcomp), fag(6, newcomp), fag(8, newcomp)  )] * rep(rep(agw[2:7], each=2 ),3) )) ) 
    # vacc_costs 
    # cqoutcome[,4,i] <- apply(dataf, 1, function(x) sum( C2v* x[c(  fag(9, newcomp)  )], C4v_v* x[c(  fag(10, newcomp)  )], C9v_v* x[c(  fag(11, newcomp)  )]  ) )
    
    cqoutcome[,4,i] <- apply(dataf, 1, function(x) sum( C2v* (as.vector(p2v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp  )
                                                        + C4v_v* (as.vector(p4v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp  )
                                                        + C9v_v * (as.vector(p9v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp  ) ))
    # AE_costs 
    # cqoutcome[,5,i] <- apply(dataf, 1, function(x) sum( (rme*Cme + rse*Cse)* x[c(  fag(9, newcomp)  )], (rme*Cme + rse*Cse)* x[c(  fag(10, newcomp)  )], (rme*Cme + rse*Cse)* x[c(  fag(11, newcomp)  )]  ) )
    
    cqoutcome[,5,i]  <- apply(dataf, 1, function(x) sum((rme*Cme + rse*Cse)* (as.vector(p2v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp )
                                                        + (rme*Cme + rse*Cse)* (as.vector(p4v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp )
                                                        + (rme*Cme + rse*Cse)* (as.vector(p9v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp ) ))
    # AE_QALY 
    # cqoutcome[,6,i] <- apply(dataf, 1, function(x) sum( (rme*Qme + rse*Qse)* x[c(  fag(9, newcomp)  )], (rme*Qme + rse*Qse)* x[c(  fag(10, newcomp)  )], (rme*Qme + rse*Qse)* x[c(  fag(11, newcomp)  )]  ) )
    
    cqoutcome[,6,i] <- apply(dataf, 1, function(x) sum( (rme*Qme + rse*Qse)*(as.vector(p2v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp ) #yp is an additional factor to adapt the number of school aged population to be vaccinated (reversed age-pyramid)
                                                        +  (rme*Qme + rse*Qse)* (as.vector(p4v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp )
                                                        +  (rme*Qme + rse*Qse)* (as.vector(p9v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp ) ))
    # TOT_QALY
    cqoutcome[,7,i] <- apply(dataf, 1, function(x) (  sum(x[c(fag(1, newcomp) ,fag(11, newcomp), fag(12, newcomp), fag(13, newcomp)) ] * rep(rep(agw[2:7], each=2 ),4) ) #sum healthy women
                                                      + sum( rep( q[1:2,y], each= length(riskgroups)*length(agegroups) )* x[c(fag(2, newcomp), fag(3, newcomp) ) ] * rep(rep(agw[2:7], each=2 ),2) ) #CIN2 and 3 qalys
                                                      + sum(q[3,y] * x[c(fag(4, newcomp), fag(5, newcomp) )] * rep(rep(agw[2:7], each=2 ),2) ) #localised cancer
                                                      + sum(q[4,y] * x[c(fag(6, newcomp), fag(7, newcomp) )] * rep(rep(agw[2:7], each=2 ),2) ) #regional
                                                      + sum(q[5,y] * x[c(fag(8, newcomp), fag(9, newcomp) )] * rep(rep(agw[2:7], each=2 ),2)) #distant
                                                      - sum( (rme*Qme + rse*Qse)* x[c(  fag(11, newcomp)  )] * rep(rep(agw[2:7], each=2 ),1), 
                                                             (rme*Qme + rse*Qse)* x[c(  fag(12, newcomp)  )]* rep(rep(agw[2:7], each=2 ),1), 
                                                             (rme*Qme + rse*Qse)* x[c(  fag(13, newcomp)  )]* rep(rep(agw[2:7], each=2 ),1)  )  # ae_qaly
    ) )
    
    # cqoutcome[,7,i] <- apply(dataf, 1, function(x) (  sum(x[c(fag(1, newcomp) ,fag(11, newcomp), fag(12, newcomp), fag(13, newcomp)) ] ) #sum healthy women
    #                                                   + sum( rep( q[1:2,y], each= length(riskgroups)*length(agegroups) )* x[c(fag(2, newcomp), fag(3, newcomp) ) ] ) #CIN2 and 3 qalys
    #                                                   + sum(q[3,y] * x[c(fag(4, newcomp), fag(5, newcomp) )]) #localised cancer
    #                                                   + sum(q[4,y] * x[c(fag(6, newcomp), fag(7, newcomp) )]) #regional
    #                                                   + sum(q[5,y] * x[c(fag(8, newcomp), fag(9, newcomp) )]) #distant
    #                                                   - sum( (rme*Qme + rse*Qse)*(as.vector(p2v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ]) ) # ae_qualy
    #                                                          +  (rme*Qme + rse*Qse)* (as.vector(p4v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ]) )
    #                                                          +  (rme*Qme + rse*Qse)* (as.vector(p9v)[1:4] *g[length(agegroups)]* as.vector(N[,,length(agegroups) ]) ) )) )
    #TOT_costs
    cqoutcome[,8,i] <-  apply(cqoutcome[,1:5,i], 1, function(x) sum(x )  )
  }
  
  
  return(cqoutcome)
}







