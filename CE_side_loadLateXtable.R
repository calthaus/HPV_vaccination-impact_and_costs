# MR, 13.07.2018 
########### Write parameter table and save it in the tex table folder


# copied from R_ce folder: CE_params_table_v2.R

# Write parameters tables
require(plyr)
# sexbehavdata <- "sbehav_ch3.csv"
source("CE_main_model_parameters.R")
source("CE_main_model_costsparameters.R")

popmod <- 8.4 #population end 2016

#directories

#where to retrieve data from
# dirl <- ""
dirl <- "tab_params/"

#where to save the data to
dirlts <- "O:/Cost_effectiveness project/Preparation_drafts/tab/"

#define date of betaparams to be used 
dobeta <- "2018_07_23"
dopr3c <- "2018_07_23"


#load scaled data
pr3c <- as.numeric(read.csv(paste(dirl, "optimised_pr3c",dopr3c, ".csv", sep="") )[,2] )
pr2c <- 0.2*pr3c
betap[1,1:13] <- as.numeric(read.csv(paste(dirl,"optimised_beta", dobeta ,".csv", sep=""))[1,2:14] )
betap[2,1:13] <- as.numeric(read.csv(paste(dirl,"optimised_beta", dobeta, ".csv", sep=""))[2,2:14] )


hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

#These tables are written in a text readble by Latex in order to directly produce a latex table.

## General parameters, valid for all types and disease outcomes: ----------------------------------------------------------
rownames_generalparams <- c( "$N_{sra}$", "$r_s$", "$r_r$", "$r_a$", "$\\Omega_a$"    , "$C_{sra}$", "$m$", 
                             "$\\epsilon_{r}$" , "$\\epsilon_{a}$" )  
tab_generalparams <- as.data.frame( matrix(ncol=2, nrow=length(rownames_generalparams ) ) ) 
row.names(tab_generalparams) <- rownames_generalparams
colnames(tab_generalparams) <- c( "Description", "Source")

tab_generalparams[ "$N_{sra}$",1] <- paste("Number of 18-64 year olds of gender 's', and activity group 'r' and age 'a', total population is estimated at: ", 
                                           popmod , "mio." , sep="")
tab_generalparams["$r_s$",1] <- paste("Gender: (proportion in each group): female (1)=", 0.5, "male (2)=", 0.5)
tab_generalparams["$r_r$",1] <- paste("Sexual activity: low (1) =", round(rgr[1],2), "high (2) =", round(rgr[2],2))
tab_generalparams["$r_a$",1] <- paste("Age groups: ", paste(paste( agegroups," = ", round(agr,2),", " , sep=""), collapse="" ), sep="" )
tab_generalparams["$\\Omega_a$" ,1] <- paste("Rate at which individuals leave their age group: ", 
                                             paste(paste( agegroups,"=", round(g,2),", " , sep=""), collapse="" ), sep="" )
tab_generalparams["$C_{sra}$",1] <- paste("Mean number of sexual contacts per gender, sexual activity group and age: ", 
                                          "female in sexual activity group 1: ",paste("\\emph{", agegroups  ,"=}" , adply(round(C,2),1:2)[1,3:8], collapse = ", " ), 
                                          ", male in sexual activity group 1: ",paste("\\emph{", agegroups  ,"=}" ,adply(round(C,2),1:2)[2,3:8], collapse = ", " ), 
                                          ", female in sexual activity group 2: ",paste("\\emph{", agegroups  ,"=}" , adply(round(C,2),1:2)[3,3:8], collapse = ", " ), 
                                          ", male in sexual activity group 2: ",paste("\\emph{", agegroups  ,"=}" , adply(round(C,2),1:2)[4,3:8], collapse = ", " ), "for age groups 1 to 6, resp.")
tab_generalparams['$m$',1] <- paste("Rate at which individuals can change activity group: ", "m = ", m, sep="" )    
tab_generalparams[ "$\\epsilon_{r}$",1] <- paste("Assortativity index between sexual activity groups: ", "$\\epsilon_{r}$"," = ",
                                                 eps_r, sep="" )   
tab_generalparams[ "$\\epsilon_{a}$",1] <- paste("Assortativity index between age groups: ", "$\\epsilon_{a}$ = ",
                                                 round(eps_a,2), sep="" )   
#source:
tab_generalparams[,2] <- c("FSO (2016)","FSO (2016)", "adapted from SHS (2012)", "data driven", "data driven", "adapted from SHS (2012)",
                           "Althaus (2015), Fingerhuth (2016)", "Althaus (2015), Fingerhuth (2016)", "adapted from Natsal'3 (2010)" )

# write.table(tab_generalparams, "O:/Cost_effectiveness project/Preparation_drafts/tab/P_general.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
# write.table(tab_generalparams, "E:/Cost_effectiveness project/Preparation_drafts/tab/P_general.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tab_generalparams, paste(dirlts,"P_general_", Sys.Date(), ".txt" ,sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)


## Model parameters, which vary by HPV types: ----------------------------------------------------------
rownames_hpvtparams <- c( "$\\beta_h$" ,  "$\\gamma_{1h}$" ,   "$\\gamma_{2h}$" , "$\\omega_h$" , "$\\eta_h$" ,  "$\\alpha_h$", "$\\nu_{oh}$",
                          "$\\phi_{12}$" , "$\\phi_{13}$", "$\\phi_{23}$" , "$\\phi_{3c}$","$\\phi_{2c}$" )  
tab_hpvtparams <- as.data.frame( matrix(ncol=( length(hpvtype )+1), nrow=length(rownames_hpvtparams ) ) ) 
row.names(tab_hpvtparams) <- rownames_hpvtparams
colnames(tab_hpvtparams) <- c("Descr.", hpvtype)

tab_hpvtparams["$\\beta_h$" ,1] <- "Transmission probability: "
tab_hpvtparams["$\\beta_h$" ,2:( length(hpvtype )+1)] <-  c(round(betap[1,],2) )

tab_hpvtparams["$\\gamma_{1h}$",1] <- c("Clearance rate (women): ") 
tab_hpvtparams["$\\gamma_{1h}$",2:( length(hpvtype )+1)] <- c(round(gamma[1,],2))

tab_hpvtparams["$\\gamma_{2h}$",1] <- c("Clearance rate (men): ") 
tab_hpvtparams["$\\gamma_{2h}$",2:( length(hpvtype )+1)] <- c(round(gamma[2,],2))

tab_hpvtparams["$\\omega_h$",1] <- c("Immunity decline rate: ") 
tab_hpvtparams["$\\omega_h$",2:( length(hpvtype )+1)] <- c(round(omega[1,],3))

tab_hpvtparams["$\\eta_h$",1] <- c("Vaccine efficacy (\\%): ") 
tab_hpvtparams["$\\eta_h$",2:( length(hpvtype )+1)] <- ve *100

tab_hpvtparams["$\\alpha_h$",1] <- c("Type-specific antibodies developement (prop.): ") 
tab_hpvtparams["$\\alpha_h$",2:( length(hpvtype )+1)] <- c(round(alpha,2))

tab_hpvtparams["$\\nu_{oh}$",1] <- c("Relative risk of re-infection (cervix): ") 
tab_hpvtparams["$\\nu_{oh}$",2:( length(hpvtype )+1)] <- c(tau[,1])

tab_hpvtparams["$\\phi_{12}$",1] <- "Rate of progression (cervix), from infection to CIN2: "
tab_hpvtparams["$\\phi_{12}$",2:( length(hpvtype )+1)] <- c(pri1[1,,1])

tab_hpvtparams["$\\phi_{13}$",1] <- "Rate of progression (cervix), from infection to CIN3: "
tab_hpvtparams["$\\phi_{13}$",2:( length(hpvtype )+1)] <- c(pri2[1,,1])

tab_hpvtparams["$\\phi_{23}$",1] <- "Rate of progression (cervix), from CIN2 to CIN3: "
tab_hpvtparams["$\\phi_{23}$",2:( length(hpvtype )+1)] <- c(round(pr23[,1],2))

tab_hpvtparams["$\\phi_{3c}$",1] <- "Rate of progression (cervix), from CIN3 to cancer: "
tab_hpvtparams["$\\phi_{3c}$",2:( length(hpvtype )+1)] <- c(round(pr3c,3))

tab_hpvtparams["$\\phi_{2c}$",1] <- "Rate of progression (cervix), from CIN2 to cancer: "
tab_hpvtparams["$\\phi_{2c}$",2:( length(hpvtype )+1)] <- c(round(pr2c,3))

# write.table(tab_hpvtparams[,2:( length(hpvtype )+1)], "O:/Cost_effectiveness project/Preparation_drafts/tab/P_perhtype.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
# write.table(tab_hpvtparams[,2:( length(hpvtype )+1)], "E:/Cost_effectiveness project/Preparation_drafts/tab/P_perhtype.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
# write.table(tab_hpvtparams[,2:( length(hpvtype )+1)], paste("O:/Cost_effectiveness project/Preparation_drafts/tab/P_perhtype_", Sys.Date(), ".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tab_hpvtparams[,c(2:11, 14)], paste(dirlts,"P_perhtype_", Sys.Date(), ".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)



#second table with description and source
#source:
tab_hpvtparams2 <- as.data.frame( matrix(ncol=2, nrow=length(rownames_hpvtparams ) ) ) 
row.names(tab_hpvtparams2) <- rownames_hpvtparams
colnames(tab_hpvtparams2) <- c("Descr.", "Source")
tab_hpvtparams2[,1] <- tab_hpvtparams[,1]
tab_hpvtparams2[,2] <- c("scaled", "Jaisamrarn (2013)", "Moreira (2014)", "Johnson (2012)", "not incl. yet",
                         "Faust (2013)", "Castellsagu\\'e (2014)", "Durham (2012)", " '' ", " Moscicki (2010)", "scaled", "scaled")

# write.table(tab_hpvtparams2, "O:/Cost_effectiveness project/Preparation_drafts/tab/P_perhtype2.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
# write.table(tab_hpvtparams2, "E:/Cost_effectiveness project/Preparation_drafts/tab/P_perhtype2.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tab_hpvtparams2, paste(dirlts,"P_perhtype2_", Sys.Date(), ".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)


## Model parameters, which varies by disease: ----------------------------------------------------------
rownames_hpvdparams <- c( "$\\sigma_{oa}$",  "$\\varsigma$" ,  "$\\tau$" ,
                          "$\\psi_{od}$"  ,  "$\\vartheta_{od}$" ,"$\\phi_{od}$" , "$\\varphi_{od}$", "$\\mu_{od}$") 

tab_hpvdparams <- as.data.frame( matrix(ncol=2, nrow=length(rownames_hpvdparams ) ) ) 
row.names(tab_hpvdparams) <- rownames_hpvdparams
# colnames(tab_hpvdparams) <-  c( "Param.","Description", "Source")

tab_hpvdparams["$\\sigma_{oa}$",1:2] <- c(paste("Age-specific screening rates for cervical cancer in women, ", "$\\sigma_{a}$",
                                                ": ", paste(paste("\\emph{",agegroups, "=}" ,round(s[1,,1],2), collapse=", "), sep="" ) ),
                                          "SHS (2012)")

tab_hpvdparams["$\\varsigma$", 1:2] <- c(paste("Sensitivity of screening (for cervical cancer, in women), ",  "$\\varsigma$", " = ",
                                               ssens, sep=""), "Bigras (2005)")
#"Durham (2012)")
tab_hpvdparams["$\\tau$" , 1:2] <- c(paste("Specificity of screening (for cervical cancer, in women), ", "$\\tau$", " = ",
                                           sspec, sep=""),"Bigras (2005)")
# "Durham (2012)")

tab_hpvdparams["$\\psi_{od}$", 1:2] <- c(paste("Probability of cure with treatement, per cancer stage, ", 
                                               "for cervical cancer: ",  paste("\\emph{", c("local", "distant", "regional"),"=}", paste(pi[,1]), collapse = ", " ), ".", sep="" ),
                                         "Durham (2012)" )

tab_hpvdparams["$\\vartheta_{od}$" , 1:2] <- c(paste("Rate of diagnostic, per cancer stage, ", 
                                                     "for cervical cancer: ",  paste("\\emph{", c("local", "distant", "regional"),"=}",  paste( round(z[,1],2) ), collapse = ", " ), ".", sep="" ),
                                               "Campos (2014)")  
# "Durham (2012)")

tab_hpvdparams["$\\phi_{od}$" , 1:2] <- c(paste("Progression rate to the next disease stage.", 
                                                # c(paste("Progression rate to the next disease stsage, from CIN2 to L (2 to 4): ", paste( round(pr2c,3), collapse = ", " ), ".", 
                                                " From L to R (4 to 6): ", paste( round(prlcrc, 3), collapse = ", " ), ".", 
                                                " From R to D (4 to 6): ", paste( round(prrcdc, 3), collapse = ", " ), "."
                                                ,sep="" ),    "Campos (2012), McCredie (2008), 2cd")
# "Durham (2012)")


tab_hpvdparams["$\\varphi_{od}$" , 1:2] <- c(paste("Regression rate to normal cytology (S), from CIN2, ", 
                                                   "for cervix in HPV16/18 and other HPV types, resp.: ", paste(round(pr2s[2:3,1],2), collapse = ", " ), ". " ,
                                                   "Regression from CIN3, for cervix: ",  paste(paste(round(pr3s[1,1],2)), collapse = ", " ), " for all HPV types", ". " , sep="" ),
                                             "Moscicki (2010) and 2cd McCredie but X")# 
# "Durham (2012)")
tab_hpvdparams["$\\mu_{od}$" , 1:2] <- c(paste("Disease outcome and disease stage specific mortality, ", 
                                               "for cervical cancer, diagnosed: ", 
                                               paste("\\emph{", c("local", "regional", "distant") ,"=}", round(dm,2), collapse = ", " ) ,
                                               " and for undiagnosed: ", 
                                               paste("\\emph{", c("local", "regional", "distant") ,"=}", round(um,2), collapse = ", " ) , ".", sep="" ),
                                         "Campos (2014)")

# write.table(tab_hpvdparams, "O:/Cost_effectiveness project/Preparation_drafts/tab/P_perdisease.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
# write.table(tab_hpvdparams, "E:/Cost_effectiveness project/Preparation_drafts/tab/P_perdisease.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tab_hpvdparams, paste(dirlts,"P_perdisease_", Sys.Date(), ".txt" , sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)



## Model parameters, for economic evaluation: COSTs and QALY: ----------------------------------------------------------
# rownames_economicev <- c( "$cv$" ,  "$ct$" , "$cns$", "$cfp$" , "$cae$",  "$qds$", "$qae$",  "$rse$" )  
# tab_economicev <- as.data.frame( matrix(ncol=2, nrow=length(rownames_economicev ) ) ) 
# row.names(tab_economicev) <- rownames_economicev
# 
# tab_economicev["$cv$" ,1:2] <- c( paste("Vaccination costs (for two doses) in CHF, ", "for the quadrivalent vaccine (V4v): ", C4v,
#                                         " and for the nonavalent vaccine (V9v): ",  C9v, sep=""),
#                                   "GDK (2008) for the V4v")
# tab_economicev["$ct$" ,1:2] <- c(paste("Treatement costs for the different disease stages, in CHF: ",
#                                        paste("\\emph{", c("CIN2", "CIN3", "local", "regional", "distant") ,"=}", Ct[,1], collapse=", "), sep=""), 
#                                  "Durham (2012), Szucs (2008)" )
# tab_economicev["$cns$" ,1:2] <- c(paste("Costs for per negative cervical screen, in CHF: ", Cpt, sep=""),
#                                   "Szucs (2008)" )
# tab_economicev["$cfp$" ,1:2] <- c(paste("Costs of follow-up for false-positive cervical screen, in CHF: ", Cfp, sep=""),
#                                   "Szucs (2008)")
# tab_economicev["$cae$" ,1:2] <- c(paste("Costs of mild and severe adverse effects due to vaccination, resp., in CHF: ", Cme, " and ", Cse, sep=""),
#                                   "Durham (2012)")
# tab_economicev["$qds$" ,1:2] <- c(paste("Life-year weight (QALY) per disease stage for cervical cancer: ",
#                                         paste("\\emph{", c("CIN2", "CIN3", "local", "regional", "distant") ,"=}",q[,1], collapse=", "), sep=""),
#                                   "Durham (2012)")
# tab_economicev["$qae$" ,1:2] <- c(paste("QALY lost per mild and severe adverse effects due to vaccination, resp. : ", Qme, " and ", Qse, sep=""),
#                                   "Durham (2012)")
# tab_economicev["$rse$" ,1:2] <- c(paste("Rates of mild and severe adverse effects due to vaccination, resp. : ", rme, " and ", rse, sep=""),
#                                   "Durham (2012)")

rownames_economicev <- c( "$f_1$" ,  "$f_2$" , "$f_{3-7}$" , "$f_{8-9}$", "$f_{10-11}$", "$q_{1-5}$",  "$q_{6-7}$", "$e_{1-2}$"  )  
tab_economicev <- as.data.frame( matrix(ncol=2, nrow=length(rownames_economicev ) ) ) 
row.names(tab_economicev) <- rownames_economicev

tab_economicev["$f_1$" ,1:2] <- c( paste("Costs of a negative pap-screen: ",  Cpt, sep=""),
                                   "Szucs (2008)")
tab_economicev["$f_2$" ,1:2] <- c(paste("Costs of follow-up for false-positive cervical screen: ",
                                        Cfp,   sep="") , "Szucs (2008)" )
tab_economicev["$f_{3-7}$" ,1:2] <- c(paste("Treatement costs for the different disease stages, in CHF: ",
                                            paste("\\emph{", c("CIN2", "CIN3", "local", "regional", "distant") ,"=}", Ct[,1], collapse=", "), sep=""), 
                                      "Durham (2012), Szucs (2008)" )
tab_economicev["$f_{8-9}$" ,1:2] <- c( paste("Vaccination costs (for two doses) in CHF, ", "for the quadrivalent vaccine (V4v): ", C4v,
                                             " and for the nonavalent vaccine (V9v): ",  C9v, sep=""),
                                       "GDK (2008) for the V4v")

tab_economicev["$f_{10-11}$" ,1:2] <- c(paste("Costs of mild and severe adverse effects due to vaccination, resp., in CHF: ", Cme, " and ", Cse, sep=""),
                                        "Durham (2012)")
tab_economicev["$q_{1-5}$" ,1:2] <- c(paste("Life-year weight (QALY) per disease stage for cervical cancer: ",
                                            paste("\\emph{", c("CIN2", "CIN3", "local", "regional", "distant") ,"=}",q[,1], collapse=", "), sep=""),
                                      "Durham (2012)")
tab_economicev["$q_{6-7}$" ,1:2] <- c(paste("QALY lost per mild and severe adverse effects due to vaccination, resp. : ", Qme, " and ", Qse, sep=""),
                                      "Durham (2012)")
tab_economicev["$e_{1-2}$" ,1:2] <- c(paste("Rates of mild and severe adverse effects due to vaccination, resp. : ", rme, " and ", rse, sep=""),
                                      "Durham (2012)")

# write.table(tab_economicev, "O:/Cost_effectiveness project/Preparation_drafts/tab/P_economic.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tab_economicev, paste(dirlts,"P_economic_", Sys.Date(), ".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)



##Cancer incidence, from NICER: ----------------------------------------------------------
rownames_cinc <- c( "Cervical cancer inc."  )  
tab_cinc <- as.data.frame( matrix(ncol= ( length(hpvtype )), nrow=length(rownames_cinc ) ) ) 
row.names(tab_cinc) <- rownames_cinc
colnames(tab_cinc) <- c(hpvtype)
#load cervical cancer incidence per type 
ccpt1  <- read.csv( paste( "tab_params/cc_inc_pertype4e6.csv", sep="") )
ccpt <- as.data.frame(ccpt1[,2])

cc_inc_pertype  <- as.numeric(ccpt$`ccpt1[, 2]`)


tab_cinc["Cervical cancer inc." ,] <- c(round( cc_inc_pertype,2) )


# write.table(tab_cinc, "O:/Cost_effectiveness project/Preparation_drafts/tab/P_c_inc.txt", quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(tab_cinc[-c(11,12)], paste(dirlts,"P_c_inc_", Sys.Date(), ".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)

rownames_prevh <- c( "HPV prev."  )  
tab_prevh <- as.data.frame( matrix(ncol= ( length(hpvtype )), nrow=length(rownames_prevh ) ) ) 
row.names(tab_prevh) <- rownames_prevh
colnames(tab_prevh) <- c(hpvtype)

#load type ratio for normal cytology
nGUAN_typerationormal <- read.csv(paste( "tab_params/nGUAN_typerationormal.csv", sep=""))[,2]
prev_pt_2024 <- nGUAN_typerationormal * 0.2435849

tab_prevh["HPV prev." ,] <- round(nGUAN_typerationormal*0.2435849*100, 2)

write.table(tab_prevh[-c(11,12)],  paste(dirlts,"hpvprev2024.txt",sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)


################### Table for age groups weights, retrieved from the ESP 2013
rownames_agw <- c("Age group weight"  )  
tab_agw  <- as.data.frame( matrix(ncol= ( length(agegroups )+2), nrow=length(rownames_agw ) ) ) 
row.names(tab_agw) <- rownames_agw
colnames(tab_agw) <- as.vector(names(agw) )

tab_agw["Age group weight" ,] <- round(agw, 2)


write.table(tab_agw, paste( dirlts,"agegroupespw" , Sys.Date(),".txt", sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)



