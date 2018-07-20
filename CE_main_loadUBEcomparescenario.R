# MR, 13.07.2018 
########### For loading results of UBELIX cluster scripts: load results from different vaccination scenarios, CE_UBE_comparescenario.R


# copied from R_ce folder: CE_load_final_v2.R

# Load Final data from models for specific outcomes and vaccination strategies but for all different hpv strains

#directories
dirl <- "tab_params/" #where to load parameter data from
plotdir <- "plots/" #where to save plots
loaddir <- "outp/m11/p" #where the UBE data was stored


### Load optimbeta from ubelix

source("CE_main_model_parameters.R")
sex <- c("f", "m")
riskgroups <- c("l", "h")
hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "other_HR")
time <- seq(0, 200, 1)
vaccstrat <- c("V4v_05","V4v_06" , "V4v_07", "V9v_05","V9v_06" , "V9v_07" )
diffv <- c(0.5,0.6,0.7)
y <-  which(dis_outcome=="cervix")

##############
dfinal <- array(0, dim=c(length(time),length(col_nam),length(vaccstrat), length(hpvtype)), dimnames = list(time,col_nam,
                                                                                                           vaccstrat, hpvtype) )  
for(i in 1:13){
  load( paste(loaddir,i,".RData", sep=""))
  dfinal[,,,i] <- dperHPVt
}


dim(dfinal)

###########

########
############################################# plot cancer incidences:
########
# in total, in women:
dim(dfinal)

totmodinc <- array(0, dim=c(length(time), length(vaccstrat), 1),  dimnames = list(time, vaccstrat , "incidence")  )
tempin <- c()
for(j in 1:length(time) ){
  for( i in 1:length(vaccstrat)){
    for(k in 0:(length(agegroups)-1) ){
      tempin[k+1] <- sum(rowSums(dfinal[j, 1+c(k*60+5,k*60+7 ,k*60+9, k*60+30+5,k*60+30+7 ,k*60+30+9) ,i,1:13] )  * rep( z[,y], 2) / (sr*agr[1+k])  * agw[2+k])
    }
    totmodinc[j,i,] <- sum(tempin  )
  }}

totmodinc[1,1,]
totmodinc[1,1,] * popmod*1e6*sr
totmodinc[1,1,] * 1e5

popmod <- 8.4

######### plotting:
col3 <- colorRampPalette(c("darkblue", "lightblue"))(3)
col4 <- colorRampPalette(c("darkgreen", "lightgreen"))(3)
col1 <- c( col3, col4)

plott <- 100
stt <- 0
denom <- 1e5
denom <- popmod*1e6*sr

pdf(paste(plotdir, "CCincidence_", denom,"_", Sys.Date(), ".pdf", sep=""), width=8, height=6)
par(mar=c(5,7,7,2))

plot(NA, xlim=c(stt,plott ), ylim=c(min(totmodinc[,,1]*denom),max(totmodinc[,,1]*denom)), 
     xlab= "years",
     ylab="HPV related \ncervical cancer incidence",
     cex.lab=2, frame.plot=FALSE)

for(i in 1:length(vaccstrat)){
  lines(time, totmodinc[,i,]*denom, col=col1[i], lwd=2) }

legend(50,max(totmodinc[,,1]*denom), legend=c( paste("V4v:", (diffv*100), "%", sep="") , paste("V9v:", (diffv*100), "%", sep="") ), 
       pch=16, col=col1, lwd=2, bty="n", title="Vaccine type and coverage" )

dev.off()

#decrease after 50 years
timt <- 101
totmodinc[1,,]
totmodinc[timt,,]

1- (totmodinc[timt,,]/totmodinc[1,,])
round(1- (totmodinc[timt,,]/totmodinc[1,,]) ,2)


########
##################################################### plot type ratios
########
spw <- array(0, dim=c( dim=length(time), length(vaccstrat), length(agegroups), length(comp), length(hpvtype)  ),
             dimnames =list( time, vaccstrat, agegroups, comp, hpvtype) )

for(k in 1: length(time)){
  for(l in 1: length(vaccstrat)){
    for(n in 1:length(hpvtype)){
      for(i in 1: length(agegroups)) {
        for(j in 1:length(comp)){
          spw[k,l,i,j,n] <- (dfinal[k, ( (i-1)* (length(comp)*4 )  + (j+1)),l,n] + dfinal[k,(i-1)*(length(comp)*4 ) + (length(comp)*2 )+(j+1),l,n] ) / (sr*agr[i])
        } }}}}



dim(spw)
spw[100, , ,,1]
spw[100, 1, ,,1]

spw[1,1,,,1]
colSums(spw[1,1,,,1]  )
#sum of all types
tt <- 1
tempd <- (spw[tt, 1, ,2,1:13])
rowSums(spw[tt, 1, ,2,c(1:10, 13)])

##################### cancer incidence per type
cincpt <- rep(0,length(hpvtype))
names(cincpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  cincpt[ht] <- ( z[1,y]*sum(spw[tt, 1, 1:6,5,ht]*agw[2:7]) + z[2,y]*sum(spw[tt, 1, 1:6,7,ht]*agw[2:7]) + z[3,y]*sum(spw[tt, 1, 1:6,9,ht]*agw[2:7] ) ) * 1e5
}

cincpt
cc_typratio <- cincpt / sum(cincpt)

ccpt1  <- read.csv(paste( dirl,"cc_inc_pertype4e6.csv", sep=""))
ccpt <- as.data.frame(ccpt1[,2])
row.names(ccpt) <- hpvtype
sum(cincpt)
sum(ccpt)
cincpt - ccpt

cc_typratio_lit <- ccpt /sum(ccpt)

#previously used screening rates:
# s[1,,1] <- c(0.28,0.56,0.66,0.68,0.61,0.38)

##################### CIN3 incidences
cin3incpt <- rep(0,length(hpvtype))
names(cin3incpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  cin3incpt[ht] <- ( pri2[1,ht,1]*sum( spw[tt, 1, 1:6,2,ht] *agw[2:7]) +  pr23[1,1]* sum( spw[tt, 1, 1:6,3,ht]  *agw[2:7])     ) 
}

##################### CIN3 treatements
cin3treatpt <- rep(0,length(hpvtype))
names(cin3treatpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  cin3treatpt[ht] <- sum( s[1,,1]*ssens* spw[tt, 1, 1:6,4,ht] *agw[2:7])       
}


##################### CIN2 incidences
cin2incpt <- rep(0,length(hpvtype))
names(cin2incpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  cin2incpt[ht] <- ( pri1[1,ht,1]*sum( spw[tt, 1, 1:6,2,ht] *agw[2:7])      ) 
}

##################### CIN treatements
cintreatpt <- rep(0,length(hpvtype))
names(cintreatpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  cintreatpt[ht] <- sum( s[1,,1]*ssens* rowSums(spw[tt, 1, 1:6,3:4,ht]) *agw[2:7])       
}

##################### annual deaths
dccincpt <- rep(0,length(hpvtype))
names(dccincpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  dccincpt[ht] <- sum( (colSums( um[1:3]* t(spw[tt, 1, 1:6,c(5,7,9),ht]) ) + colSums( dm[1:3]* t(spw[tt, 1, 1:6,c(6,8,10),ht]) ) )*agw[2:7] )     
}



sum( cin3incpt * 1e5 )

sum( cin3incpt) * (popmod/2)*1e6 #incidence
sum( cin2incpt) * (popmod/2)*1e6 #incidence
sum( dccincpt) * (popmod/2)*1e6 #mortality per year
sum(cintreatpt) * (popmod/2)*1e6 #total CIN cases treated
sum( cincpt) * (popmod/2)*1e1 #cancer incidence
sum(cin3treatpt) * (popmod/2)*1e6

sum(cincpt) * (popmod/2)*1e1

##### hpv types ratio for CIN3 (to compare it to CIN3+plus data)
tt <- 1
comp

#compare model ratio with data from the litteratues (Guan 2012, CIN3+plus study and Clifford 2003)
nGUAN_typerationormal <- read.csv(paste(dirl,"nGUAN_typerationormal.csv", sep=""))[,2] 
nGUAN_typeratiocin3 <- read.csv(paste(dirl,"nGUAN_typeratiocin3.csv", sep=""))[,2] 
nSwissratiocin3 <- read.csv(paste(dirl,"nSwissratiocin3.csv", sep=""))[,2] 

names(nGUAN_typeratiocin3) %in% hpvtype
hpvtype %in% names(nGUAN_typeratiocin3)

#for hpv infection prevalence
i_prevpt <- rep(0,length(hpvtype))
names(i_prevpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  i_prevpt[ht] <- ( sum(spw[tt, 1, 1:6,2,ht]*agw[2:7])  )
  # i_prevpt[ht] <- ( sum(spw[tt, 1, 1:6,2,ht]*agr)  )
  i_prevpt[ht] <- ( sum(spw[tt, 1, 2,2,ht])  )
}
i_typratio <- i_prevpt /sum(i_prevpt)

#for CIN2 prevalence
cin2_prevpt <- rep(0,length(hpvtype))
names(cin2_prevpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  cin2_prevpt[ht] <- ( sum(spw[tt, 1, 1:6,3,ht]*agw[2:7])  ) 
}
cin2_typratio <- cin2_prevpt /sum(cin2_prevpt)


#for CIN3 prevalence
cin3_prevpt <- rep(0,length(hpvtype))
names(cin3_prevpt) <- hpvtype
for(i in 1:11){
  ht <- c(1:10, 13)[i]
  cin3_prevpt[ht] <- ( sum(spw[tt, 1, 1:6,4,ht]*agw[2:7])  ) 
}
cin3_typratio <- cin3_prevpt /sum(cin3_prevpt)

##plot

pdf(paste(plotdir, "hpvtyperatio_viv2", Sys.Date(), ".pdf", sep=""), width=6, height=6)
par(mfrow=c(3,1))
par(oma = c(4, 4, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 4, 3, 2)) # make the plots be closer together
################################### First normal cytology
# par(mar=c(5,5,3,2))
plot(NA, xlim=c(1,11 ), ylim=c(0,0.8) , 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, xaxt="n", main="normal cytology \n(fitted)")
axis(side=1, at=1:11, labels=FALSE)
text(seq(1, 11, by=1), par("usr")[3] - 0.08, labels = hpvtype[-1*c(11,12)], srt = 45, pos = 1, xpd = TRUE)
points(1:11 ,i_typratio[-1*c(11,12)], pch=17)
points(1:11 ,nGUAN_typerationormal[-1*c(11,12)], col="lightgreen", pch=17)
lines(1:11 ,i_typratio[-1*c(11,12)], pch=17)
lines(1:11 ,nGUAN_typerationormal[-1*c(11,12)], col="lightgreen", pch=17)
legend(9,0.8, legend=c( "model based", "Guan (2012)"), 
       pch=17, col=c("black", "lightgreen"), lwd=1, cex=0.9, bty="n" )
################################### Third Cervical Cancer 
# par(mar=c(5,5,3,2))
plot(NA, xlim=c(1,11 ), ylim=c(0,0.8) , 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, xaxt="n", main="cevical cancer \n(fitted)")
axis(side=1, at=1:11, labels=FALSE)
text(seq(1, 11, by=1), par("usr")[3] - 0.08, labels = hpvtype[-1*c(11,12)], srt = 45, pos = 1, xpd = TRUE)
points( (1:11) ,cc_typratio[-1*c(11,12)], pch=19)
points((1:11) ,cc_typratio_lit$`ccpt1[, 2]`[-1*c(11,12)], col="lightblue", pch=19)
lines( (1:11) ,cc_typratio[-1*c(11,12)], pch=19)
lines((1:11) ,cc_typratio_lit$`ccpt1[, 2]`[-1*c(11,12)], col="lightblue", pch=19)
legend(9,0.8, legend=c( "model based", "Clifford (2003)"), 
       pch=19, col=c("black", "lightblue"), lwd=1, cex=0.9, bty="n" )
################## Overall x and y labels
mtext('HPV types', side = 1, outer = TRUE, line = 2)
mtext('Type-attributable ratio', side = 2, outer = TRUE, line = 2)
################################### Second CIN3 
# par(mar=c(5,5,3,2))
plot(NA, xlim=c(1,11 ), ylim=c(0,0.8) , 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, xaxt="n", main="CIN3")
axis(side=1, at=1:11, labels=FALSE)
text(seq(1, 11, by=1), par("usr")[3] - 0.08, labels = hpvtype[-1*c(11,12)], srt = 45, pos = 1, xpd = TRUE)
points((1:11) ,cin3_typratio[-1*c(11,12)], pch=15)
points((1:11) , nSwissratiocin3[-1*c(11,12)], col="orange", pch=15)
points((1:11) , nGUAN_typeratiocin3[-1*c(11,12)], col="red", pch=15)
lines((1:11) ,cin3_typratio[-1*c(11,12)], pch=15)
lines((1:11) , nSwissratiocin3[-1*c(11,12)], col="orange", pch=15)
lines((1:11) , nGUAN_typeratiocin3[-1*c(11,12)], col="red", pch=15)
legend(8.5,0.8, legend=c( "model based", "CIN3+plus study", "Guan (2012)"), 
       pch=15, col=c("black", "orange", "red"), lwd=1, cex=0.9, bty="n" )

dev.off()

########### -----------------------------------------    plot cancer incidences and hpv prev
prev_pt_2024 <- nGUAN_typerationormal * 0.2435849
pdf(paste(plotdir, "hpvprevandccinc_vi3", Sys.Date(), ".pdf", sep=""), width=6, height=6)
par(mfrow=c(3,1))
par(oma = c(3, 2, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together
################################### First normal cytology
par(mar=c(5,5,3,2))
plot(NA, xlim=c(1,11 ), ylim=c(0,0.06) , 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, xaxt="n", main="prevalence in women (20-24 y.o) with normal cytology\n(fitted)")
axis(side=1, at=1:11, labels=FALSE)
text(seq(1, 11, by=1), par("usr")[3] - 0.08, labels = hpvtype[-1*c(11,12)], srt = 45, pos = 1, xpd = TRUE)
points(1:11 ,i_prevpt[-1*c(11,12)], pch=17)
points(1:11 ,prev_pt_2024[-1*c(11,12)], col="lightgreen", pch=17)
lines(1:11 ,i_prevpt[-1*c(11,12)], pch=17)
lines(1:11 ,prev_pt_2024[-1*c(11,12)], col="lightgreen", pch=17)
legend(9,0.06, legend=c( "model based", "Guan (2012)"), 
       pch=17, col=c("black", "lightgreen"), lwd=1, cex=0.9, bty="n" )
################################### Third Cervical Cancer 
par(mar=c(5,5,3,2))
plot(NA, xlim=c(1,11 ), ylim=c(0,5) , 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, xaxt="n", main="cevical cancer incidence per 10,000 women \n(fitted)")
axis(side=1, at=1:11, labels=FALSE)
text(seq(1, 11, by=1), par("usr")[3] - 0.08, labels = hpvtype[-1*c(11,12)], srt = 45, pos = 1, xpd = TRUE)
points( (1:11) ,cincpt[-1*c(11,12)], pch=19)
points((1:11) ,ccpt$`ccpt1[, 2]`[-1*c(11,12)], col="lightblue", pch=19)
lines( (1:11) ,cincpt[-1*c(11,12)], pch=19)
lines((1:11) ,ccpt$`ccpt1[, 2]`[-1*c(11,12)], col="lightblue", pch=19)
legend(9,5, legend=c( "model based", "Clifford (2003)"), 
       pch=19, col=c("black", "lightblue"), lwd=1, cex=0.9, bty="n" )
################## Overall x and y labels

################################### Second CIN3 type ratio 
par(mar=c(5,5,3,2))
plot(NA, xlim=c(1,11 ), ylim=c(0,0.7) , 
     xlab= "",
     ylab="",
     cex.lab=2, frame.plot=FALSE, xaxt="n", main="Type attributable ratio for CIN3")
axis(side=1, at=1:11, labels=FALSE)
text(seq(1, 11, by=1), par("usr")[3] - 0.08, labels = hpvtype[-1*c(11,12)], srt = 45, pos = 1, xpd = TRUE)
points((1:11) ,cin3_typratio[-1*c(11,12)], pch=15)
points((1:11) , nSwissratiocin3[-1*c(11,12)], col="orange", pch=15)
points((1:11) , nGUAN_typeratiocin3[-1*c(11,12)], col="red", pch=15)
lines((1:11) ,cin3_typratio[-1*c(11,12)], pch=15)
lines((1:11) , nSwissratiocin3[-1*c(11,12)], col="orange", pch=15)
lines((1:11) , nGUAN_typeratiocin3[-1*c(11,12)], col="red", pch=15)
legend(9,0.8, legend=c( "model based", "CIN3+plus study", "Guan (2012)"), 
       pch=15, col=c("black", "orange", "red"), lwd=1, cex=0.9, bty="n" )

mtext('HPV types', side = 1, outer = TRUE, line = -1)
mtext('Type-attributable ratio', side = 2, outer = TRUE, line = -1)

dev.off()


############## Plot mortality reduction and reduction in CIN treatements
col3 <- colorRampPalette(c("darkblue", "lightblue"))(3)
col4 <- colorRampPalette(c("darkgreen", "lightgreen"))(3)
col1 <- c( col3, col4)

plott <-101

outcomoi <- c("CIN treated", "CC_cases",  "mortality")

totcasesinSwitz <- array(0, dim=c(length(time), length(vaccstrat), length(outcomoi)),  dimnames = list(time, vaccstrat ,outcomoi)  )

for(j in 1:length(vaccstrat)){
  for(tt in 1:length(time)){
    cintreatpt <- rep(0,length(hpvtype))
    cincpt <- rep(0,length(hpvtype))
    dccincpt <- rep(0,length(hpvtype))
    
    for(i in 1:11){
      ht <- c(1:10, 13)[i]
      cintreatpt[ht] <- sum( s[1,,1]*ssens* rowSums(spw[tt, j, 1:6,3:4,ht]) *agw[2:7]) 
      cincpt[ht] <- ( z[1,y]*sum(spw[tt, j, 1:6,5,ht]*agw[2:7]) + z[2,y]*sum(spw[tt, j, 1:6,7,ht]*agw[2:7]) + z[3,y]*sum(spw[tt, j, 1:6,9,ht]*agw[2:7] ) ) 
      dccincpt[ht] <- sum( (colSums( um[1:3]* t(spw[tt, j, 1:6,c(5,7,9),ht]) ) + colSums( dm[1:3]* t(spw[tt, j, 1:6,c(6,8,10),ht]) ) )*agw[2:7] )     
    }
    totcasesinSwitz[tt,j,1] <- sum(cintreatpt)* (popmod/2) *1e6
    totcasesinSwitz[tt,j,2] <- sum(cincpt)* (popmod/2) *1e6
    totcasesinSwitz[tt,j,3] <- sum(dccincpt)* (popmod/2) *1e6
  }}



pdf(paste(plotdir, "MortalityandCINtreatements", Sys.Date(), ".pdf", sep=""), width=8, height=6)
par(mfrow=c(1,2))
par(oma = c(3, 3, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together

for(j in c(1,3) ){
  plot(NA, xlim=c(0,length(1:plott) ), ylim=c(min(totcasesinSwitz[,,j]),max(totcasesinSwitz[,,j])) , 
       ylab= "",
       xlab="",
       frame.plot=FALSE,  main=outcomoi[j])
  
  for(i in 1:length(vaccstrat)){
    lines(time, totcasesinSwitz[,i,j], col=col1[i], lwd=2)
  }
  if(j==2){
    legend(60, max(totcasesinSwitz[,,j], na.rm=TRUE) , legend=c( paste("V4v:", (diffv*100), "%", c("","", ""), sep="") ,
                                                                 paste("V9v:", (diffv*100), "%", sep="") ), 
           pch=16, col=col1, lwd=2, bty="n", title="Vaccine type \nand coverage" )
  }
}
mtext('years after vaccination onset', side = 1, outer = TRUE, line = -1, cex=1.5)
mtext('number of cases', side = 2, outer = TRUE, line = 1, cex=1.5)
dev.off()

#differences in number of cin treatements and deaths between strategies and baseline
vaccstrat2 <- c(  "V4v: 0.5" ,"V4v: 0.6", "V4v: 0.7", "V9v: 0.5" ,"V9v: 0.6", "V9v: 0.7")

casesaverted <- array(0, dim=c(length(vaccstrat2), length(outcomoi)),  dimnames = list( vaccstrat2 ,outcomoi)  )

for( i in 1:length(vaccstrat2)){
  casesaverted[i,1] <- sum(totcasesinSwitz[2:102,1,1]) - sum(totcasesinSwitz[2:102,i,1]) 
  casesaverted[i,2] <-  sum(totcasesinSwitz[2:102,1,2]) - sum(totcasesinSwitz[2:102,i,2]) 
  casesaverted[i,3] <-  sum(totcasesinSwitz[2:102,1,3]) - sum(totcasesinSwitz[2:102,i,3])
}

write.table(round(casesaverted,0), paste("O:/Cost_effectiveness project/Preparation_drafts/tab/cases_avertedcompB", Sys.Date(), ".txt" ,sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)

################ Plot HPV prevalence per age group
#for hpv infection prevalence
dim(spw)

defendtime <- 70 

t2 <- matrix(NA, ncol=length(agegroups), nrow = length(1:defendtime))

prevperage_OT <- array(0, dim=c( dim=length(time), length(vaccstrat), length(agegroups)  ),  dimnames =list( time, vaccstrat, agegroups ) )

for(i in 1:length(vaccstrat)){
  prevperage_OT[,i,] <- t( apply(spw[,i,,2,], 1, function(x) rowSums(x) ) )
}

ptitles <- c("50% quadrivalent", "60% quadrivalent", "70% quadrivalent", "50% nonavalent", "60% nonavalent", "70% nonavalent")


pdf(paste(plotdir, "prevperage_OT", Sys.Date(), ".pdf", sep=""), width=8, height=6)
par(mfrow=c(2,3))
par(oma = c(3, 3, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together

col4 <- colorRampPalette(c("darkblue", "lightblue"))(length(agegroups))

for(j in 1:length(vaccstrat)){
  
  plot(NA, xlim=c(1,length(1:defendtime) ), ylim=c(0,max(prevperage_OT)+0.05) , 
       ylab= "",
       xlab="",
       frame.plot=FALSE,  main=ptitles[j])
  
  for(i in 1:length(agegroups)){
    lines(time, prevperage_OT[,j,i], col=col4[i], lwd=2)
  }
  if(j ==6){
    legend(defendtime-20,max(prevperage_OT)+0.06, legend=agegroups, 
           pch=19, col=col4, lwd=1, cex=0.9, bty="n", title="age groups")
  }
  ################## Overall x and y labels
  mtext('HPV prevalence', side = 1, outer = TRUE, line = 1, cex=1)
  mtext('years after vaccination onset', side = 2, outer = TRUE, line = 1, cex=1)
}
dev.off()




############# other plot: age group after 0,10,20,30,40, 50 y
tavta <- c(0, 10, 20, 50) #time after vaccination to look at 
prevperage_OT2 <- prevperage_OT[ tavta+1,,]

col5 <- colorRampPalette(c( "darkblue", "lightblue"))(length(tavta))


ptitles <- c("50% quadrivalent", "60% quadrivalent", "70% quadrivalent", "50% nonavalent", "60% nonavalent", "70% nonavalent")


pdf(paste(plotdir,"prevperage_OT2", Sys.Date(), ".pdf", sep=""), width=8, height=6)
par(mfrow=c(2,3))
par(oma = c(3, 3, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(5, 2, 1, 1)) # make the plots be closer together

for(j in 1:length(vaccstrat)){
  plot(NA, xlim=c(1,length(1:length(agegroups)) ), ylim=c(0,max(prevperage_OT2)+0.05) , 
       ylab= "",
       xlab="",
       frame.plot=FALSE, xaxt="n"  ,main=ptitles[j])
  
  axis(side=1, at=1:length(agegroups), labels=FALSE)
  text(seq(1, length(agegroups), by=1), par("usr")[3] - 0.005, labels = agegroups, srt = 0, pos = 1, xpd = TRUE)
  # mtext('age groups', side = 1, outer = TRUE, line = -2)
  # mtext('HPV prevalence', side = 2, outer = TRUE, line = 0)
  
  for(i in 1:length(tavta)){
    points(1:length(agegroups) , prevperage_OT2[i,j,], pch=19, col=col5[i]) 
    lines(1:length(agegroups) , prevperage_OT2[i,j,], col=col5[i]) 
  }
  if(j== 6){
    legend(length(agegroups)-1.5,max(prevperage_OT2)+0.04, legend=tavta, 
           pch=19, col=col5, lwd=1, cex=0.9, bty="n" ,title="years after \nvaccination start")
  }
  mtext('HPV prevalence', side = 2, outer = TRUE, line = 1, cex=1)
  mtext('age groups', side = 1, outer = TRUE, line = 1, cex=1)
}

dev.off()





######### plotting:

plott <- 100
stt <- 0

pdf("CCincidence260418.pdf", width=8, height=6)
par(mar=c(5,7,7,2))

plot(NA, xlim=c(stt,plott ), ylim=c(min(totmodinc[,,1]),max(totmodinc[,,1])), 
     xlab= "years",
     ylab="HPV related \ncervical cancer incidence",
     cex.lab=2, frame.plot=FALSE)

for(i in 1:length(vaccstrat)){
  lines(time, totmodinc[,i,], col=col1[i], lwd=2)
}

legend(50,max(totmodinc[,,1]), legend=c( paste("V4v:", (diffv*100), "%", sep="") , paste("V9v:", (diffv*100), "%", sep="") ), 
       pch=16, col=col1, lwd=2, bty="n", title="Vaccine type and coverage" )

dev.off()



