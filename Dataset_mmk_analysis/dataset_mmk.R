# This script implements diffussion controled mass transfer of substrate into MM kinetics to generate 
# moisture sensitivity prediction of microbial respiration 
rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(GGally)

mpdata <- read.csv(file = 'datasets.csv') # dataset use a dummy colomn to scale all ra measurements
#==subsetting unique sample ====
rate<-select(filter(mpdata,ra_s1==1), c(id, rs, mvol, mwp, porosity, clay, bd,org))

#===calculating accessible soc an oxy at steady-state====
opt<-function(po,mvol,bd,org){
  soc<-org/12*bd*1e-3  #assum ss 1% as DOC  (mol/cm3) and 10% of these DOC are accessible
  oxy<-(po-mvol)*0.21*1.3e-6  #(mol/cm3)Henry's law constant (Caq/p=1.3e-3 mol/L/atm)   
  return (list(soc, oxy))}
po<-rate$porosity
mvol<-rate$mvol
org<-rate$org
bd<-rate$bd
ss_opt<-opt(po,mvol,bd,org)
rate$soc_ss<-ss_opt[[1]] ##approximately 100mol/cm3
rate$oxy_ss<-ss_opt[[2]]

#=====preparing dynamic enzyme-substrate binding (as Km)======================================
Npsite=3000 # number of transporter per cell
k2_p=100  #transporter specific substrate uptake rate, unit:s-1
rc=1.e-6 #microbial cell radius unit:m
rp=1.e-9  #transporter radius unit:m
Na=6.e23  #Avogadro number
Ratm=50 #atmospheric resistance, 50 s/m

temp=25+273.15#define temperature

#=calculating gaseous and aqueous tortuosity (Original paper based on Morldrup, 2003)
tau<-function (mvol, po, clay){
  wpo<-mvol  #water filled porosity
  gpo<-po-mvol #gas filled porosity
  b<-2.91+0.195*clay #shape parameter 
  taug<- (po-mvol)*((po-mvol)/po)^3/b
  tauw<- mvol*(mvol/po)^(3/b-1)
  return (list(wpo=wpo, gpo=gpo,taug=taug,tauw=tauw))
}
mvol<-rate$mvol
po<-rate$porosity
clay<-rate$clay
tortuo<-tau(mvol, po,clay)

#==calculating O2 and C substrate diffusivity
Diffu<-function (gpo, wpo, taug, tauw){
  Dw_o2=1.4e-9*temp/298.0  #aqueous tracer diffusivity at 25
  Dg_o2=1.8e-5*(temp/273.0)**1.82 #oxygen diffusivity in gas phase
  henry_o2=3.2e-2*exp(-1500.*(1/temp-1/298.15))
  Dwo2= 0.5*(Dg_o2*taug*gpo*henry_o2+Dw_o2*tauw*wpo) #bulk aqueous molecular diffusivity as a colume weighted average between aquesous and gaseous phases
  Dw_s=6e-9 #oxygen diffusivity in water
  Dws=Dw_s*tauw*wpo#bulk substrate diffusivity (between the soil matrix and microsite)
  return(list(Dwo2=Dwo2, Dws=Dws))
}

gpo<-tortuo$gpo
wpo<-tortuo$wpo
taug<-tortuo$taug
tauw<-tortuo$tauw
Dw<-Diffu(gpo, wpo, taug, tauw)

#==calculating Km for OC (mol m-3)========
Kaff<-function(Ncell, gpo, wpo, taug, tauw, Dws,Dwg,mp){
  Bdens<-Ncell/Na #free microsite microbial abundance mol/m3
  Rm<-rc*(40*Ncell)^(1/3)  #microsite radius
  vm<-pi*4/3*Rm^3  #microsite volume 
  kw2<-k2_p*Npsite #(Npsite*rp+pi*rc) #cell specific uptake rate for OX, unit s-1
  Dw_s0 = 6e-9  #for claculating reference affinity
  Dw_g0 = 1.4e-9
  fin<-Npsite*rp/ (Npsite*rp+pi*rc)#interception probability (number of molecules that 1 mol cell will encounter and be able to intercept)
  ksw1<-4*pi*Dw_s0*rc*fin*Na   #substrate delivery parameter unit m3 mo-1 s-1
  kow1<-4*pi*Dw_g0*rc*fin*Na
  ksw0<-kw2/ksw1#reference affinity (Km used in MM kinetics in a well-mixed solution)
  kow0<-kw2/kow1#reference affinity (Km used in MM kinetics in a well-mixed solution)
  film<-exp(-13.65-0.857*log(mp/1000))    #calculting water film thickness
  ks_con<- (film/(Rm*Dw_s0*(Rm+film)) + 1/(Dws*(Rm+film)))*vm/(4*pi)#conductance coefficient
  ks_aff<-ksw0*(1+ks_con*ksw1*Ncell/Na/vm)
  ko_con<- (film/(Rm*Dw_g0*(Rm+film)) + 1/(Dwg*(Rm+film)))*vm/(4*pi)#conductance coefficient
  ko_aff<-kow0*(1+ko_con*kow1*Ncell/Na/vm)
  return(list(ksw0=ksw0, kow0=kow0,KaffOC=ks_aff, Kaffoxy=ko_aff))
}

mp<-rate$mwp
Dws<-Dw$Dws
Dwg<-Dw$Dwo2
Ncell<-2.68e-10*rate$soc_ss*1e-6*6.02e23 
test<-Kaff(Ncell, gpo, wpo, taug, tauw, Dws,Dwg,mp)

rate$Keffs<-test$KaffOC*1e-3
rate$Keffo<-test$Kaffoxy*1e-3

#====preparing dataframe for ss and nonss modeling============================================
userate<-select(rate, c(id, soc_ss, oxy_ss,Keffs, Keffo))
usesoil<-select(mpdata, c(id,rs,ra, mvol, mwp,  bd, porosity, clay, org))
newdata<-merge(usesoil, userate, by="id")


#==========Diffusion based moisture-respiration relationship I==============================================================================
fh<-function(mvol,po,bd){
  Ds0=2e-9 #aqueous tracer diffusivity at 25 (m2/s)
  Dg0=2.1e-5 #oxygen diffusivity in gas
  fDg<-((po-mvol)/(po))^0.5#(po)^1.5*((po-mvol)/(po))^2.5 #gas phase relative diffusivity
  fDs<-((mvol)/(po))^0.5#(po)^1.5*((mvol)/(po))^2.5 #aqueous phase relative diffusivity
  H_o2<-1.3e-6  #mol/cm3/atm
  hs<-6/(mvol+bd*10)*Ds0*fDs/(0.00002^2) #DOM delivery (mass transfer rate in d-1)
  hg<-6/(mvol+bd*1)*Dg0*fDg*H_o2/(0.00002^2)  #DO delivery (mass transfer rate)
  return (list(fDs,fDg, hs, hg))}
mvol<-newdata$mvol
po<-newdata$porosity
bd<-newdata$bd
fhout<-fh(mvol, po, bd)
hs<-fhout[[3]]
hg<-fhout[[4]]

DiffMM<-function (hs,hg, soc, oxy, kmc, kmg,vmax){
  fm<-newdata$rs^(1/2)  #microbial hydrological sensitivity
  ac<-kmc/(soc)
  bc<-vmax/hs/(soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)
  css<-F1c*kmc/(1-F1c)
  ag<-kmg/(oxy)
  bg<-vmax/hg/(oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  oss<-F1g*kmg/(1-F1g)
  Ft<-fm*F1c*F1g
  return(list(Ft))
}

soc<-newdata$soc_ss
oxy<-newdata$oxy_ss
vmax<-newdata$Vmax
kmo<-newdata$Keffo
kms<-newdata$Keffs

out1<-DiffMM(hs,hg, soc,oxy,1e-7,1e-7,1e-4) #Km_ref
out2<-DiffMM(hs,hg, soc,oxy,kms,kmo,1e-4) #Km_eff


newdata$out1<-out1[[1]]
newdata$out2<-out2[[1]]
newdata$out3<-newdata$rs^(1/2) *fhout[[1]]*fhout[[2]] #linear
newdata$out4<-3.11*newdata$rs-2.42*newdata$rs^2 #empirical 


testplot<-data.frame(rs=newdata$rs)
testplot$Obv<-newdata$ra
testplot$Km<-newdata$out1_s
testplot$Km_eff<-newdata$out2_s
testplot$Linear<-newdata$out3_s
testplot$Empirical<-newdata$out4

#plot(testplot)

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, breaks=20,plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}


pairs(testplot, pch = 19,  cex = 0.5, col=c("#8073ac"),diag.panel=panel.hist,
      upper.panel=NULL)

#=================================result statistics and visulization===================================
ptest<-data.frame(rs=newdata$rs) 
ptest$res<-newdata$out3_s-mpdata$ra

ave<-mean(mpdata$ra)
ssres<-sum(ptest$res^2)
sstot<-sum((mpdata$ra-ave)^2)  
RR<-1-ssres/sstot

simple.fit = lm(res~rs, data=ptest)
summary(simple.fit)

#residual plot
plot1<-ggplot(ptest, aes(x = rs,y = res)) +
  scale_x_continuous(name = expression(paste("Relative saturation")),limits = c(0,1)) +
  scale_y_continuous(name = 'Residual',limits = c(-1,1))+geom_point(col="#642ba6", pch=1 , cex=1.6)+geom_hline(yintercept=0, color="black", linetype="dashed",size=0.5)+
  geom_smooth(method=lm, color="tomato1", size=1, fill="tomato1", se=TRUE)

plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)

pdf("Diff_MMKphy.pdf", width=3.3, height=3)
plot2
dev.off() 

#====residual plot for two data series====
  diff<-data.frame(rs=newdata$rs)
  diff$mmk1<-newdata$out2_s-mpdata$ra
  diff$mmk2<-newdata$out1_s-mpdata$ra
  pdiff<-reshape2::melt(diff, id.vars=("rs"))
  
  plot1<-ggplot(pdiff, aes(x=rs, y=value))+geom_point(aes(color=variable),cex=1,shape=1, alpha=0.9)+
    geom_smooth(aes(color=variable, fill=variable),method = lm, se = TRUE, size=1)+
    geom_hline(yintercept=0, color="black",linetype="dashed",size=0.6)+
    scale_color_manual(values = c("#8073ac","#b35806","#662788"))+
    scale_fill_manual(values = c("#8073ac","#b35806","#662788"))+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(-1,1))+
    theme(text = element_text(size=16)) 
  plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
  print(plot2)
pdf("trait2phy.pdf", width=4.6, height=3)
plot2
dev.off() 



#=============================MM with alternative diffusivity=================
SimMM<-function(mvol, po, soc, oxy, kmc, kmg,vmax){
  fm1<-newdata$rs^(1/2.2)
  fm2<-exp(-0.0002*newdata$mwp) 
  Ds0=6e-9 #aqueous tracer diffusivity at 25 (m2/s)
  Dg0=2.1e-5 #oxygen diffusivity in gas
  fDg<-((po-mvol)/(po))^0.5 #gas phase relative diffusivity
  fDs<-((mvol)/(po))^0.5 #aqueous phase relative diffusivity
  #fDg<-(po)^1.5*((po-mvol)/(po))^2.5 #gas phase relative diffusivity
  #fDs<-(po)^1.5*((mvol)/(po))^2.5 #aqueous phase relative diffusivity
  H_o2<-1.3e-6  #mol/cm3/atm
  hs<-6/(mvol+bd*10)*Ds0*fDs/(0.00002^2) #DOM delivery (mass transfer rate in d-1)
  hg<-6/(mvol+bd*1)*Dg0*fDg*H_o2/(0.00002^2)  #DO delivery (mass transfer rate)
  ac<-kmc/(soc)
  bc<-vmax/hs/(soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)
  ag<-kmg/(oxy)
  bg<-vmax/hg/(oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  ft1<-fm1*vmax*F1c*F1g
  ft2<-fm2*vmax*F1c*F1g
  return (list(ft1,ft2))}
soc<-newdata$soc_ss
oxy<-newdata$oxy_ss
mvol<-newdata$mvol
bd<-newdata$bd
po<-newdata$porosity
out1<-SimMM(mvol, po, soc,oxy,1e-7,1e-7, 1e-4) 
newdata$out1<-out1[[1]]
newdata$out2<-out1[[2]]



#========== Sensitivity analysis (delta) ==================
fh<-function(mvol,po,bd, delta){
  Ds0=6e-9 #aqueous tracer diffusivity at 25 (m2/s)
  Dg0=2.1e-5 #oxygen diffusivity in gas
  fDg<-((po-mvol)/(po))^0.5#(po)^1.5*((po-mvol)/(po))^2.5 #gas phase relative diffusivity
  fDs<-((mvol)/(po))^0.5#(po)^1.5*((mvol)/(po))^2.5 #aqueous phase relative diffusivity
  H_o2<-1.3e-6  #mol/cm3/atm
  hs<-6/(mvol+bd*10)*Ds0*fDs/(delta^2) #DOM delivery (mass transfer rate in d-1)
  hg<-6/(mvol+bd*1)*Dg0*fDg*H_o2/(delta^2)  #DO delivery (mass transfer rate)
  return (list(fDs,fDg, hs, hg))}
mvol<-newdata$mvol
po<-newdata$porosity
bd<-newdata$bd


fhout1<-fh(mvol, po, bd, 1e-5)
fhout2<-fh(mvol, po, bd, 2e-5)
fhout3<-fh(mvol, po, bd, 4e-5)
fhout4<-fh(mvol, po, bd, 8e-5)
fhout5<-fh(mvol, po, bd, 16e-5)

DiffMM<-function (hs,hg, soc, oxy, kmc, kmg,vmax){
  #fm<-exp(-0.0003*newdata$mwp) 
  fm<-newdata$rs^(1/2.5)
  ac<-kmc/(soc)
  bc<-vmax/hs/(soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)
  css<-F1c*kmc/(1-F1c)
  ag<-kmg/(oxy)
  bg<-vmax/hg/(oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  oss<-F1g*kmg/(1-F1g)
  Ft<-fm*vmax*F1c*F1g
  Fl<-fm*vmax*hs*soc/(hs*soc+kmc)*hg*oxy/(hg*oxy+kmg)
  return(list(Ft,Fl))
}

#=====
soc<-newdata$soc_ss
oxy<-newdata$oxy_ss
vmax<-newdata$Vmax

out1<-DiffMM(fhout1[[3]],fhout1[[4]], soc,oxy,1e-7,1e-7,1e-6) #kmo =10uM vmax=0.216 d-1(Allion), 0.1-10 d-1 (MEND)
out2<-DiffMM(fhout2[[3]],fhout2[[4]], soc,oxy,1e-7,1e-7,1e-6) 
out3<-DiffMM(fhout3[[3]],fhout3[[4]], soc,oxy,1e-7,1e-7,1e-6)
out4<-DiffMM(fhout4[[3]],fhout4[[4]], soc,oxy,1e-7,1e-7,1e-6) 
out5<-DiffMM(fhout5[[3]],fhout5[[4]], soc,oxy,1e-7,1e-7,1e-6) 

#=====Sensitivity analysis (Km and Vmax) ====================
hs<-fhout4[[3]]
hg<-fhout4[[4]]
fDs<-fhout4[[1]]
fDg<-fhout4[[2]]

DiffMMt<-function (hs,hg, fDs, fDg,soc, oxy, kmc, kmg,vmax){
  #fm<-exp(-0.0003*newdata$mwp) 
  fm<-newdata$rs^(1/2.5)
  ac<-kmc/(soc)
  bc<-vmax/hs/(soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)
  css<-F1c*kmc/(1-F1c)
  ag<-kmg/(oxy)
  bg<-vmax/hg/(oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  oss<-F1g*kmg/(1-F1g)
  Ft<-fm*vmax*F1c*F1g
  Fl<-fm*vmax*fDs*soc/(fDs*soc+kmc)*fDg*oxy/(fDg*oxy+kmg)
  return(list(Ft,Fl))
}

out1<-DiffMMt(hs,hg,fDs, fDg,  soc,oxy,1e-3,1e-3,1e-6) #kmo =10uM vmax=0.216 d-1(Allion), 0.1-10 d-1 (MEND)
out2<-DiffMMt(hs,hg, fDs, fDg, soc,oxy,1e-5,1e-5,1e-6)
out3<-DiffMMt(hs,hg,fDs, fDg,  soc,oxy,1e-7,1e-7,1e-6) 
out4<-DiffMMt(hs,hg, fDs, fDg, soc,oxy,1e-9,1e-9,1e-6) 
out5<-DiffMMt(hs,hg, fDs, fDg, soc,oxy,1e-11,1e-11,1e-6) 


out1<-DiffMM(hs,hg, soc,oxy,1e-7,1e-7,1e-8) #kmo =10uM vmax=0.216 d-1(Allion), 0.1-10 d-1 (MEND)
out2<-DiffMM(hs,hg, soc,oxy,1e-7,1e-7,1e-7) 
out3<-DiffMM(hs,hg, soc,oxy,1e-7,1e-7,1e-6)
out4<-DiffMM(hs,hg, soc,oxy,1e-7,1e-7,1e-5) 
out5<-DiffMM(hs,hg, soc,oxy,1e-7,1e-7,1e-4) 


newdata$out1<-out1[[2]]
newdata$out2<-out2[[2]]
newdata$out3<-out3[[2]]
newdata$out4<-out4[[2]]
newdata$out5<-out5[[2]]
newdata$out6<-newdata$rs^(1/2.5) *fhout1[[1]]*fhout1[[2]]

testplot<-data.frame(RS=newdata$out6_s)
testplot$Linear<-newdata$out6_s
testplot$Vmax1<-newdata$out1_s
testplot$Vmax2<-newdata$out2_s
testplot$Vmax3<-newdata$out3_s
testplot$Vmax4<-newdata$out4_s
testplot$Vmax5<-newdata$out5_s
plot(testplot)


#=====
lowerfun <- function(data,mapping){
  ggplot(testplot, mapping= mapping)+
    geom_point(color="#542788", cex=0.3,pch=1,stroke=0.3)+
    #geom_hline(yintercept=0, type="dashed", color="black", size=0.2)+
    scale_x_continuous(limits = c(0,1),breaks = seq(0, 1, by = 1))+
    scale_y_continuous(limits = c(0,1),breaks = seq(0, 1, by = 1))+theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}  

plot1<-ggpairs(testplot,
        lower = list(continuous = wrap(lowerfun)),
        upper = NULL,
        diag = NULL)
plot2<-plot1+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=10))
print(plot2)

pdf("Delta_lin.pdf", width=6, height=5)
plot2
dev.off()


