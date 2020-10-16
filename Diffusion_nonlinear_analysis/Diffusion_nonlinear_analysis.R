#This script predicts moisture sensitivity from the combined effect of mass transfer and michealis menten kinetics
#within a hypothetical microbial-substrate system=====
rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(scales)
library(MASS)

#===========create dataset====================
rs<-seq(0.1,0.9,0.01)
po<-0.6
mvol<-po*rs
bd<-(1-po)*2.65
soc<-1e-7 #0.05/12*bd*1e-3 #mol/cm3
oxy<-1e-7 # po*0.21*1.3e-6 #mol/cm3
theta<-mvol/0.8#thaR=0,thaS=0.7
t1<-theta^(-1/0.33)-1#n=1.5, alpha=0.03, m=0.33
t2<-t1^(1/1.5)/0.03
data<-data.frame(rs=rs)
data$mvol<-mvol
data$po<-po
data$bd<-bd
data$mp<-t2
data$clay<-0.25
data$soc<-soc
data$oxy<-oxy


#====Mass transfer and MM kinetics====
fh<-function(mvol,po,bd){
  Ds0=1.4e-9  #aqueous tracer diffusivity at 25 (m2/s)
  Dg0=2.1e-5 #oxygen diffusivity in gas
  fDg<-(po)^1.5*((po-mvol)/(po))^2.5 #gas phase relative diffusivity
  fDs<-(po)^1.5*((mvol)/(po))^2.5 #aqueous phase relative diffusivity
  H_o2<-1.3e-6  #mol/cm3/atm
  hs<-6/(mvol+bd*10)*Ds0*fDs/(0.00001^2) #DOM delivery (mass transfer rate in d-1)
  hg<-6/(mvol+bd*1)*Dg0*fDg*H_o2/(0.00001^2)  #DO delivery (mass transfer rate)
  return (list(fDs,fDg, hs, hg))}

fhout<-fh(mvol, po, bd)
hs<-fhout[[3]]
hg<-fhout[[4]]

#===steady state mass transfer uptake balance solution Ca*====
Lee<-function (hs,hg, soc, oxy, kmc, kmg,vmax){
  ac<-kmc/soc
  bc<-vmax/(hs*soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)#uptake flux
  css<-F1c*kmc/(1-F1c)
  ag<-kmg/oxy
  bg<-vmax/(hg*oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  oss<-F1g*kmg/(1-F1g)
  Ft<-F1c*F1g
  return(list(bc, bg, Ft,F1c,F1g, css, oss))
}

# subplot a--High microbial conversion/High affinity (lower km values)
out1<-Lee(hs,hg, soc,oxy,1e-13,1e-13,1.16e-4) #vmax=10 d-1
out2<-Lee(hs,hg, soc,oxy,1e-12,1e-12,1.16e-4) 
out3<-Lee(hs,hg, soc,oxy,1e-11,1e-11,1.16e-4) 
out4<-Lee(hs,hg, soc,oxy,1e-10,1e-10,1.16e-4) 
out5<-Lee(hs,hg, soc,oxy,1e-9,1e-9,1.16e-4) 
out6<-Lee(hs,hg, soc,oxy,1e-8,1e-8,1.16e-4) 

# subplot c--Low microbial conversion/High affinity (lower km values)
out1<-Lee(hs,hg, soc,oxy,1e-13,1e-13,1.16e-9) #vmax=0.0001 d-1
out2<-Lee(hs,hg, soc,oxy,1e-12,1e-12,1.16e-9) 
out3<-Lee(hs,hg, soc,oxy,1e-11,1e-11,1.16e-9) 
out4<-Lee(hs,hg, soc,oxy,1e-10,1e-10,1.16e-9) 
out5<-Lee(hs,hg, soc,oxy,1e-9,1e-9,1.16e-9) 
out6<-Lee(hs,hg, soc,oxy,1e-8,1e-8,1.16e-9) 

#subplot b--High microbial conversion/Lowaffinity 
out1<-Lee(hs,hg, soc,oxy,1e-6,1e-6,1.16e-4) # vmax=10 d-1
out2<-Lee(hs,hg, soc,oxy,1e-5,1e-5,1.16e-4) 
out3<-Lee(hs,hg, soc,oxy,1e-4,1e-4,1.16e-4) 
out4<-Lee(hs,hg, soc,oxy,1e-3,1e-3,1.16e-4) 
out5<-Lee(hs,hg, soc,oxy,1e-2,1e-2,1.16e-4) 
out6<-Lee(hs,hg, soc,oxy,1e-1,1e-1,1.16e-4) 

#s ubplot d--Low microbial conversion/Low affinity 
out1<-Lee(hs,hg, soc,oxy,1e-6,1e-6,1.16e-9) #vmax=0.0001 d-1
out2<-Lee(hs,hg, soc,oxy,1e-5,1e-5,1.16e-9) 
out3<-Lee(hs,hg, soc,oxy,1e-4,1e-4,1.16e-9) 
out4<-Lee(hs,hg, soc,oxy,1e-3,1e-3,1.16e-9) 
out5<-Lee(hs,hg, soc,oxy,1e-2,1e-2,1.16e-9) 
out6<-Lee(hs,hg, soc,oxy,1e-1,1e-1,1.16e-8) 

sdata<-data.frame(rs=rs)
sdata$km1<-out1[[3]]/max(out1[[3]])
sdata$km2<-out2[[3]]/max(out2[[3]])
sdata$km3<-out3[[3]]/max(out3[[3]])
sdata$km4<-out4[[3]]/max(out4[[3]])
sdata$km5<-out5[[3]]/max(out5[[3]])
sdata$km6<-out6[[3]]/max(out6[[3]])
ptest<-reshape2::melt(sdata, id.vars=("rs"))

#==== "#ffeebf","#ffd19b", "#fdb863", "#fba238", "#e08214", "#c76618", "#b35806", "#904c20", "#7c3b08"==
#===="#d8daeb","#c2c2d9","#b2abd2","#a28fbd","#8073ac","#6a4d9a","#542788","#411464","#2d004b"
#organe figure
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Moisture sensitivity",limits = c(0,1))+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("#ffeebf", "#fdb863",  "#e08214",  "#b35806", "#904c20", "#7c3b08"))+
  scale_fill_manual(values = c("#ffeebf", "#fdb863",  "#e08214",  "#b35806", "#904c20", "#7c3b08"))+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid"))

#purple figure
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Moisture sensitivity",limits = c(0,1))+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("#d8daeb","#b2abd2","#8073ac","#6a4d9a","#542788","#2d004b"))+
  scale_fill_manual(values = c("#d8daeb","#b2abd2","#8073ac","#6a4d9a","#542788","#2d004b"))+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid"))

plot2<-plot1+theme_linedraw()+theme(panel.background = element_rect(fill = "transparent",colour = NA),
                                    plot.background = element_rect(fill = "transparent",colour = NA),
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=16))
print(plot2)

pdf("Figure2c.pdf", width=4.4, height=3)
plot2
dev.off()

#explain Figure 2c
out1<-Lee(hs,hg, soc,oxy,1e-7,1e-7,1.16e-12) 
out2<-Lee(hs,hg, soc,oxy,1e-7,1e-7,1.16e-11)
out3<-Lee(hs,hg, soc,oxy,1e-7,1e-7,1.16e-10) 
out4<-Lee(hs,hg, soc,oxy,1e-7,1e-7,1.16e-9) 
out5<-Lee(hs,hg, soc,oxy,1e-7,1e-7,1.16e-8) 
sdata<-data.frame(rs=rs)
sdata$km1<-log(out1[[2]]*out1[[1]])
sdata$km2<-log(out2[[2]]*out2[[1]])
sdata$km3<-log(out3[[2]]*out3[[1]])
sdata$km4<-log(out4[[2]]*out4[[1]])
sdata$km5<-log(out5[[2]]*out5[[1]])
ptest<-reshape2::melt(sdata, id.vars=("rs"))

plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Moisture sensitivity")+
  geom_line(aes(color=variable,linetype=variable), size=1)+geom_hline(yintercept=0, linetype="dashed", col="gray")+
  scale_color_manual(values = c("#ffeebf", "#fdb863",  "#e08214",  "#b35806", "#904c20", "#7c3b08"))+
  scale_fill_manual(values = c("#ffeebf", "#fdb863",  "#e08214",  "#b35806", "#904c20", "#7c3b08"))+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid"))
plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)

pdf("FigureS2_2c.pdf", width=4.4, height=3)
plot2
dev.off()
#============================================end================================================================



#====Mass transfer distance====
fhrange<-function(mvol,po,bd,delta){
  Ds0=1.4e-9 #(glucose 10^-11) #aqueous tracer diffusivity at 25 (m2/s)
  Dg0=2.1e-5#oxygen diffusivity in water
  fDg<-(po)^1.5*((po-mvol)/(po))^2.5 #gas phase relative diffusivity
  fDs<-(po)^1.5*((mvol)/(po))^2.5 #aqueous phase relative diffusivity
  H_o2<-1.3e-6  #mol/cm3/atm
  hs<-6/(mvol+bd*10)*Ds0*fDs/(delta^2) #DOM delivery (mass transfer rate in d-1)
  hg<-6/(mvol+bd*1)*Dg0*fDg*H_o2/(delta^2)  #DO delivery (mass transfer rate)
  return (list(fDs,fDg, hs, hg))}

fhout1<-fhrange(mvol, po, bd, 1e-5) #10
fhout2<-fhrange(mvol, po, bd, 5e-5) #50
fhout3<-fhrange(mvol, po, bd, 1e-4) #100
fhout4<-fhrange(mvol, po, bd, 5e-4) #500

hs1<-fhout1[[3]]
hg1<-fhout1[[4]]
hs2<-fhout2[[3]]
hg2<-fhout2[[4]]
hs3<-fhout3[[3]]
hg3<-fhout3[[4]]
hs4<-fhout4[[3]]
hg4<-fhout4[[4]]
#===use ss solution ====
Lee<-function (hs,hg, soc, oxy, kmc, kmg,vmax){
  ac<-kmc/soc
  bc<-vmax/(hs*soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)#uptake flux
  css<-F1c*kmc/(1-F1c)
  ag<-kmg/oxy
  bg<-vmax/(hg*oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  oss<-F1g*kmg/(1-F1g)
  Ft<-F1c*F1g
  return(list(ac, bc, Ft))
}

#varying vmax
out1<-Lee(hs4,hg4, soc,oxy,1e-6,1e-6,1.16e-9) #vmax=0.0001 d-1
out2<-Lee(hs4,hg4, soc,oxy,1e-6,1e-6,1.16e-8) #vmax=0.001 d-1
out3<-Lee(hs4,hg4, soc,oxy,1e-6,1e-6,1.16e-7) #vmax=0.01 d-1
out4<-Lee(hs4,hg4, soc,oxy,1e-6,1e-6,1.16e-6) #vmax=0.1 d-1
out5<-Lee(hs4,hg4, soc,oxy,1e-6,1e-6,1.16e-5) #vmax=1 d-1 
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[3]]/max(out1[[3]])
sdata$km2<-out2[[3]]/max(out2[[3]])
sdata$km3<-out3[[3]]/max(out3[[3]])
sdata$km4<-out4[[3]]/max(out4[[3]])
sdata$km5<-out5[[3]]/max(out5[[3]])
ptest<-reshape2::melt(sdata, id.vars=("rs"))
#orange lines
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Moisture sensitivity",limits=c(0,1))+#,trans = log_trans(), breaks = base_breaks())+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("#7c3b08","#b35806","#e08214", "#fdb863","#fee0b6"))+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid"))
plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)
pdf("fw_vmax_delta4.pdf", width=4.4, height=3)
plot2
dev.off()


#varying km
out1<-Lee(hs2,hg2, soc,oxy,1e-8,1e-8,1.16e-7) #km=1uM 
out2<-Lee(hs2,hg2, soc,oxy,1e-7,1e-7,1.16e-7) #km=10uM
out3<-Lee(hs2,hg2, soc,oxy,1e-6,1e-6,1.16e-7) #km=100uM
out4<-Lee(hs2,hg2, soc,oxy,1e-5,1e-5,1.16e-7) #km=1000uM
out5<-Lee(hs2,hg2, soc,oxy,1e-4,1e-4,1.16e-7) #km=10000uM
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[3]]/max(out1[[3]])
sdata$km2<-out2[[3]]/max(out2[[3]])
sdata$km3<-out3[[3]]/max(out3[[3]])
sdata$km4<-out4[[3]]/max(out4[[3]])
sdata$km5<-out5[[3]]/max(out5[[3]])
ptest<-reshape2::melt(sdata, id.vars=("rs"))
#purple lines
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Moisture sensitivity",limits = c(0,1))+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))+
  scale_fill_manual(values = c("#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid"))
plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)
pdf("fw_km_delta2.pdf", width=4.4, height=3)
plot2
dev.off()
#============================================end================================================================

#===========================varying C depleting or supply=======================================================
fh<-function(mvol,po,bd){
  Ds0=6e-10 #aqueous tracer diffusivity at 25 (m2/s)
  Dg0=2.1e-5 #oxygen diffusivity in gas
  fDg<-(po)^1.5*((po-mvol)/(po))^2.5 #gas phase relative diffusivity
  fDs<-(po)^1.5*((mvol)/(po))^2.5 #aqueous phase relative diffusivity
  H_o2<-1.3e-6  #mol/cm3/atm
  hs<-6/(mvol+bd*10)*Ds0*fDs/(0.00001^2) #DOM delivery (mass transfer rate in d-1)
  hg<-6/(mvol+bd*1)*Dg0*fDg*H_o2/(0.00001^2)  #DO delivery (mass transfer rate)
  return (list(fDs,fDg, hs, hg))}

fhout<-fh(mvol, po, bd)
hs<-fhout[[3]]
hg<-fhout[[4]]

Leedep<-function (hs,hg, soc, oxy, kmc, kmg,vmax,delta){
  fw<-1*data$rs
  ac<-kmc/(soc-delta*fw*soc)
  bc<-vmax/hs/(soc-delta*fw*soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)
  css<-F1c*kmc/(1-F1c)
  ag<-kmg/(oxy)
  bg<-vmax/hg/(oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  oss<-F1g*kmg/(1-F1g)
  Ft<-F1c*F1g
  Fcss<-hs*(soc-delta*fw*soc-css)
  mass<-hs*(soc-delta*fw*soc)
  cratio<-css/soc
  return(list(bc,Fcss, Ft, mass,cratio))
}

Leesup<-function (hs,hg, soc, oxy, kmc, kmg,vmax,delta){
  fw<-1*data$rs
  ac<-kmc/(soc+delta*fw*soc)
  bc<-vmax/hs/(soc+delta*fw*soc)
  t1c<-(1-4*bc/(1+ac+bc)^2)^0.5
  F1c<-(1+ac+bc)/2/bc*(1-t1c)
  css<-F1c*kmc/(1-F1c)
  ag<-kmg/(oxy)
  bg<-vmax/hg/(oxy)
  t1g<-(1-4*bg/(1+ag+bg)^2)^0.5
  F1g<-(1+ag+bg)/2/bg*(1-t1g)
  oss<-F1g*kmg/(1-F1g)
  Ft<-F1c*F1g
  Fcss<-hs*(soc+delta*fw*soc-css)
  mass<-hs*(soc+delta*fw*soc)
  cratio<-css/soc
  return(list(bc,Fcss, Ft,mass, cratio))
}


#c delpetiong linear 
out1<-Leedep(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0) #kmo =10uM vmax=0.216 d-1(Allion), 0.1-10 d-1 (MEND)
out2<-Leedep(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.2) 
out3<-Leedep(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.4) 
out4<-Leedep(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.6) 
out5<-Leedep(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.8) 
out6<-Leesup(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0) #kmo =10uM vmax=0.216 d-1(Allion), 0.1-10 d-1 (MEND)
out7<-Leesup(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.2) 
out8<-Leesup(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.4) 
out9<-Leesup(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.6) 
out10<-Leesup(hs,hg, soc,oxy,1e-5,1e-5,1.16e-7,0.8) 

out1<-Leedep(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0) #kmo =10uM vmax=0.216 d-1(Allion), 0.1-10 d-1 (MEND)
out2<-Leedep(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.2) 
out3<-Leedep(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.4) 
out4<-Leedep(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.6) 
out5<-Leedep(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.8) 
out6<-Leesup(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0) #kmo =10uM vmax=0.216 d-1(Allion), 0.1-10 d-1 (MEND)
out7<-Leesup(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.2) 
out8<-Leesup(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.4) 
out9<-Leesup(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.6) 
out10<-Leesup(hs,hg, soc,oxy,1e-9,1e-9,1.16e-7,0.8)

#====plotting against reference run
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[3]]/max(out1[[3]])
sdata$km2<-out2[[3]]/max(out1[[3]])
sdata$km3<-out3[[3]]/max(out1[[3]])
sdata$km4<-out4[[3]]/max(out1[[3]])
sdata$km5<-out5[[3]]/max(out1[[3]])
sdata$km7<-out7[[3]]/max(out1[[3]])
sdata$km8<-out8[[3]]/max(out1[[3]])
sdata$km9<-out9[[3]]/max(out1[[3]])
sdata$km10<-out10[[3]]/max(out1[[3]])
ptest<-reshape2::melt(sdata, id.vars=("rs"))
#====plotting rescaled pattern
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[3]]/max(out1[[3]])
sdata$km2<-out2[[3]]/max(out2[[3]])
sdata$km3<-out3[[3]]/max(out3[[3]])
sdata$km4<-out4[[3]]/max(out4[[3]])
sdata$km5<-out5[[3]]/max(out5[[3]])
sdata$km7<-out7[[3]]/max(out7[[3]])
sdata$km8<-out8[[3]]/max(out8[[3]])
sdata$km9<-out9[[3]]/max(out9[[3]])
sdata$km10<-out10[[3]]/max(out10[[3]])
ptest<-reshape2::melt(sdata, id.vars=("rs"))

plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Relative activity",limits=c(0,2))+#,trans = log_trans(), breaks = base_breaks())+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("gray30", "#b2abd2","#8073ac","#542788","#411464","#fdb863",  "#e08214", "#b35806", "#904c20","#b2abd2","#8073ac","#542788","#411464"))+
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid","solid","solid","dashed","dashed"))

plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)
pdf("lowAff.pdf", width=4.6, height=3)
plot2
dev.off()

#=====relative steady state concentration
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[5]]
sdata$km2<-out2[[5]]
sdata$km3<-out3[[5]]
sdata$km4<-out4[[5]]
sdata$km5<-out5[[5]]
sdata$km7<-out7[[5]]
sdata$km8<-out8[[5]]
sdata$km9<-out9[[5]]
sdata$km10<-out10[[5]]

ptest<-reshape2::melt(sdata, id.vars=("rs"))
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Ca*/Csoil")+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("gray30", "#b2abd2","#8073ac","#542788","#411464","#fdb863",  "#e08214", "#b35806", "#904c20","#b2abd2","#8073ac","#542788","#411464"))+
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid","solid","solid","dashed","dashed"))
plot3<-plot1+annotation_logticks()
plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)
#====C transport with diffusion control only
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[4]]
sdata$km2<-out2[[4]]
sdata$km3<-out3[[4]]
sdata$km4<-out4[[4]]
sdata$km5<-out5[[4]]
sdata$km7<-out7[[4]]
sdata$km8<-out8[[4]]
sdata$km9<-out9[[4]]
sdata$km10<-out10[[4]]

ptest<-reshape2::melt(sdata, id.vars=("rs"))
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Fmass",trans = log_trans(), breaks =c(1e-9,1e-8,1e-7),limits = c(1e-10,1e-6))+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("gray30", "#b2abd2","#8073ac","#542788","#411464","#fdb863",  "#e08214", "#b35806", "#904c20","#b2abd2","#8073ac","#542788","#411464"))+
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid","solid","solid","dashed","dashed"))
plot3<-plot1+annotation_logticks()
plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)
#====C transport with conbined action of diffusion and nonlinear kinetics
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[2]]
sdata$km2<-out2[[2]]
sdata$km3<-out3[[2]]
sdata$km4<-out4[[2]]
sdata$km5<-out5[[2]]
sdata$km7<-out7[[2]]
sdata$km8<-out8[[2]]
sdata$km9<-out9[[2]]
sdata$km10<-out10[[2]]

ptest<-reshape2::melt(sdata, id.vars=("rs"))
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Fmass",trans = log_trans(), breaks =c(1e-9,1e-8,1e-7),limits = c(1e-10,1e-6))+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  scale_color_manual(values = c("gray30", "#b2abd2","#8073ac","#542788","#411464","#fdb863",  "#e08214", "#b35806", "#904c20","#b2abd2","#8073ac","#542788","#411464"))+
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid","solid","solid","dashed","dashed"))
plot3<-plot1+annotation_logticks()
plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)


#===Microbial conversion capacity
sdata<-data.frame(rs=rs)
sdata$km1<-out1[[1]]
sdata$km2<-out2[[1]]
sdata$km3<-out3[[1]]
sdata$km4<-out4[[1]]
sdata$km5<-out5[[1]]
sdata$km7<-out7[[1]]
sdata$km8<-out8[[1]]
sdata$km9<-out9[[1]]
sdata$km10<-out10[[1]]

ptest<-reshape2::melt(sdata, id.vars=("rs"))
plot1<-ggplot(ptest, aes(x = rs,y=value, color=variable)) + 
  scale_x_continuous(name = "Relative saturation",limits = c(0,1)) +
  scale_y_continuous(name = "Microbial conversion capacity",trans = log_trans(), breaks = c(1,10,100),labels = scales::trans_format("log10", function(x) 10^x))+
  geom_line(aes(color=variable,linetype=variable), size=1)+
  geom_hline(yintercept=1, linetype="dashed", color="black", size=0.6)+
  scale_color_manual(values = c("gray30", "#b2abd2","#8073ac","#542788","#411464","#fdb863",  "#e08214", "#b35806", "#904c20","#b2abd2","#8073ac","#542788","#411464"))+
  scale_linetype_manual(values=c("dashed","solid","solid","solid","solid","solid","solid","solid","solid","dashed","dashed"))
plot3<-plot1+annotation_logticks()
plot2<-plot1+theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(text = element_text(size=16))
print(plot2)

