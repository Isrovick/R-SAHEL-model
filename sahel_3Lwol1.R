rm(list=ls())
library("xlsx")
library("ggplot2")

SRC=file.choose()
DIR= file.path(dirname(SRC))
data  <- read.xlsx(SRC, sheetIndex = 1)


# CALCULO DE EL ULTIMO DIA SIN LLUVIA(SUBSTITUTO DE ORIGINAL DSLR)
DSLR_ <- function(TIME,RAINT,DATEB){
   d=0
        DATEP  = AMOD(DATEB + TIME - d + 364 , 365 ) +1  #Experimental tesista
   while(RAINT[DATEP]<0.5){
       d=d+1
       
        DATEP  = AMOD(DATEB + TIME - d + 364 , 365 ) + 1  #Experimental tesista
   }

   return(d+1)
}

#####################  CSMP FUNCTIONS  (partial)  ###############################

# Input switch (RELAY)
INSW <- function(X1,X2,X3){
    
    if( X1 < 0 ){
      return(X2)
    }
    else{
      return(X3)
    }
}

#LIMITER
LIMIT <-function(P1,P2,X){
    
    if(X < P1){
        return(P1)
    } else if( X > P2){
        return(P2)
    } else {
        return(X)
    }

}

INTGRL <- function(IC,X,T,C){

  f <- function(x) X
  F <-Vectorize(f)
  Y <- integrate(F , 0, T)

  return(Y$value + IC)
}

# modulo function
AMOD <-function(X,P){
  
  Y<- X
  n<- 0;
  while (!(Y<P && Y>=0)) {
    n<- n+1
    Y<- X-n*P
  }
  return(Y)
}

#####################  APPENDIX B, T12  partial  ###############################

# checks whether X is between limits.  in FUNCTIONS,  SUBROUTINES. 
SUERRM <- function(MNR,X,XMIN,XMAX,NUNIT){

  if(((X < XMIN*0.99) && (XMIN != -99.)) || ((X > XMAX*1.01) && (XMAX != -99.))){
    sprintf("fatal error in variable or parameter value message number, value,  minimum and maximum: %f %f %f %f",MNR,X,XMIN,XMAX)
    stop("fatal error suerrm...")
  } 
}

# astronomical  standard computations.  used  in LID,  L1Q.
# SUBROUTINE SUASTC (DATE,LAT,RDTM,RDTC,FRDIF,COSLD, SINLD,DSINBE,SOLC,DLA)
SUASTC <- function(DATE,LAT,RDTM){
  
  PI <- 3.1415926  
  RAD <- 0.0174533

  DEC =-asin(sin(23.45*RAD)*cos(2.*PI*(DATE+10.)/365.))
  COSLD   =cos(DEC)*cos(LAT*RAD)
  SINLD   =sin(DEC)*sin(LAT*RAD)
  AOB     =SINLD/COSLD
    SUERRM(2.1,DATE,0.,365.,6.)
    SUERRM(2.2,AOB,-1.0,1.0,6.)
  DLA     =12.*(1.+2.0*asin(AOB)/PI)
  DSINBE  =3600.*(DLA*(SINLD+0.4*(SINLD*SINLD+COSLD*COSLD*0.5))+12.0*COSLD*(2.0+3.0*0.4*SINLD)*sqrt(1.-AOB*AOB)/PI) 
  DSINB   =3600.*(SINLD*DLA+24./PI*COSLD*sqrt(1.-AOB**2))
  SOLC    =1370.*(1.0+0.033*cos(2.*PI*DATE/365.))
  RDTC    =SOLC*DSINB
    SUERRM(2.3,RDTM,0.,RDTC,6.)
  ATMTR  =RDTM/RDTC
  if(ATMTR > 0.75)  FRDIF  =0.23
  if(ATMTR <= 0.75 && ATMTR > 0.35)   FRDIF  =1.33-1.46*ATMTR
  if(ATMTR <= 0.35 && ATMTR > 0.07)   FRDIF  =1.-2.3*(ATMTR-0.07)**2
  if(ATMTR <= 0.07)  FRDIF  =1.00

  result=list("RDTM" = RDTM,"RDTC" =RDTC,"FRDIF" =FRDIF,"COSLD" =COSLD,"SINLD" =SINLD,"DSINBE" =DSINBE,"SOLC" =SOLC,"DLA" =DLA)

  return(result)
}

# computes  daylength,  daily total radiation clear.  in LlD,  L1Q.
SUASTR <- function(DATE,LAT,RDTM){

  INSP<- -4.0   
  PI<-3.1415926  
  RAD<-0.0174533
  
  suastc = SUASTC(DATE,LAT,RDTM)
  SINLD = suastc$SINLD
  COSLD = suastc$COSLD
  DLP   =12.*(PI+2.*asin((-sin(INSP*RAD)+SINLD)/COSLD))/PI

  result=list("RDTC"= suastc$RDTC,"DLA"= suastc$DLA, "DLP"= DLP)

  return(result)
}

# vapour pressure  (kPa)  relation to temperature.  in LlQ, 
FUVP <-function(TP){
  fuvp    =0.100*6.11*exp(17.47*TP/(TP+239.))
  return(fuvp)
}

#  potential  evapotranspiration rates  crop,  soil.  used in L2C.
SUEVTR <-function(RDTC,RDTM,RF,FRD,TPAD,VPA,RSL,RSB,RST,HUAA){

  SUERRM(3.1,FRD,0.,1.,6.)
  #VPAS  =FUVP(TPAD)
  VPAS   = 0.01 * HUAA * FUVP(TPAD)  #HUAA in relative percent 6.1.5 Table  29. 
  SUERRM(3.2,VPAS,0.0,12.55,6.)
  SUERRM(3.3,VPA,0.0,VPAS,6.)
  SLOPE  =4158.6*10.*VPAS/(TPAD+239.)**2 
  APSCH  =0.67*(RSB+RST+RSL)/(RSB/0.93+RST)
  RLWI   =4.8972E-3*(TPAD+273.)**4*(0.618+0.0365*sqrt(10.*VPA)) 
  RLWO   =4.8972E-3*1.00*(TPAD+273.)**4 
  RDTN   =RDTM*(1.-RF)-(RLWO-RLWI)*(RDTM/(0.75*RDTC))*FRD 
  EVPR   =0.001*RDTN*SLOPE/((SLOPE+APSCH)*2390.)
  DRYP   =(VPAS-VPA)*10.*1200./(RSB+RST) *FRD 
  EVPD   =86400.*0.001*DRYP/((SLOPE+APSCH)*2390.)
  #browser()
  result <-list("EVPR" = EVPR, "EVPD" = EVPD)
  return(result)

}

#calculates  canopy resistance  upper layers. in L2C.
FURSC <-function(WDS,ALV,PLHT,ZREF){

  ZR    = max(c(ZREF,PLHT+1.))
  D     = max(c(0.1,0.63*PLHT))
  ZNOT  = max(c(0.05,0.1*PLHT)) 
  ALVX  = max(c(1.,ALV)) 
  WDSX  = max(c(0.2,WDS))
  FURSC_  =0.74*(log((ZR-D)/ZNOT))**2/(0.16*WDSX)*ALVX

  return(FURSC_)
}

# calculates windspeed near soil  surface.  in L2C. 
FUWRED <-function(WDLV,ALV,PLHT,WDS){

  PLHTX   = max(c(0.05,PLHT)) 
  ALVX    = max(c(0.01,ALV))
  MIXL    = sqrt(1.2732*max(c(0.005,WDLV))/(ALVX/PLHTX)) 
  A       = sqrt(0.2*ALVX*PLHTX/(2.*MIXL*0.5)) 
  FUWRED_  = max(c(0.2,WDS))*exp(-A*(1.0-0.05/PLHTX)) 

  return(FUWRED_)
}

# Computes  reduction of water uptake,  used in L2SU,  L2SS. 
FUWS <-function(TRC,ALV,WCL,WSSC,WFSC,WCWP,WCFC,WCST){

  A <-0.76
  B <-0.15
  ALVMAX <-2
  FUWSX<-0
  if(WCL  <=  WCFC){
    SDPF   =1./(A+B*ALVMAX*TRC/(ALV+1.E-10))-(1.-WSSC)*0.4 
    if(WSSC < 0.6){
      SDPF =SDPF+0.025*min(c(0.,ALVMAX*TRC/(ALV+1.E-10)-6.))/(1.+5.*WSSC+4.*WSSC*WSSC)
    }
    WCX    =WCWP+(WCFC-WCWP)*(1.00-min(c(1.,max(0.,SDPF)))) 
    FUWSX   =(WCL-WCWP)/(WCX-WCWP+1.E-10)
  }
  else{
    FUWSX  =1.-(1.-WFSC)*(WCL-WCFC)/(WCST-WCFC+1.E-10)
  }

  fuws = min(c(1.,max(c(0.,FUWSX))))

  return(fuws)
}

#check on soil water balance. used in L2SU, L2SS. 
FUWCHK <-function(CKWFL,CKWIN,TIME,IDATE){

  FUWCHK_=2.0*(CKWIN-CKWFL)/(CKWIN+CKWFL+1.E-10)
  if((abs(FUWCHK_) > 0.01)  && abs(CKWIN) > 0.2){
    
    print(c(error="Error  in water balance,  please  check: CKWRD (CKWIN, CKWFL, AT TIME)",IDATE=IDATE,ckwin=CKWIN,cwfl=CKWFL,time=TIME))
    #sprintf("Error  in water balance,  please  check: CKWRD CKWIN=%f CKWFL=%f AT TIME %f",CKWFL,CKWIN,TIME)
    
    return(1)
  } 
  return (0)
}

#####################  INPUT DATA  ###############################

SOILWATER=data$Soilwater #Experimental Tesista
STNNS1=data$stoniness1[1] #Experimental Tesista
STNNS2=data$stoniness2[1] #Experimental Tesista
STNNS3=data$stoniness3[1] #Experimental Tesista
STNNS4=data$stoniness4[1] #Experimental Tesista
SDATES=data$Fecha         #Experimental Tesista

# # # # # # # # #  Listing  .3  INPUT DATA (Partial)  # # # # # # # # #

##RUN CONTROL AND OUTPUT
#METHOD RECT
DELT  = data$DELT[1]
TIME  = data$TIME[1]
FINTIM= data$FINTIM[1]
PRDEG = data$PRDEG[1]
OUTDEG= data$OUTDEG[1]

# # # # # # # # #  Listing  .5  INPUT DATA (Partial)  # # # # # # # # #

##WATER  RELATIONS  AND ROOT  GROWTH; TABLES  22,24,25 
WSSC  = data$WSSC[1]
WFSC  = data$WFSC[1]
##PHENOLOGICAL  DEVELOPMENT
SLC = data$SLC[1]
SSC = data$SSC[1]
WDLV  = data$WDLV[1]
##INITIALIZATION
DATEB = data$DATEB[1]
FIN_DS= data$FIN_DS[1]
FIN_CELVN= data$FIN_CELVN[1]
FIN_TAVP= data$FIN_TAVP[1]

# # # # # # # # #  Listing  .7  INPUT DATA (Partial)  # # # # # # # # #
# # # # # # # # #  Listing  .8  INPUT DATA (Partial)  # # # # # # # # #
# # # # # # # # #  Listing  .10 INPUT DATA (Partial)  # # # # # # # # #

#DATA  FOR  MODULE  L2SU
TKL1  = (data$tkllitter[1])*0.001 *(1-STNNS1)
TKL2  = (data$tklmin1[1])  *0.001 *(1-STNNS2)
TKL3  = (data$tklmin2[1])  *0.001 *(1-STNNS3)
TKL4  = (data$tklmin3[1])  *0.001 *(1-STNNS4)

WCFC1 = data$WCFC1[1]
WCFC2 = data$WCFC2[1]
WCFC3 = data$WCFC3[1]
WCFC4 = data$WCFC4[1]

WCWP1 = data$WCWP1[1]
WCWP2 = data$WCWP2[1]
WCWP3 = data$WCWP3[1]
WCWP4 = data$WCWP4[1]

WCAD1 = data$WCAD1[1]
WCAD2 = data$WCAD2[1]
WCAD3 = data$WCAD3[1]
WCAD4 = data$WCAD4[1]

WCST1 = data$WCST1[1]
WCST2 = data$WCST2[1]
WCST3 = data$WCST3[1]
WCST4 = data$WCST4[1]

WCLI1 = data$initialmoisturel1[1] #INCON
WCLI2 = data$initialmoisturel2[1] #INCON
WCLI3 = data$initialmoisturel3[1] #INCON
WCLI4 = data$initialmoisturel4[1] #INCON

##SURFACE  AND  OTHER  SOIL  CHARACTERISTICS 
FRNOF = data$FRNOF[1]
RFSD  = data$RFSD[1]
WDCL  = data$WDCL[1]
EES   = data$EES[1]

# # # # # # # # #  Listing  .11 INPUT DATA (Partial)  # # # # # # # # #

LAT   = data$lat[1]
RDUCF = data$RDUCF[1]
RDTMT = as.double(data$RDTMT)
TPHT  = as.double(data$tmmx)
TPLT  = as.double(data$tmmn)
RAINT = as.double(data$Rain)
HUAAT = as.double(data$HUAAT)
WDST  = as.double(data$WDST)

#####################  INITAL  ###############################

# # # # # # # # #  Listing  .3  INITIAL (Partial)  # # # # # # # # #
# # # # # # # # #  Listing  .5  INITIAL (Partial)  # # # # # # # # #
# # # # # # # # #  Listing  .7  INITIAL (Partial)  # # # # # # # # #

IDATE = DATEB
TPSI  = (TPLT[IDATE]+TPHT[IDATE])/2.

# # # # # # # # #  Listing  .8  INITIAL (Partial)  # # # # # # # # #


WL1I  = WCLI1*TKL1*1.E4 
WL2I  = WCLI2*TKL2*1.E4 
WL3I  = WCLI3*TKL3*1.E4
WL4I  = WCLI4*TKL4*1.E4 
TKLT  = TKL1+TKL2+TKL3+TKL4

# # # # # # # # #  Listing  .10 INITIAL (Partial)  # # # # # # # # #
# # # # # # # # #  Listing  .11 INITIAL (Partial)  # # # # # # # # #

#################   EXPERIMENTAL TESISTA  ###############################

TIMES <- c(1:365)
WCL1T <-c()
WCL2T <-c()
WCL3T <-c()
WCL4T <-c()
WCLER <-c()
TPAVT <-c()
RAINT_ <-c()
EXPV <-c()
EVSWT <- c()
TRWTA <- c()
cont=0
#####################  DYNAMIC  ###############################

for(TIME in  TIMES){

    WL2=WL2I
    WL3=WL3I
    WL4=WL4I

    WCL2=WCLI2
    WCL3=WCLI3
    WCL4=WCLI4

    # # # # # # # # #  Listing  .3 DYNAMIC (Partial)  # # # # # # # # #

    ##LEAF AREA 
    #Explanation  in  section  3.3
    #ALV    =INTGRL(ALVI,GLA-LLA+GSA)
    ALV  = 0.0 # For fallow soil
    
    ##WEATHER DATA AND TIME 
    ##Explanation  in chapter  6  and  section  3.4 
    
    DATE  = AMOD(DATEB + TIME + 364 , 365 ) + 1
    IDATE = DATE
    
    RDTM  = RDTMT[IDATE]*RDUCF
    suastr= SUASTR(IDATE,LAT,RDTM)
    
    RDTC  = suastr$RDTC
    DLA   = suastr$DLA
    DLP   = suastr$DL

    TPAV  = (TPLT[IDATE]+TPHT[IDATE])/2.
    TPAD  = (TPHT[IDATE]+TPAV)/2.
    

    TPAVT[IDATE]<-TPAV               #Experimental tesista
    EXPV[IDATE]<-SOILWATER[IDATE]  #Experimental tesista

    # # # # # # # # #  Listing  .5 DYNAMIC (Partial)  # # # # # # # # # 
    # # # # # # # # #  Listing  .7 DYNAMIC (Partial)  # # # # # # # # #

    ##ROOTED  DEPTH  AND  CROP  HEIGHT 
    ##Explanation in section 4.2
    #ZRT   =INTGRL(ZRTI,GZRT*AND(ZRTM-ZRT,l.0-DS))
    ZRT  = 0.0 # For fallow soil
    #PLHT  =AFGEN(PLHTT,DS)
    PLHT = 0.0 # For fallow soil
    
    ##POTENTIAL TRANSPIRATION AND  DIFFUSION RESISTANCES  CANOPY 
    ##Explanation in sections 4.1,  4.4 
    #TRC    =TRCPR*(l.-EXP(-0.5*ALV))+TRCPD*AMIN1(2.5,ALV)
    TRC  = 0.0 # For fallow soil
    TRRM  = TRC/(ZRT+1.E-10)


    ##EXTRA  WEATHER  DATA 
    ##Explanation in section  6.1,  5.1
    WDSAV = max(c(0.2,WDST[IDATE]))
    #VPA   = min(c(FUVP(TPAD),HUAAT[IDATE])) #HUAA in Kpa
    VPA   = 0.01 * HUAAT[IDATE] * FUVP(TPAV)  #HUAA in relative percent 6.1.5 Table  29. 
    RAIN  = RAINT[IDATE]
    RAINT_[IDATE]=RAIN
    #DSLR  = INTGRL(1.,INSW(RAINT[IDATE+1]-0.5,1.,1.00001-DSLR)/DELT,TIME,"DSLR")
    DSLR = DSLR_(TIME-1,RAINT,DATEB)#Experimental Tesista
   
    ##POTENTIAL  EVAPORATION  SOIL 
    ##Explanation in section 5.1 
    WDSS  = FUWRED(WDLV,ALV,PLHT,WDSAV)
    RSBS  = 172.*sqrt(WDCL/WDSS) 
    RSTS  = FURSC(WDSAV,1.,0.1*PLHT,0.63*PLHT)
    RFS   = RFSD*(1.-0.5*WCL2/WCST2)
    suevtr= SUEVTR(RDTC,RDTM,RFS,1.00,TPAV,VPA,0.00,RSBS,RSTS,HUAAT[IDATE])
    EVSPR = suevtr$EVPR
    EVSPD = suevtr$EVPD
    EVSC  = EVSPR*exp(-0.5*ALV)+EVSPD

    # # # # # # # # #  Listing  .8 DYNAMIC  (Complete)  # # # # # # # # #

    ##ACTUAL  TRANSPIRATION  (WATER  UPTAKE)
    ##Explanation  in  section  5.2 

    WSE2  = FUWS(TRC,ALV,WCL2,WSSC,WFSC,WCWP2,WCFC2,WCST2) 
    WSE3  = FUWS(TRC,ALV,WCL3,WSSC,WFSC,WCWP3,WCFC3,WCST3) 
    WSE4  = FUWS(TRC,ALV,WCL4,WSSC,WFSC,WCWP4,WCFC4,WCST4) 

    ZRT2  = LIMIT(0.,TKL2,ZRT-TKL1) 
    ZRT3  = LIMIT(0.,TKL3,ZRT-TKL1-TKL2)
    ZRT4  = LIMIT(0.,TKL4,ZRT-TKL1-TKL2-TKL3)

    TRWL2 = TRRM*WSE2*ZRT2 
    TRWL3 = TRRM*WSE3*ZRT3
    TRWL4 = TRRM*WSE4*ZRT4
    TRW   = TRWL2+TRWL3+TRWL4
    TRWT  = TRW
    TRWTA[IDATE] = TRW

    ##GROWTH  ROOTED  DEPTH
    ##Explanation  in  subsection  4.2.3
    WSERT = INSW(ZRT-TKL2,WSE2,INSW(ZRT-TKL2-TKL3,WSE3,WSE4))
    #TERT   =AFGEN(PLMTT,TPS)

    ##EVAPORATION
    ##Explanation  in  section  5.2
    
    WLFL2 = RAIN*(1.0-FRNOF)#sorted
    
    EVSH  = min(EVSC, (WL2*0.0001-WCAD2*TKL2)*1000./DELT+WLFL2)
    #EVSD  = min(EVSC, 0.6*EVSC*(sqrt(DSLR)-sqrt(DSLR-1))+WLFL1)
    EVSD  = min(EVSC, EVSC*(sqrt(DSLR)-sqrt(DSLR-1))+WLFL2)#Experimental tesista
    EVSW  = INSW(DSLR-1.1,EVSH,EVSD)
    EVSWT[IDATE]= EVSW #Experimental tesista

    FEVL2 <- max(WL2-WCAD2*TKL2*1.E4,0.)*exp(-EES*(TKL1+(0.25*TKL2)))
    FEVL3 <- max(WL3-WCAD3*TKL3*1.E4,0.)*exp(-EES*(TKL1+TKL2+(0.25*TKL3)))
    FEVL4 <- max(WL4-WCAD4*TKL4*1.E4,0.)*exp(-EES*(TKL1+TKL2+TKL3+(0.25*TKL4)))

    FEVLT = FEVL2+FEVL3+FEVL4 
    if(FEVLT <= 0 ) FEVLT=1 #Experimental tesista
    

    EVSW2 = EVSW*(FEVL2/FEVLT) 
    EVSW3 = EVSW*(FEVL3/FEVLT)
    EVSW4 = EVSW*(FEVL4/FEVLT)

    ##AVAILABLE  AND  TOTAL  SOIL  WATER
    ##Explanation  in  subsection  4.2.3
    WLFL3 = max(0.,WLFL2-(WCFC2*TKL2*1000.-WL2*0.10)/DELT)
    WLFL4 = max(0.,WLFL3-(WCFC3*TKL3*1000.-WL3*0.10)/DELT)
    WLFL5 = max(0.,WLFL4-(WCFC4*TKL4*1000.-WL4*0.10)/DELT)
    

    
    WL2   = WL2I+(WLFL2-WLFL3-EVSW2-TRWL2)*10.0
    WL3   = WL3I+(WLFL3-WLFL4-EVSW3-TRWL3)*10.0
    WL4   = WL4I+(WLFL4-WLFL5-EVSW4-TRWL4)*10.0
    
    WCUM  = (WL2+WL3+WL4)/10.0 

    WCL2  = WL2/(TKL2*1.E4) 
    WCL3  = WL3/(TKL3*1.E4)
    WCL4  = WL4/(TKL4*1.E4)
  
        ##WATER  BALANCE  CHECK 
    ##Explanation  in  section  5.4
    
    CKWFL =(WLFL2-EVSW-TRW-WLFL5)*10. 
    CKWIN = WL2-WL2I+WL3-WL3I+WL4-WL4I
    CKWRD = FUWCHK(CKWFL,CKWIN,TIME,IDATE)
    if(CKWRD) cont=cont+1
    # # # # # # # # #  Listing  .10 DYNAMIC (Partial)  # # # # # # # # #
    # # # # # # # # #  Listing  .11 DYNAMIC (Partial)  # # # # # # # # #   

    WCL2T[IDATE]<-WCL2
    WCL3T[IDATE]<-WCL3
    WCL4T[IDATE]<-WCL4
    WCLER[IDATE]<-CKWRD
    

    WL2I=WL2
    WL3I=WL3
    WL4I=WL4

    WCLI2=WCL2
    WCLI3=WCL3
    WCLI4=WCL4

}
    # # # # # # # # #  PLOTTING AND AESTHETICS  # # # # # # # # #
porc_= cont/length(WCL2T)*100   
porc<-paste(c("porcentaje de error: ",porc_,"%"), collapse = "")
print(porc)

sec=seq(from=1, to=365, by=4)

sahel.data <- data.frame(TIMES,
                         
                         RAINT = RAINT_,
                         TPAVT = TPAVT,
                         RAINTS = RAINT_/100,
                         TPAVTS = TPAVT/100,
                         WCL2T, 
                         WCL3T,
                         WCL4T,
                         WCLER,
                         EXPV)



sahel.graph2 <- ggplot(sahel.data,aes(TIMES,WCLER)) + 
  
  scale_y_continuous("RELATIVE WATER CONTENT LAYER 2", sec.axis=sec_axis(~.*100, name= "RAIN AND TEMPERATURE"))+
  
  geom_col(aes(x=TIMES , y  =WCLER), color="black", size=.001,alpha=0.2)+
  
  geom_line(aes(x=TIMES , y = TPAVTS), color="red",size=1)+ 
 
  annotate("text",  x =200, y =1.01, label = "TPAVT (Cº)", color="red")+
  
  geom_col(aes(x=TIMES , y = RAINTS),  color="lightblue4" ,size=0.2, alpha=.5)+
  annotate("text",  x =240, y =1.01, label = "RAIN (mm)", color="lightblue4")+
  
  geom_hline(aes(yintercept  =WCST2), color="red", linetype=5,size=.5) + 
  annotate("text", x =-10, y =WCST2+0.02, label = "WCST", color="red")+
  geom_hline(aes(yintercept  =WCFC2), color="lawngreen", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCFC2+0.02, label = "WCFC", color="lawngreen")+  
  geom_hline(aes(yintercept  =WCWP2), color="yellow4", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCWP2+0.02, label = "WCWP", color="yellow4")+
  geom_hline(aes(yintercept  =WCAD2), color="lightgoldenrod4", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCAD2+0.02, label = "WCAD", color="lightgoldenrod4")+

  geom_line(aes(x=TIMES , y  =WCL2T), color="dodgerblue4",size=1)+ 
  annotate("text", x =60, y =1.01, label = "WCL2", color="dodgerblue4")+
  
  geom_point(aes(x=TIMES , y  =EXPV), color="gold", size=3)+
  annotate("text", x =120, y =1.01, label = "EXPERIMENTAL VALUE", color="gold")+
  
  theme(legend.position = "bottom")+
  theme_classic()+
  labs(    x = "DAYS",  
           y = "RELATIVE WATER CONTENT LAYER 2",
           title ="L2SU SOILS (SAHEL MODEL)", 
           subtitle = "SEMI-HARID SOILS THAT EASILY LEAK",
           caption = SRC,
           colour="variables"
  ) + theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

  for (val in sec) {
    sahel.graph2 <- sahel.graph2 + annotate("text",  x =val, y =-0.07, label =SDATES[val], color="grey", angle= 315, size=3)
  }


sahel.graph3 <- ggplot(sahel.data,aes(TIMES,WCLER)) + 
  
  scale_y_continuous("RELATIVE WATER CONTENT LAYER 3", sec.axis=sec_axis(~.*100, name= "RAIN AND TEMPERATURE"))+
  
  geom_line(aes(x=TIMES , y = TPAVTS), color="red",size=1.5)+ 
  annotate("text",  x =200, y =1.01, label = "TPAVT (Cº)", color="red")+
  
  geom_col(aes(x=TIMES , y = RAINTS),  color="lightblue4" ,size=0.5, alpha=.5)+
  annotate("text",  x =240, y =1.01, label = "RAIN (mm)", color="lightblue4")+
  
  geom_hline(aes(yintercept  =WCST3), color="red", linetype=5,size=.5) + 
  annotate("text", x =-10, y =WCST3+0.02, label = "WCST", color="red")+
  geom_hline(aes(yintercept  =WCFC3), color="lawngreen", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCFC3+0.02, label = "WCFC", color="lawngreen")+  
  geom_hline(aes(yintercept  =WCWP3), color="yellow4", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCWP3+0.02, label = "WCWP", color="yellow4")+
  geom_hline(aes(yintercept  =WCAD3), color="lightgoldenrod4", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCAD3+0.02, label = "WCAD", color="lightgoldenrod4")+
  
  geom_line(aes(x=TIMES , y  =WCL3T), color="darkslateblue",size=1.5)+
  annotate("text", x =80, y =1.01, label = "WCL3", color="darkslateblue")+
  
  theme(legend.position = "bottom")+
  theme_classic()+
  labs(    x = "DAYS",  
           y = "RELATIVE WATER CONTENT LAYER 3",
           title ="L2SU SOILS (SAHEL MODEL)", 
           subtitle = "SEMI-HARID SOILS THAT EASILY LEAK",
           caption = SRC,
           colour="variables"
  ) + theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) 
  for (val in sec) {
    sahel.graph3 <- sahel.graph3 + annotate("text",  x =val, y =-0.07, label =SDATES[val], color="grey", angle= 315, size=3)
  }


sahel.graph4 <- ggplot(sahel.data,aes(TIMES,WCLER)) + 
  
  scale_y_continuous("RELATIVE WATER CONTENT LAYER 4", sec.axis=sec_axis(~.*100, name= "RAIN AND TEMPERATURE"))+
  
  geom_line(aes(x=TIMES , y = TPAVTS), color="red",size=1.5)+ 
  annotate("text",  x =200, y =1.01, label = "TPAVT (Cº)", color="red")+
  
  geom_col(aes(x=TIMES , y = RAINTS),  color="lightblue4" ,size=0.5, alpha=.5)+
  annotate("text",  x =240, y =1.01, label = "RAIN (mm)", color="lightblue4")+
  
  geom_hline(aes(yintercept  =WCST4), color="red", linetype=5,size=.5) + 
  annotate("text", x =-10, y =WCST4+0.02, label = "WCST", color="red")+
  geom_hline(aes(yintercept  =WCFC4), color="lawngreen", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCFC4+0.02, label = "WCFC", color="lawngreen")+  
  geom_hline(aes(yintercept  =WCWP4), color="yellow4", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCWP4+0.02, label = "WCWP", color="yellow4")+
  geom_hline(aes(yintercept  =WCAD4), color="lightgoldenrod4", linetype=5,size=.5) +
  annotate("text", x =-10, y =WCAD4+0.02, label = "WCAD", color="lightgoldenrod4")+
  
  geom_line(aes(x=TIMES , y  =WCL4T), color="blue",size=1.5)+
  annotate("text", x =100, y =1.01, label = "WCL4", color="blue")+
  
  theme(legend.position = "bottom")+
  theme_classic()+
  labs(    x = "DAYS",  
           y = "RELATIVE WATER CONTENT LAYER 4",
           title ="L2SU SOILS (SAHEL MODEL)", 
           subtitle = "SEMI-HARID SOILS THAT EASILY LEAK",
           caption = SRC,
           colour="variables"
  ) + theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

  for (val in sec) {
    sahel.graph4 <- sahel.graph4 + annotate("text",  x =val, y =-0.07, label =SDATES[val], color="grey", angle= 315, size=3)
  }

x11(width=16)
print(sahel.graph4,print.eval=TRUE)
x11(width=16)
print(sahel.graph3,print.eval=TRUE)
x11(width=16)
print(sahel.graph2,print.eval=TRUE)

sahel.dataToExcel <- data.frame(
  
  RAIN = RAINT_,
  TPAV = TPAVT,
  EVSW= EVSWT,
  TRW=TRWTA,
  
  EXPV,
  WCL2T,
   
  WCL3T,
  WCL4T,
  WCLER)

name <-paste(c(DIR,"/sahel3L",format(Sys.time(), "%a %b %d %X %Y"),".xls"), collapse = "")

write.xlsx(sahel.dataToExcel, file=name,
           sheetName="data Sahel", append=FALSE)




