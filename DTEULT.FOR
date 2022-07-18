      SUBROUTINE DTEULT(D,BO,B1,DO1,B2,DO2,DBO,B3,D3,B4,D4,H1,KODE,RATO,
     1PERGIV,PERMIN,PERMAX,ENI,EMXI,EMYI,NART,DY,FE1U,FE2U,CAPA,E1U,EEU,
     2SCU,MODU,NAX1,ENP,EMXP,EMYP)
      COMMON/A/KMX,KFEX,    BR,BS,EB1,EBMAX,EEMAX,EBAX,EEAX,EOB,EOE,AN,A
     1NO,KMXX,LC,AAG,BBG,AKA,AKB,AKZ,DBR,NBET,LCX,ESH,TSH,ZED,BBZ,IPROP
      COMMON/TEST/ITEST,IER
      COMMON/PAR/RA,IRA,RO,GO,GOF,ETOD,YNMAXO,RAC,RL,ICICL
      COMMON/PAGE/LINES,LINEX,NPAGE
      DIMENSION XI(22),YI(22),XP(22),YP(22),XN(22),YN(22),XNR(22),YNR(22
     1),XW(22),YW(22),XU(22),YU(22),XPF(22),YPF(22),XIF(22),YIF(22),XNF(
     222),YNF(22),XAF(22),YAF(22),FEGIV(22),FEMIN(22),FEMAX(22),FE(22),F
     3ER(22),IPRI(22)
      common/conx/conx1,conx2
      COMMON/Q/QMAX,QMIN,TORSX,QMAX2,QMIN2,eng,wpf,wmf,f
      DIMENSION ENI(1),EMXI(1),EMYI(1),NART(1),DY(1),FE1U(1),FE2U(1),CAP
     1A(1),E1U(1),EEU(1),SCU(1),NAX1(1),ENP(1),EMXP(1),EMYP(1)
      IF(ITEST.EQ.2)WRITE(*,555)
      CALL CORDSE(D,BO,B1,DO1,B2,DO2,DBO,B3,D3,B4,D4,H1,XI,YI,XIF,YIF)
      IPRI(1)=1
      IPRI(2)=1
      CALL CGFJ(XI,YI,XS,YS,F,PHI,AJ1,AJ2,1)
c     print*,' *DTEULT ys... ',ys,xs,f,aj1,aj2,phi
      wpf=aj1/(f*ys)
      wmf=aj1/(f*(d-ys))
      DO 513 I=1,LCX
      IF(DY(I).EQ.0.)DY(I)=YS
      EMXI(I)=EMXI(I)+ENI(I)*DY(I)
  513 EMYI(I)=0.
      CALL TRANSP(XI,YI,XP,YP,XIF,YIF,XPF,YPF,ENI,EMXI,EMYI,ENP,EMXP,EMY
     1P,XS,YS,PHI,AKX,AKY,F,AJ1,AJ2)
      DIV=1.+RATO
      FEMI=F*PERMIN/(DIV*100.)
      FEMA=F*PERMAX/(DIV*100.)
      FEGI=F*PERGIV/(DIV*100.)
      FEMIN(2)=FEMI
      FEMAX(2)=FEMA
      FEGIV(2)=FEGI
      FEMIN(1)=FEMI*RATO
      FEMAX(1)=FEMA*RATO
      FEGIV(1)=FEGI*RATO
      INEQ=0
      IF(MODU.EQ.2.or.modu.eq.4)INEQ=1
      IF(MODU.EQ.2.or.modu.eq.4)MODE=1
      conx1=0.
      conx2=0.
      DO 1000 LC=1,LCX
      ycon1=0.
      ycon2=0.
c     print*,' *DTEULT1 eni... ',lc,eni(lc),emxi(lc),emyi(lc),ineq
c     print*,' *dteult2 enp... ',lc,enp(lc),emxp(lc),emyp(lc),mode
      IF(INEQ.EQ.0)GO TO 40
      CALL INITAt(XP,YP,XN,YN,XPF,YPF,XNF,YNF,ENP,EMXP,EMYP,XO,YO,GAM,DE
     1L,BET,E1,E2,YNMAX,YPAR,YE,D,H,H1,F,AJ1,AJ2,NAX,EN,0,
     2 XNR,YNR,XW,YW,XU,YU,AKX,AKY)
c     print*,' *Dteult3 ynmax...',ynmax,ye,gam,del,bet,nax
      NAX1(LC)=NAX
      EE=YE*E1/YNMAX
      IF(EN.GT.0..AND.nax.eq.3)GO TO 10
      IF(NAX.EQ.1)GOTO 400
      IF(E1.LE.0..AND.Ee.LE.0.)GO TO 40
c     IF(ETOD.LT.0..AND.ETOD.GT.-.40)GOTO 40
      EEX=3.
      CALL DOUBFE(EEX,XN,YN,XNR,YNR,XW,YW,XU,YU,YNMAX,XNF,YNF,FE,FEGIV,G
     1AM,XO,YO,YE,EE,NAX,XPF,YPF,MODE,FE1,FE2,ID1,ID2,E1,E2,SC,F,ENP,EMX
     2P,D,H)
      FEC=FE1
      IF(ID1.EQ.2)FEC=FE2
      FET=FE2
      IF(ID2.EQ.1)FET=FE1
      CAP=1.
      IF(FE1.GE.0..AND.FE2.GE.0.)then
        if(fet.lt.fec)then
        fe(1)=fec
        fe(2)=fec
        endif
c     if(-e1.gt.1.5)then
c       if(emxi(lc).lt.0.)then
c       ycon1=ynmax*(1.+1.5/e1)
c       else
c       ycon2=ynmax*(1.+1.5/e1)
c       endif
c     endif
c     print*,' dteult1',e1,ynmax,ycon1,ycon2,nax
      goto 1004
      endif
c     if(-e1.gt.1.5)then
c       if(emxi(lc).lt.0.)then
c       ycon1=ynmax*(1.+1.5/e1)
c       else
c       ycon2=ynmax*(1.+1.5/e1)
c       endif
c     endif
c     print*,' dteult2',e1,ynmax,ycon1,ycon2,nax
      IF(FET.LT.0.)GOTO 25
      CALL SMPLFE(XN,YN,XNR,YNR,XW,YW,XU,YU,YNMAX,XNF,YNF,FE,F
     1EGIV,GAM,XO,YO,YE,EE,NAX,XPF,YPF,MODE,FE1,FE2,ID1,ID2,E1,E2,SC,F,
     2  ENP,EMXP,D,H,H1)
c     if(-e1.gt.1.5)then
c       if(emxi(lc).lt.0.)then
c       ycon1=ynmax*(1.+1.5/e1)
c       else
c       ycon2=ynmax*(1.+1.5/e1)
c       endif
c     endif
c     print*,' dteult3',e1,ynmax,ycon1,ycon2,nax
      GOTO 1005
   25 CALL FDESGN(XP,YP,XN,YN,XNR,YNR,XW,YW,XU,YU,XPF,YPF,XNF,YNF,XAF,YA
     1F,FEMIN,FEMAX,FE,IPRI,ENP,EMXP,EMYP,NART,AJ1,AJ2,F,AKX,AKY,PHI,GAM
     2,DEL,E1,E2,EE,SC,YNMAX,RER,IEND,NAX,BET,MODE,FEGIV,ID1,0,XO,YO,D)
      NAX1(LC)=NAX
c     print*,' dteult4',e1,ynmax,ycon1,ycon2,nax
      IF(NAX.EQ.2)GO TO 1000
      CAP=1.
      GO TO 1005
   10 CALL TENSML(ENP,EMXP,YPF,FE)
      CAP=1.
      SC=AN
      GO TO 1004
   40 IF(NART(LC).EQ.1)GO TO 90
      DO 5 I=1,KFEX
    5 FE(I)=FEMIN(I)
      RO=PERMIN/100.
      JU=2
      GO TO 1001
  120 IF(CAP.GE.1)GO TO 1005
      DO 15 I=1,KFEX
   15 FE(I)=FEMAX(I)
      RO1=PERMIN/(CAP*100.)
      IF(RO1.LE.PERMAX*.5)THEN
      RO=RO1
      GOTO 25
      ENDIF
      RO=PERMAX
      JU=3
      GO TO 1001
   30 IF(CAP.LE.1)GO TO 1005
      RO=(RO1+PERMAX/CAP)*.5
      GO TO 25
   90 DO 95 I=1,KFEX
   95 FE(I)=FEGIV(I)
      RO=PERGIV
      JU=1
 1001 CALL CAPFAC(XP,YP,XN,YN,XNR,YNR,XW,YW,XU,YU,XPF,YPF,XNF,YNF,XAF,YA
     1F,FE,ENP,EMXP,EMYP,AJ1,AJ2,F,AKX,AKY,CAP,E1,E2,EE,SC,YNMAX,GAM,DEL
     2,NAX,RER,IEND,MODE,0)
      NAX1(LC)=NAX
      IF(NAX.EQ.1)GO TO 400
      IF(NAX.EQ.2)GO TO 1000
      GO TO (1005,120,30),JU
 1004 RER=0.
      IEND=1
 1005 FE1U(LC)=FE(1)*10000.
      FE2U(LC)=FE(2)*10000.
      CAPA(LC)=CAP
      E1U(LC)=E1
      EEU(LC)=EE
      SCU(LC)=SC
c     print*,' As ',fe2u(lc),fe1u(lc),lc,cap,e1,ee,sc
      if(-e1/cap.gt.1.5)then
        if(emxi(lc).lt.0.)then
        ycon1=ynmax*(1.+1.5*cap/e1)
        else
        ycon2=ynmax*(1.+1.5*cap/e1)
        endif
      endif
      GO TO 999
  400 IF(ENP(LC).GT.0.)THEN
      E1=EEAX
      EE=EEAX
      GOTO 10
      ENDIF
      FEA=(-ANO*ENP(LC)-BR*F)/BS
      FEMI=PERMIN*F/100.
      IF(FEA.LT.FEMI)FEA=FEMI
      CAP=(-BR*F-FEA*BS)/(ANO*ENP(LC))
      FEA=FEA*10000.
      FE1U(LC)=FEA*(YS-H1)/(D-2.*H1)
      FE2U(LC)=FEA-FE1U(LC)
      CAPA(LC)=CAP
      E1U(LC)=EBAX
      EEU(LC)=EBAX
      SCU(LC)=ANO
      e1=ebax
      if(-e1/cap.gt.1.5)then
      ycon1=d/2.
      ycon2=d/2.
      endif
  999 if(conx1.lt.ycon1)conx1=ycon1
      if(conx2.lt.ycon2)conx2=ycon2
c     print*,'*dteult',e1,ynmax,ycon1,ycon2,conx1,conx2,lc
 1000 CONTINUE
      RETURN
  555 FORMAT(5X,'*dteult')
      END
