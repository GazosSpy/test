      Program retain
      character*3 gru
      character*5 free
      character*4 name(20)
      real ka,kp,ng,ngam,nc,ic,iq,igam,kh,kv
      CHARACTER*3 BETON
      CHARACTER*2 UNI
      DIMENSION ENI(50),EMXI(50),EMYI(50),NART(50),DY(50),FE1U(50),FE2U(
     150),CAPA(50),E1U(50),EEU(50),SCU(50),NAX1(50),ENP(50),EMXP(50)
     2 ,EMYP(50)
      COMMON/A/KMX,KFEX,    BR,BS,EB1,EBMAX,EEMAX,EBAX,EEAX,EOB,EOE,AN,A
     1NO,KMXX,LT,AAG,BBG,AKA,AKB,AKZ,DBR,NBET,LCX,ESH,TSH,ZED,BBZ,IPROP
      COMMON /A1/PERMIN,PERMAX,TO1,TO2,TO3,BP,NFER,KFER,MSHEAR
      COMMON/XXH/   IDEP,MIA,MIAX,IC,ICX,KODU,MEBU,NO,IDJO,NBJJ,MAM,ETA,
     1 MSEC,IQ1,ETEO,KIND,PI,NTRAU,MID,KIC
      COMMON /PAGE/LINES,LINEX,NPAGE
      COMMON /LAT/ LAT
      COMMON/PAR/RA,IRA,RO,GO,GOF,ETOD,YNMAXO,RAC,RL,ICICL
      COMMON/JFE/JFER
      COMMON/PERCEN/pminb,pmaxb,pminc,pmaxc
      COMMON/UNI/UNI
C READ/WRITE INPUT DATA
      read(*,10)gru
      if(gru.ne.'ret')then
      write(*,*)'*Incorrect program name  ',gru
      stop
      endif
      read(*,11)jse,kh ,kv,nn,linex,npage,free,ggc,cov,iopt,ipasiv,jopt,
     * bincr,eo,int,red1,reds1,mon
      eps=kh
      if(nn.eq.0)nn=5
      dx1=1./nn
      if(ggc.eq.0.)ggc=25
      gc=ggc
      if(cov.eq.0.)cov=.05
      if(bincr.eq.0.)bincr=.05
      if(eo.eq.0.)eo=.30e+8
      if(red1.eq.0.)red1=1.20
      if(reds1.eq.0.)reds1=1.5
      lat=1
      if(free.eq.'greek')lat=0
      if(linex.eq.0)linex=43
      if(npage.eq.0)npage=1
      read(*,20)name
      write(*,30)name
    1 read(*,10)gru
      if(gru.eq.'   ')goto 1
    2 if(gru.eq.'geo')then
      read(*,40)hw,top,bot,slop,slopi,dcmin,dcmax,b1min,b1max,ipri1,
     * b2min,b2max,ipri2
      bot1=top+(slop+slopi)*hw
      bot=max(bot,bot1)
      write(*,50)hw,slop,slopi
      if(b1max.eq.0.)b1max=10.*b1min
      if(b2max.eq.0.)b2max=10.*b2min
      if(dcmax.eq.0.)dcmax=10.*dcmin
      write(*,51)top,bot,dcmin,dcmax,b1min,b1max,ipri1
     * b2min,b2max,ipri2
      elseif(gru.eq.'soi')then
      read (*,60)qall,gg1,g2,ph1,ph2,c1,c2,ht,bet,saf1,safs,over,slide,
     *slides,del
      if(gg1.eq.0.)gg1=18.
      if(g2.eq.0.)g2=gg1
      g1=gg1
      if(ph1.eq.0.)then
      write(*,*)'*Internal friction angle of backfill must be inputed!'
      goto 444
      endif
      if(ph2.eq.0.)ph2=ph1
      if(safs.eq.0.)safs=1.5
      if(slides.eq.0.)slides=1.10
      qalls=qall
      if(del.eq.0.)del=bet
      if(g2.eq.0.)g2=g1
      if(over.eq.0.)over=1.5
      if(slide.eq.0.)slide=1.5
      write(*,70)g1,g2,ph1,ph2,c2
      write(*,71)ht,bet
      if(ipasiv.eq.0)then
      write(*,72)over,slide
      else
      write(*,73)over,slide
      endif
        if(qall.gt.0.)then
        write(*,80)qall
        else
        write(*,85)
        endif
      elseif(gru.eq.'qua')then
      call qualit(eo)
      elseif(gru.eq.'sur')then
      read(*,90)sur,px,py,pm,pxs
      write(*,100)sur,px,pxs,py,pm
      endif
    3 read(*,10)gru
      if(gru.eq.'   ')goto 3
      if(gru.ne.'las')goto 2
C COMPUTE EARTH PRESSURE COEFFICIENTS
      pi=3.14159265
      pic=pi/180.
      ph1=ph1*pic
      ph2=ph2*pic
      bet=bet*pic
      del=del*pic
      cbet=cos(bet)
      cdel=cos(del)
      if(mon.eq.0)then
      ka=cos(ph1)**2/cos(del)/(1.+sqrt(sin(ph1+del)*sin(ph1-bet)/(cos(d
     *el)*cos(bet))))**2
      kp=cos(ph2)**2/(1.-sin(ph2))**2
      if(ipasiv.eq.2)ka=1-sin(ph1)
      write(*,115)ka,kp
      elseif(mon.eq.2)then
      ka=tan(pi/4.-ph1/2.)
      ka=ka*ka
      if(ipasiv.eq.2)ka=1.-sin(ph1)
      kp=tan(pi/4.+ph2/2.)
      kp=kp*kp
      write(*,110)ka,kp
      else
      a=cbet
      b=a*a
      c=cos(ph1)
      c=c*c
      cc=cos(ph2)
      cc=cc*cc
      ka=a*(a-sqrt(b-c))/(a+sqrt(b-c))
      if(ipasiv.eq.2)ka=1.-sin(ph1)
      kp=(1.+sqrt(1.-cc))/(1.-sqrt(1.-cc))
      write(*,120)ka,kp
      endif
C STEM THICKNESS REVISION
  400 if(jse.eq.1)then
        if(eps.gt.0.)then
        red=red1
        reds=reds1
        jopt=1
          if(mon.eq.0)then
          akh=1.0*kh
          if(ipasiv.eq.2)akh=1.0*kh
          thet=atan(akh/(1.-kv))
          betlim=ph1-thet
          if(bet.gt.betlim)then
          write(*,*)'*Solution not possible for earthquake specified!'
          stop
          endif
          ka=cos(ph1-thet)**2/(cos(thet)*cos(thet+del)*
     *       (1.+sqrt(sin(ph1+del)*sin(ph1-thet-bet)/(cos(thet+del)*
     *       cos(bet))))**2)
          gc=ggc*(1.-kv)
          g1=gg1*(1.-kv)
          se1=eps/(1.-kv)
          se2=1.
          else
          se1=eps
          se2=1.+eps*2.
          endif
        else
        write(*,*)'*Earthquake acceleration is zero - eartquake analysis
     * not possible!'
        stop
        endif
        write(*,785)kh,kv,ka
      else
        write(*,786)
        se1=0.
        se2=1.
        red=1.
        reds=1.
      endif
      tto=to1*red*1000.
      qq=(sur+g1*hw/2.)*hw*ka*cdel*se2
      if(ipasiv.eq.1)qq=qq-g2*ht*ht*kp/2.
      botq=.85*(bot-cov)
      tt=qq/botq
      bot1=bot
      if(tt.gt.tto)then
      bot=qq/(tto*.85)+cov
      if(iopt.eq.1)top=top+(bot-bot1)
      endif
      tt=qq/((bot-cov)*.85)
      bor=hw*slop
      bol=bot-top-bor
      write(*,140)top,bot
C STEM SECTION FORCES 
      write(*,150)
      x1=-dx1
  606 x1=x1+dx1
      x=x1*hw
      bxr=x1*bor
      bxl=x1*bol
      h5=sur*x*ka*cdel*se2
      h6=g1*x*x*ka*cdel/2.*se2
      w1=gc*top*x
      w2=.5*gc*bxr*x
      w3=.5*gc*bxl*x
      w4=.5*g1*bxl*x
      w7=sur*bxl
      h1=w1*se1
      h2=w2*se1
      h3=w3*se1
      h4=w4*se1
      h7=w7*se1
      y=0.
      if(ipasiv.eq.1)then
      y=x-hw+ht
      if(y.lt.0.)y=0
      endif
      h8=g2*y*y*kp/2
      qxf=-h8
      qxus=h1+h2+h3+h4+h7+pxs
      qxu=px+h5+h6+qxus
      qyf=py+w1+w2+w3+w4+w7
      qyu=0.
      qmf=pm-w1*top/2.+w2*bxr/3.-w3*(bxl/3.+top)-w4*(2.*bxl/3.+top)
     * -w7*(bxl/2.+top)-h8*y/3.
      qmus=h1*x/2.+(h2+h3)*x/3.+h4*x*2./3.+h7*x+pxs*x
      qmu=h5*x/2.+h6*x/3.+px*x+qmus
      qxx=qxf+qxu
      qyy=qyf+qyu
      qmm=qmf+qmu
      bx=top+bxl+bxr
      ex=bx/2.-bxr
      qmmc=qmm+qyy*ex
      lcx=1
      eni(1)=-qyy*.1/red
      emxi(1)=qmmc*.1/red
      emyi(1)=0.
      dy(1)=bx/2.
      call dteult(bx,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,cov,1,1.,0.,pminb,
     * pmaxb,eni,emxi,emyi,0,dy,fe1u,fe2u,capa,e1u,eeu,scu,2,nax1,
     2 enp,emxp,emyp)
      write(*,145)x,qmmc,qyy,qxx,bx,fe1u(1),fe2u(1),e1u(1),eeu(1)
      if(abs(x-hw).gt..01)goto 606
      b1=b1min
      b2=b2min
      dc=dcmin
      goto 250
  500 if(b2.lt.b2max.and.ipri.le.ipri1)then
      b2=b2+bincr
      if(b2.gt.b2max)b2=b2max
      elseif(b1.lt.b1max)then
      b1=b1+bincr
      if(b1.gt.b1max)b1=b1max
      else
      write(*,*)' *Solution not possible'
      goto 2000
      endif        
C OVERTURNING-SLIDING-SOIL PRESSURE
  250 toe=b1-bor
      heel=b2-(bot-bor)
      w9=toe*ht*g2
      w10=heel*hw*g1
      w11=sur*heel
      h9=w9*se1
      h10=w10*se1
      h11=w11*se1
      htot=hw+dc+(b2-top)*tan(bet)
      h12=sur*htot*ka*cos(del)*se2
      h13=g1*htot*htot*ka/2.*cos(del)*se2
      pah=h12+h13
      h14=0.
      h15=0.
      if(ipasiv.eq.1)then
      h14=g2*ht*dc*kp
      h15=g2*dc*dc*kp/2
      endif
      w16=gc*(b1+b2)*dc
      w17=g1*(b2-top)**2*tan(bet)/2
      h16=w16*se1
      h17=w17*se1
      w18=pah*tan(del)
      sxf=qxf-h14-h15
      sxu=px+h12+h13+qxus+h9+h10+h11+h16+h17
      syf=qyf+w9+w10+w11+w16+w17+w18
      syu=qyu
      smf=qmf+qxf*dc-h14*dc/2.-h15*dc/3.+w9*(bor+toe/2.)-(w10+w11)*
     * (bot-bor+heel/2.)-w16*((b1+b2)*.5-b1)-w17*(b2-(b2-top)/3.)-w18*b2
      smu=px*(hw+dc)+h12*htot/2.+h13*htot/3.+qmus+qxus*dc+h9*(ht/2.+dc)
     * +h10*(hw/2.+dc)+h11*(dc+hw+htot)/2.+h16*dc/2.+h17*(htot-
     2 (b2-top)*tan(bet)*2./3.)
      sxx=sxf+sxu
      syy=syf+syu
      smm=smf+smu
      bb=b1+b2
      cx=bb/2.-b1
      smmc=smm+syy*cx
      s1=syy/bb+smmc*6./(bb*bb)
      s2=syy/bb-smmc*6./(bb*bb)
      ex=smmc/syy
      ss=s1
      bbpri=bb
      if(s2.lt.0.)then
      xx=bb/2.-ex
      ss=2.*syy/(3.*xx)
      bbpri=3.*xx
      endif
      exmax=bb/6
      if(jopt.eq.1)goto 410
      if(ex.gt.exmax)goto 500
  410 if(ex.gt.bb/3.)goto 500
      if(int.eq.1)write(*,160)ex,exmax
      smft=smf-syf*b1
      smut=smu-syu*b1
      sfov=-smft/smut
      fs=tan(ph2)*syy+.67*c2*bbpri
      sfsl=(fs-sxf)/sxu
      if(int.eq.1)write(*,170)sfov,sfsl
      if(jse.eq.0.)then
        if(sfsl.lt.slide)goto 500
      else
        if(sfsl.lt.slides)goto 500
      endif
      if(sfov.lt.over.and.jse.ne.1)goto 500
C COMPUTE BEARING PRESSURE
      if(qalls.ne.0.)then
      qall=qalls*reds
      goto 275
      endif
      if(jse.eq.1)then
      saf=safs
      elseif(saf1.ne.0)then
      saf=saf1
      elseif(c2.eq.0.)then
      saf=2.
      else
      saf=3.
      endif
      d=dc+ht
      if(ph2.gt.0.)then
      nq=exp(pi*tan(ph2))*kp
      nc=(nq-1.)/tan(ph2)
      ngam=1.5*(nq-1.)*tan(ph2)
      else
      nc=5.14
      nq=1.
      ngam=0.
      endif
      bprim=bb-2.*ex
      q=g2*d
      dcc=1.+.4*d/bprim
      if(ph2.eq.0.)dcc=.4*d/bprim
      dq=1.+2.*tan(ph2)*(1.-sin(ph2))**2*atan(d/bprim)
C     iq=(1.-.5*sxu/(syy+bprim/tan(ph2)*c2))**5
      iq=(1.-.7*sxu*saf/(syy*saf+bprim*c2/tan(ph2)))**3
      ic=iq-(1.-iq)/(nq-1.)
C     if(ph2.eq.0.)ic=.5-.5*sqrt(1.-sxu/(bprim*c2))
      if(ph2.eq.0.)then
      awx=1.-sxu*saf/(bprim*c2)
        if(awx.ge.0.)then
        ic=.5+.5*sqrt(awx)
        else
        write(*,*)'*Negative sqrt argument in ic-factor calculation!'
        stop
        endif
      endif
C     igam=(1.-.7*sxu/(syy+bprim*c2/tan(ph2)))**5
      igam=(1.-sxu*saf/(syy*saf+bprim*c2/tan(ph2)))**3
      if(ph2.gt.0.)then
      qult=c2*nc*dcc*ic+d*g2*nq*dq*iq+.5*g2*bprim*ngam*igam
      else
      qult=c2*nc*(1.+dcc-ic)+d*g2*nq*dq*iq
      endif
      qall=qult/saf
      if(int.eq.1)write(*,180)nc,nq,ngam,ic,iq,igan,dcc,dq,qult,qall
  275 continue
      if(int.eq.1)write(*,190)s1,ss,s2
      if(int.eq.1)write(*,195)b1,b2,bb
      if(ss.gt.qall)goto 500
C TOE FORCES
      bt=b1-bor
      call toef(bt,s1,s2,ss,bb,bbpri,dc,gc,vtoe,vmtoe)
      dcq=.85*(dc-cov)
      tt=vtoe/dcq
      dct=dc
      if(tt.gt.tto.and.dc.lt.dcmax)then
      dct=vtoe/(tto*.85)+cov
      idct=dct*100.+1.
      dct=idct*.01
      dc=dct
      goto 250
      endif
C HEEL FORCES
      bh=b2+bor-bot
      call heelf(bh,s1,s2,ss,bb,bbpri,dc,hw,top,b2,bet,sur,pah,gc,g2,
     * vheel,vmheel,del)
      dcq=.85*(dc-cov)
      tt=vheel/dcq
      dch=dc
      if(tt.gt.tto.and.dc.lt.dcmax)then
      dch=vheel/(tto*.85)+cov
      idch=dch*100.+1
      dch=idch*.01
      endif
      dcx=max(dct,dch)
      if(dcx.gt.dc)then
      dc=dcx
      if(dc.gt.dcmax)dc=dcmax
      goto 250
      endif
 2000 write(*,160)ex,exmax
      write(*,170)sfov,sfsl
C     write(*,180)nc,nq,ngam,ic,iq,igam,dcc,dq,qult,qall
      write(*,181)qult,qall
      write(*,190)s1,ss,s2
      write(*,195)b1,b2,bb
      write(*,227)dc
      write(*,1001)
      x1=-dx1
  707 x1=x1+dx1
      x=x1*bh
      call heelf(x,s1,s2,ss,bb,bbpri,dc,hw,top,b2,bet,sur,pah,gc,g1,
     * vheel,vmheel,del)
      lcx=1
      eni(1)=0.
      emxi(1)=vmheel*.1/red
      emyi(1)=0.
      call dteult(dc,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,cov,1,1.,0.,pminb,
     1 pmaxb,eni,emxi,emyi,0,dy,fe1u,fe2u,capa,e1u,eeu,scu,2,nax1,
     2 enp,emxp,emyp)
      write(*,237)x,vmheel,vheel,fe1u(1),fe2u(1),e1u(1),eeu(1)
      if(abs(x-bh).gt..01)goto 707
      write(8,1002)
      x1=-dx1
  808 x1=x1+dx1
      x=x1*bt
      call toef(x,s1,s2,ss,bb,bbpri,dc,gc,vtoe,vmtoe)
      lcx=1
      eni(1)=0.
      emxi(1)=vmtoe*.1/red
      emyi(1)=0.
      call dteult(dc,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,cov,1,1.,0.,pminb,              
     1 pmaxb,eni,emxi,emyi,0,dy,fe1u,fe2u,capa,e1u,eeu,scu,2,nax1,
     2 enp,emxp,emyp)
      write(*,237)x,vmtoe,vtoe,fe2u(1),fe1u(1),e1u(1),eeu(1)
      if(abs(x-bt).gt..01)goto 808
      if(jse.eq.2)then
      jse=1
      goto 400
      endif
  444 stop
   10 FORMAT(A3)
   11 FORMAT( 5X,I5,2F5.2,3I5,A5,2F5.2,3I5,F5.2,E10.0,I5,2F5.2,i5)
   20 FORMAT(20A4)
   40 FORMAT(9F5.2,I5,2F5.2,I5)
   30 FORMAT( //T5,' R E T A I N I N G   W A L L   D E S I G N'//20A4)
   50 FORMAT(//T5,'WALL HEIGHT ABOVE FOOTING =',F8.2,' M'
     1        /T5,'WALL BATTER:   FRONT FACE =',F8.3,' %'
     2        /T5,'                BACK FACE =',F8.3,' %')
   51 FORMAT( /T5,'TRIAL WALL THICKNESS:     TOP =',F8.2,' M'
     1        /T5,'                       BOTTOM =',F8.2,' M'
     2        /T5,'BASE SLAB THICKNESS :   Dc =',F8.2,'  /',F8.2,' M'
     3        /T5,'BASE SLAB DIMENSIONS:   B1 =',F8.2,'  /',F8.2,' M'
     *          ,I10
     4        /T5,'                        B2 =',F8.2,'  /',F8.2,' M'
     *          ,I10)
   60 FORMAT( 15F5.2)
   70 FORMAT( /T5,'SOIL DATA:  WEIGHT OF BACKFILL =',F10.3,' KN/M3'
     1        /T5,'            WEIGHT OF BASE SOIL=',F10.3,' KN/M3'
     2        /T5,'            INTERNAL FRICTION OF BACKFILL =',
     3          F10.3,' DEG'
     4        /T5,'            INTERNAL FRICTION OF BASE SOIL=',
     5          F10.3,' DEG'
     6        /T5,'            COHESION OF BASE SOIL =        ',
     7          F10.3,' KN/M2')
   71 FORMAT( /T5,'BACKFILL HEIGHT OVER TOE=',F8.2,' M'
     1        /T5,'SLOPE OF WALL BACKFILL  =',F8.2,' DEG')
   72 FORMAT( /T5,'SAFETY FACTORS:   OVERTURNING =',F7.2,
     1        /T5,'                  SLIDING     =',F7.2,'  PASSIVE RESI
     2STANCE NEGLECTED')
   73 FORMAT( /T5,'SAFETY FACTORS:   OVERTURNING =',F7.2,
     1        /T5,'                  SLIDING     =',F7.2,'  PASSIVE RESI
     2STANCE CONSIDERED')
   80 FORMAT( /T5,'ALLOWABLE SOIL PRESSURE=',F10.2,' KN/M2')
   85 FORMAT( /T5,'ALLOWABLE SOIL PRESSURE TO BE COMPUTED USING GIVEN DA
     1TA')
   90 FORMAT( 5X,F5.2,6F10.2)
  100 FORMAT( /T5,'SURCHARGE  = ',F7.2,' KN/M2'
     1        /T5,'FORCES AT STEM TOP:  HORIZONTAL =',F10.2,' KN/M',5X,
     * 'HORIZ.EARTHQ =',F10.2,' KN/M'
     2        /T5,'                     VERTICAL   =',F10.2,' KN/M'
     3        /T5,'                     MOMENT     =',F10.2,' KNM/M')
  120 FORMAT( /T5,'RANKINE EARTH PRESSURE COEFFICIENTS - BACKFILL: Ka ='
     1         ,F10.5
     2        /T5,'                                     BASE SOIL: Kp ='
     3         ,F10.5)
  115 FORMAT( /T5,'MONONOBE-OKABE EARTH PRESSURE COEFFICIENTS - BACKFILL
     1: Ka ='  ,F10.5
     2        /T5,'                                            BASE SOIL
     3: Kp ='  ,F10.5)
  110 FORMAT( /T5,'COULOMB EARTH PRESSURE COEFFICIENTS - BACKFILL: Ka ='
     1         ,F10.5
     2        /T5,'                                     BASE SOIL: Kp ='
     3         ,F10.5)
  140 FORMAT( /T5,'CALCULATED STEM THICKNESS:     TOP =',F8.2,' M'
     1        /T5,'                            BOTTOM =',F8.2,' M')
  150 FORMAT( /T5,'S T E M   S E C T I O N   F O R C E S   A N D   R E I
     * N F O R C E M E N T'
     1        /T5,'    X        M         N         V       d        As1
     2     As2     e1      ee')
  145 FORMAT( T5,F7.2,3F10.2,F8.3,2F8.2,2F8.2)
  160 FORMAT( /T5,'ECCENTRICITY =',F8.3,'  -  SECTION CORE AT B/6=',
     1          F8.3)  
  170 FORMAT( /T5,'SAFETY FACTORS:  AGAINST OVERTURNING=',F8.3
     1        /T5,'                 AGAINST SLIDING    =',F8.3)
  190 FORMAT( /T5,'SOIL PRESSURES:     q_toe =',F10.2,' KN/M2   /  ',
     1         'q_toemax =', F10.2,' KN/M2'
     2        /T5,'                    q_heel=',F10.2,' KN/M2')
  195 FORMAT( /T5,'BASE SLAB DIMENSIONS:  B_toe  =',F8.2,' M'
     1        /T5,'                       B_heel =',F8.2,' M'
     2        /T5,'                       B_total=',F8.2,' M')
  180 FORMAT( /T5,'B E A R I N G   C A P A C I T Y   F A C T O R S'
     1        /T5,'Nc =',F8.3,'     Nq =',F8.3,'     Ng =',F8.3
     2        /T5,'ic =',F8.3,'     iq =',F8.3,'     ig =',F8.3
     3        /T5,'dc =',F8.3,'     dq =',F8.3
     4        /T5,'ULTIMATE SOIL PRESSURE  : q_ult =',F10.2,' KN/M2'
     5        /T5,'ALLOWABLE SOIL PRESSURE : q_all =',F10.2,' KN/M2')
  181 FORMAT( /T5,'ULTIMATE SOIL PRESSURE  : q_ult =',F10.2,' KN/M2'
     *        /T5,'ALLOWABLE SOIL PRESSURE : q_all =',F10.2,' KN/M2')
  227 FORMAT( /T5,'BASE SLAB THICKNESS :  Dc     =',F8.2,' M')
  237 FORMAT(  T5,F7.2,2F10.2,2F8.2,2F8.2)
 1001 FORMAT( /T5,'H E E L   F O R C E S   A N D   R E I N F O R C E M E
     1 N T'   /T5,'    X        M         V       As1      As2     e1
     2   ee')
 1002 FORMAT(//T5,'T O E   F O R C E S   A N D   R E I N F O R C E M E N
     1 T'     /T5,'    X        M         V        As1     As2     e1
     2   ee')
  785 FORMAT( /////T5,'D E S I G N   W I T H   E A R T H Q U A K E   I N
     * C L U D E D  -  kh=',F6.2,'  kv=',F6.2,'  ka=',f8.5)
  786 FORMAT( /////T5,'D E S I G N   W I T H O U T   E A R T H Q U A K E
     * ')
      END         




      
                 
