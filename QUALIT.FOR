      SUBROUTINE QUALIT(EO)
      CHARACTER*3 QUALI,quali1*2,kkod,duc*1,squalo*4,xx*2,a*5,x1*1,x2*1
      CHARACTER*2 UNI,ekos*9,kkod1*3,duc1*1
      character*3 cub
      common/brbs/bbr,bbs
      COMMON/A/KMX,KFEX,    BR,BS,EB1,EBMAX,EEMAX,EBAX,EEAX,EOB,EOE,AN,A
     1NO,KMXX,LT,AAG,BBG,AKA,AKB,AKZ,DBR,NBET,LCX,ESH,TSH,ZED,BBZ,IPROP
      COMMON /A1/PERMIN,PERMAX,TO1,TO2,TO3,BP,NFER,KFER,MSHEAR
      common/A2/to11
      COMMON/XXH/   IDEP,MIA,MIAX,IC,ICX,KODU,MEBU,NO,IDJO,NBJJ,MAM,ETA,
     1 MSEC,IQ1,ETEO,KIND,PI,NTRAU,MID,KIC,msecy
      COMMON /PAGE/LINES,LINEX,NPAGE
      COMMON /LAT/ LAT
      COMMON/PAR/RA,IRA,RO,GO,GOF,ETOD,YNMAXO,RAC,RL,ICICL,dra,rax,ran
     * ,ici
      COMMON/PERCEN/pminb,pmaxb,pminc,pmaxc
      COMMON/UNI/UNI
      common/quali/quali,kkod,duc,squalo
      common/push1/prox,acc
      COMMON/JFE/JFER,jpts,jspe,jcap,fck,fyk,ispir,nstr
      common/nbp/nbp,nbp1,trd,kshe,kerth(10),vvx(4,2,6),fcd,fywd,crit1,
     * crit2,fbd,zak,redis,gamc,gams,romin,fyd,zndx,gamce,gamse
      common/irt/irt
      ill=0
      READ(irt,900)QUALI,xx,acc,BS,UNI,kkod,duc,PMINc,PMAXc,EOE,eb1
     * ,ebmax,fywd,jpts,zed,gamc,gams,esx,nbet,dra,rax,ran,ici,ispir
      br=0.
      if(kkod.eq.'   ')kkod='GRK'
      if(acc.eq.0.)then
       if(kkod.eq.'EC2')then
       acc=1.
       elseif(kkod.eq.'GRK')then
       acc=.85
       elseif(kkod.eq.'din')then
       acc=.85
       endif
      endif
c     print*,'qualit beg+',dra,rax,ran
      if(dra.eq.0.)dra=.01
      if (rax.eq.0..and.ran.eq.0.)rax=1.31
      eb1=-eb1
      ebmax=-ebmax
      a=quali//xx
      squalo=a(3:)
      x1=xx(:1)
      x2=xx(2:)
      if(x1.eq.'0')then
      k1=0
      elseif(x1.eq.'1')then
      k1=1
      elseif(x1.eq.'2')then
      k1=2
      elseif(x1.eq.'3')then
      k1=3
      elseif(x1.eq.'4')then
      k1=4
      elseif(x1.eq.'5')then
      k1=5
      elseif(x1.eq.'6')then
      k1=6
      elseif(x1.eq.'7')then
      k1=7
      elseif(x1.eq.'8')then
      k1=8
      elseif(x1.eq.'9')then
      k1=9
      else
c     write(0,*)'*Illegal concrete quality'
      ill=1
      endif
      if(x2.eq.'0')then
      k2=0
      elseif(x2.eq.'1')then
      k2=1
      elseif(x2.eq.'2')then
      k2=2
      elseif(x2.eq.'3')then
      k2=3
      elseif(x2.eq.'4')then
      k2=4
      elseif(x2.eq.'5')then
      k2=5
      elseif(x2.eq.'6')then
      k2=6
      elseif(x2.eq.'7')then
      k2=7
      elseif(x2.eq.'8')then
      k2=8
      elseif(x2.eq.'9')then
      k2=9
      else
c     write(0,*)'*Illegal concrete quality'
      ill=1
      endif
      if(ill.eq.1)goto 4000
      nbp=10*k1+k2
      bp=nbp
c      print 34,squalo,nbp,bp
c  34 format('qualit:-------',a4,i5,f5.0)
      if(jpts.eq.0)jpts=1
      if(kkod.eq.'   ')kkod='GRK'
      if(duc.eq.' ')duc='M'
      if(quali.eq.'  C'.or.quali.eq.'  B')then
      if(bs.eq.0)bs=500.
      fck=bp
      fyk=bs
      romin=.1725*fck**(2./3.)/fyk
      bbz=.30*fck**(2./3.)
      if(kkod.eq.'EC2')then
       if(duc.eq.'H')then
       cc=.35
       zndx=.55
       elseif(duc.eq.'M')then
       cc=.65
       zndx=.65
       else
       cc=.75
       zndx=.75
       endif
      elseif(kkod.eq.'GRK')then
      cc=.65
      zndx=.65
      elseif(kkod.eq.'CYP')then
      cc=.65
      zndx=.50
      else
      cc=.65
      zndx=1.5
      endif
c     print*,'qualit',zak,fck,fyk,cc
      ibs=int(bs)
      else
      if(bs.eq.0)bs=2.4
       if(abs(bs-2.4).lt..0001)then
       fyk=420.
       elseif(abs(bs-2.8).lt..0001)then
       fyk=500.
       else
       fyk=220.
       endif
      if(br.ne.0.)bp=br
      an=1.75
      ano=2.10
      if(eb1.eq.0.)EB1=-2.
      if(ebmax.eq.0.)EBMAX=-3.5
      EEMAX=5.
      if(esx.ne.0.)eemax=abs(esx)
      EBAX=-2.
      EEAX=3.
      EOB=EO*.000001
      IF(EOE.EQ.0.)EOE=200.
      endif
      if(fywd.eq.0.)fywd=bs
      ibsw=int(fywd)
      JFER=1
      if(uni.eq.'  ')uni='KN'
      if(uni.eq.'kN')uni='KN'
      if(uni.eq.'‘œ')uni='TO'
      if(uni.eq.' Õ')uni='KN'
      if(uni.eq.'ÃÕ')uni='MN'
      if(gamc.eq.0.)gamc=1.5
      if(quali.eq.' S2'.or.quali.eq.' S3'.or.quali.eq.' ST'.or.quali
     * .eq.' St')then
      if(gams.eq.0.)gams=1.10
      else
      if(gams.eq.0.)gams=1.15
      endif
      gamce=gamc
      gamse=gams
      if(pminc.eq.0)pminc=.01
      if(pmaxc.eq.0)pmaxc=.04
      permin=pminc
      permax=pmaxc
      if(quali.eq.'  C')zak=100.*cc*fck/(fyk*gamc/gams)
      if(quali.eq.' B3')then
      if(br.eq.0.)bp=100.
      br=24.0
      to11=.50
      to1=.8
      to2=2.0
      to3=2.0
      quali1='00'
      goto 1000
      elseif(quali.eq.' B2')then
      if(br.eq.0.)bp=80.
      br=19.5
      to11=.40
      to1=.7
      to2=1.8
      to3=1.8
      quali1='25'
      goto 1000
      elseif(quali.eq.' B1')then
      if(br.eq.0.)bp=60.
      br=14.4
      to11=.30
      to1=.6
      to2=1.6
      to3=1.6
      quali1='60'
      goto 1000
      endif
      if(quali.eq.' S2'.or.quali.eq.' S3'.or.quali.eq.' ST'.or.quali
     * .eq.' St')then
       if(bs.eq.0.)then
        if(nbp.eq.35)fyk=235.
        if(nbp.eq.75)fyk=275.
        if(nbp.eq.55)fyk=355.
       endif
        if(lines+4.ge.linex)then
        call ffeed(npage)
        lines=0
        endif
        if(lat.eq.0)then
        write(*,200)quali,nbp,uni
        else
        write(*,2000)quali,nbp,uni
        endif
        lines=lines+4
      goto 4000
      endif
      if(eb1.eq.0.)EB1=-2.
      if(ebmax.eq.0.)EBMAX=-3.5
      EEMAX=5.
      if(esx.ne.0)then
      eemax=abs(esx)
      else
      if(quali.eq.'  C')eemax=20.
      if(kkod.eq.'din')eemax=25.
      if(kkod.eq.'CYP')eemax=10.
      endif
      EBAX=-2.
      EEAX=3.
      EOB=EO*.000001
      IF(EOE.EQ.0.)EOE=200.
      AN=1.75
      ANO=2.10
      if(quali.eq.'  C')then
      an=1.
      ano=1.
      endif
      if(quali.eq.'  C')then
      bs=bs/gams
      fywd=fywd/gams
      br=acc*bp/gamc
      trd=.035*bp**(2./3.)
      fcd=bp/gamc
      fyd=fyk/gams
c     print*,'*qualit fcd=',fcd
      fbd=.3315*bp**(2./3.)
      goto 10
      endif
      if(nbp.eq.25)then
      IF(BR.EQ.0.)BR=17.5
      to11=.50
      TO1=.75
      TO2=1.8
      TO3=3.0
      fbd=1.8*1.75
      elseif(nbp.eq.15)then
      IF(BR.EQ.0.)BR=10.5
      to11=.35
      TO1=.50
      TO2=1.2
      TO3=2.0
      fbd=1.4*1.75
      elseif(nbp.eq.35)then
      if(br.eq.0.)br=23.0
      to11=.60
      to1=1.0
      to2=2.4
      to3=4.0
      fbd=2.2*1.75
      elseif(nbp.eq.45)then
      if(br.eq.0.)br=27.0
      to11=.70
      to1=1.1
      to2=2.7
      to3=4.5
      fbd=2.6*1.75
      elseif(nbp.eq.55)then
      if(br.eq.0.)br=30.0
      to11=.80
      to1=1.25
      to2=3.0
      to3=5.0
      fbd=3.0*1.75
      else
c     write(*,*)'*Illegal concrete quality'
      ier=ier+1
      goto 4000
      endif
   10 PERMI=PERMIN*100.
      PERMA=PERMAX*100.
      IRA=0
      if(lines+17.ge.linex)then
      call ffeed(npage)
      lines=0
      endif
      eemaxu=eemax
      if(esx.lt.0.)then
      if(quali.eq.'  C')eemaxu=20.
      if(kkod.eq.'din')eemaxu=25.
      if(kkod.eq.'CYP')eemaxu=10.
      endif
      ekos='         '
      if(kkod.eq.'GRK'.and.quali.eq.'  C')ekos='≈ Ÿ” 2000'
      kkod1=kkod
      duc1=duc
      if(kkod.eq.'din'.or.kkod.eq.'DIN')then
      kkod1='   '
      duc1=' '
      endif
      if(kkod.eq.'din')ekos='din 102  '
      if(kkod.eq.'DIN')ekos='DIN 4227 '
      cub='   '
      if(quali.eq.'  C')call cube(nbp,cub)
      IF(LAT.EQ.0)GO TO 358
      if(nbet.eq.0)then
      WRITE(*,1977)QUALI,nbp,cub,ibs,kkod1,duc1,ekos,ibsw,BR,BS,EB1,
     * ebmax,ebax,eemaxu,     EOE,AN, ANO,gamc,gams,UNI
      else
      WRITE(*,1988)QUALI,nbp,cub,ibs,kkod1,duc1,ekos,ibsw,BR,BS,EB1,
     * ebmax,ebax,eemaxu,     EOE,AN, ANO,gamc,gams,UNI
      endif
      GO TO 359
  358 if(nbet.eq.0)then
      WRITE(*, 977)QUALI,nbp,cub,ibs,kkod1,duc1,ekos,ibsw,BR,BS,EB1,
     * ebmax,ebax,eemaxu,     EOE,AN,ANO,gamc,gams,UNI
      else
      WRITE(*, 988)QUALI,nbp,cub,ibs,kkod1,duc1,ekos,ibsw,BR,BS,EB1,
     * ebmax,ebax,eemaxu,     EOE,AN,ANO,gamc,gams,UNI
      endif
  359 LINES=LINES+17
  500 BR=BR*100.
      BS=BS*100.
      EOB=EOB*1000000.
      EOE=EOE*1000000.
c     print*,'qualit eoe',eoe,eob
      goto 4000
 1000 fbd=to2*1.75
      if(lat.eq.0)then
      write(*,2977)quali,quali1,bs,fywd,uni
      else
      write(*,3977)quali,quali1,bs,fywd,uni
      endif
      lines=lines+5
      fywd=fywd*175.
c4000 print*,'qualit',fck,fcd,br,fyd,bs,eob,eoe,squalo
c     print*,eb1,ebmax,eemax,ebax,eeax,permin,permax
 4000 bbr=br
      bbs=bs
c     print*,'qualit end',dra,rax,ran,ici
      return
  900 FORMAT(A3, a2, 2F5.0,3X,A2,2x,a3,4x,a1,6F5.0,i5,
     * 4f5.0,i5,3f5.0,2i5)
 1977 FORMAT(/10X,'CONCRETE QUALITY',A3,I2,a3,'   STEEL QUALITY  S',i3,
     *12x,a3,2x,a,3x,a9  /31x,'                     S',i3,' STIRRUPS'
     *             /10X,                        'DESIGN STRENGTH OF CONC
     1RETE fcd=',F6.2,2X,'MN/M2 '/10X,'DESIGN STRENGTH OF STEEL   fyd=',
     *F6.1,2X,'MN/M2 '
     2 //10X,'STRESS-STRAIN CURVE FOR CONCRETE PARABOLIC UP TO Âc1=' ,
     3 F6.1,1X,'0/00'/10X,'ULTIMATE CONCRETE STRAIN IN FLEXURE Âcu   =',
     4F6.1,1X,'0/00'/10X,'ALLOWED CONCRETE STRAIN IN AXIAL COMPRESSION Â
     5cu   =',F6.1,' 0/00'/10X,'ULTIMATE ELONGATION OF STEEL FOR SECTION
     6 IN FLEXURE Âsu   =',F6.1,1X,'0/00'
     8//10X,'ELASTICITY MODULUS OF STEEL   Es=',F7.0,2X,'GN/M2'
     */10X,'SAFETY COEFFICIENT N=',F5.2,'/',F5.2
     */10x,'SAFETY COEFFICIENTS „M: „c/„s = ',f5.2,'/',f5.2
     9//10x,'UNITS: ',A2,' ,M')
 1988 FORMAT(/10X,'CONCRETE QUALITY',A3,I2,a3,'   STEEL QUALITY  S',i3,
     *12x,a3,2x,a,3x,a9  /31x,'                     S',i3,' STIRRUPS'
     *             /10X,                        'DESIGN STRENGTH OF CONC
     1RETE fcd=',F6.1,2X,'MN/M2 '/10X,'DESIGN STRENGTH OF STEEL   fyd=',
     *F6.1,2X,'MN/M2 '
     2 //10X,'STRESS-STRAIN CURVE FOR CONCRETE LINEAR UP TO ec1=' ,
     3 F6.2,1X,'0/00'/10X,'ULTIMATE CONCRETE STRAIN IN FLEXURE Âcu   =',
     4F6.1,1X,'0/00'/10X,'ALLOWED CONCRETE STRAIN IN AXIAL COMPRESSION Â
     5cu   =',F6.1,' 0/00'/10X,'ULTIMATE ELONGATION OF STEEL FOR SECTION
     6 IN FLEXURE Âsu   =',F6.1,1X,'0/00'
     8//10X,'ELASTICITY MODULUS OF STEEL   Es=',F7.0,2X,'GN/M2'
     */10X,'SAFETY COEFFICIENT N=',F5.2,'/',F5.2
     */10x,'SAFETY COEFFICIENTS „M: „c/„s = ',f5.2,'/',f5.2
     9//10x,'UNITS: ',A2,' ,M')
C
  977 FORMAT(/10X,'–OIOTH” ”KYPOƒEMATO”',A3,i2,a3,'    ¡‘«√œ—…¡ ◊¡À’¬¡',
     *'  S',i3,12x,a3,2x,a,3x,a9 /45x,'              S',i3,' ”’Õƒ≈‘«—≈”'
     *                                        /10X,'Y–OÀO√I”TIKH ANTOXH
     1” YPOƒEMATO” fcd=',F6.2,2X,'MN/M2'/10X,'’–œÀœ√ ¡Õ‘œ◊« XAÀYBO” O–ÀI ”
     2”MOY fyd=',F6.1,2X,'MN/M2'//10X,     'ƒIA√PAMMA TA”EŸN ¬—¡◊YN”EŸN
     3” ’—œƒ≈MATO” –A—ABOÀIKO MEXPI Âc1=',F5.1,1X,'0/00'/10X,  'ME√I”TH
     4¬—¡◊’Õ”«',1X,'”KYPOƒEMATO” ”E KAMÿH Âcu   =',F6.1,1X,'0/00'
     5                                               /10X,'ME√I”TH ¬—¡◊’
     6Õ”’ ” ’—œƒ. ”≈  ≈Õ‘—… « »À…ÿH Âcu   =',F6.1,' 0/00'/10X,'Ã≈√…”‘« M
     7« ’Õ”« œ–À…”Ãœ’ ƒ…¡‘œÃ«” ”≈  ¡Ãÿ«   Âsu=',F4.1,1X,'0/00'
     * //10X,
     8'METPO EÀA”TIKOTHTO” XAÀYBO”  Es=',F7.0,2X,'GN/M2'/10X,'”YNTEÀE”TH
     9” A”÷AÀEIA” Ì=',F5.2,'/',F5.2/10X,'”’Õ‘≈À≈”‘≈” ¡”÷¡À≈…¡” ’À… ŸÕ „Ã
     *: „c/„s = ',f5.2,'/',f5.2
     *//10X,'MONAƒE”: ',A2,' ,M')
  988 FORMAT(/10X,'–OIOTH” ”KYPOƒEMATO”',A3,i2,a3,'    ¡‘«√œ—…¡ ◊¡À’¬¡',
     *'  S',i3,12x,a3,2x,a,3x,a9 /45x,'              S',i3,' ”’Õƒ≈‘«—≈”'
     *                                        /10X,'Y–OÀO√I”TIKH ANTOXH
     1” YPOƒEMATO” fcd=',F6.1,2X,'MN/M2'/10X,'’–œÀœ√ ¡Õ‘œ◊« XAÀYBO” O–ÀI
     2”MOY fyd=',F6.1,2X,'MN/M2'//10X,     'ƒIA√PAMMA TA”EŸN ¬—¡◊YN”EŸN
     3” ’—œƒ≈MATO” √—¡ÃÃ… œ MEXPI Âc1=',F5.2,1X,'0/00'/10X,'ME√I”TH
     4¬—¡◊’Õ”«',1X,'”KYPOƒEMATO” ”E KAMÿH Âcu   =',F6.1,1X,'0/00'
     5                                               /10X,'ME√I”TH ¬—¡◊’
     6Õ”’ ” ’—œƒ. ”≈  ≈Õ‘—… « »À…ÿH Âcu   =',F6.1,' 0/00'/10X,'Ã≈√…”‘« M
     7« ’Õ”« œ–À…”Ãœ’ ƒ…¡‘œÃ«” ”≈  ¡Ãÿ«   Âsu   =',F4.1,1X,'0/00'
     * //10X,
     8'METPO EÀA”TIKOTHTO” XAÀYBO”  Es=',F7.0,2X,'GN/M2'/10X,'”YNTEÀE”TH
     9” A”÷AÀEIA” Ì=',F5.2,'/',F5.2/10X,'”’Õ‘≈À≈”‘≈” ¡”÷¡À≈…¡” ’À… ŸÕ „Ã
     *: „c/„s = ',f5.2,'/',f5.2
     *//10X,'MONAƒE”: ',A2,' ,M')
  200 format( /10x,'–œ…œ‘«” ◊¡À’¬œ”',a3,i2//10x,'ÃœÕ¡ƒ≈”: ',a2,' ,M')
 2000 format( /10x,'STEEL QUALITY',a3,i2//10x,'UNITS: ',a2,' ,M')
 2977 FORMAT(/10X,'–OIOTH” ”KYPOƒEMATO”',A3,a2/10X,'≈–…‘—≈–œÃ≈Õ« ‘¡”« ◊¡
     *À’¬¡     Ûe=',F6.2,2X,'T/CM2'/10X,'≈–…‘—≈–œÃ≈Õ« ‘¡”« ”’Õƒ≈‘«—ŸÕ Ûe
     *=',F6.2,2x,'T/CM2'//10x,'UNITS: ',A2,' ,M')
 3977 FORMAT(/10X,'CONCRETE QUALITY ',A3,a2/10X,'ALLOWABLE STRESS OF STE
     *EL    Ûe=',F6.2,2X,'T/CM2 '/10X,'ALLOWABLE STRESS OF STIRRUPS Se='
     *,F6.2,2x,'T/CM2'//10x,'UNITS: ',A2,' ,M')
      END
