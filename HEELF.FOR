      subroutine heelf(x,qtoe,qheel,qtoex,bb,bbpri,dc,hw,top,bheel,
     * bet,sur,pah,gc,g1,vheel,vmheel,del)
      hwx=hw+(bheel-top-x)*tan(bet)
      if(qheel.gt.0.)then
      qx=qheel+x*(qtoe-qheel)/bb
      else
      bxa=x-bb+bbpri
      if(bxa.lt.0.)bxa=0.
      qx=bxa*qtoex/bbpri
      endif
      v4=g1*x*hwx
      v5=g1*x*x*tan(bet)/2.
      v6=x*sur
      v7=gc*x*dc
      v8=x*qheel
      if(v8.lt.0.)v8=0.
      v9=x*(qx-qheel)/2.
      if(qheel.lt.0.)v9=bxa*qx/2.
      v10=pah*tan(del)
      bxpri=x
      if(qheel.lt.0.)bxpri=bxa
      vheel=v4+v5+v6+v7-v8-v9+v10
      vmheel=v4*x/2.+v5*x*2./3.+v6*x/2.+v7*x/2.-v8*x/2.-v9*bxpri/3.
     * +v10*x
      return
      end
