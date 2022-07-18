      subroutine toef(x,qtoe,qheel,qtoex,bb,bbpri,dc,gc,vtoe,vmtoe)
      if(qheel.gt.0.)then
      qx=qheel+(bb-x)*(qtoe-qheel)/bb
      v2=(qtoe-qx)*x/2.
      else
      qx=(bbpri-x)*qtoex/bbpri
      v2=(qtoex-qx)*x/2.
      endif
      v1=qx*x
      v3=gc*x*dc
      vtoe=v1+v2-v3
      vmtoe=v1*x/2.+2.*v2*x/3.-v3*x/2.
      return
      end
