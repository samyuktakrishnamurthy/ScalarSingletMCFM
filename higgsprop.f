      function higgsprop(s)
      implicit none
      include 'types.f'
      complex(dp):: higgsprop
c--- computes Higgs propagator for Higgs boson four-momentum squared s
c--- if CPscheme = .true. then it is computed in the complex pole
c---   scheme (Goria, Passarino, Rosco, arXiv:1112.5517,
c---           and c.f. Eq. (2.11) of arXiv:1206.4803)
c--- otherwise it takes the usual Breit-Wigner form
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'cpscheme.f'
      include 'first.f'

      real(dp):: s,mhbarsq,mhbar,gammahbar,vevsq
      complex(dp):: cfac

c---  Defining all the random parts that need to be added to run the code
      complex(dp):: modHiggsSE, delZn

      complex(dp):: p1_1,p1_2,p1_3
      complex(dp):: p2_1,p2_2,p2_3

      complex(dp):: p3_1,p3_2,p2
c--- Define the mass of the scalar ms
      real(dp):: lambdaS, msc
      complex(dp):: muH

      save mhbarsq,cfac,muH
      lambdaS = 0.1_dp
      msc = 120_dp
c---  In this code, all of the newly defined variables are msc - mass of the scalar singlet, lambdaS - the coupling betweent the scalar and the higgs
c---  Define Complex higgs mass squared

      muH = SQRT(cplx2(hmass**2,-hmass*hwidth))
c-- Defining the renormalized SE in parts to avoid errors
      p1_1 = sqrt(s**2 - 4*s*msc**2)
      p1_2 = LOG((-s +2*msc**2+p1_1)/(2*msc**2))
      p1_3 = p1_2*p1_1/s

      p2_1 = sqrt(muH**4 - 4*muH**2*msc**2)
      p2_2 = LOG((-muH**2 +2*msc**2 +p2_1))/(2*msc**2)
      p2_3 = p2_2*P2_1/muH**2

      p3_1 = ((2*p2_2*msc**2)/p2_1) - 1
      p3_2 = (s-muH**2)*p3_1/muH**2

      p2 = p2_3 - p3_2
      modHiggsSE = vevsq*lambdaS**2*(p1_3 - p2)

      delZn = vevsq*lambdaS**2*(p3_1)/(2*muH**2)

      if (CPscheme) then
c--- complex pole scheme propagator      
        if (first) then
          mhbarsq=hmass**2+hwidth**2
          mhbar=sqrt(mhbarsq)
          gammahbar=mhbar/hmass*hwidth
          cfac=cplx2(one,gammahbar/mhbar)
          first=.false.
        write(6,*)
        write(6,*)'****************************************************'
        write(6,*)'*  Using complex pole scheme for Higgs propagator  *'
        write(6,99) mhbar,gammahbar
        write(6,*)'****************************************************'
        write(6,*)
        endif
        higgsprop=cfac/(s*cfac-cplx2(mhbarsq,zip)+modHiggsSE)
      else
c--- Breit Wigner propagator      
        higgsprop=one/cplx2(s-hmass**2,hmass*hwidth)
      endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end
      
      
      
