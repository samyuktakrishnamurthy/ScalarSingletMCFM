      function modifiedHiggsSE(s)
      implicit none
      include 'types.f'
      complex(dp):: modifiedHiggsSE

      complex(dp):: MHSE_1_1
      complex(dp):: MHSE_1_2
      complex(dp):: MHSE_1

      complex(dp):: modifiedHiggsSE
      complex(dp):: modifiedHiggsSE_1
      complex(dp):: modifiedHiggsSE_2
      complex(dp):: modifiedHiggsSE_3
c--- computes Higgs renormalized SE for Higgs boson four-momentum squared s
c--- otherwise it takes the usual Breit-Wigner form
c--- 

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'cpscheme.f'
      include 'first.f'
      real(dp):: s,mhbarsq,mhbar,gammahbar
      complex(dp):: cfac
      save mhbarsq,cfac
      
c--- Define the Coupling between Higgs and Scalar - constant for now
      real(dp):: lambdaS

c--- Define the mass of the scalar ms 
      real(dp):: mS
c-- Find where VEV is defined in the code / what it is called
      real(dp)::Vev
c-- Define Complex higgs mass squared
      real(dp)::muH_sq
      muH_sq = cplx2(hmass**2,hmass*hwidth)
c-- Defining the renormalized SE in parts to avoid errors
      MHSE_1_1 = SQRT(s(1,2)**2 - 4*s(1,2)*mS**2)
      MHSE_1_2 = -s(1,2) + 2*mS**2
      MHSE_1 = (LOG((MHSE_1_2 + MHSE_1_1)/2*mS**2)*MHSE_1_1)/s(1,2)

      MHSE_2_1 = SQRT(muH_sq**2 - 4*muH_sq*mS**2)
      MHSE_2_2 = -muH_sq + 2*mS**2
      MHSE_2 = (LOG(MHSE_2 + MHSE_2_1))

      
