MODULE B85
!-----------------------------------------------------------------------
!
!    File:         B85.90
!    Module:       B85
!    Type:         Module w/ Functions
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      3/31/16
!    WeakLib ver:  
!
!    Purpose:
!      Provides the function need considering the physics in Bruenn 85
!
!    CONTAINS:
!
!      Function totalECapEm: 
!                   gives opacity with np as approximation 
!
!      Function totalElasticScatteringKernel: 
!                   gives isoenergetic scattering kernel with nMoment
!                   l = 0 or 1
!                            
!
!    Modules used:
!      wlKindModule
!      wlExtPhysicalConstantsModule
!      ( function fexp is called )
!-----------------------------------------------------------------------
  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: &
    h, kMeV, therm1, therm2, dmnp, me, mbG, mp, mn, cvel_inv, cvel, ergmev,&
      cv_p, cv_n, ca_p, ca_n, gf, hbarc, cv, ca
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero 

  implicit none
   
  PUBLIC totalECapEm, &
         totalElasticScatteringKernel, &
         NESKern

CONTAINS

!========================Function=============================

  REAL(dp) FUNCTION totalECapEm &
    ( energy, rho, T, Z, A, chem_e, chem_n, chem_p, xheavy, xn, xp )
!------------------------------------------------------------------------------
! Purpose:
!   To compute the neutrino absorptivity.
!   (1) Absorptivity = emissivity + inverse of mean path
!------------------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), INTENT(in) :: energy, rho, T, Z, A, chem_e, chem_n, chem_p, &
                            xheavy, xn, xp

    REAL(dp) :: TMeV, n, qpri, nhn, npz, etapn, jnucleon, jnuclear, midFe, &
                            midE, chem_v, mpG, mnG, fexp, &
                            rop, ron, midFexpp, midFep, midEp, midCons
    REAL(dp) :: emitnp, absornp, emitni, absorni

    TMeV   = T * kMeV                         ! kmev = 8.61733d-11 [MeV K^{-1}]
      N    = A - Z
    qpri   = chem_n - chem_p + 3.0_dp + dmnp  ! [MeV] 3 = energy of the 1f5/2 level 
    chem_v = chem_e + chem_p - chem_n - dmnp  ! neutrino chemical potential
    
   
     mpG   = mp * ergmev * cvel_inv * cvel_inv ! proton mass [g]
     mnG   = mn * ergmev * cvel_inv * cvel_inv ! neutron mass [g]

    if(n.le.34.0)               nhn = 6.0_dp
    if(n.gt.34.0.and.n.le.40.0) nhn = 40.0_dp - N
    if(n.gt.40.0)               nhn = 0.0_dp

    if(z.le.20.0)               npz = 0.0
    if(z.gt.20.0.and.z.le.28.0) npz = z - 20.0
    if(z.gt.28.0)               npz = 8.0
    
    
!    etapn = rho * ( xn - xp ) / ( mbG * ( FEXP( (chem_n-chem_p)/TMeV ) - 1.0_dp ) )
     rop   = rho * xp / mpG                   ! Approxiation in the nondegenerate regime
     ron   = rho * xn / mnG 

!----------------------------------------------------------------------------
! Stop the function if any of etapn/rop/ron is negative
!--------------------------------------------------------------------------- 
    IF ( ( rop < 0.0_dp ) .or. ( ron < 0.0_dp ) ) THEN
!    IF ( (etapn < 0.0_dp ) .or. ( rop < 0.0_dp ) .or. ( ron < 0.0_dp ) ) THEN
      WRITE(*,*)'xn - xp is ', xn - xp
      WRITE(*,*)'fexp term is ', FEXP( (chem_n-chem_p)/TMeV ) - 1.0_dp 
      WRITE(*,*)'chem_n - chem_p is ', chem_n - chem_p 
      STOP 
    END IF

!-----------------------------------------------------------------------------
!   j_nuclear(emitni) and chi_nuclear(absorni) 
!-----------------------------------------------------------------------------
    IF ( xheavy * npz * nhn == 0.0_dp ) THEN
       emitni  = 0.0_dp
      absorni  = 0.0_dp
    ELSE
       midEp   = (energy+qpri)**2 * SQRT( & 
                    MAX( 1.0_dp - ( me / (energy+qpri) )**2, 0.0_dp ) )     
      midFexpp = FEXP( (energy+qpri-chem_e) / TMeV )
       midFep  = 1.0_dp / ( midFexpp + 1.0_dp )
      midCons  = therm2 * rho * xheavy * npz * nhn * midEp / (mbG * A)

       emitni  = midCons * midFep 
      absorni  = midCons * FEXP( (chem_n + dmnp - chem_p - qpri) /TMeV ) &
               * ( 1.0_dp - midFep) 
    END IF

!-----------------------------------------------------------------------------
!   j_nucleon(emitnp) and chi_nucleon(absornp)
!-----------------------------------------------------------------------------
!------------------------------------------------------------------
!  Set emitnp + absornp = zero and return if both xn and xp are zero
!------------------------------------------------------------------
    IF ( xn == 0.0 .and. xp == 0.0 ) THEN
      totalECapEm = emitni + absorni
      WRITE(*,*) 'xn and xp = 0'
      RETURN
    END IF

      midFe   = 1.0_dp / ( FEXP( (energy+dmnp-chem_e) / TMeV ) + 1.0_dp )
       midE   = (energy+dmnp)**2 &
                 * SQRT( 1.0_dp - ( me / (energy+dmnp) )**2 )

     jnucleon = therm1 * rop * midE * midFe
     absornp  = therm1 * ron * midE * ( 1.0_dp - midFe )
       emitnp = jnucleon

    totalECapEm = ( emitni + absorni ) + ( emitnp + absornp )

    IF ( ISNAN(totalECapEm) ) THEN
      WRITE(*,*) "totalECapEm is NAN! "
      STOP
    END IF

    RETURN
  END FUNCTION totalECapEm


  REAL(dp) FUNCTION totalElasticScatteringKernel&
     ( energy, rho, T, xh, A, Z, xn, xp, l )

!------------------------------------------------------------------------------
! Purpose:
!   To compute the zero and first legendre coefs for the neutrino-electron 
!   elastic scattering kernel. 
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------
!   Input Variables
!-----------------------------------------------------------------------
    REAL(dp), INTENT(in) :: energy, rho, T, xh, A, Z, xn, xp
    INTEGER, INTENT(in)  :: l
!-----------------------------------------------------------------------
!   Physical Constants 
!-----------------------------------------------------------------------
    REAL(dp)             :: N, nucleiTP, & ! 'TP' for thermal parameter
                            nucleonTP, nucleiExp, Cv0, Cv1,&
                            etann, etapp, Npara

!-----------------------------------------------------------------------
!   Local Variables
!-----------------------------------------------------------------------   

    REAL(dp) :: ESNucleiKernel_0, ESNucleonKernel_0, &
                ESNucleiKernel_1, ESNucleonKernel_1, &
                tempC0, TempC1
    INTEGER  :: nquad = 20

        N     = A - Z
       Cv0    = half * ( cv_p + cv_n) 
       Cv1    = cv_p - cv_n
     Npara    = cvel_inv**4.0 * energy**2.0 / h**3.0
    
    nucleiExp = 4.0_dp * 4.8_dp * 10**(-6.0_dp) * &
                A**(2.0_dp/3.0_dp) * energy**2.0     
    nucleiExp = MAX( nucleiExp, SQRT( TINY( 1.0_dp ) ) )

    nucleiTP  = ( (twpi*gf)**2 / h ) * ( rho*xh/mbG ) * &
                A * ( Cv0 - ( (N-Z)*Cv1 )/(2.0_dp*A) )**2

    nucleonTP = ( twpi * gf )**2 / h

!------------------------------
!  scattering on nuclei
!------------------------------

    tempC0 = IntegralESNuclei( nucleiExp, 0, nquad )
    tempC1 = IntegralESNuclei( nucleiExp, 1, nquad )

    ESNucleiKernel_0 = nucleiTP * tempC0 / 2.0_dp

    ESNucleiKernel_1 = nucleiTP * tempC1 * 3.0_dp / 2.0_dp

    IF ( ISNAN(ESNucleiKernel_0) .or. ISNAN(ESNucleiKernel_1)  ) THEN
     WRITE(*,*) "ERROR AT B85.f90 MARK 1003 !"
     WRITE(*,*) "nucleiExp is ", nucleiExp
     WRITE(*,*) "ESNucleiKernel_0 ", ESNucleiKernel_0
     WRITE(*,*) "ESNucleiKernel_1 ", ESNucleiKernel_1
     STOP
    END IF

!--------------------------------------
!   Scattering on Nucleons
!-------------------------------------

    CALL etaxx( rho, T, xn, xp, etann, etapp )

    ESNucleonKernel_0 = (0.5_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 + 3.0_dp * ca_n**2) + &
                          etapp * ( cv_p**2 + 3.0_dp * ca_p**2) )
 
    ESNucleonKernel_1 = (1.5_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 - ca_n**2) + &
                          etapp * ( cv_p**2 - ca_p**2) )

    IF ( ISNAN(ESNucleonKernel_0) .or. ISNAN(ESNucleonKernel_1)  ) THEN
     WRITE(*,*) "ERROR AT B85.f90 MARK 1004 !"
     STOP
    END IF 

    IF ( l == 0 ) THEN
    
     totalElasticScatteringKernel = Npara * ( ESNucleiKernel_0 &
                                            + ESNucleonKernel_0 )

    ELSE IF ( l == 1) THEN

     totalElasticScatteringKernel = Npara * ( ESNucleiKernel_1 &
                                            + ESNucleonKernel_1 )

    ELSE

     WRITE(*,*) "ERROR: Unable to provide Legendre Moment with &
                        l other than 0 and 1 "
    END IF
    
    IF ( ISNAN(totalElasticScatteringKernel) ) THEN
      WRITE(*,*) "totalElasticScatteringKernel is NAN! "
      WRITE(*,*) "l is", l
      WRITE(*,*) "ESNucleiKernel_0 + ESNucleonKernel_0 ", ESNucleiKernel_0+ESNucleonKernel_0 
      WRITE(*,*) "ESNucleiKernel_1 + ESNucleonKernel_1 ", ESNucleiKernel_1+ESNucleonKernel_1 
      STOP
    END IF

    RETURN

  END FUNCTION totalElasticScatteringKernel

  
  SUBROUTINE NESKern( energygrid, omega, T, chem_e, nesktab )
!----------------------------------------------------------------------
! Purpose:
!    To compute the neutrino-electron scattering (OUT) kernel 
!    (1) e /= ep
!    R_out = cons * ( 1 / e / ep) * 
!                            ( beta1 * I1 + beta2 * I2 + beta3 * I3 )
!    (2) e == ep
!    new form of I1,2,3
!----------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in) :: energygrid, omega
    REAL(dp), INTENT(in) :: T, chem_e      
    REAL(dp), DIMENSION(:,:,:), INTENT(out) :: nesktab
    REAL(dp)             :: tsq, tinv, eta, FEXP, &
                            x1, x2, x2inv, &
                            FA_, FA0, FA1, &
                            y, tpiet, &
                            COMBO1, COMBO2, &
                            G0, G1, G2, buffer1, buffer2
    REAL(dp)             :: beta1, beta2, beta3, log10_sig0, log10_cons
    INTEGER              :: D_e, D_ome, i_e1, i_e2, i_ome
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ediff, etap, fgamm, esq
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: delta, y0, I1, I2, I3,&
                                               A, B, C

       tsq = T*T
      tinv = 1.0/T
       eta = chem_e * tinv
     tpiet = twpi * T                      ! 2*pi*T
     beta1 = ( cv + ca )**2 
     beta2 = ( cv - ca )**2
     beta3 = ca**2 - cv**2
       D_e = SIZE( energygrid )
     D_ome = SIZE( omega )
     log10_sig0  = LOG10(1.764) - 44  ! log10(sig0)            
     log10_cons  = LOG10(half * pi * cvel / ( twpi**3 * me**2 )) + log10_sig0

    ALLOCATE(  ediff( D_e, D_e )        )  ! e1 - e2
    ALLOCATE(   etap( D_e, D_e )        )  ! eta prim
    ALLOCATE(   esq ( D_e, D_e )        )  ! e*e
    ALLOCATE(  fgamm( D_e, D_e )        )  ! gamma function
    ALLOCATE( delta ( D_ome, D_e, D_e ) )
    ALLOCATE(    y0 ( D_ome, D_e, D_e ) )
    ALLOCATE(     A ( D_ome, D_e, D_e ) )
    ALLOCATE(     B ( D_ome, D_e, D_e ) )
    ALLOCATE(     C ( D_ome, D_e, D_e ) )
    ALLOCATE(     I1( D_ome, D_e, D_e ) )
    ALLOCATE(     I2( D_ome, D_e, D_e ) )
    ALLOCATE(     I3( D_ome, D_e, D_e ) )

    DO i_e1 = 1, D_e
     DO i_e2 = 1, D_e

       ediff(i_e1,i_e2) = energygrid( i_e1 ) - energygrid( i_e2 )
        etap(i_e1,i_e2) = eta - ediff(i_e1,i_e2)*tinv
         esq(i_e1,i_e2) = energygrid( i_e1 ) * energygrid( i_e2 )
         
     IF ( i_e1 .ne. i_e2 ) THEN
       fgamm(i_e1,i_e2) = 1.0/(FEXP(-ediff(i_e1,i_e2)*tinv)-1.0)
     ELSE
     END IF     

     END DO ! i_e2
    END DO ! i_e1

    DO i_ome = 1, D_ome
     DO  i_e1 = 1, D_e
      DO  i_e2 = 1, D_e

       delta(i_ome,i_e1,i_e2) =                                          &
                             SQRT( esq(i_e1,i_e1) + esq(i_e2,i_e2)       &
                             - 2.0*esq(i_e1,i_e2)*omega(i_ome)  )
       y0(i_ome,i_e1,i_e2) =                                             &
                             tinv*(-0.5*ediff(i_e1,i_e2)                 &
                              +0.5*delta( i_ome, i_e1, i_e2 )            &
                                *sqrt(1.0+2.0*me*me                      &
                                /(esq(i_e1,i_e2)*(1.0-omega(i_ome)))     & 
                                     )                                   &
                                  )
       A(i_ome,i_e1,i_e2) = esq(i_e1,i_e1) + esq(i_e2,i_e2)              &
                          + esq(i_e1,i_e2)*(3.0 + omega(i_ome))

       B(i_ome,i_e1,i_e2) = energygrid(i_e1) *                           &
                            ( 2.0*esq(i_e1,i_e1)                         &
                            + esq(i_e1,i_e2)*(3.0 - omega(i_ome))        &
                            - esq(i_e2,i_e2)*(1.0 + 3.0*omega(i_ome) ) ) 
       C(i_ome,i_e1,i_e2) = esq(i_e1,i_e1) *                             &
                            ( (energygrid(i_e1) -                        &
                                  energygrid(i_e2)*omega(i_ome) )**2     &
                             - half*esq(i_e2,i_e2)*(1.0-omega(i_ome)**2) &
                             - half*( 1.0+omega(i_ome) )*me*me           &
                                         *(delta(i_ome,i_e1,i_e2)**2)    &
                               /( (1.0-omega(i_ome)) * esq(i_e1,i_e1) ) )  

      END DO ! i_e2
     END DO ! i_e1
    END DO ! i_ome  

!-------------------------------
!   energy_in == energy_out
!-------------------------------
   DO i_ome = 1, D_ome
    DO  i_e1 = 1, D_e
        
      CALL NESKernFsame( eta, y0(i_ome,i_e1,i_e1), FA_, FA0, FA1 )

      I1(i_ome,i_e1,i_e1)= &
         ( tpiet*esq(i_e1,i_e1)*esq(i_e1,i_e1)                      &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                &
              /delta(i_ome,i_e1,i_e1)**5  )                         &
              *(A(i_ome,i_e1,i_e1)*tsq                              &
                            *(2.0*FA1                               &
                             +2.0*y0(i_ome,i_e1,i_e1)*FA0           &
                             +y0(i_ome,i_e1,i_e1)**2 *FA_           &
                             )                                      &
               +B(i_ome,i_e1,i_e1)*T                                &
                            *(FA0                                   &
                             +y0(i_ome,i_e1,i_e1)*FA_               &
                             )                                      &
               +C(i_ome,i_e1,i_e1)*FA_                              &
               )      

      I2(i_ome,i_e1,i_e1) = I1(i_ome,i_e1,i_e1)                    

      I3(i_ome,i_e1,i_e1) =                                         &
         tpiet*esq(i_e1,i_e1)                                       &
              *(1.0-omega(i_ome))                                   &
              *me*me*FA_/delta(i_ome,i_e1,i_e1)            
    END DO ! i_e1
   END DO ! i_ome

!-------------------------------
!   energy_in \= energy_out
!-------------------------------
    DO i_ome = 1, D_ome
     DO  i_e1 = 1, D_e
      DO  i_e2 = 1, D_e

      IF( i_e1.ne.i_e2 ) THEN      

      CALL NESKernGvalue( etap(i_e1,i_e2), &
                          eta, y0(i_ome,i_e1,i_e2), G0, G1, G2 )

      COMBO1=G2+2.0*y0(i_ome,i_e1,i_e2)*G1&
             +y0(i_ome,i_e1,i_e2)*y0(i_ome,i_e1,i_e2)*G0                                              
      COMBO2=G1+y0(i_ome,i_e1,i_e2)*G0                                                           

      I1(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e1)*esq(i_e2,i_e2)                            &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                     &
              /delta(i_ome,i_e1,i_e2)**5                                 &
              *fgamm(i_e1,i_e2)                                          &
              *(A(i_ome,i_e1,i_e2)*tsq*COMBO1                            &
               +B(i_ome,i_e1,i_e2)*T*COMBO2                              &
               +C(i_ome,i_e1,i_e2)*G0                                    &
               )                                                         
                                  
      I2(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e1)*esq(i_e2,i_e2)                            &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                     &
              /delta(i_ome,i_e1,i_e2)**5                                 &
              *fgamm(i_e1,i_e2)                                          &
              *(A(i_ome,i_e1,i_e2)*tsq*COMBO1                            &
               -B(i_ome,i_e2,i_e1)*T*COMBO2                              &
               +C(i_ome,i_e2,i_e1)*G0                                    &
               )                                                         
                                                                                
      I3(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e2)*(1.0-omega(i_ome))                        &
             *me*me*fgamm(i_e1,i_e2)*G0/delta(i_ome,i_e1,i_e2)             

      ELSE
      END IF

      buffer1                   = ( beta1 * I1(i_ome,i_e1,i_e2)    &
                                        + beta2 * I2(i_ome,i_e1,i_e2)    &
                                        + beta3 * I3(i_ome,i_e1,i_e2) )  &
                                      /  esq(i_e1,i_e2)
      buffer2 = LOG10(buffer1) + log10_cons 

      nesktab(i_ome,i_e1,i_e2) = 10**buffer2

      END DO ! i_e2
     END DO ! i_e1
    END DO ! i_ome

   END SUBROUTINE NESKern

   SUBROUTINE NESKernFsame( eta, y0, FA_, FA0, FA1) 

     REAL(dp), INTENT(in)    :: eta, y0
     REAL(dp), INTENT(out)   :: FA_, FA0, FA1
     REAL(dp)                 :: FEXP
     REAL(dp)                 :: x1, x2, x2inv
     REAL(dp)                 :: y, sumFA
     REAL(dp), DIMENSION(27)  :: Tarr, CH1arr
     INTEGER                  :: i_t

     FA_ = 1.0/( FEXP( -(eta - y0) ) + 1.0) 
  
     CH1arr(1)  = 1.935064300869969
     CH1arr(2)  = 0.166073032927855
     CH1arr(3)  = 0.024879329924228
     CH1arr(4)  = 0.004686361959447
     CH1arr(5)  = 0.001001627496164
     CH1arr(6)  = 0.000232002196094
     CH1arr(7)  = 0.000056817822718
     CH1arr(8)  = 0.000014496300557
     CH1arr(9)  = 0.000003816329463
     CH1arr(10) = 0.000001029904264
     CH1arr(11) = 0.000000283575385
     CH1arr(12) = 0.000000079387055
     CH1arr(13) = 0.000000022536705
     CH1arr(14) = 0.000000006474338
     CH1arr(15) = 0.000000001879117
     CH1arr(16) = 0.000000000550291
     CH1arr(17) = 0.000000000162421
     CH1arr(18) = 0.000000000048274
     CH1arr(19) = 0.000000000014437
     CH1arr(20) = 0.000000000004342
     CH1arr(21) = 0.000000000001312
     CH1arr(22) = 0.000000000000398
     CH1arr(23) = 0.000000000000121
     CH1arr(24) = 0.000000000000037
     CH1arr(25) = 0.000000000000011
     CH1arr(26) = 0.000000000000004
     CH1arr(27) = 0.000000000000001

     x1 = eta -y0
     x2 = -FEXP(x1)
     x2inv = 1.0/x2
      
     FA0 = LOG(1.0 - x2)

     IF(x2 .lt. -1.0) x2 = x2inv ! x1 .lt. 0
     y = (4.0*x2+1.0)/3.0

     Tarr(1) = 1.0  ! T0
     Tarr(2) = y    ! T1
     DO i_t = 3, SIZE(Tarr)
       Tarr(i_t) = 2.0*y*Tarr(i_t-1)-Tarr(i_t-2)
     END DO
    
     sumFA = 0.0
     DO i_t = 2, SIZE(Tarr)
       sumFA = CH1arr(i_t) * Tarr(i_t) + sumFA
     END DO
     FA1 = &
           -x2*( sumFA + 0.5*CH1arr(1)*Tarr(1) )

     IF( (1.0/x2inv) .lt. -1.0) FA1 = &  ! x2 .lt. -1
                               -FA1 + 0.5*x1**2 + pi*pi/6.0

   END SUBROUTINE NESKernFsame

   SUBROUTINE NESKernFdiff&
              ( etap, eta, y0, FA00, FA01, FA10, FA11, FA20, FA21)

     REAL(dp), INTENT(in)     :: etap, eta, y0
     REAL(dp), INTENT(out)    :: FA00, FA01, FA10, FA11, FA20, FA21

     REAL(dp)                 :: x1, x2, x2inv, FEXP
     REAL(dp)                 :: y1, sumFA
     REAL(dp), DIMENSION(27)  :: Tarr, CH1arr
     REAL(dp), DIMENSION(24)  :: CH2arr
     INTEGER                  :: i_t

     CH1arr(1)  = 1.935064300869969
     CH1arr(2)  = 0.166073032927855
     CH1arr(3)  = 0.024879329924228
     CH1arr(4)  = 0.004686361959447
     CH1arr(5)  = 0.001001627496164
     CH1arr(6)  = 0.000232002196094
     CH1arr(7)  = 0.000056817822718
     CH1arr(8)  = 0.000014496300557
     CH1arr(9)  = 0.000003816329463
     CH1arr(10) = 0.000001029904264
     CH1arr(11) = 0.000000283575385
     CH1arr(12) = 0.000000079387055
     CH1arr(13) = 0.000000022536705
     CH1arr(14) = 0.000000006474338
     CH1arr(15) = 0.000000001879117
     CH1arr(16) = 0.000000000550291
     CH1arr(17) = 0.000000000162421
     CH1arr(18) = 0.000000000048274
     CH1arr(19) = 0.000000000014437
     CH1arr(20) = 0.000000000004342
     CH1arr(21) = 0.000000000001312
     CH1arr(22) = 0.000000000000398
     CH1arr(23) = 0.000000000000121
     CH1arr(24) = 0.000000000000037
     CH1arr(25) = 0.000000000000011
     CH1arr(26) = 0.000000000000004
     CH1arr(27) = 0.000000000000001

     CH2arr(1)  = 1.958417213383495
     CH2arr(2)  = 0.085188131486831
     CH2arr(3)  = 0.008559852220133
     CH2arr(4)  = 0.001211772144129
     CH2arr(5)  = 0.000207227685308
     CH2arr(6)  = 0.000039969586914 
     CH2arr(7)  = 0.000008380640658
     CH2arr(8)  = 0.000001868489447
     CH2arr(9)  = 0.000000436660867
     CH2arr(10) = 0.000000105917334
     CH2arr(11) = 0.000000026478920 
     CH2arr(12) = 0.000000006787000
     CH2arr(13) = 0.000000001776536 
     CH2arr(14) = 0.000000000473417
     CH2arr(15) = 0.000000000128121
     CH2arr(16) = 0.000000000035143
     CH2arr(17) = 0.000000000009754   
     CH2arr(18) = 0.000000000002736   
     CH2arr(19) = 0.000000000000775 
     CH2arr(20) = 0.000000000000221 
     CH2arr(21) = 0.000000000000064  
     CH2arr(22) = 0.000000000000018 
     CH2arr(23) = 0.000000000000005
     CH2arr(24) = 0.000000000000002

     x1 = eta - y0
     x2 = - FEXP(x1)
     x2inv = 1.0/x2

     FA00 = LOG(1.0-x2)
    
     IF( x2 .lt. -1.0) x2 = x2inv  ! x1 .lt. 0
     y1 = (4.0*x2+1.0)/3.0

     Tarr(1) = 1.0
     Tarr(2) = y1
     DO i_t = 3, SIZE(Tarr)
       Tarr(i_t) = 2.0*y1*Tarr(i_t-1)-Tarr(i_t-2)
     END DO

     sumFA = 0.0
     DO i_t = 1, SIZE(CH2arr)
       sumFA = sumFA + CH2arr(i_t)*Tarr(i_t) 
     END DO

     FA20 = -2.0*x2*( sumFA - 0.5*CH2arr(1)*Tarr(1) )
    
     sumFA = 0.0 
     DO i_t = 1, SIZE(CH1arr)
       sumFA = sumFA + CH1Arr(i_t)*Tarr(i_t)
     END DO 

     FA10 = -x2*( sumFA - 0.5*CH1arr(1)*Tarr(1) )

     IF( (1.0/x2inv).lt.-1.0) THEN  ! x2 .lt. -1
       FA10 = -FA10 + 0.5*x1*x1+pi*pi/6.0
       FA20 =  FA20 + pi*pi*x1/3.0 + x1*x1*x1/3.0
     END IF

     x1 = etap - y0 
     x2 = - FEXP(x1)
     x2inv = 1.0/x2

     FA01 = LOG(1.0-x2)
     IF( x2 .lt. -1.0) x2 = x2inv
     y1 = (4.0*x2+1.0)/3.0

     Tarr(1) = 1.0
     Tarr(2) = y1
     DO i_t = 3, SIZE(Tarr)
       Tarr(i_t) = 2.0*y1*Tarr(i_t-1)-Tarr(i_t-2)
     END DO

     sumFA = 0.0
     DO i_t = 1, SIZE(CH2arr)
       sumFA = sumFA + CH2arr(i_t)*Tarr(i_t)
     END DO

     FA21 = -2.0*x2*( sumFA - 0.5*CH2arr(1)*Tarr(1) )

     sumFA = 0.0
     DO i_t = 1, SIZE(CH1arr)
       sumFA = sumFA + CH1Arr(i_t)*Tarr(i_t)
     END DO

     FA11 = -x2*( sumFA - 0.5*CH1arr(1)*Tarr(1) )

     IF( (1.0/x2inv).lt.-1.0) THEN
       FA11 = -FA11 + 0.5*x1*x1+pi*pi/6.0
       FA21 =  FA21 + pi*pi*x1/3.0 + x1*x1*x1/3.0
     END IF

   END SUBROUTINE NESKernFdiff

   SUBROUTINE NESKernGvalue( etap, eta, y0, G0, G1, G2 )

    REAL(dp), INTENT(in)      :: etap, eta, y0
    REAL(dp), INTENT(out)     :: G0, G1, G2
    REAL(dp)                 :: FA00, FA01, FA10, FA11, FA20, FA21

    CALL NESKernFdiff( etap, eta, y0, FA00, FA01, FA10, FA11, FA20, FA21)

    G0 = FA01 - FA00
    G1 = FA11 - FA10
    G2 = FA21 - FA20

  END SUBROUTINE NESKernGvalue


  SUBROUTINE etaxx( rho, T, xn, xp, etann, etapp )
  

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

  REAL(dp), INTENT(in)    :: rho           ! density (g/cm3)
  REAL(dp), INTENT(in)    :: T             ! temperature [K]
  REAL(dp), INTENT(in)    :: xn            ! neutron mass fraction
  REAL(dp), INTENT(in)    :: xp            ! proton mass fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

  REAL(dp), INTENT(out)   :: etann         ! neutron number corrected for blocking
  REAL(dp), INTENT(out)   :: etapp         ! proton number corrected for blocking

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

  REAL(dp)                :: nn            ! neutron number uncorrected for blocking (cm^{-3})
  REAL(dp)                :: np            ! proton number uncorrected for blocking (cm^{-3})
  REAL(dp)                :: d_n            ! neutron number uncorrected for blocking (fm^{-3})
  REAL(dp)                :: d_p            ! proton number uncorrected for blocking (fm^{-3})
  REAL(dp)                :: efn           ! degenerate expression
  REAL(dp)                :: efp           ! degenerate expression
  REAL(dp)                :: etanndgnt     ! nondegenerate expression
  REAL(dp)                :: etappdgnt     ! nondegenerate expression
  REAL(dp)                :: mpG, mnG
  REAL(dp), PARAMETER     :: tthird = 2.d0/3.d0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  nn, np
!-----------------------------------------------------------------------
  
  mpG   = mp * ergmev * cvel_inv * cvel_inv ! proton mass [g]
  mnG   = mn * ergmev * cvel_inv * cvel_inv ! neutron mass [g]
  nn                 = xn * rho/mpG
  np                 = xp * rho/mnG

  IF ( nn <= zero  .or.  np <= zero ) THEN
    WRITE(*,*) "ERROR! nn or np less than zero."
  END IF

!-----------------------------------------------------------------------
!  etann, etanp (analytic approximation)
!-----------------------------------------------------------------------

  d_n          = nn * 1.d-39
  d_p          = np * 1.d-39
  efn          = ( hbarc**2/( 2.d+00 * mn ) ) &
                    * ( 3.d+00 * pi**2 * d_n )**tthird
  efp          = ( hbarc**2/( 2.d+00 * mp ) ) &
                    * ( 3.d+00 * pi**2 * d_p )**tthird
  etanndgnt    = 1.5d+00 * ( kMev * T/efn )
  etappdgnt    = 1.5d+00 * ( kMev * T/efp )
  etann        = nn * etanndgnt/DSQRT( 1.d+00 + etanndgnt**2 )
  etapp        = np * etappdgnt/DSQRT( 1.d+00 + etappdgnt**2 )

  END SUBROUTINE etaxx


  FUNCTION IntegralESNuclei( a, l, nquad )

  REAL(dp), INTENT(in)       :: a
  INTEGER,  INTENT(in)       :: l, nquad

  REAL(dp), DIMENSION(nquad) :: roots, weights
  REAL(dp)                   :: IntegralESNuclei, buffer, func
  INTEGER                    :: ii

  CALL gaquad( nquad, roots, weights, -1.0_dp , 1.0_dp )

  buffer = 0.0

  IF ( l == 0 ) THEN
 
   DO ii = 1, nquad
     func = ( 1.0 + roots(ii) )* EXP( a * ( roots(ii) - 1.0 ) )
     buffer = buffer + weights(ii) * func 
   END DO

  ELSE IF ( l == 1 ) THEN

   DO ii = 1, nquad
     func = roots(ii) * ( roots(ii) + 1.0 ) * EXP ( a * ( roots(ii) - 1.0 ) )
     buffer = buffer + weights(ii) * func
   END DO

  ELSE 
    WRITE(*,*) "ERROR when calling IntegralESNuclei function."
  END IF

  IntegralESNuclei = buffer

  END FUNCTION IntegralESNuclei

END MODULE B85
