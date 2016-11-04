MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  implicit none

  PUBLIC LogInterpolateSingleVariable
  PUBLIC LogInterpolateAllVariables
  PUBLIC LogInterpolateDifferentiateSingleVariable
  PUBLIC LogInterpolateDifferentiateAllVariables
  PUBLIC MonotonicityCheck
  PUBLIC GetGamma1
  PUBLIC ComputeTempFromIntEnergy
  PUBLIC ComputeTempFromEntropy
  PUBLIC EOSTableQuery
  !PUBLIC InterpolateEOS

  REAL(dp), PARAMETER :: ln10 = LOG(10.d0)

  INTERFACE LogInterpolateSingleVariable
    MODULE PROCEDURE LogInterpolateSingleVariable_3D
    MODULE PROCEDURE LogInterpolateSingleVariable_3D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_4D
    MODULE PROCEDURE LogInterpolateSingleVariable_4D_Custom
  END INTERFACE LogInterpolateSingleVariable

  INTERFACE LogInterpolateAllVariables
    MODULE PROCEDURE LogInterpolateAllVariables_3D
    MODULE PROCEDURE LogInterpolateAllVariables_3D_Custom
  END INTERFACE LogInterpolateAllVariables

  INTERFACE LogInterpolateDifferentiateSingleVariable
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_3D
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_4D
  END INTERFACE LogInterpolateDifferentiateSingleVariable

CONTAINS


  SUBROUTINE locate( xx, n, x, j )

    INTEGER, INTENT(in)      :: n
    INTEGER, INTENT(out)     :: j
    REAL(dp), INTENT(in)     :: x,xx(n)
    INTEGER                  :: jl,jm,ju

    jl = 0
    ju = n+1
    DO WHILE ( ju - jl > 1 )
      jm = (ju+jl)/2
      IF ((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) THEN
        jl = jm
      ELSE
        ju = jm
      END IF
    END DO

    IF (x.eq.xx(1)) THEN
      j = 1
    ELSEIF (x.eq.xx(n)) THEN
      j = n-1  
    ELSE
      j = jl
    END IF

  END SUBROUTINE locate


  PURE INTEGER FUNCTION Index1D( x, xx, n )

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    INTEGER :: il, im, iu

    il = 0
    iu = n+1
    DO WHILE ( iu - il > 1 )
      im = (iu+il)/2
      IF ((xx(n).ge.xx(1)).eqv.(x.ge.xx(im))) THEN
        il = im
      ELSE
        iu = im
      END IF
    END DO

    IF (x.eq.xx(1)) THEN
      Index1D = 1
    ELSEIF (x.eq.xx(n)) THEN
      Index1D = n-1
    ELSE
      Index1D = il
    END IF

    RETURN
  END FUNCTION Index1D


  PURE REAL(dp) FUNCTION TriLinear &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2, dX3 )

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX2, dX3

    REAL(dp) :: ddX1, ddX2, ddX3

    ddX1 = 1.0_dp - dX1
    ddX2 = 1.0_dp - dX2
    ddX3 = 1.0_dp - dX3

    TriLinear                                        &
      = ddX3                                         &
         * (   ddX2 * ( ddX1 * p000 + dX1 * p100 )   &
             +  dX2 * ( ddX1 * p010 + dX1 * p110 ) ) &
      +  dX3                                         &
         * (   ddX2 * ( ddX1 * p001 + dX1 * p101 )   &
             +  dX2 * ( ddX1 * p011 + dX1 * p111 ) )

    RETURN
  END FUNCTION TriLinear


  PURE REAL(dp) FUNCTION TetraLinear &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3, dX4 )

    REAL(dp), INTENT(in) :: &
      p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3, dX4

    REAL(dp) :: ddX1, ddX2, ddX3, ddX4

    ddX1 = 1.0_dp - dX1
    ddX2 = 1.0_dp - dX2
    ddX3 = 1.0_dp - dX3
    ddX4 = 1.0_dp - dX4

    TetraLinear                                                    &
      = ddX4                                                       &
         * (   ddX3 * (  ddX2 * ( ddX1 * p0000 + dX1 * p1000 )     &
                        + dX2 * ( ddX1 * p0100 + dX1 * p1100 ) )   &
             +  dX3 * (  ddX2 * ( ddX1 * p0010 + dX1 * p1010 )     &
                        + dX2 * ( ddX1 * p0110 + dX1 * p1110 ) ) ) &
      +  dX4                                                       &
         * (   ddX3 * (  ddX2 * ( ddX1 * p0001 + dX1 * p1001 )     &
                        + dX2 * ( ddX1 * p0101 + dX1 * p1101 ) )   &
             +  dX3 * (  ddX2 * ( ddX1 * p0011 + dX1 * p1011 )     &
                        + dX2 * ( ddX1 * p0111 + dX1 * p1111 ) ) )

    RETURN
  END FUNCTION TetraLinear


  SUBROUTINE LogInterpolateSingleVariable_3D &
               ( x1, x2, x3, Coordinate1, Coordinate2, Coordinate3, &
                 LogInterp, Offset, Table, Interpolant, MaskVar )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp 
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    REAL(dp), INTENT(in)                :: Offset
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize
    REAL(dp), DIMENSION(:), INTENT(out) :: Interpolant 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: i, j, k, il1, il2, il3
  
    epsilon = 1.d-200

    ! Check the usage of x2
    ! x2 is used to size the run of the do loop because the input is expected 
    ! to be a series of 3-tuple rho, T, Ye points; 
    !     so SIZE(x1) = SIZE(x2) = SIZE(x3) 
    ! If this subroutine is being used to make an array, make sure the other 
    ! inputs are the same size (i.e. x1(1) = x(2) = x(3) etc.)
    !-------------------------------

    Masksize = SIZE( x2 )
    ALLOCATE( work_mask( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, SIZE( x2 )
      
      IF ( .not.work_mask(i) ) CYCLE

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )

      p000 = ( Table( il1  , il2  , il3   ) )
      p100 = ( Table( il1+1, il2  , il3   ) )
      p010 = ( Table( il1  , il2+1, il3   ) )
      p110 = ( Table( il1+1, il2+1, il3   ) )
      p001 = ( Table( il1  , il2  , il3+1 ) )
      p101 = ( Table( il1+1, il2  , il3+1 ) )
      p011 = ( Table( il1  , il2+1, il3+1 ) )
      p111 = ( Table( il1+1, il2+1, il3+1 ) )

      IF ( LogInterp(1) == 1 ) THEN 
        delta(1) = LOG10( x1(i) / Coordinate1(il1) ) &
                     / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        delta(1) = ( x1(i) - Coordinate1(il1) ) &
                     / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN 
        delta(2) = LOG10( x2(i) / Coordinate2(il2) ) &
                     / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        delta(2) = ( x2(i) - Coordinate2(il2) ) &
                     / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN 
        delta(3) = LOG10( x3(i) / Coordinate3(il3) ) &
                     / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        delta(3) = ( x3(i) - Coordinate3(il3) ) &
                     / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      Interpolant(i) &
        = 10.d0**( &
            (1.0_dp - delta(3)) &
              * (   (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &
                  +           delta(1)  * (1.0_dp - delta(2)) * p100   &
                  + (1.0_dp - delta(1)) *           delta(2)  * p010   &
                  +           delta(1)  *           delta(2)  * p110 ) &
            +         delta(3) &
              * (   (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                  +           delta(1)  * (1.0_dp - delta(2)) * p101   &
                  + (1.0_dp - delta(1)) *           delta(2)  * p011   &
                  +           delta(1)  *           delta(2)  * p111 ) ) &
          - Offset
   
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_3D


  SUBROUTINE LogInterpolateSingleVariable_3D_Custom &
               ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: D,  T,  Y
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Ds, Ts, Ys
    REAL(dp),                   INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),     INTENT(out) :: Interpolant

    INTEGER  :: &
      iP, iD, iT, iY
    REAL(dp) :: &
      dD, dT, dY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111

    IF( .NOT. ALL( [ SIZE(T), SIZE(Y) ] == SIZE(D) ) )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'LogInterpolateSingleVariable_3D_Custom'
      WRITE(*,'(A4,A)') &
        '', 'ERROR: arrays of interpolation points have different sizes'
      WRITE(*,*)
      RETURN
    END IF

    DO iP = 1, SIZE( D )

      iD = Index1D( D(iP), Ds, SIZE( Ds ) )
      iT = Index1D( T(iP), Ts, SIZE( Ts ) )
      iY = Index1D( Y(iP), Ys, SIZE( Ys ) )

      dD = LOG10( D(iP) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
      dT = LOG10( T(iP) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
      dY = ( Y(iP) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      p000 = ( Table( iD  , iT  , iY   ) )
      p100 = ( Table( iD+1, iT  , iY   ) )
      p010 = ( Table( iD  , iT+1, iY   ) )
      p110 = ( Table( iD+1, iT+1, iY   ) )
      p001 = ( Table( iD  , iT  , iY+1 ) )
      p101 = ( Table( iD+1, iT  , iY+1 ) )
      p011 = ( Table( iD  , iT+1, iY+1 ) )
      p111 = ( Table( iD+1, iT+1, iY+1 ) )

      Interpolant(iP) &
        = 10.0d0**( &
            TriLinear &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, dD, dT, dY ) ) - OS

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_3D_Custom


  SUBROUTINE LogInterpolateSingleVariable_4D &
               ( x1, x2, x3, x4, Coordinate1, Coordinate2, Coordinate3, &
                 Coordinate4, LogInterp, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: x1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x4
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate4
    INTEGER,  DIMENSION(4),       INTENT(in)  :: LogInterp
    REAL(dp),                     INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),       INTENT(out) :: Interpolant

    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111,&
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
    REAL(dp), DIMENSION(4) :: alpha, delta
    INTEGER :: i, il1, il2, il3, il4

    IF ( ( SIZE(x1) .NE. SIZE(x2) ) .OR. ( SIZE(x2) .NE. SIZE(x3) ) &
         .OR. ( SIZE(x3) .NE. SIZE(x4) ) .OR. ( SIZE(x4) .NE. SIZE (x1) ) ) &
    THEN
      WRITE(*,*) &
        'ERROR: describe arrays (of interpolation point) have diff size.'
      RETURN
    END IF
    
    DO i = 1, SIZE( x1 )

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )
      CALL locate( Coordinate4, SIZE( Coordinate4 ), x4(i), il4 )

      p0000 = ( Table( il1  , il2  , il3  , il4   ) )
      p0001 = ( Table( il1  , il2  , il3  , il4+1 ) )
      p0010 = ( Table( il1  , il2  , il3+1, il4   ) )
      p0011 = ( Table( il1  , il2  , il3+1, il4+1 ) )
      p0100 = ( Table( il1  , il2+1, il3  , il4   ) )
      p0101 = ( Table( il1  , il2+1, il3  , il4+1 ) )
      p0110 = ( Table( il1  , il2+1, il3+1, il4   ) )
      p0111 = ( Table( il1  , il2+1, il3+1, il4+1 ) )
      p1000 = ( Table( il1+1, il2  , il3  , il4   ) )
      p1001 = ( Table( il1+1, il2  , il3  , il4+1 ) )
      p1010 = ( Table( il1+1, il2  , il3+1, il4   ) )
      p1011 = ( Table( il1+1, il2  , il3+1, il4+1 ) )
      p1100 = ( Table( il1+1, il2+1, il3  , il4   ) )
      p1101 = ( Table( il1+1, il2+1, il3  , il4+1 ) )
      p1110 = ( Table( il1+1, il2+1, il3+1, il4   ) )
      p1111 = ( Table( il1+1, il2+1, il3+1, il4+1 ) )

      IF ( LogInterp(1) == 1 ) THEN
        delta(1) = LOG10( x1(i) / Coordinate1(il1) ) &
                    / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        delta(1) = ( x1(i) - Coordinate1(il1) ) &
                     / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        delta(2) = LOG10( x2(i) / Coordinate2(il2) ) &
                     / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        delta(2) = ( x2(i) - Coordinate2(il2) ) &
                     / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        delta(3) = LOG10( x3(i) / Coordinate3(il3) ) &
                     / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        delta(3) = ( x3(i) - Coordinate3(il3) ) &
                     / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      IF ( LogInterp(4) == 1 ) THEN
        delta(4) = LOG10( x4(i) / Coordinate4(il4) ) &
                     / LOG10( Coordinate4(il4+1) / Coordinate4(il4) )
      ELSE
        delta(4) = ( x4(i) - Coordinate4(il4) ) &
                     / ( Coordinate4(il4+1) - Coordinate4(il4) )
      END IF

      Interpolant(i) &
      = 10.d0**( &
          (1.0_dp - delta(4)) &
            * (   (1.0_dp - delta(3)) * &
                        (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0000   &
                + (1.0_dp - delta(3)) * &
                        (1.0_dp - delta(2)) *           delta(1)  * p1000   &
                + (1.0_dp - delta(3)) * &
                                  delta(2)  * (1.0_dp - delta(1)) * p0100   &
                + (1.0_dp - delta(3)) * &
                                  delta(2)  *           delta(1)  * p1100   &
                +           delta(3)  * &
                        (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0010   &
                +           delta(3)  * &                    
                        (1.0_dp - delta(2)) *           delta(1)  * p1010   &
                +           delta(3)  * &
                                  delta(2)  * (1.0_dp - delta(1)) * p0110   &
                +           delta(3)  * &
                                  delta(2)  *           delta(1)  * p1110 ) &
          +         delta(4)  &
            * (   (1.0_dp - delta(3)) * &                
                        (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0001   &
                + (1.0_dp - delta(3)) * &                    
                        (1.0_dp - delta(2)) *           delta(1)  * p1001   &
                + (1.0_dp - delta(3)) * &
                                  delta(2)  * (1.0_dp - delta(1)) * p0101   &
                + (1.0_dp - delta(3)) * &
                                  delta(2)  *           delta(1)  * p1101   &
                +           delta(3)  * &
                        (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0011   &
                +           delta(3)  * &
                        (1.0_dp - delta(2)) *           delta(1)  * p1011   &
                +           delta(3)  * &
                                  delta(2)  * (1.0_dp - delta(1)) * p0111   &
                +           delta(3)  * &
                                  delta(2)  *           delta(1)  * p1111 ) ) &
        - Offset

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_4D


  SUBROUTINE LogInterpolateSingleVariable_4D_Custom &
               ( E, D, T, Y, Es, Ds, Ts, Ys, OS, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: E,  D,  T,  Y
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Es, Ds, Ts, Ys
    REAL(dp),                     INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),       INTENT(out) :: Interpolant

    INTEGER  :: &
      iP, iE, iD, iT, iY
    REAL(dp) :: &
      p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111, &
      dE, dD, dT, dY

    IF( .NOT. ALL( [ SIZE(D), SIZE(T), SIZE(Y) ] == SIZE(E) ) )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'LogInterpolateSingleVariable_4D_Custom'
      WRITE(*,'(A4,A)') &
        '', 'ERROR: arrays of interpolation points have different sizes'
      WRITE(*,*)
      RETURN
    END IF

    DO iP = 1, SIZE( E )

      iE = Index1D( E(iP), Es, SIZE( Es ) )
      iD = Index1D( D(iP), Ds, SIZE( Ds ) )
      iT = Index1D( T(iP), Ts, SIZE( Ts ) )
      iY = Index1D( Y(iP), Ys, SIZE( Ys ) )

      p0000 = Table( iE  , iD  , iT  , iY   )
      p1000 = Table( iE+1, iD  , iT  , iY   )
      p0100 = Table( iE  , iD+1, iT  , iY   )
      p1100 = Table( iE+1, iD+1, iT  , iY   )
      p0010 = Table( iE  , iD  , iT+1, iY   )
      p1010 = Table( iE+1, iD  , iT+1, iY   )
      p0110 = Table( iE  , iD+1, iT+1, iY   )
      p1110 = Table( iE+1, iD+1, iT+1, iY   )
      p0001 = Table( iE  , iD  , iT  , iY+1 )
      p1001 = Table( iE+1, iD  , iT  , iY+1 )
      p0101 = Table( iE  , iD+1, iT  , iY+1 )
      p1101 = Table( iE+1, iD+1, iT  , iY+1 )
      p0011 = Table( iE  , iD  , iT+1, iY+1 )
      p1011 = Table( iE+1, iD  , iT+1, iY+1 )
      p0111 = Table( iE  , iD+1, iT+1, iY+1 )
      p1111 = Table( iE+1, iD+1, iT+1, iY+1 )

      dE = LOG10( E(iP) / Es(iE) ) / LOG10( Es(iE+1) / Es(iE) )
      dD = LOG10( D(iP) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
      dT = LOG10( T(iP) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
      dY = ( Y(iP) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      Interpolant(iP) &
        = 10.0d0**( &
            TetraLinear &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dE, dD, dT, dY ) ) - OS

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_4D_Custom


  SUBROUTINE LogInterpolateAllVariables_3D &
               ( x1, x2, x3, LogInterp, TS, DV, Interpolants, MaskVar )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar

    REAL(dp), DIMENSION(:,:), INTENT(out) :: Interpolants 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: i, j, il1, il2, il3
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize
  
    epsilon = 1.d-200

    Masksize = SIZE( x2 )
    ALLOCATE( work_mask( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, SIZE( x2 )

      IF ( .not.work_mask(i) ) CYCLE

      ASSOCIATE( Coordinate1 => TS % States(1) % Values, &    
                 Coordinate2 => TS % States(2) % Values, &    
                 Coordinate3 => TS % States(3) % Values )   

      CALL locate( Coordinate1, SIZE(Coordinate1), x1(i), il1 )
      CALL locate( Coordinate2, SIZE(Coordinate2), x2(i), il2 )
      CALL locate( Coordinate3, SIZE(Coordinate3), x3(i), il3 )


      IF ( LogInterp(1) == 1 ) THEN
        delta(1) = LOG10( x1(i) / Coordinate1(il1) ) / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        delta(1) = ( x1(i) - Coordinate1(il1) ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        delta(2) = LOG10( x2(i) / Coordinate2(il2) ) / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        delta(2) = ( x2(i) - Coordinate2(il2) ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        delta(3) = LOG10( x3(i) / Coordinate3(il3) ) / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        delta(3) = ( x3(i) - Coordinate3(il3) ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      END ASSOCIATE

      DO j = 1, DV % nVariables

      
        ASSOCIATE( Table => DV % Variables(j) % Values(:,:,:), &
                   Offset => DV % Offsets(j) )
      
        p000 = ( Table( il1  , il2  , il3   ) )
        p100 = ( Table( il1+1, il2  , il3   ) )
        p010 = ( Table( il1  , il2+1, il3   ) )
        p110 = ( Table( il1+1, il2+1, il3   ) )
        p001 = ( Table( il1  , il2  , il3+1 ) )
        p101 = ( Table( il1+1, il2  , il3+1 ) )
        p011 = ( Table( il1  , il2+1, il3+1 ) )
        p111 = ( Table( il1+1, il2+1, il3+1 ) )
  
        Interpolants(i,j) &
          = 10.d0**( &
                (1.0_dp - delta(3)) * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &
                                     +            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                     + ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                     +            delta(1)  *           delta(2)  * p110 ) &
                        + delta(3)  * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                     +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                     +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                     +            delta(1)  *           delta(2)  * p111 ) &
  
                   ) - Offset
        
        END ASSOCIATE
    
      END DO
    END DO

  END SUBROUTINE LogInterpolateAllVariables_3D


  SUBROUTINE LogInterpolateAllVariables_3D_Custom &
               ( D, T, Y, Ds, Ts, Ys, DV, Interpolants )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: D,  T,  Y
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Ds, Ts, Ys
    TYPE(DependentVariablesType), INTENT(in)  :: DV
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Interpolants

    INTEGER :: &
      iP, iV, iD, iT, iY
    REAL(dp) :: &
      dD, dT, dY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111

    IF( .NOT. ALL( [ SIZE(T), SIZE(Y) ] == SIZE(D) ) )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'LogInterpolateAllVariables_3D_Custom'
      WRITE(*,'(A4,A)') &
        '', 'ERROR: arrays of interpolation points have different sizes'
      WRITE(*,*)
      RETURN
    END IF

    DO iP = 1, SIZE( D )

      iD = Index1D( D(iP), Ds, SIZE( Ds ) )
      iT = Index1D( T(iP), Ts, SIZE( Ts ) )
      iY = Index1D( Y(iP), Ys, SIZE( Ys ) )

      dD = LOG10( D(iP) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
      dT = LOG10( T(iP) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
      dY = ( Y(iP) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      DO iV = 1, DV % nVariables

        ASSOCIATE &
          ( Table => DV % Variables(iV) % Values(:,:,:), &
            OS    => DV % Offsets  (iV) )

        p000 = Table( iD  , iT  , iY   )
        p100 = Table( iD+1, iT  , iY   )
        p010 = Table( iD  , iT+1, iY   )
        p110 = Table( iD+1, iT+1, iY   )
        p001 = Table( iD  , iT  , iY+1 )
        p101 = Table( iD+1, iT  , iY+1 )
        p011 = Table( iD  , iT+1, iY+1 )
        p111 = Table( iD+1, iT+1, iY+1 )

        Interpolants(iV, iP) &
          = 10.0d0**( &
              TriLinear &
                ( p000, p100, p010, p110, &
                  p001, p101, p011, p111, dD, dT, dY ) ) - OS

        END ASSOCIATE ! Table, etc.

      END DO

    END DO

  END SUBROUTINE LogInterpolateAllVariables_3D_Custom


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D                    &
               ( x1, x2, x3, Coordinate1, Coordinate2, Coordinate3,          &
                 LogInterp, Offset, Table, Interpolant, Derivative, MaskVar)     

    REAL(dp), DIMENSION(:), INTENT(in)     :: x1
    REAL(dp), DIMENSION(:), INTENT(in)     :: x2
    REAL(dp), DIMENSION(:), INTENT(in)     :: x3
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate3
    INTEGER, DIMENSION(3), INTENT(in)      :: LogInterp 
    REAL(dp), INTENT(in)                   :: Offset
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    REAL(dp), DIMENSION(:), INTENT(out)    :: Interpolant
    REAL(dp), DIMENSION(:,:), INTENT(out)  :: Derivative 
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize
    
    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: alpha, delta
    INTEGER :: i, j, k, il1, il2, il3

    epsilon = 1.d-200

    Masksize = SIZE( x2 )
    ALLOCATE( work_mask( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, SIZE( x2 )

      IF ( .not.work_mask(i) ) CYCLE

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )

      p000 = ( Table( il1  , il2  , il3   ) )
      p100 = ( Table( il1+1, il2  , il3   ) )
      p010 = ( Table( il1  , il2+1, il3   ) )
      p110 = ( Table( il1+1, il2+1, il3   ) )
      p001 = ( Table( il1  , il2  , il3+1 ) )
      p101 = ( Table( il1+1, il2  , il3+1 ) )
      p011 = ( Table( il1  , il2+1, il3+1 ) )
      p111 = ( Table( il1+1, il2+1, il3+1 ) )

      IF ( LogInterp(1) == 1 ) THEN
        alpha(1) &
          = 1.0d0 / ( x1(i) * LOG10( Coordinate1(il1+1) / Coordinate1(il1) ) )
        delta(1) &
          = LOG10( x1(i) / Coordinate1(il1) ) &
              / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        alpha(1) &
          = ln10 / ( Coordinate1(il1+1) - Coordinate1(il1) )
        delta(1) &
          = ( x1(i) - Coordinate1(il1) ) &
              / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        alpha(2) &
          = 1.0d0 / ( x2(i) * LOG10( Coordinate2(il2+1) / Coordinate2(il2) ) )
        delta(2) &
          = LOG10( x2(i) / Coordinate2(il2) ) &
              / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        alpha(2) &
          = ln10 / ( Coordinate2(il2+1) - Coordinate2(il2) )
        delta(2) &
          = ( x2(i) - Coordinate2(il2) ) &
              / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        alpha(3) &
          = 1.0d0 / ( x3(i) * LOG10( Coordinate3(il3+1) / Coordinate3(il3) ) )
        delta(3) &
          = LOG10( x3(i) / Coordinate3(il3) ) &
              / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        alpha(3) &
          = ln10 / ( Coordinate3(il3+1) - Coordinate3(il3) )
        delta(3) &
          = ( x3(i) - Coordinate3(il3) ) &
              / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      Interpolant(i) &
        = 10.d0**( &
            (1.0_dp - delta(3)) &
              * (   (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &
                  +           delta(1)  * (1.0_dp - delta(2)) * p100   &
                  + (1.0_dp - delta(1)) *           delta(2)  * p010   &
                  +           delta(1)  *           delta(2)  * p110 ) &
            +         delta(3) &
              * (   (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                  +           delta(1)  * (1.0_dp - delta(2)) * p101   &
                  + (1.0_dp - delta(1)) *           delta(2)  * p011   &
                  +           delta(1)  *           delta(2)  * p111 ) ) &
          - Offset

      Derivative(i,1) &
        = ( (Interpolant(i) ) * alpha(1) & 
            * ( (1.0_dp - delta(3)) * ( (delta(2) - 1.0_dp) * p000   &
                                    +  ( 1.0_dp - delta(2)) * p100   &
                                    -             delta(2)  * p010   &
                                    +             delta(2)  * p110 ) &
                         + delta(3) * ( (delta(2) - 1.0_dp) * p001   &
                                    +  ( 1.0_dp - delta(2)) * p101   &
                                    -             delta(2)  * p011   &
                                    +             delta(2)  * p111 ) ) )

      Derivative(i,2) &
        = ( ( Interpolant(i) ) * alpha(2) &
              * ( (1.0_dp - delta(3) ) * ( (delta(1) - 1.0_dp) * p000   &
                                       -             delta(1)  * p100   &
                                       +  ( 1.0_dp - delta(1)) * p010   &
                                       +             delta(1)  * p110 ) &
                            + delta(3) * ( (delta(1) - 1.0_dp) * p001   &
                                       -             delta(1)  * p101   &
                                       +   (1.0_dp - delta(1)) * p011   &
                                       +             delta(1)  * p111 ) ) )

      Derivative(i,3) &
        = ( ( Interpolant(i) ) * alpha(3) &
            * ( ( (delta(1) - 1.0_dp)) * (1.0_dp - delta(2)) * p000   &
                -            delta(1)  * (1.0_dp - delta(2)) * p100   &
                - ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                -            delta(1)  *           delta(2)  * p110   &
                +  (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                +            delta(1)  *           delta(2)  * p111 ) )

    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D

  SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D &
               ( x1, x2, x3, x4, Coordinate1, Coordinate2, Coordinate3, &
                 Coordinate4, LogInterp, Offset, Table, Interpolant, &
                 Derivative, debug )     

    REAL(dp), DIMENSION(:), INTENT(in)     :: x1
    REAL(dp), DIMENSION(:), INTENT(in)     :: x2
    REAL(dp), DIMENSION(:), INTENT(in)     :: x3
    REAL(dp), DIMENSION(:), INTENT(in)     :: x4
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate3
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate4
    INTEGER, DIMENSION(4), INTENT(in)      :: LogInterp
    REAL(dp), INTENT(in)                   :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in) :: Table
    REAL(dp), DIMENSION(:), INTENT(out)    :: Interpolant
    REAL(dp), DIMENSION(:,:), INTENT(out)  :: Derivative
    LOGICAL                                :: debug

    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111,&
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
    REAL(dp), DIMENSION(4) :: alpha, delta
    INTEGER :: i, j, k, l, il1, il2, il3, il4
    
    IF ( ( SIZE(x1) .NE. SIZE(x2) ) .OR. ( SIZE(x2) .NE. SIZE(x3) ) &
         .OR. ( SIZE(x3) .NE. SIZE(x4) ) .OR. ( SIZE(x4) .NE. SIZE (x1) ) ) &
    THEN
      WRITE(*,*) &
        'ERROR: describe arrays (of interpolation point) have diff size.'
      RETURN
    END IF

    IF (debug) THEN
      WRITE(*,*) ' Now we are in 4D differentiate routine'
      WRITE(*,*) ' SIZE(x1) (number of loop) is ', SIZE( x1 )
    END IF

    DO i = 1, SIZE( x1 )

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )
      CALL locate( Coordinate4, SIZE( Coordinate4 ), x4(i), il4 )

      IF (debug) THEN
        WRITE(*,*) ' The located postion is ', il1, il2, il3, il4
      END IF
 
      p0000 = ( Table( il1  , il2  , il3  , il4   ) )
      p0001 = ( Table( il1  , il2  , il3  , il4+1 ) )
      p0010 = ( Table( il1  , il2  , il3+1, il4   ) )
      p0011 = ( Table( il1  , il2  , il3+1, il4+1 ) )
      p0100 = ( Table( il1  , il2+1, il3  , il4   ) )
      p0101 = ( Table( il1  , il2+1, il3  , il4+1 ) )
      p0110 = ( Table( il1  , il2+1, il3+1, il4   ) )
      p0111 = ( Table( il1  , il2+1, il3+1, il4+1 ) )
      p1000 = ( Table( il1+1, il2  , il3  , il4   ) )
      p1001 = ( Table( il1+1, il2  , il3  , il4+1 ) )
      p1010 = ( Table( il1+1, il2  , il3+1, il4   ) )
      p1011 = ( Table( il1+1, il2  , il3+1, il4+1 ) )
      p1100 = ( Table( il1+1, il2+1, il3  , il4   ) )
      p1101 = ( Table( il1+1, il2+1, il3  , il4+1 ) )
      p1110 = ( Table( il1+1, il2+1, il3+1, il4   ) )
      p1111 = ( Table( il1+1, il2+1, il3+1, il4+1 ) )

      IF (debug) THEN
        WRITE(*,*) ' 16 corners are loaded. Then are:'
        WRITE(*,*) p0000, p0001, &
                   p0010, p0011, p0100, p0101, p0110, p0111, p1000, &
                   p1001, p1010, p1011, p1100, p1101, p1110, p1111
      END IF

      IF ( LogInterp(1) == 1 ) THEN
        alpha(1) &
          = 1.0_dp / ( x1(i) * LOG10( Coordinate1(il1+1) / Coordinate1(il1) ) )
        delta(1) &
          = LOG10( x1(i) / Coordinate1(il1) ) &
              / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        alpha(1) &
          = ln10 / ( Coordinate1(il1+1) - Coordinate1(il1) )
        delta(1) &
          = ( x1(i) - Coordinate1(il1) ) &
              / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        alpha(2) &
          = 1.0_dp / ( x2(i) * LOG10( Coordinate2(il2+1) / Coordinate2(il2) ) )
        delta(2) &
          = LOG10( x2(i) / Coordinate2(il2) ) &
              / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        alpha(2) &
          = ln10 / ( Coordinate2(il2+1) - Coordinate2(il2) )
        delta(2) &
          = ( x2(i) - Coordinate2(il2) ) &
              / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        alpha(3) &
          = 1.0_dp / ( x3(i) * LOG10( Coordinate3(il3+1) / Coordinate3(il3) ) )
        delta(3) &
          = LOG10( x3(i) / Coordinate3(il3) ) &
              / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        alpha(3) &
          = ln10 / ( Coordinate3(il3+1) - Coordinate3(il3) )
        delta(3) &
          = ( x3(i) - Coordinate3(il3) ) &
              / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      IF ( LogInterp(4) == 1 ) THEN
        alpha(4) &
          = 1.0_dp / ( x4(i) * LOG10( Coordinate4(il4+1) / Coordinate4(il4) ) )
        delta(4) &
          = LOG10( x4(i) / Coordinate4(il4) ) &
              / LOG10( Coordinate4(il4+1) / Coordinate4(il4) )
      ELSE
        alpha(4) &
          = ln10 / ( Coordinate4(il4+1) - Coordinate4(il4) )
        delta(4) &
          = ( x4(i) - Coordinate4(il4) ) &
              / ( Coordinate4(il4+1) - Coordinate4(il4) )
      END IF
     
      IF (debug) THEN
        WRITE(*,*) '  Alpha and delta are calculated.'
        WRITE(*,*) ' delta = ', delta
        WRITE(*,*) ' alpha = ', alpha
      END IF

      Interpolant(i) &
        = 10.d0**( &
            (1.0_dp - delta(4)) &
              * (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0000   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1000   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0100   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1100   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0010   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1010   &
                  +           delta(3)  * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0110   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1110 ) &
            +         delta(4)  &
              * (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0001   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1001   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0101   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1101   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0011   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1011   &
                  +           delta(3)  * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0111   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1111 ) ) &
          - Offset

      IF (debug) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'Interpolant is calculated =', Interpolant(i)
        WRITE(*,*) 'LOG Interpolant is calculated =', LOG10( Interpolant(i) )
        WRITE(*,*) ''
      END IF

      Derivative(i,1) &  ! E
        = ( Interpolant(i) ) * alpha(1) &
          * ( (1.0_dp - delta(4)) * ( (1.0_dp - delta(3) ) * &
                                      ( ( delta(2) - 1.0_dp )  * p0000   &
                                      + ( 1.0_dp   - delta(2)) * p1000   &
                                                   - delta(2)  * p0100   &
                                                   + delta(2)  * p1100 ) &
                                    +           delta(3)   * &
                                      ( ( delta(2) - 1.0_dp )  * p0010   &
                                      - ( delta(2) - 1.0_dp )  * p1010   &
                                                   - delta(2)  * p0110   &
                                        +            delta(2)  * p1110 ))&
                      + delta(4) * ( (1.0_dp - delta(3) )  * &
                                     ( ( delta(2) - 1.0_dp )   * p0001   &
                                     + ( 1.0_dp   - delta(2))  * p1001   &
                                                  - delta(2)   * p0101   &
                                                  + delta(2)   * p1101 ) &
                                     +           delta(3)  * &
                                      ( ( delta(2) - 1.0_dp )  * p0011   &
                                      + ( 1.0_dp - delta(2) )  * p1011   &
                                                  - delta(2)   * p0111   &
                                                  + delta(2)   * p1111)) )
     
      Derivative(i,2) &  ! rho
        = ( Interpolant(i) ) * alpha(2) &
          * ( &
            (1.0_dp - delta(4)) &
              * (   (1.0_dp - delta(3)) * &
                                    ( -1.0_dp + delta(1)) * p0000   &
                  + (1.0_dp - delta(3)) * &
                                    (         - delta(1)) * p1000   &
                  + (1.0_dp - delta(3)) * &
                                    (  1.0_dp - delta(1)) * p0100   &
                  + (1.0_dp - delta(3)) * &
                                                delta(1)  * p1100   &
                  +           delta(3)  * &
                                    ( -1.0_dp + delta(1)) * p0010   &
                  +           delta(3)  * &
                                    (         - delta(1)) * p1010   &
                  +           delta(3)  * &
                                    (  1.0_dp - delta(1)) * p0110   &
                  +           delta(3)  * &
                                                delta(1)  * p1110 ) &
            +         delta(4)  &
              * (   (1.0_dp - delta(3)) * &
                                    ( -1.0_dp + delta(1)) * p0001   &
                  + (1.0_dp - delta(3)) * &
                                    (         - delta(1)) * p1001   &
                  + (1.0_dp - delta(3)) * &
                                    (  1.0_dp - delta(1)) * p0101   &
                  + (1.0_dp - delta(3)) * &
                                                delta(1)  * p1101   &
                  +           delta(3)  * &
                                    ( -1.0_dp - delta(1)) * p0011   &
                  +           delta(3)  * &
                                    (         - delta(1)) * p1011   &
                  +           delta(3)  * &
                                    (  1.0_dp - delta(1)) * p0111   &
                  +           delta(3)  * &
                                                delta(1)  * p1111 ) )

      Derivative(i,3) &  ! T
        = ( Interpolant(i) )  * alpha(3) &
          * ( (1.0_dp - delta(4)) &
              * (      (  (delta(2) - 1.0_dp) * (1.0_dp - delta(1)) * p0000   &
                  +       (delta(2) - 1.0_dp) *           delta(1)  * p1000   &
                  +                 delta(2)  * (delta(1) - 1.0_dp) * p0100   &
                  -                 delta(2)  *           delta(1)  * p1100 ) &
                  +    (  (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0010   &
                  +       (1.0_dp - delta(2)) *           delta(1)  * p1010   &
                  +                 delta(2)  * (1.0_dp - delta(1)) * p0110   &
                  +                 delta(2)  *           delta(1)  * p1110)) &
            +           delta(4)  &
              * (      (  (delta(2) - 1.0_dp) * (1.0_dp - delta(1)) * p0001   &
                  +       (delta(2) - 1.0_dp) *           delta(1)  * p1001   &
                  +                 delta(2)  * (delta(1) - 1.0_dp) * p0101   &
                  -                 delta(2)  *           delta(1)  * p1101 ) &
                  +    (  (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0011   &
                  +       (1.0_dp - delta(2)) *           delta(1)  * p1011   &
                  +                 delta(2)  * (1.0_dp - delta(1)) * p0111   &
                  +                 delta(2)  *           delta(1)  * p1111)) ) 

      Derivative(i,4) &  ! Ye
        = ( Interpolant(i) ) * alpha(4) &
          * ( - (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0000   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1000   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0100   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1100   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0010   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1010   &
                  +           delta(3)  * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0110   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1110 ) &
            +   (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0001   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1001   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0101   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1101   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0011   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1011   &
                  +           delta(3)  * &
                                   delta(2)  * (1.0_dp - delta(1)) * p0111   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1111 ) )
    
    IF (debug) THEN
      WRITE(*,*) ' '
      WRITE(*,*) 'Derivative is calculated =', Derivative(i,:)
      WRITE(*,*) ''
    END IF

    END DO

    IF (debug) THEN
      WRITE(*,*) 'End of differentiate routine'
      WRITE(*,*) ''
    END IF

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D_Custom &
       ( E, D, T, Y, Es, Ds, Ts, Ys, OS, Table, Interpolant, Derivatives )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: E,  D,  T,  Y
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Es, Ds, Ts, Ys
    REAL(dp),                     INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),       INTENT(out) :: Interpolant
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Derivatives

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D_Custom


  SUBROUTINE LogInterpolateDifferentiateAllVariables &
               ( x1, x2, x3, LogInterp, TS, DV, Interpolants, Derivatives, MaskVar )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar

    REAL(dp), DIMENSION(:,:), INTENT(out) :: Interpolants 

    REAL(dp), DIMENSION(:,:,:), INTENT(out) :: Derivatives 
    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: alpha, delta
    INTEGER :: i, j, il1, il2, il3
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize

    epsilon = 1.d-200

    Masksize = SIZE( x2 )
    ALLOCATE( work_mask( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, SIZE(x2)

      IF ( .not.work_mask(i) ) CYCLE

      ASSOCIATE( Coordinate1 => TS % States(1) % Values, &
                 Coordinate2 => TS % States(2) % Values, &
                 Coordinate3 => TS % States(3) % Values )

      CALL locate( Coordinate1, SIZE(Coordinate1), x1(i), il1 )
      CALL locate( Coordinate2, SIZE(Coordinate2), x2(i), il2 )
      CALL locate( Coordinate3, SIZE(Coordinate3), x3(i), il3 )

      IF ( LogInterp(1) == 1 ) THEN
      alpha(1) = ( 1.0d0 ) / ( x1(i) * LOG10( Coordinate1(il1+1) / Coordinate1(il1) ) )
      delta(1) = LOG10( x1(i) / Coordinate1(il1) ) / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
      alpha(1) = ( ln10 ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      delta(1) = ( x1(i) - Coordinate1(il1) ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
      alpha(2) = ( 1.0d0 ) / ( x2(i) * LOG10( Coordinate2(il2+1) / Coordinate2(il2) ) )
      delta(2) = LOG10( x2(i) / Coordinate2(il2) ) / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
      alpha(2) = ( ln10 ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      delta(2) = ( x2(i) - Coordinate2(il2) ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
      alpha(3) = ( 1.0d0 ) / ( x3(i) * LOG10( Coordinate3(il3+1) / Coordinate3(il3) ) )
      delta(3) = LOG10( x3(i) / Coordinate3(il3) ) / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
      alpha(3) = ( ln10 ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      delta(3) = ( x3(i) - Coordinate3(il3) ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      END ASSOCIATE

      DO j = 1, DV % nVariables

        ASSOCIATE( Table => DV % Variables(j) % Values(:,:,:), &
                   Offset => DV % Offsets(j) )

        p000 = ( Table( il1  , il2  , il3   ) )
        p100 = ( Table( il1+1, il2  , il3   ) )
        p010 = ( Table( il1  , il2+1, il3   ) )
        p110 = ( Table( il1+1, il2+1, il3   ) )
        p001 = ( Table( il1  , il2  , il3+1 ) )
        p101 = ( Table( il1+1, il2  , il3+1 ) )
        p011 = ( Table( il1  , il2+1, il3+1 ) )
        p111 = ( Table( il1+1, il2+1, il3+1 ) )

        Interpolants(i,j) &
          = 10.d0**( &
                (1.0_dp - delta(3)) * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &
                                     +            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                     + ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                     +            delta(1)  *           delta(2)  * p110 ) &
                        + delta(3)  * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                     +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                     +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                     +            delta(1)  *           delta(2)  * p111 ) &

                   ) - Offset

        Derivatives(i,1,j) &
          = ( (Interpolants(i,j) ) * alpha(1) &
              * ( (1.0_dp - delta(3)) * ( (delta(2) - 1.0_dp) * p000   &
                                      +  ( 1.0_dp - delta(2)) * p100   &
                                      -             delta(2)  * p010   &
                                      +             delta(2)  * p110 ) &
                           + delta(3) * ( (delta(2) - 1.0_dp) * p001   &
                                      +  ( 1.0_dp - delta(2)) * p101   &
                                      -             delta(2)  * p011   &
                                      +             delta(2)  * p111 ) ) )
  
        Derivatives(i,2,j) &
          = ( ( Interpolants(i,j) ) * alpha(2) &
              * ( (1.0_dp - delta(3) ) * ( (delta(1) - 1.0_dp) * p000   &
                                       -             delta(1)  * p100   & 
                                       +  ( 1.0_dp - delta(1)) * p010   & 
                                       +             delta(1)  * p110 ) & 
                            + delta(3) * ( (delta(1) - 1.0_dp) * p001   & 
                                       -             delta(1)  * p101   & 
                                       +   (1.0_dp - delta(1)) * p011   & 
                                       +             delta(1)  * p111 ) ) )

        Derivatives(i,3,j) &
          = ( ( Interpolants(i,j) ) * alpha(3) &
                                     * ( ( (delta(1) - 1.0_dp)) * (1.0_dp - delta(2)) * p000   &
                                         -            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                         - ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                         -            delta(1)  *           delta(2)  * p110   &
                                         +  (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                         +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                         +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                         +            delta(1)  *           delta(2)  * p111 ) )
  
  
        END ASSOCIATE

      END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateAllVariables 

  SUBROUTINE GetGamma1( x1, x2, x3, Coordinate1, Coordinate2, &
                Coordinate3, LogInterp, TS, DV, Gamma1 ) 

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV

    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: Gamma1

    INTEGER :: i
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Interpolant 
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative 

    ALLOCATE( Interpolant( SIZE(x1) ) )
    ALLOCATE( Gamma1( SIZE(x1) ) )
    ALLOCATE( Derivative( SIZE(x1), 3 ) )

      CALL LogInterpolateDifferentiateSingleVariable( x1, x2, x3,                 &
                                    TS % States(1) % Values(:),        &
                                    TS % States(2) % Values(:),        &
                                    TS % States(3) % Values(:),        &
                                    LogInterp, DV % Offsets(1),        &
                                    DV % Variables(1) % Values(:,:,:), &
                                    Interpolant(:), Derivative(:,:) )
       
      DO i = 1, SIZE(x1) 
        Gamma1(i) =  ( x1(i)/Interpolant(i) ) * Derivative(i, 1 ) 
      END DO

  END SUBROUTINE GetGamma1 
  
  SUBROUTINE MonotonicityCheck ( Table, Nrho, NT, NYe, Axis, Repaired )

    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    INTEGER, DIMENSION(:,:,:), INTENT(in) :: Repaired
    INTEGER, INTENT(in) :: Nrho 
    INTEGER, INTENT(in) :: NT
    INTEGER, INTENT(in) :: NYe
    INTEGER, INTENT(in) :: Axis

    INTEGER :: i, j, k, count

    97 FORMAT ("Table not monotonic in rho at (Nrho, NT, NYe) = ", 3(1x,i4) )
    98 FORMAT ("Table not monotonic in T at (Nrho, NT, NYe) = ", 3(1x,i4) )
    99 FORMAT ("Table not monotonic in Ye at (Nrho, NT, NYe) = ", 3(1x,i4) )
 
    count = 0
    
    SELECT CASE ( Axis ) 
     
    CASE( 1 )
      DO k = 1, NYe
        DO j = 1, NT  
          DO i = 2, Nrho - 1

            IF ( ( ( Table(i+1, j, k) - Table(i, j, k) ) * &
                 ( Table(i, j, k) - Table(i-1, j, k) ) ) < 0. ) THEN
              WRITE (*,97) i, j, k
              WRITE (*,*) "Repaired =", Repaired(i,j,k), Repaired(i+1,j,k), Repaired(i-1,j,k), &
                Repaired(i,j+1,k), Repaired(i,j-1,k), Repaired(i,j,k+1), Repaired(i,j,k-1)
              count = count + 1
            END IF
          END DO
        END DO
      END DO

    CASE( 2 )
      DO k = 1, NYe
        DO j = 2, NT - 1 
          DO i = 1, Nrho

            IF ( ( ( Table(i, j+1, k) - Table(i, j, k) ) * &
                 ( Table(i, j, k) - Table(i, j-1, k) ) ) < 0.) THEN 
              WRITE (*,98) i, j, k
              WRITE (*,*) "Repaired =", Repaired(i,j,k), Repaired(i+1,j,k), Repaired(i-1,j,k), &
                Repaired(i,j+1,k), Repaired(i,j-1,k), Repaired(i,j,k+1), Repaired(i,j,k-1)
              count = count + 1
            END IF
          END DO
        END DO
      END DO

   CASE( 3 )
      DO k = 2, NYe - 1
        DO j = 1, NT
          DO i = 1, Nrho

            IF ( ( ( Table(i, j, k+1) - Table(i, j, k) ) * &
                 ( Table(i, j, k) - Table(i, j, k-1) ) ) < 0. ) &
              WRITE (*, 99) i, j, k
              WRITE (*,*) "Repaired =", Repaired(i,j,k), Repaired(i+1,j,k), Repaired(i-1,j,k), &
                Repaired(i,j+1,k), Repaired(i,j-1,k), Repaired(i,j,k+1), Repaired(i,j,k-1)
              count = count + 1
          END DO
        END DO
      END DO

    CASE DEFAULT
      WRITE (*,*) "Invalid Axis", Axis
      STOP

    END SELECT
    WRITE (*,*) count, " Non-monotonic out of " , NYe*NT*Nrho

  END SUBROUTINE MonotonicityCheck 


  SUBROUTINE ComputeTempFromIntEnergy &
               ( rho, e_int, ye, density_table, temp_table, ye_table, &
                 LogInterp, energy_table, Offset, Temperature )

    REAL(dp), INTENT(in)                    :: rho
    REAL(dp), INTENT(in)                    :: e_int
    REAL(dp), INTENT(in)                    :: ye
    REAL(dp), DIMENSION(:), INTENT(in)      :: density_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: ye_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: temp_table
    INTEGER, DIMENSION(3), INTENT(in)       :: LogInterp 
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: energy_table
    REAL(dp), INTENT(in) :: Offset

    REAL(dp), DIMENSION(1)                  :: eibuff ! internal energy buffer
    INTEGER                                 :: nPoints
    INTEGER                                 :: i, j
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rhobuff 
    REAL(dp), DIMENSION(:), ALLOCATABLE :: energy_array
    REAL(dp), DIMENSION(:), ALLOCATABLE :: yebuff

    REAL(dp), DIMENSION(1), INTENT(out)     :: Temperature

  nPoints = SIZE(temp_table)

  ALLOCATE( energy_array( nPoints ), rhobuff( nPoints ), yebuff( nPoints) )

  rhobuff(1:nPoints) = rho
  eibuff(1) = e_int
  yebuff(1:nPoints) = ye

  CALL LogInterpolateSingleVariable                                     &
         ( rhobuff, temp_table, yebuff,                                 &
           density_table,                                               &
           temp_table,                                                  &
           ye_table,                                                    &
           LogInterp,                                                   &
           Offset,                                                      &
           energy_table(:,:,:), energy_array )

    DO j = 1, SIZE( eibuff )

      CALL locate( energy_array, SIZE( energy_array ), eibuff(j), i )
      IF ( i == SIZE(energy_array) ) THEN
        STOP 'energy too high'
      END IF

      IF ( i == 0 ) THEN
        Temperature(j) = 0.d0
        CYCLE
      END IF
      Temperature(j) = 10.d0**( &
             LOG10( temp_table(i) ) + LOG10( temp_table(i+1) / temp_table(i) ) &
                        * LOG10( ( eibuff(j) / energy_array(i) ) )             &
                        / LOG10( energy_array(i+1) / energy_array(i) ) )
    END DO

    DEALLOCATE( energy_array, rhobuff, yebuff )

  END SUBROUTINE ComputeTempFromIntEnergy


  SUBROUTINE ComputeTempFromEntropy &
               ( rho, s, ye, density_table, temp_table, ye_table, &
                 LogInterp, entropy_table, Offset, Temperature )

    REAL(dp), INTENT(in)                    :: rho
    REAL(dp), INTENT(in)                    :: s
    REAL(dp), INTENT(in)                    :: ye
    REAL(dp), DIMENSION(:), INTENT(in)      :: density_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: ye_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: temp_table
    INTEGER, DIMENSION(3), INTENT(in)       :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: entropy_table
    REAL(dp), INTENT(in) :: Offset

    REAL(dp), DIMENSION(1)                  :: sbuff ! entropy buffer
    INTEGER                                 :: nPoints
    INTEGER                                 :: i, j
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rhobuff
    REAL(dp), DIMENSION(:), ALLOCATABLE :: entropy_array
    REAL(dp), DIMENSION(:), ALLOCATABLE :: yebuff

    REAL(dp), DIMENSION(1), INTENT(out)     :: Temperature

  nPoints = SIZE(temp_table)

  ALLOCATE( entropy_array( nPoints ), rhobuff( nPoints ), yebuff( nPoints) )

  rhobuff(1:nPoints) = rho
  sbuff(1) = s
  yebuff(1:nPoints) = ye

  CALL LogInterpolateSingleVariable                                     &
         ( rhobuff, temp_table, yebuff,                                 &
           density_table,                                               &
           temp_table,                                                  &
           ye_table,                                                    &
           LogInterp,                                                   &
           Offset,                                                      &
           entropy_table(:,:,:), entropy_array )

    DO j = 1, SIZE( sbuff )

      CALL locate( entropy_array, SIZE( entropy_array ), sbuff(j), i )
      IF ( i == SIZE(entropy_array) ) THEN
        STOP
      END IF

      IF ( i == 0 ) THEN
        Temperature(j) = 0.d0
        CYCLE
      END IF
      Temperature(j) = 10.d0**( &
             LOG10( temp_table(i) ) + LOG10( temp_table(i+1) / temp_table(i) ) &
                        * LOG10( ( sbuff(j) / entropy_array(i) ) )             &
                        / LOG10( entropy_array(i+1) / entropy_array(i) ) )
    END DO

  DEALLOCATE( entropy_array, rhobuff, yebuff )

  END SUBROUTINE ComputeTempFromEntropy

  SUBROUTINE ComputeTempFromPressure &
               ( rho, p, ye, density_table, temp_table, ye_table, &
                 LogInterp, pressure_table, Offset, Temperature )

    REAL(dp), INTENT(in)                    :: rho
    REAL(dp), INTENT(in)                    :: p
    REAL(dp), INTENT(in)                    :: ye
    REAL(dp), DIMENSION(:), INTENT(in)      :: density_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: ye_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: temp_table
    INTEGER, DIMENSION(3), INTENT(in)       :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: pressure_table
    REAL(dp), INTENT(in) :: Offset

    REAL(dp), DIMENSION(1)                  :: pbuff ! pressure buffer
    INTEGER                                 :: nPoints
    INTEGER                                 :: i, j
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rhobuff
    REAL(dp), DIMENSION(:), ALLOCATABLE :: pressure_array
    REAL(dp), DIMENSION(:), ALLOCATABLE :: yebuff

    REAL(dp), DIMENSION(1), INTENT(out)     :: Temperature

  nPoints = SIZE(temp_table)

  ALLOCATE( pressure_array( nPoints ), rhobuff( nPoints ), yebuff( nPoints) )

  rhobuff(1:nPoints) = rho
  pbuff(1) = p
  yebuff(1:nPoints) = ye

  CALL LogInterpolateSingleVariable                                     &
         ( rhobuff, temp_table, yebuff,                                 &
           density_table,                                               &
           temp_table,                                                  &
           ye_table,                                                    &
           LogInterp,                                                   &
           Offset,                                                      &
           pressure_table(:,:,:), pressure_array )

    DO j = 1, SIZE( pbuff )

      CALL locate( pressure_array, SIZE( pressure_array ), pbuff(j), i )
      IF ( i == SIZE(pressure_array) ) THEN
        STOP
      END IF

      IF ( i == 0 ) THEN
        Temperature(j) = 0.d0
        CYCLE
      END IF
      Temperature(j) = 10.d0**( &
             LOG10( temp_table(i) ) + LOG10( temp_table(i+1) / temp_table(i) ) &
                        * LOG10( ( pbuff(j) / pressure_array(i) ) )             &
                        / LOG10( pressure_array(i+1) / pressure_array(i) ) )
    END DO

  DEALLOCATE( pressure_array, rhobuff, yebuff )

  END SUBROUTINE ComputeTempFromPressure


  SUBROUTINE EOSTableQuery &
               ( rho, T, Ye, LogInterp, TS, DV, Interpolants )

    INTEGER                                    :: i, j
    REAL(dp), DIMENSION(:), INTENT(in)         :: rho
    REAL(dp), DIMENSION(:), INTENT(in)         :: T
    REAL(dp), DIMENSION(:), INTENT(in)         :: Ye
    INTEGER, DIMENSION(3), INTENT(in)          :: LogInterp
    TYPE(ThermoStateType), INTENT(in)          :: TS
    TYPE(DependentVariablesType), INTENT(in)   :: DV
    REAL(dp), DIMENSION(:,:), INTENT(out)      :: Interpolants

!    CALL InitializeHDF( )

!    CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

    CALL LogInterpolateAllVariables( rho, T, Ye, LogInterp, &
                                     TS, DV, Interpolants )
    DO i = 1, SIZE(rho)
      WRITE(*,*) 'Rho=', rho(i), 'T=', T(i), 'Ye=', Ye(i)
      DO j = 1, DV % nVariables
        WRITE(*,*) DV % Names(j), Interpolants(i,j)
      END DO
    END DO

  END SUBROUTINE EOSTableQuery

  !SUBROUTINE InterpolateSingleVariableEOS( rho, T, Ye, maskvar, interpolants )
  !           ( x1, x2, x3, Coordinate1, Coordinate2, Coordinate3, &
  !               LogInterp, Offset, Table, Interpolant )
!
!    REAL(dp), DIMENSION(:), INTENT(in) :: x1
!    REAL(dp), DIMENSION(:), INTENT(in) :: x2
!    REAL(dp), DIMENSION(:), INTENT(in) :: x3
!    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
!    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
!    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
!    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp
!    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
!    REAL(dp), INTENT(in) :: Offset
!    REAL(dp), DIMENSION(:), INTENT(out) :: Interpolant
!
!    REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(in) :: rho, T, Ye
!    LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(in) :: maskvar
!    LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(in) :: in_mask
!
!    maskvar(:) = ( nse_c(:) == 1 )
!
!    IF ( PRESENT(in_mask) ) THEN
!      work_mask = in_mask
!    ELSE
!      work_mask = true
!
!    DO i = 1, SIZE(rho)

!      IF ( .not. maskvar(i) ) CYCLE
!
!      interpolants(:,i) = ...

!    END DO

!  END SUBROUTINE InterpolateEOS

END MODULE wlInterpolationModule
