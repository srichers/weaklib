PROGRAM wlOpacityInterpolationTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateDifferentiateSingleVariable
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid

  IMPLICIT NONE

!--------- parameters for creating energy grid 
  INTEGER, PARAMETER     :: Inte_nPointE = 30
  REAL(dp)               :: Inte_Emin = 2.0d00
  REAL(dp)               :: Inte_Emax = 2.0d02
  TYPE(GridType)   :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r, Inte_rho, Inte_T, Inte_Ye, e_int,&
                                         database
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,h
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, datasize
  REAL(dp)                            :: Offset

!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: Interpolant
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative

!----------------------------------------
!   interpolated energy 
!----------------------------------------
 
  Format1 = "(5A12)"
  Format2 = "(5ES12.3)"
  Format3 = "(8A12)"
  Format4 = "(8ES12.3)"

  OPEN(1, FILE = "Output0ms.d", FORM = "formatted", ACTION = 'read')
  datasize = 292

! OPEN(1, FILE = "Output100ms.d", FORM = "formatted", ACTION = 'read')
! datasize = 217
  
  CALL AllocateGrid( Inte_E, Inte_nPointE )

  Inte_E % Unit = 'MeV                  '
  Inte_E % Name = 'Intepolated Energy   '
  Inte_E % MinValue = Inte_Emin
  Inte_E % MaxValue = Inte_Emax
  Inte_E % LogInterp = 1
  Inte_E % nPoints = Inte_nPointE
  LogInterp(1) = 1                 ! EnergyGrid is LogGrid

  CALL MakeLogGrid &
          ( Inte_E % MinValue, Inte_E % MaxValue, &
            Inte_E % nPoints, Inte_E % Values )

  CALL DescribeGrid( Inte_E ) 

!---------------------------------------
!    interpolated rho, T, Ye
!---------------------------------------

  ALLOCATE( database( datasize * 5) )
  ALLOCATE( r( datasize ) )
  ALLOCATE( e_int( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( Interpolant( Inte_nPointE ) )
  ALLOCATE( Derivative( Inte_nPointE, 4 ) )

  READ( 1, Format1 ) a,b,c,d,e
  READ( 1, Format2 ) database

  CLOSE( 1, STATUS = 'keep')  

  DO i = 1, datasize  
    r(i) = database(i*5-4)
    Inte_rho(i) = database(i*5-3)
    Inte_T(i) = database(i*5-2)
    Inte_Ye(i) = database(i*5-1)
  END DO 

  LogInterp(2:4) = (/1, 1, 0/)     ! rho and T is LogGrid, Ye is linear
 
!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  Offset = OpacityTable % thermEmAb % Offset
!--------------------------------------
!   do interpolation
!--------------------------------------
  e = ('    energy  ')
  f = ('    Inter  ')
  g = ('    deriv T ')
  h = ('   deriv Ye')  

  OPEN( 10, FILE = "IntOutput0ms.d", FORM = "formatted", ACTION = 'write')
!  OPEN( 10, FILE = "IntOutput100ms.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) a,b,c,d,e,f,g,h

  ASSOCIATE( Table  => OpacityTable % thermEmAb % Absorptivity(1) % Values,&
             Energy => Inte_E % Values )

  DO i = 1, datasize

    buffer1(:) = Inte_rho(i)
    buffer2(:) = Inte_T(i)
    buffer3(:) = Inte_Ye(i)

    CALL LogInterpolateDifferentiateSingleVariable & 
           ( Energy, buffer1, buffer2, buffer3, & 
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset, Table, Interpolant, Derivative, .FALSE. )
  
    DO ii = 1, Inte_nPointE
      WRITE(10, Format4) r(i), buffer1(ii), buffer2(ii), buffer3(ii), &
                         Inte_E % Values(ii), Interpolant(ii), Derivative(ii,3),&
                         Derivative(ii,4)
    END DO ! ii

  END DO ! i

  END ASSOCIATE ! Table

  CLOSE( 10, STATUS = 'keep')  

  WRITE(*,*) 'File IntOutput*ms.d was written/rewrtited.'

END PROGRAM wlOpacityInterpolationTest
