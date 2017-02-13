PROGRAM wlInterpolation

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateDifferentiateSingleVariable, &
    LogInterpolateSingleVariable
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
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_rho, Inte_T, Inte_Ye, database
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,h
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, datasize
  REAL(dp)                            :: Offset_EOS_cmpe
  REAL(dp)                            :: Offset_EOS_cmpp
  REAL(dp)                            :: Offset_EOS_cmpn
  REAL(dp)                            :: Offset_EOS_cmpnu
  REAL(dp)                            :: Offset_Em

!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: Interpolant
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: bufferO1, bufferO2, bufferO3, bufferO4
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative

!----------------------------------------
!   interpolated energy 
!----------------------------------------
 
  Format1 = "(3A12)"
  Format2 = "(3ES12.3)"
  Format3 = "(7A12)"
  Format4 = "(7ES12.3)"

  OPEN(1, FILE = "Inputfile_stand.d", FORM = "formatted", ACTION = 'read')
  datasize = 2  ! Vary from different input files

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

  ALLOCATE( database( datasize * 3) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( bufferO1( datasize ) )
  ALLOCATE( bufferO2( datasize ) )
  ALLOCATE( bufferO3( datasize ) )
  ALLOCATE( bufferO4( datasize ) )
  ALLOCATE( Interpolant( Inte_nPointE ) )
  ALLOCATE( Derivative( Inte_nPointE, 4 ) )

  READ( 1, Format1 ) a,b,c
  READ( 1, Format2 ) database

  CLOSE( 1, STATUS = 'keep')  

  DO i = 1, datasize  
    Inte_rho(i) = database(i*3-2)
    Inte_T(i) = database(i*3-1)
    Inte_Ye(i) = database(i*3)
  END DO 

  LogInterp(2:4) = (/1, 1, 0/)     ! rho and T is LogGrid, Ye is linear
 
!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  Offset_Em = OpacityTable % thermEmAb % Offset
  Offset_EOS_cmpe = OpacityTable % EOSTable % DV % Offsets(4)
  Offset_EOS_cmpp = OpacityTable % EOSTable % DV % Offsets(5)
  Offset_EOS_cmpn = OpacityTable % EOSTable % DV % Offsets(6)

!--------------------------------------
!   do interpolation
!--------------------------------------
  e = ('    energy  ')
  f = (' thermEmAb  ')
  g = ('Em_deriv T ')
  h = ('Em_deriv Ye')  

  OPEN( 10, FILE = "Output.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) a,b,c,e,f,g,h

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
             LogInterp, Offset_Em, Table, Interpolant, Derivative, .FALSE. )
  
    DO ii = 1, Inte_nPointE
      WRITE(10, Format4) buffer1(ii), buffer2(ii), buffer3(ii), &
                         Inte_E % Values(ii), Interpolant(ii), Derivative(ii,3),&
                         Derivative(ii,4)
    END DO ! ii

  END DO ! i

  END ASSOCIATE ! Table

  CLOSE( 10, STATUS = 'keep')  

  CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,0/), Offset_EOS_cmpe, &
             OpacityTable % EOSTable % DV % Variables(4) % Values, &
             bufferO1  )

  WRITE(*,*) 'electron chemical potential = ', bufferO1

  CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,0/), Offset_EOS_cmpp, &
             OpacityTable % EOSTable % DV % Variables(5) % Values, &
             bufferO2  )

  WRITE(*,*) 'proton chemical potential = ', bufferO2

  CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,0/), Offset_EOS_cmpn, &
             OpacityTable % EOSTable % DV % Variables(6) % Values, &
             bufferO3  )

  WRITE(*,*) 'neutron chemical potential = ', bufferO3

  bufferO4 = bufferO1 + bufferO2 - bufferO3 - 1.29333d+00
            ! chem_e + chem_p - chem_n - dmnp

  WRITE(*,*) 'neutrino chemical potential = ', bufferO4

  WRITE(*,*) 'File Output.d was written/rewrtited.'

END PROGRAM wlInterpolation
