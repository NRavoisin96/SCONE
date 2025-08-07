module analysisDistribution

  use numPrecision,                  only : shortInt
  use genericProcedures,             only : fatalError
  use, intrinsic :: iso_fortran_env, only : int32, int64

  implicit none
  private

  public :: getIsInitialisedANALYSIS, initDistributionANALYSIS, updateDistributionTypeANALYSIS, &
            updateDistributionValenceANALYSIS, exportDistributionANALYSIS, updateCountANALYSIS

  integer(int32), dimension(:,:), allocatable    :: distributionMatANALYSIS
  integer(int32), dimension(:),   allocatable    :: distributionANALYSIS
  integer(shortInt)                              :: maxValenceANALYSIS, layerNoANALYSIS
  integer(int64)                                 :: countANALYSIS, newCountANALYSIS

contains

  !!
  !!
  !!
  function getIsInitialisedANALYSIS() result(isInitialised)
    logical                     :: isInitialised

    isInitialised = .TRUE.
    if (.NOT. allocated(distributionANALYSIS)) isInitialised = .FALSE.

  end function getIsInitialisedANALYSIS

  !!
  !!
  !!
  subroutine initDistributionANALYSIS(max)
    integer(shortInt), intent(in)            :: max

    allocate(distributionANALYSIS(max+7))
    allocate(distributionMatANALYSIS(max+7, 50))
    distributionANALYSIS = 0
    maxValenceANALYSIS = max
    countANALYSIS     = 0_int64
    layerNoANALYSIS   = 1
    newCountANALYSIS  = 0_int64

  end subroutine initDistributionANALYSIS

  !!
  !!
  !!
  subroutine updateDistributionTypeANALYSIS(type)
    integer(shortInt), intent(in)                   :: type

    distributionANALYSIS(maxValenceANALYSIS + type) = distributionANALYSIS(maxValenceANALYSIS + type) + 1

  end subroutine updateDistributionTypeANALYSIS

  !!
  !!
  !!
  subroutine updateDistributionValenceANALYSIS(valence)
    integer(shortInt), intent(in)                   :: valence

    distributionANALYSIS(valence) = distributionANALYSIS(valence) + 1

  end subroutine updateDistributionValenceANALYSIS

  !!
  !!
  !!
  subroutine exportDistributionANALYSIS(filename)
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none

    character(len=*), intent(in)  :: filename
    integer(int64), parameter     :: UNIT=99_int64
    integer                       :: i, j, ios
    integer                       :: nrows, ncols

    ! Determine matrix dimensions
    nrows = size(distributionMatANALYSIS, 1)
    ncols = size(distributionMatANALYSIS, 2)

    ! Open the file for writing (replace if it exists)
    open(unit=UNIT, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Could not open ', trim(filename), ' (iostat=', ios, ')'
      return
    end if

    ! Write each column as one CSV line
    do j = 1, ncols
      do i = 1, nrows
        if (i < nrows) then
          ! Write value plus comma, no newline
          write(UNIT,'(I0,",")', advance='no') distributionMatANALYSIS(i,j)
        else
          ! Last value in row: no trailing comma
          write(UNIT,'(I0)',    advance='no') distributionMatANALYSIS(i,j)
        end if
      end do
      ! Now write the newline
      write(UNIT,*)
    end do

    close(UNIT)
    write(*,*) 'Wrote ', nrows, '×', ncols, ' matrix to ', trim(filename)

    call fatalError("analysisDistribution.f90", "CSV file for distribution has been exported.")
  
  end subroutine exportDistributionANALYSIS

  !!
  !!
  !!
  subroutine updateCountANALYSIS()

    countANALYSIS     = countANALYSIS    + 1_int64
    newCountANALYSIS  = newCountANALYSIS + 1_int64

    if (newCountANALYSIS == 2500000) then

      distributionMatANALYSIS(:,layerNoANALYSIS) = distributionANALYSIS
      layerNoANALYSIS = layerNoANALYSIS + 1
      newCountANALYSIS = 0_int64

      if (countANALYSIS == 2500000*50) then
        call exportDistributionANALYSIS('distributionMatANALYSIS.csv')
      end if

    end if

  end subroutine updateCountANALYSIS

end module analysisDistribution
