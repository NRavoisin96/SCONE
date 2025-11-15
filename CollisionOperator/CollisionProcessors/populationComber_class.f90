module populationComber_class

  use dictionary_class,            only : dictionary
  use endfConstants,               only : N_N_SPLIT
  use errors_mod,                  only : fatalError
  use geometryReg_mod,             only : fieldIdx, fieldPtr
  use numPrecision
  use particleDungeon_class,       only : particleDungeon
  use physicalParticle_inter,      only : physicalParticle
  use physicalParticleState_class, only : castPhysicalParticleStatePtr, physicalParticleState
  use tallyAdmin_class,            only : tallyAdmin
  use transportObjectState_class,  only : transportObjectState
  use universalVariables,          only : nameWW
  use weightWindowsField_class,    only : weightWindowsField, weightWindowsField_TptrCast

  implicit none
  private

  ! Parameters.
  integer(shortInt), parameter :: DEFAULT_N_MAX_SPLITS = 1000
  real(defReal), parameter     :: DEFAULT_AVERAGE_WEIGHT = 0.5_defReal, DEFAULT_MAX_WEIGHT = 1.25_defReal, &
                                  DEFAULT_MIN_WEIGHT = 0.25_defReal

  !!
  !!
  !!
  type, public :: populationComber
    private
    integer(shortInt)                 :: nMaxSplits = 0
    logical(defBool)                  :: usesRussianRoulette = .false., usesSplitting = .false., usesWeightWindows = .false.
    real(defReal)                     :: averageWeight = ZERO, maxWeight = ZERO, minWeight = ZERO
    type(weightWindowsField), pointer :: weightWindowsMap
  contains
    procedure          :: cutoffs
    procedure          :: getUsesRussianRoulette
    procedure          :: getUsesWeightWindows
    procedure          :: init
    procedure          :: kill
    procedure, private :: russianRoulette
    procedure, private :: split
  end type populationComber

contains
  !!
  !!
  !!
  subroutine cutoffs(self, p, dungeon, tally)
    class(populationComber), intent(in)    :: self
    class(physicalParticle), intent(inout) :: p
    type(particleDungeon), intent(inout)   :: dungeon
    type(tallyAdmin), intent(inout)        :: tally
    real(defReal)                          :: maxWeight, weight
    real(defReal), dimension(3)            :: weightWindowsValues ! 1 = minWeight, 2 = maxWeight, 3 = averageWeight

    ! Retrieve current particle weight.
    weight = p % getWeight()

    if (self % usesWeightWindows) then
      ! Weight windows treatment.
      weightWindowsValues = self % weightWindowsMap % at(p)

      ! If a particle is outside the WW map and all the weight limits
      ! are zero nothing happens. NOTE: this holds for positive weights only
      maxWeight = weightWindowsValues(2)
      if (maxWeight < weight .and. maxWeight /= ZERO .and. p % getSplitsNumber() < self % nMaxSplits) then
        call self % split(maxWeight, p, dungeon, tally)
        return

      end if

      if (weight < weightWindowsValues(1)) then
        call self % russianRoulette(weightWindowsValues(3), p)

      end if

    elseif (self % usesSplitting .and. self % maxWeight < weight) then
      ! Splitting with fixed threshold.
      call self % split(self % maxWeight, p, dungeon, tally)

    elseif (self % usesRussianRoulette .and. weight < self % minWeight) then
      ! Russian roulette with fixed threshold and survival weight.
      call self % russianRoulette(self % averageWeight, p)

    end if

  end subroutine cutoffs

  !!
  !!
  !!
  elemental function getUsesRussianRoulette(self) result(usesRussianRoulette)
    class(populationComber), intent(in) :: self
    logical(defBool)                    :: usesRussianRoulette

    usesRussianRoulette = self % usesRussianRoulette

  end function getUsesRussianRoulette

  !!
  !!
  !!
  elemental function getUsesWeightWindows(self) result(usesWeightWindows)
    class(populationComber), intent(in) :: self
    logical(defBool)                    :: usesWeightWindows

    usesWeightWindows = self % usesWeightWindows

  end function getUsesWeightWindows

  !!
  !!
  !!
  subroutine init(self, dict)
    class(populationComber), intent(inout) :: self
    class(dictionary), intent(in)          :: dict
    character(*), parameter                :: HERE = 'init (populationComber_class.f90)'

    ! Obtain settings for variance reduction.
    call dict % getOrDefault(self % nMaxSplits, 'maxSplit', DEFAULT_N_MAX_SPLITS)
    
    call dict % getOrDefault(self % usesRussianRoulette, 'roulette', .false.)
    if (self % usesRussianRoulette) then
      call dict % getOrDefault(self % averageWeight, 'avWgt', DEFAULT_AVERAGE_WEIGHT)
      call dict % getOrDefault(self % minWeight, 'minWgt', DEFAULT_MIN_WEIGHT)

    end if
    
    call dict % getOrDefault(self % usesSplitting, 'split', .false.)
    if (self % usesSplitting) then
      call dict % getOrDefault(self % maxWeight, 'maxWgt', DEFAULT_MAX_WEIGHT)
      if (self % maxWeight < TWO * self % minWeight) &
      call fatalError(HERE, 'Upper weight bound must greater than twice the lower weight bound.')

    end if

    ! Sets up the weight windows field if requested.
    call dict % getOrDefault(self % usesWeightWindows, 'weightWindows', .false.)
    if (self % usesWeightWindows) &
    self % weightWindowsMap => weightWindowsField_TptrCast(fieldPtr(fieldIdx(nameWW)))

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(populationComber), intent(inout) :: self

    ! Local.
    self % nMaxSplits = 0
    self % usesRussianRoulette = .false.
    self % usesSplitting = .false.
    self % usesWeightWindows = .false.
    self % averageWeight = ZERO
    self % maxWeight = ZERO
    self % minWeight = ZERO
    call self % weightWindowsMap % kill()

  end subroutine kill

  !!
  !!
  !!
  subroutine russianRoulette(self, averageWeight, p)
    class(populationComber), intent(in)    :: self
    real(defReal), intent(in)              :: averageWeight
    class(physicalParticle), intent(inout) :: p
    real(defReal)                          :: randomNumber

    call p % generateRandomNumber(randomNumber)
    if (randomNumber < ONE - p % getWeight() / averageWeight) then
      call p % setIsDead(.true.)

    else
      call p % setWeight(averageWeight)

    end if

  end subroutine russianRoulette

  !!
  !!
  !!
  subroutine split(self, maxWeight, p, dungeon, tally)
    class(populationComber), intent(in)      :: self
    real(defReal), intent(in)                :: maxWeight
    class(physicalParticle), intent(inout)   :: p
    type(particleDungeon), intent(inout)     :: dungeon
    type(tallyAdmin), intent(inout)          :: tally
    class(physicalParticleState), pointer    :: physicalParticleStatePtr
    class(transportObjectState), allocatable :: stateCopy
    integer(shortInt)                        :: i, mult, nSplits
    real(defReal)                            :: newWeight, weight

    ! Compute multiplier. This value must be at least two. Cap it to maximum number of splits.
    nSplits = p % getSplitsNumber()
    weight = p % getWeight()
    mult = min(ceiling(weight / maxWeight), self % nMaxSplits - nSplits + 1)

    ! Compute new weight, copy the particle's current state and update its weight.
    newWeight = weight / mult
    call p % copyCurrentState(stateCopy)
    physicalParticleStatePtr => castPhysicalParticleStatePtr(stateCopy)
    call physicalParticleStatePtr % setWeight(newWeight)

    ! Add new state to the dungeon.
    do i = 1, mult - 1
      call dungeon % detain(physicalParticleStatePtr)
      call tally % reportSpawn(N_N_SPLIT, p, physicalParticleStatePtr)

    end do

    ! Update the number of splits and weight for the particle.
    call p % setSplitsNumber(nSplits + mult)
    call p % setWeight(newWeight)

  end subroutine split

end module populationComber_class