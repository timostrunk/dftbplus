!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines implementing on Effective Core Potential (ECP) environment ecpenv.
module ecpenv
  use accuracy
  use assert
  use commontypes, only: TOrbitals
  use typegeometry, only: TGeometry
  implicit none
  private

  public :: TECPEnvInp, TECPEnv, ECPEnv_init


  !> Input for the ECPEnv module
  type TECPEnvInp

    !> Geometry of the environment
    type(TGeometry), pointer :: envGeo

    !> ECP Parameters
    real(dp), allocatable :: param(:,:)

    !> Inner molecule information
    integer :: nAtom, nSpecies

  end type TECPEnvInp


  !> Internal status of the ECPEnv module.
  type TECPEnv
    integer :: nAt, nSp  ! Number of inner atoms and different species
    integer :: nAtEnv, nSpEnv  ! Number of environment atoms and different species
    integer, allocatable :: speciesEnv(:)
    real(dp), allocatable :: param(:,:)
    real(dp), allocatable :: coordsEnv(:,:)
    real(dp), allocatable :: potential(:)
  contains
    procedure :: updateCoords
    procedure :: addGradientDC
    procedure :: getEnergyPerAtom
  end type TECPEnv

contains


  !> Initializes instance.
  subroutine ECPEnv_init(this, inp)

    !> Instance.
    type(TECPEnv), intent(out) :: this

    !> Input data.
    type(TECPEnvInp), intent(in) :: inp

    this%nAt = inp%nAtom  ! Number of inner molecule atoms
    this%nSp = inp%nSpecies  ! Number of inner molecule elements
    this%nAtEnv = inp%envGeo%nAtom   ! Number of environment atoms
    this%nSpEnv = inp%envGeo%nSpecies   ! Number of environment elements
    allocate(this%speciesEnv(this%nAtEnv))
    this%speciesEnv = inp%envGeo%species  ! Environment species identifier, shape: [nAtEnv]
    allocate(this%coordsEnv(3, this%nAtEnv))
    this%coordsEnv(:,:) = inp%envGeo%coords(:,:)  ! Environment Coordinates, shape: [3, nAtEnv]
    allocate(this%param(2, this%nSpEnv))
    this%param(:,:) = inp%param(:,:)   ! ECP Parameters per species, shape: [2, nSpEnv]

    allocate(this%potential(this%nAt))
    this%potential(:) = 0.0_dp

  end subroutine ECPEnv_init


  !> Updates data structures if there are changed coordinates for the instance.
  subroutine updateCoords(this, coords, species)

    !> Instance.
    class(TECPEnv), intent(inout) :: this

    !> Inner molecule coordinates
    real(dp) :: coords(:,:)

    !> Species for all atoms, shape: [nAllAtom].
    integer, intent(in) :: species(:)

    integer :: iAt, iAtEnv, iSpEnv
    real(dp) :: dist

    do iAt = 1, this%nAt
      do iAtEnv = 1, this%nAtEnv
        dist = sqrt(sum((coords(:, iAt) - this%coordsEnv(:, iAtEnv))**2))
        iSpEnv = this%speciesEnv(iAtEnv)
        this%potential(iAt) = this%potential(iAt) + getECP(this%param(:,iSpEnv), dist)
      end do
    end do

  end subroutine updateCoords


  !> Updates data structures if there are changed coordinates for the instance.
  subroutine addGradientDC(this, coords, species, derivs)

    !> Instance.
    class(TECPEnv), intent(inout) :: this

    !> Inner molecule coordinates
    real(dp) :: coords(:,:)

    !> Species for all atoms, shape: [nAllAtom].
    integer, intent(in) :: species(:)

    !> Gradient on exit.
    real(dp), intent(inout) :: derivs(:,:)

    integer :: iAt, iAtEnv, iSpEnv
    real(dp) :: dist, tmp

    do iAt = 1, this%nAt
      do iAtEnv = 1, this%nAtEnv
        dist = sqrt(sum((coords(:, iAt) - this%coordsEnv(:, iAtEnv))**2))
        iSpEnv = this%speciesEnv(iAtEnv)
        tmp = getECPDeriv(this%param(:,iSpEnv), dist)
        derivs(:, iAt) = derivs(:, iAt) + tmp * (coords(:, iAt) - this%coordsEnv(:, iAtEnv))
      end do
    end do

  end subroutine addGradientDC


  !> Returns energy per atom.
  subroutine getEnergyPerAtom(this, energyPerAtom)
    class(TECPEnv), intent(inout) :: this
    real(dp), intent(out) :: energyPerAtom(:)

    @:ASSERT(size(energyPerAtom) == this%nAt)

    energyPerAtom(:) = this%potential(:)
    print *, "ECPEnv Potential:"
    print *, this%potential(:)

  end subroutine getEnergyPerAtom


! Private routines


  !> Gets a effective core potential value
  function getECP(param, dist) result(res)
    real(dp), intent(in) :: param(:)
    real(dp), intent(in) :: dist
    real(dp) :: res

    @:ASSERT(size(param, dim=1) == 2)

    res = param(1) * exp(-param(2) * dist)

  end function getECP


  !> Gets a effective core potential value
  function getECPDeriv(param, dist) result(res)
    real(dp), intent(in) :: param(:)
    real(dp), intent(in) :: dist
    real(dp) :: res

    @:ASSERT(size(param, dim=1) == 2)

    res = param(1) * param(2) / dist * exp(-param(2) * dist)

  end function getECPDeriv

end module ecpenv

