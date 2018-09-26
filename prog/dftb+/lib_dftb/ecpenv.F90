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
    allocate(this%param(3, this%nSpEnv))
    this%param(:,:) = inp%param(:,:)   ! ECP Parameters per species, shape: [2, nSpEnv]

    allocate(this%potential(this%nAt))

  end subroutine ECPEnv_init


  !> Updates data structures if there are changed coordinates for the instance.
  subroutine updateCoords(this, coords, species)

    !> Instance.
    class(TECPEnv), intent(inout) :: this

    !> Inner molecule coordinates
    real(dp) :: coords(:,:)

    !> Species for all atoms, shape: [nAllAtom].
    integer, intent(in) :: species(:)

    integer :: iAt, iAtEnv, iSp, iSpEnv
    real(dp) :: rr0, alpha, epsilon, dist

    this%potential(:) = 0.0_dp
    do iAt = 1, this%nAt
      do iAtEnv = 1, this%nAtEnv
        dist = sqrt(sum((coords(:, iAt) - this%coordsEnv(:, iAtEnv))**2))
        iSp = species(iAt)
        iSpEnv = this%speciesEnv(iAtEnv)
        epsilon = sqrt(this%param(1, iSp) * this%param(1, iSpEnv))
        alpha = 0.5_dp * (this%param(2, iSp) + this%param(2, iSpEnv))
        rr0 = 0.5_dp * (this%param(3, iSp) + this%param(3, iSpEnv))
        this%potential(iAt) = this%potential(iAt) + getECP(epsilon, alpha, rr0, dist)
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

    integer :: iAt, iAtEnv, iSp, iSpEnv, ii
    real(dp) :: rr0, alpha, epsilon, tmpR1, dist
    real(dp) :: eperatom(this%nAt), tmpR2, tmpgrad(3,this%nAt)
    real(dp) :: tmpVec1(3), tmpVec2(3), tmpCoord(3,this%nAt)

    tmpgrad(:, :) = 0.0_dp
    do iAt = 1, this%nAt
      do iAtEnv = 1, this%nAtEnv
        dist = sqrt(sum((coords(:, iAt) - this%coordsEnv(:, iAtEnv))**2))
        iSp = species(iAt)
        iSpEnv = this%speciesEnv(iAtEnv)
        epsilon = sqrt(this%param(1, iSp) * this%param(1, iSpEnv))
        alpha = 0.5_dp * (this%param(2, iSp) + this%param(2, iSpEnv))
        rr0 = 0.5_dp * (this%param(3, iSp) + this%param(3, iSpEnv))
        tmpR1 = getECPDeriv(epsilon, alpha, rr0, dist)
        tmpgrad(:, iAt) = tmpgrad(:, iAt) - tmpR1 * (coords(:, iAt) - this%coordsEnv(:, iAtEnv))
      end do
    end do

    derivs(:, :) = derivs(:, :) - tmpgrad(:, :)

  end subroutine addGradientDC


  !> Returns energy per atom.
  subroutine getEnergyPerAtom(this, energyPerAtom)
    class(TECPEnv), intent(inout) :: this
    real(dp), intent(out) :: energyPerAtom(:)

    @:ASSERT(size(energyPerAtom) == this%nAt)

    energyPerAtom(:) = this%potential(:)

  end subroutine getEnergyPerAtom


! Private routines


  !> Gets a effective core potential value
  function getECP(epsilon, alpha, rr0, dist) result(res)
    real(dp), intent(in) :: epsilon
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: rr0
    real(dp), intent(in) :: dist
    real(dp) :: tmpR1, tmpR2
    real(dp) :: res

    if (rr0 <= tolSameDist) then
      res = 0.0_dp
    else if (alpha <= tolSameDist) then
      res = -epsilon
    else
      tmpR1 = epsilon / (alpha - 6.0_dp)
      tmpR2 = alpha - alpha * dist / rr0
      res = tmpR1 * (6.0_dp * exp(tmpR2) - alpha * (rr0 / dist)**6)
    end if

  end function getECP


  !> Gets a effective core potential value
  function getECPDeriv(epsilon, alpha, rr0, dist) result(res)
    real(dp), intent(in) :: epsilon
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: rr0
    real(dp), intent(in) :: dist
    real(dp) :: tmpR1, tmpR2
    real(dp) :: res

    tmpR1 = 6.0_dp * epsilon * alpha / (alpha - 6.0_dp)
    tmpR2 = alpha - alpha * dist / rr0
    res = tmpR1 * (rr0**6 / dist**8 - exp(tmpR2) / rr0 / dist)

  end function getECPDeriv

end module ecpenv

