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
    real(dp), allocatable :: paramEnv(:,:)
    real(dp), allocatable :: paramInner(:,:)

    !> Inner molecule information
    integer :: nAtom, nSpecies

  end type TECPEnvInp


  !> Internal status of the ECPEnv module.
  type TECPEnv
    integer :: nAt, nSp  ! Number of inner atoms and different species
    integer :: nAtEnv, nSpEnv  ! Number of environment atoms and different species
    integer, allocatable :: speciesEnv(:)
    real(dp), allocatable :: paramEnv(:,:)
    real(dp), allocatable :: paramInner(:,:)
    real(dp), allocatable :: paramLJ(:,:,:)
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
    allocate(this%paramEnv(4, this%nSpEnv))
    this%paramEnv(:,:) = inp%paramEnv(:,:)   ! ECP Parameters per species, shape: [4, nSpEnv]
    allocate(this%paramInner(4, this%nSp))
    this%paramInner(:,:) = inp%paramInner(:,:)   ! ECP Parameters per species, shape: [4, nSp]
    allocate(this%paramLJ(3, this%nSp, this%nSpEnv))  ! Switch parameters, shape [3, nSp, nSpEnv]
    call initLennardJones(this)
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
    real(dp) :: dist, rr0, alpha, epsilon, rIP, lj1, lj2, tmpR1

    this%potential(:) = 0.0_dp
    do iAt = 1, this%nAt
      do iAtEnv = 1, this%nAtEnv
        dist = sqrt(sum((coords(:, iAt) - this%coordsEnv(:, iAtEnv))**2))
        iSp = species(iAt)
        iSpEnv = this%speciesEnv(iAtEnv)
        epsilon = sqrt(this%paramInner(1, iSp) * this%paramEnv(1, iSpEnv))
        alpha = 0.5_dp * (this%paramInner(2, iSp) + this%paramEnv(2, iSpEnv))
        rr0 = 0.5_dp * (this%paramInner(3, iSp) + this%paramEnv(3, iSpEnv))
        rIP = this%paramLJ(1, iSp, iSpEnv)
        lj1 = this%paramLJ(2, iSp, iSpEnv)
        lj2 = this%paramLJ(3, iSp, iSpEnv)
        tmpR1 = getECP(dist, epsilon, alpha, rr0, rIP, lj1, lj2)
        this%potential(iAt) = this%potential(iAt) + tmpR1
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
    real(dp) :: dist, rr0, alpha, epsilon, rIP, lj1, tmpR1
    real(dp) :: eperatom(this%nAt), tmpgrad(3,this%nAt)

    tmpgrad(:, :) = 0.0_dp
    do iAt = 1, this%nAt
      do iAtEnv = 1, this%nAtEnv
        dist = sqrt(sum((coords(:, iAt) - this%coordsEnv(:, iAtEnv))**2))
        iSp = species(iAt)
        iSpEnv = this%speciesEnv(iAtEnv)
        epsilon = sqrt(this%paramInner(1, iSp) * this%paramEnv(1, iSpEnv))
        alpha = 0.5_dp * (this%paramInner(2, iSp) + this%paramEnv(2, iSpEnv))
        rr0 = 0.5_dp * (this%paramInner(3, iSp) + this%paramEnv(3, iSpEnv))
        rIP = this%paramLJ(1, iSp, iSpEnv)
        lj1 = this%paramLJ(2, iSp, iSpEnv)
        tmpR1 = getECPDeriv(dist, epsilon, alpha, rr0, rIP, lj1)
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

  !> Creates parameters for inflection point LJ switch
  subroutine initLennardJones(this)

    !> Instance.
    class(TECPEnv), intent(inout) :: this

    integer :: iSp, iSpEnv
    real(dp) :: rIP  ! Inflection point
    real(dp) :: epsilon, alpha, rr0
    real(dp) :: potIP, potIPDeriv  ! Buckingham potential (and deriv) at rIP

    do iSp = 1, this%nSp
      do iSpEnv = 1, this%nSpEnv
        if (this%paramInner(4, iSp) >= this%paramEnv(4, iSpEnv)) then
          rIP = this%paramInner(4, iSp)
        else
          rIP = this%paramEnv(4, iSp)
        end if
        epsilon = sqrt(this%paramInner(1, iSp) * this%paramEnv(1, iSpEnv))
        alpha = 0.5_dp * (this%paramInner(2, iSp) + this%paramEnv(2, iSpEnv))
        rr0 = 0.5_dp * (this%paramInner(3, iSp) + this%paramEnv(3, iSpEnv))
        potIP = getECP(rIP, epsilon, alpha, rr0, 0.0_dp, 0.0_dp, 0.0_dp)
        potIPDeriv = getECPDeriv(rIP, epsilon, alpha, rr0, 0.0_dp, 0.0_dp)
        this%paramLJ(1, iSp, iSpEnv) = rIP
        this%paramLJ(2, iSp, iSpEnv) = - rIP**13 * potIPDeriv / 12.0_dp
        this%paramLJ(3, iSp, iSpEnv) = potIP + rIP * potIPDeriv / 12.0_dp
      end do
    end do

  end subroutine initLennardJones


  !> Gets a effective core potential value
  function getECP(dist, epsilon, alpha, rr0, rIP, lj1, lj2) result(res)
    real(dp), intent(in) :: dist
    real(dp), intent(in) :: epsilon
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: rr0
    real(dp), intent(in) :: rIP
    real(dp), intent(in) :: lj1
    real(dp), intent(in) :: lj2
    real(dp) :: tmpR1, tmpR2
    real(dp) :: res

    if (dist < rIP) then  ! Use Lennard-Jones potential
      res = lj1 / dist**12 + lj2
    else
      if (rr0 <= tolSameDist) then
        res = 0.0_dp
      else if (alpha <= tolSameDist) then
        res = -epsilon
      else  ! Use Buckingham potential
        tmpR1 = epsilon / (alpha - 6.0_dp)
        tmpR2 = alpha - alpha * dist / rr0
        res = tmpR1 * (6.0_dp * exp(tmpR2) - alpha * (rr0 / dist)**6)
      end if
    end if

  end function getECP


  !> Gets the derivative of the  effective core potential
  function getECPDeriv(dist, epsilon, alpha, rr0, rIP, lj1) result(res)
    real(dp), intent(in) :: dist
    real(dp), intent(in) :: epsilon
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: rr0
    real(dp), intent(in) :: rIP
    real(dp), intent(in) :: lj1
    real(dp) :: tmpR1, tmpR2
    real(dp) :: res

    if (dist < rIP) then  ! Use Lennard-Jones derivative
      res = -12.0_dp * lj1 / dist**13
    else  ! Use Buckingham derivative
      tmpR1 = 6.0_dp * epsilon * alpha / (alpha - 6.0_dp)
      tmpR2 = alpha - alpha * dist / rr0
      res = tmpR1 * (rr0**6 / dist**8 - exp(tmpR2) / rr0 / dist)
    end if

  end function getECPDeriv

end module ecpenv

