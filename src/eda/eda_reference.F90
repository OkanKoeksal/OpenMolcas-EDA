!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module EDA_Reference

use Definitions, only: wp, iwp

implicit none
private

character(len=8), parameter :: FragRunName(10) = [ &
  'EDAFRG1 ', 'EDAFRG2 ', 'EDAFRG3 ', 'EDAFRG4 ', 'EDAFRG5 ', &
  'EDAFRG6 ', 'EDAFRG7 ', 'EDAFRG8 ', 'EDAFRG9 ', 'EDAFR10 ' ]

type, public :: EDA_Reference_Type
  integer(kind=iwp) :: nBas = 0
  integer(kind=iwp) :: nOrb = 0
  integer(kind=iwp) :: nOcc = 0
  integer(kind=iwp) :: nAtoms = 0
  real(kind=wp) :: EHF = 0.0_wp
  real(kind=wp) :: ETotal = 0.0_wp
  logical :: HasMP2 = .false.
  character(len=8) :: Method = '        '
  real(kind=wp), allocatable :: COcc(:,:)
end type EDA_Reference_Type

public :: EDA_Read_Super, EDA_Read_Fragment, EDA_Free_Reference

contains

subroutine EDA_Read_Super(Ref,ierr,Message)

  type(EDA_Reference_Type), intent(inout) :: Ref
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  call Read_Current_RunFile(Ref,ierr,Message)

end subroutine EDA_Read_Super

subroutine EDA_Read_Fragment(iFrag,Ref,ierr,Message)

  integer(kind=iwp), intent(in) :: iFrag
  type(EDA_Reference_Type), intent(inout) :: Ref
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  logical(kind=iwp) :: Exists

  ierr = 1
  Message = ''

  if ((iFrag < 1) .or. (iFrag > size(FragRunName))) then
    Message = 'EDA fragment index is outside the supported range.'
    return
  end if

  call f_Inquire(trim(FragRunName(iFrag)),Exists)
  if (.not. Exists) then
    write(Message,'(A,I0,A)') 'Missing EDA fragment RunFile ',iFrag,'.'
    return
  end if

  call NameRun(trim(FragRunName(iFrag)))
  call Read_Current_RunFile(Ref,ierr,Message)
  call NameRun('#Pop')

end subroutine EDA_Read_Fragment

subroutine Read_Current_RunFile(Ref,ierr,Message)

  type(EDA_Reference_Type), intent(inout) :: Ref
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  integer(kind=iwp) :: nSym, nData
  integer(kind=iwp) :: nBasA(1), nOrbA(1), nOccA(1)
  logical(kind=iwp) :: Found
  real(kind=wp), allocatable :: CFull(:,:)

  ierr = 1
  Message = 'Unable to read an RHF/MBPT2 reference from RunFile.'

  call EDA_Free_Reference(Ref)

  call Get_iScalar('nSym',nSym)
  if (nSym /= 1) then
    Message = 'EDA stage 1 requires C1 (NoSym) references.'
    return
  end if

  call Get_iArray('nBas',nBasA,1)
  call Get_iArray('nOrb',nOrbA,1)
  call Get_iArray('nIsh',nOccA,1)
  call Get_iScalar('Unique atoms',Ref%nAtoms)
  call Get_cArray('Relax Method',Ref%Method,8)

  Ref%nBas = nBasA(1)
  Ref%nOrb = nOrbA(1)
  Ref%nOcc = nOccA(1)

  if ((Ref%nBas < 1) .or. (Ref%nOrb < Ref%nOcc) .or. (Ref%nOcc < 1)) then
    Message = 'Invalid basis/orbital dimensions in EDA reference RunFile.'
    return
  end if

  if ((Ref%Method /= 'RHF-SCF ') .and. (Ref%Method /= 'MBPT2   ')) then
    Message = 'EDA stage 1 accepts only RHF-SCF or closed-shell MBPT2 references.'
    return
  end if

  call qpg_dArray('SCF orbitals',Found,nData)
  if (.not. Found) then
    Message = 'SCF orbitals are missing from EDA reference RunFile.'
    return
  end if
  if (nData /= Ref%nBas*Ref%nOrb) then
    Message = 'SCF orbital array size is inconsistent with nBas*nOrb.'
    return
  end if

  allocate(CFull(Ref%nBas,Ref%nOrb))
  allocate(Ref%COcc(Ref%nBas,Ref%nOcc))
  call Get_dArray('SCF orbitals',CFull,nData)
  Ref%COcc(:,:) = CFull(:,1:Ref%nOcc)
  deallocate(CFull)

  call Get_dScalar('SCF energy',Ref%EHF)
  call Get_dScalar('Last energy',Ref%ETotal)
  Ref%HasMP2 = (Ref%Method == 'MBPT2   ')
  if (.not. Ref%HasMP2) Ref%ETotal = Ref%EHF

  ierr = 0
  Message = ''

end subroutine Read_Current_RunFile

subroutine EDA_Free_Reference(Ref)

  type(EDA_Reference_Type), intent(inout) :: Ref

  if (allocated(Ref%COcc)) deallocate(Ref%COcc)
  Ref%nBas = 0
  Ref%nOrb = 0
  Ref%nOcc = 0
  Ref%nAtoms = 0
  Ref%EHF = 0.0_wp
  Ref%ETotal = 0.0_wp
  Ref%HasMP2 = .false.
  Ref%Method = '        '

end subroutine EDA_Free_Reference

end module EDA_Reference
