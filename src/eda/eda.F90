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

subroutine EDA(iReturn)

use Definitions, only: wp, iwp, u6
use Constants, only: auTokcalmol
use Para_Info, only: nProcs
use EDA_Data, only: EDA_nFrag, EDA_nAtom, EDA_RefMethod, EDA_FragBasis
use EDA_Reference, only: EDA_Reference_Type, EDA_Read_Super, EDA_Read_Fragment, EDA_Free_Reference
use EDA_Energy, only: EDA_Results_Type, EDA_Compute_ClosedShell

implicit none
#include "warnings.h"

integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: ierr, iFrag, nAtomsExpected
character(len=256) :: Message
type(EDA_Reference_Type) :: Super
type(EDA_Reference_Type), allocatable :: Frag(:)
type(EDA_Results_Type) :: Results

call EDA_Read_Input(ierr)
if (ierr /= _RC_ALL_IS_WELL_) then
  iReturn = ierr
  return
end if

if (nProcs /= 1) then
  write(u6,'(6X,A)') 'EDA stage 1 is serial only. Run with MOLCAS_NPROCS=1.'
  iReturn = _RC_INPUT_ERROR_
  return
end if

allocate(Frag(EDA_nFrag))

call EDA_Read_Super(Super,ierr,Message)
if (ierr /= 0) goto 900

nAtomsExpected = sum(EDA_nAtom(1:EDA_nFrag))
if (Super%nAtoms /= nAtomsExpected) then
  write(Message,'(A,I0,A,I0,A)') 'EDA NATOM sum is ',nAtomsExpected, &
                                  ' but the supermolecule RunFile contains ',Super%nAtoms,' atoms.'
  ierr = 1
  goto 900
end if

if ((trim(EDA_RefMethod) == 'RHF') .and. (Super%Method /= 'RHF-SCF ')) then
  Message = 'REFERENCE=RHF requires the current supermolecule RunFile to end at RHF-SCF.'
  ierr = 1
  goto 900
end if
if ((trim(EDA_RefMethod) == 'MBPT2') .and. (.not. Super%HasMP2)) then
  Message = 'REFERENCE=MBPT2 requires the current supermolecule RunFile to contain an MBPT2 result.'
  ierr = 1
  goto 900
end if

fragment_loop: do iFrag=1,EDA_nFrag
  call EDA_Read_Fragment(iFrag,Frag(iFrag),ierr,Message)
  if (ierr /= 0) exit fragment_loop
  if ((trim(EDA_RefMethod) == 'RHF') .and. (Frag(iFrag)%Method /= 'RHF-SCF ')) then
    write(Message,'(A,I0,A)') 'REFERENCE=RHF but fragment ',iFrag,' is not an RHF-SCF RunFile.'
    ierr = 1
    exit fragment_loop
  end if
  if ((trim(EDA_RefMethod) == 'MBPT2') .and. (.not. Frag(iFrag)%HasMP2)) then
    write(Message,'(A,I0,A)') 'REFERENCE=MBPT2 but fragment ',iFrag,' has no MBPT2 energy.'
    ierr = 1
    exit fragment_loop
  end if
end do fragment_loop
if (ierr /= 0) goto 900

call EDA_Compute_ClosedShell(Super,Frag,EDA_nFrag,EDA_FragBasis,Results,ierr,Message)
if (ierr /= 0) goto 900

call PrintResults(Super,Frag,EDA_nFrag,Results)

iReturn = _RC_ALL_IS_WELL_
goto 990

900 continue
write(u6,'(/,6X,2A)') 'EDA failed: ',trim(Message)
iReturn = _RC_GENERAL_ERROR_

990 continue
call EDA_Free_Reference(Super)
do iFrag=1,size(Frag)
  call EDA_Free_Reference(Frag(iFrag))
end do
if (allocated(Frag)) deallocate(Frag)

contains

subroutine PrintResults(SuperRef,FragRef,nFrag,Res)

  type(EDA_Reference_Type), intent(in) :: SuperRef
  type(EDA_Reference_Type), intent(in) :: FragRef(:)
  integer(kind=iwp), intent(in) :: nFrag
  type(EDA_Results_Type), intent(in) :: Res
  integer(kind=iwp) :: i
  real(kind=wp) :: SumEHF

  SumEHF = 0.0_wp
  do i=1,nFrag
    SumEHF = SumEHF + FragRef(i)%EHF
  end do

  write(u6,'(/,6X,A)') '=============================================================='
  write(u6,'(6X,A)')   '  Su-Li closed-shell CMOEDA prototype (OpenMolcas stage 1)'
  write(u6,'(6X,A)')   '=============================================================='
  write(u6,'(6X,A)')   'This is the Su-Li electrostatic/exchange/repulsion/polarization'
  write(u6,'(6X,A)')   'decomposition, not the legacy Kitaura-Morokuma PL/CT/MIX split.'
  write(u6,'(6X,A,1X,A)') 'Fragment basis convention:',trim(EDA_FragBasis)
  write(u6,'(6X,A,I0)') 'Supermolecule basis functions:',SuperRef%nBas
  write(u6,'(6X,A,I0)') 'Supermolecule occupied RHF orbitals:',SuperRef%nOcc
  write(u6,'(6X,A,F20.10,A)') 'Sum of fragment RHF energies:',SumEHF,' a.u.'
  write(u6,'(6X,A,F20.10,A)') 'Supermolecule RHF energy:   ',SuperRef%EHF,' a.u.'
  write(u6,'(/,6X,A)') 'Intermediate frozen-fragment energies'
  write(u6,'(8X,A,T54,F20.10)') 'E(0): non-antisymmetric, non-orthogonal',Res%E0
  write(u6,'(8X,A,T54,F20.10)') 'E(1): antisymmetric, non-orthogonal    ',Res%E1
  write(u6,'(8X,A,T54,F20.10)') 'E(2): antisymmetric, orthogonalized    ',Res%E2

  write(u6,'(/,6X,A)') 'CMOEDA interaction-energy decomposition'
  write(u6,'(6X,A,T43,A,T61,A)') 'Component','Hartree','kcal/mol'
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Electrostatic',Res%Electrostatic,Res%Electrostatic*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Exchange',Res%Exchange,Res%Exchange*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Repulsion',Res%Repulsion,Res%Repulsion*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'Polarization',Res%Polarization,Res%Polarization*auTokcalmol
  write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'HF interaction',Res%HFInteraction,Res%HFInteraction*auTokcalmol
  if (Res%HasMP2) then
    write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'MP2 dispersion',Res%Dispersion,Res%Dispersion*auTokcalmol
    write(u6,'(6X,A,T43,F16.8,T61,F16.4)') 'MP2 interaction',Res%MP2Interaction,Res%MP2Interaction*auTokcalmol
  end if
  write(u6,'(6X,A,T43,ES16.6)') 'HF closure residual',Res%Closure
  write(u6,'(6X,A/)') '=============================================================='

end subroutine PrintResults

end subroutine EDA
