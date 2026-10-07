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

module EDA_Energy

use Definitions, only: wp, iwp
use EDA_Reference, only: EDA_Reference_Type
use EDA_Integrals, only: EDA_Read_OneElectron, EDA_TwoElectron_Energies

implicit none
private

public :: EDA_Results_Type, EDA_Compute_ClosedShell

type :: EDA_Results_Type
  real(kind=wp) :: E0 = 0.0_wp
  real(kind=wp) :: E1 = 0.0_wp
  real(kind=wp) :: E2 = 0.0_wp
  real(kind=wp) :: Electrostatic = 0.0_wp
  real(kind=wp) :: Exchange = 0.0_wp
  real(kind=wp) :: Repulsion = 0.0_wp
  real(kind=wp) :: Polarization = 0.0_wp
  real(kind=wp) :: HFInteraction = 0.0_wp
  real(kind=wp) :: Dispersion = 0.0_wp
  real(kind=wp) :: MP2Interaction = 0.0_wp
  real(kind=wp) :: Closure = 0.0_wp
  logical :: HasMP2 = .false.
end type EDA_Results_Type

contains

subroutine EDA_Compute_ClosedShell(Super,Frag,nFrag,FragBasis,Results,ierr,Message)

  type(EDA_Reference_Type), intent(in) :: Super
  type(EDA_Reference_Type), intent(in) :: Frag(:)
  integer(kind=iwp), intent(in) :: nFrag
  character(len=*), intent(in) :: FragBasis
  type(EDA_Results_Type), intent(out) :: Results
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  integer(kind=iwp) :: nBas,nOcc,nOccSum,nBasSum
  integer(kind=iwp) :: iFrag, i0, i1, mu, nu, p
  integer(kind=iwp) :: info, lWork
  real(kind=wp) :: Enuc,OneRaw,OneOrth,Jraw,Kraw,Jorth,Korth
  real(kind=wp) :: SumEHF,SumEMP2,OverlapMin
  real(kind=wp), allocatable :: H(:,:),S(:,:),C(:,:),Tmp(:,:),OccS(:,:),Eig(:),Work(:),InvOccS(:,:)
  real(kind=wp), allocatable :: Praw(:,:),Porth(:,:),Pfrag(:,:,:),Kfrag(:)
  logical :: AllMP2

  ierr = 1
  Message = ''
  Results = EDA_Results_Type()

  nBas = Super%nBas
  nOcc = Super%nOcc
  if ((nBas < 1) .or. (nOcc < 1)) then
    Message = 'Invalid supermolecule RHF dimensions.'
    return
  end if

  nOccSum = 0
  nBasSum = 0
  do iFrag=1,nFrag
    nOccSum = nOccSum + Frag(iFrag)%nOcc
    nBasSum = nBasSum + Frag(iFrag)%nBas
  end do
  if (nOccSum /= nOcc) then
    Message = 'Sum of fragment occupied orbitals differs from the supermolecule RHF occupation.'
    return
  end if

  if (trim(FragBasis) == 'OWN') then
    if (nBasSum /= nBas) then
      Message = 'FBASIS=OWN requires the sum of fragment AO dimensions to equal the supermolecule AO dimension.'
      return
    end if
  else if (trim(FragBasis) == 'SUPER') then
    do iFrag=1,nFrag
      if (Frag(iFrag)%nBas /= nBas) then
        Message = 'FBASIS=SUPER requires every fragment RunFile to use the full supermolecule AO basis.'
        return
      end if
    end do
  else
    Message = 'Internal EDA error: unknown fragment-basis convention.'
    return
  end if

  allocate(H(nBas,nBas),S(nBas,nBas),C(nBas,nOcc),Tmp(nBas,nOcc))
  allocate(OccS(nOcc,nOcc),Eig(nOcc),InvOccS(nOcc,nOcc))
  allocate(Praw(nBas,nBas),Porth(nBas,nBas),Pfrag(nBas,nBas,nFrag),Kfrag(nFrag))

  call EDA_Read_OneElectron(nBas,H,S,Enuc,ierr,Message)
  if (ierr /= 0) goto 900

  C = 0.0_wp
  Pfrag = 0.0_wp
  i0 = 1
  i1 = 0

  if (trim(FragBasis) == 'OWN') then
    p = 0
    do iFrag=1,nFrag
      i1 = i0 + Frag(iFrag)%nOcc - 1
      do nu=1,Frag(iFrag)%nOcc
        do mu=1,Frag(iFrag)%nBas
          C(p+mu,i0+nu-1) = Frag(iFrag)%COcc(mu,nu)
        end do
      end do
      call DGEMM_('N','T',nBas,nBas,Frag(iFrag)%nOcc,2.0_wp, &
                  C(1,i0),nBas,C(1,i0),nBas,0.0_wp,Pfrag(1,1,iFrag),nBas)
      p = p + Frag(iFrag)%nBas
      i0 = i1 + 1
    end do
  else
    do iFrag=1,nFrag
      i1 = i0 + Frag(iFrag)%nOcc - 1
      C(:,i0:i1) = Frag(iFrag)%COcc(:,:)
      call DGEMM_('N','T',nBas,nBas,Frag(iFrag)%nOcc,2.0_wp, &
                  C(1,i0),nBas,C(1,i0),nBas,0.0_wp,Pfrag(1,1,iFrag),nBas)
      i0 = i1 + 1
    end do
  end if

  Praw = 0.0_wp
  do iFrag=1,nFrag
    Praw(:,:) = Praw(:,:) + Pfrag(:,:,iFrag)
  end do

  ! Occupied-orbital overlap: O = C^T S C.
  call DGEMM_('N','N',nBas,nOcc,nBas,1.0_wp,S,nBas,C,nBas,0.0_wp,Tmp,nBas)
  call DGEMM_('T','N',nOcc,nOcc,nBas,1.0_wp,C,nBas,Tmp,nBas,0.0_wp,OccS,nOcc)

  ! Symmetric eigendecomposition O = U diag(s) U^T, then O^-1.
  allocate(Work(1))
  lWork = -1
  call DSYEV_('V','U',nOcc,OccS,nOcc,Eig,Work,lWork,info)
  if (info /= 0) then
    Message = 'Failed to query the occupied-overlap diagonalization workspace.'
    goto 900
  end if
  lWork = max(1_iwp,int(Work(1),kind=iwp))
  deallocate(Work)
  allocate(Work(lWork))
  call DSYEV_('V','U',nOcc,OccS,nOcc,Eig,Work,lWork,info)
  if (info /= 0) then
    Message = 'Failed to diagonalize the occupied fragment-orbital overlap matrix.'
    goto 900
  end if

  OverlapMin = minval(Eig)
  if (OverlapMin <= 1.0e-10_wp) then
    Message = 'Fragment occupied orbitals are linearly dependent in the supermolecule overlap metric.'
    goto 900
  end if

  InvOccS = 0.0_wp
  do p=1,nOcc
    do nu=1,nOcc
      do mu=1,nOcc
        InvOccS(mu,nu) = InvOccS(mu,nu) + OccS(mu,p)*OccS(nu,p)/Eig(p)
      end do
    end do
  end do

  ! Porth = 2 C O^-1 C^T.
  call DGEMM_('N','N',nBas,nOcc,nOcc,1.0_wp,C,nBas,InvOccS,nOcc,0.0_wp,Tmp,nBas)
  call DGEMM_('N','T',nBas,nBas,nOcc,2.0_wp,Tmp,nBas,C,nBas,0.0_wp,Porth,nBas)

  OneRaw = sum(Praw*H)
  OneOrth = sum(Porth*H)

  call EDA_TwoElectron_Energies(nBas,nFrag,Praw,Porth,Pfrag,Jraw,Kraw,Jorth,Korth,Kfrag,ierr,Message)
  if (ierr /= 0) goto 900

  Results%E0 = OneRaw + 0.5_wp*Jraw - 0.25_wp*sum(Kfrag) + Enuc
  Results%E1 = OneRaw + 0.5_wp*Jraw - 0.25_wp*Kraw + Enuc
  Results%E2 = OneOrth + 0.5_wp*Jorth - 0.25_wp*Korth + Enuc

  SumEHF = 0.0_wp
  SumEMP2 = 0.0_wp
  AllMP2 = Super%HasMP2
  do iFrag=1,nFrag
    SumEHF = SumEHF + Frag(iFrag)%EHF
    SumEMP2 = SumEMP2 + Frag(iFrag)%ETotal
    AllMP2 = AllMP2 .and. Frag(iFrag)%HasMP2
  end do

  Results%Electrostatic = Results%E0 - SumEHF
  Results%Exchange = Results%E1 - Results%E0
  Results%Repulsion = Results%E2 - Results%E1
  Results%Polarization = Super%EHF - Results%E2
  Results%HFInteraction = Super%EHF - SumEHF
  Results%Closure = Results%HFInteraction - (Results%Electrostatic + Results%Exchange + &
                    Results%Repulsion + Results%Polarization)

  Results%HasMP2 = AllMP2
  if (AllMP2) then
    Results%MP2Interaction = Super%ETotal - SumEMP2
    Results%Dispersion = Results%MP2Interaction - Results%HFInteraction
  end if

  ierr = 0
  Message = ''

900 continue
  if (allocated(H)) deallocate(H)
  if (allocated(S)) deallocate(S)
  if (allocated(C)) deallocate(C)
  if (allocated(Tmp)) deallocate(Tmp)
  if (allocated(OccS)) deallocate(OccS)
  if (allocated(Eig)) deallocate(Eig)
  if (allocated(Work)) deallocate(Work)
  if (allocated(InvOccS)) deallocate(InvOccS)
  if (allocated(Praw)) deallocate(Praw)
  if (allocated(Porth)) deallocate(Porth)
  if (allocated(Pfrag)) deallocate(Pfrag)
  if (allocated(Kfrag)) deallocate(Kfrag)

end subroutine EDA_Compute_ClosedShell

end module EDA_Energy
