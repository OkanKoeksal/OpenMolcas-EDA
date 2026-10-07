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

module EDA_Integrals

use Definitions, only: wp, iwp

implicit none
private
public :: EDA_Read_OneElectron, EDA_TwoElectron_Energies

contains

subroutine EDA_Read_OneElectron(nBas,H,S,Enuc,ierr,Message)

  use OneDat, only: sNoNuc, sNoOri

  integer(kind=iwp), intent(in) :: nBas
  real(kind=wp), intent(out) :: H(nBas,nBas), S(nBas,nBas)
  real(kind=wp), intent(out) :: Enuc
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  integer(kind=iwp) :: nPair, rc, iOpt, iComp, iSym, LuOne
  logical :: Opened
  character(len=8) :: Label
  real(kind=wp), allocatable :: Buf(:)

  ierr = 1
  Message = 'Unable to read ONEINT for EDA.'
  H = 0.0_wp
  S = 0.0_wp
  Enuc = 0.0_wp
  Opened = .false.
  LuOne = 41
  nPair = nBas*(nBas+1)/2
  allocate(Buf(nPair+4))

  call OpnOne(rc,0,'ONEINT',LuOne)
  if (rc /= 0) goto 900
  Opened = .true.

  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  iSym = 1

  Label = 'OneHam  '
  call RdOne(rc,iOpt,Label,iComp,Buf,iSym)
  if (rc /= 0) goto 900
  call UnpackSymmetric(Buf,H,nBas)

  Label = 'Mltpl  0'
  call RdOne(rc,iOpt,Label,iComp,Buf,iSym)
  if (rc /= 0) goto 900
  call UnpackSymmetric(Buf,S,nBas)

  call Get_dScalar('PotNuc',Enuc)

  ierr = 0
  Message = ''

900 continue
  if (Opened) call ClsOne(rc,0)
  if (allocated(Buf)) deallocate(Buf)

end subroutine EDA_Read_OneElectron

subroutine EDA_TwoElectron_Energies(nBas,nFrag,Praw,Porth,Pfrag, &
                                     Jraw,Kraw,Jorth,Korth,Kfrag,ierr,Message)

  integer(kind=iwp), intent(in) :: nBas, nFrag
  real(kind=wp), intent(in) :: Praw(nBas,nBas), Porth(nBas,nBas)
  real(kind=wp), intent(in) :: Pfrag(nBas,nBas,nFrag)
  real(kind=wp), intent(out) :: Jraw,Kraw,Jorth,Korth
  real(kind=wp), intent(out) :: Kfrag(nFrag)
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message

  integer(kind=iwp) :: nPair, rc, CloseRC, LuOrd
  integer(kind=iwp) :: FileNSym, FileBas(8), FileSkip(8)
  integer(kind=iwp) :: i,j,k,l,q,nMat,iOpt,iFrag
  integer(kind=iwp) :: mij,mkl
  logical(kind=iwp) :: Exists, Square
  logical :: Opened
  real(kind=wp) :: v, JFactor
  real(kind=wp), allocatable :: Buf(:)

  ierr = 1
  Message = 'Unable to read conventional two-electron integrals for EDA.'
  Jraw = 0.0_wp
  Kraw = 0.0_wp
  Jorth = 0.0_wp
  Korth = 0.0_wp
  Kfrag = 0.0_wp
  Opened = .false.
  LuOrd = 43
  nPair = nBas*(nBas+1)/2

  call f_Inquire('ORDINT',Exists)
  if (.not. Exists) then
    Message = 'EDA stage 1 requires stored ORDINT integrals: use &SEWARD with NoCholesky.'
    return
  end if

  allocate(Buf(nPair+1))
  call OpnOrd(rc,0,'ORDINT',LuOrd)
  if (rc /= 0) goto 900
  Opened = .true.

  call GetOrd(rc,Square,FileNSym,FileBas,FileSkip)
  if (rc /= 0) goto 900
  if ((FileNSym /= 1) .or. (FileBas(1) /= nBas) .or. (FileSkip(1) /= 0)) then
    Message = 'ORDINT dimensions do not match the C1 EDA supermolecule reference.'
    goto 900
  end if

  iOpt = 1
  do i=1,nBas
    do j=1,i
      Buf = 0.0_wp
      call RdOrd_(rc,iOpt,1,1,1,1,Buf,nPair+1,nMat)
      iOpt = 2
      if ((rc /= 0) .or. (nMat /= 1)) goto 900

      if (i == j) then
        mij = 1
      else
        mij = 2
      end if

      q = 0
      do k=1,nBas
        do l=1,k
          q = q+1
          v = Buf(q)
          if (k == l) then
            mkl = 1
          else
            mkl = 2
          end if
          JFactor = real(mij*mkl,wp)*v

          Jraw = Jraw + JFactor*Praw(i,j)*Praw(k,l)
          Jorth = Jorth + JFactor*Porth(i,j)*Porth(k,l)

          Kraw = Kraw + v*ExchangeMultiplicity(Praw,i,j,k,l)
          Korth = Korth + v*ExchangeMultiplicity(Porth,i,j,k,l)
          do iFrag=1,nFrag
            Kfrag(iFrag) = Kfrag(iFrag) + v*ExchangeMultiplicity(Pfrag(:,:,iFrag),i,j,k,l)
          end do
        end do
      end do
    end do
  end do

  ierr = 0
  Message = ''

900 continue
  if (Opened) then
    call ClsOrd(CloseRC)
    if ((CloseRC /= 0) .and. (ierr == 0)) then
      ierr = 1
      Message = 'Unable to close ORDINT after EDA.'
    end if
  end if
  if (allocated(Buf)) deallocate(Buf)

contains

  real(kind=wp) function ExchangeMultiplicity(P,a,b,c,d)

    real(kind=wp), intent(in) :: P(nBas,nBas)
    integer(kind=iwp), intent(in) :: a,b,c,d

    ExchangeMultiplicity = P(a,c)*P(b,d)
    if (a /= b) ExchangeMultiplicity = ExchangeMultiplicity + P(b,c)*P(a,d)
    if (c /= d) ExchangeMultiplicity = ExchangeMultiplicity + P(a,d)*P(b,c)
    if ((a /= b) .and. (c /= d)) ExchangeMultiplicity = ExchangeMultiplicity + P(b,d)*P(a,c)

  end function ExchangeMultiplicity

end subroutine EDA_TwoElectron_Energies

subroutine UnpackSymmetric(P,A,n)

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(in) :: P(:)
  real(kind=wp), intent(out) :: A(n,n)
  integer(kind=iwp) :: i,j,q

  q = 0
  do i=1,n
    do j=1,i
      q = q+1
      A(i,j) = P(q)
      A(j,i) = P(q)
    end do
  end do

end subroutine UnpackSymmetric

end module EDA_Integrals
