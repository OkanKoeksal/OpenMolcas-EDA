!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module EDA_Parallel

use Definitions, only: iwp, u6
use Para_Info, only: nProcs, MyRank, Is_Real_Par

implicit none
private

public :: EDA_Is_Parallel, EDA_Is_Root, EDA_Any, EDA_Sync_Error

contains

logical function EDA_Is_Parallel()

  EDA_Is_Parallel = (nProcs > 1) .and. Is_Real_Par()

end function EDA_Is_Parallel

logical function EDA_Is_Root()

  EDA_Is_Root = (.not. EDA_Is_Parallel()) .or. (MyRank == 0)

end function EDA_Is_Root

logical function EDA_Any(Local)

  logical, intent(in) :: Local
  integer(kind=iwp) :: Flag

  Flag = 0
  if (Local) Flag = 1
  if (EDA_Is_Parallel()) call gaIgOP_SCAL(Flag,'max')
  EDA_Any = (Flag /= 0)

end function EDA_Any

subroutine EDA_Sync_Error(Code,Message)

  integer(kind=iwp), intent(inout) :: Code
  character(len=*), intent(inout) :: Message
  logical :: FailedHere

  FailedHere = (Code /= 0)
  if (EDA_Any(FailedHere)) then
    if (FailedHere) write(u6,'(A,I6,2A)') ' CMOEDA rank ',MyRank,': ',trim(Message)
    Code = 1
    if (.not. FailedHere) Message = 'CMOEDA failed on another MPI rank: inspect rank diagnostics.'
  end if

end subroutine EDA_Sync_Error

end module EDA_Parallel
