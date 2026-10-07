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

module EDA_Data

use Definitions, only: iwp

implicit none
private

integer, parameter, public :: EDA_MaxFrag = 10

integer(kind=iwp), public :: EDA_nFrag = 0
integer(kind=iwp), public :: EDA_nAtom(EDA_MaxFrag) = 0
integer(kind=iwp), public :: EDA_Charge(EDA_MaxFrag) = 0
integer(kind=iwp), public :: EDA_Mult(EDA_MaxFrag) = 1
character(len=8), public :: EDA_RefMethod = 'AUTO    '
character(len=8), public :: EDA_FragBasis = 'OWN     '
logical, public :: EDA_Debug = .false.

public :: EDA_Reset_Input

contains

subroutine EDA_Reset_Input()

  EDA_nFrag = 0
  EDA_nAtom = 0
  EDA_Charge = 0
  EDA_Mult = 1
  EDA_RefMethod = 'AUTO    '
  EDA_FragBasis = 'OWN     '
  EDA_Debug = .false.

end subroutine EDA_Reset_Input

end module EDA_Data
