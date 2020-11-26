module var_prec_extra_mod

!    dq_coeff is a program to calculate the finite difference coefficients
!    for an arbitrary set of points using differential quadrature
!    Copyright (C) 2020 Jason Roberts
!
!    This file is part of dq_ceoff.
!
!    dq_coeff is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    dq_coeff is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with dq_ceoff.  If not, see <https://www.gnu.org/licenses/>.

use data_types_mod

implicit none

contains

subroutine split_half(a,a_high,a_low,n)

integer :: n
integer (kind=int_kind), intent(in) :: a(1:n)
integer (kind=int_kind), intent(out) :: a_high(1:n)
integer (kind=int_kind), intent(out) :: a_low(1:n)
integer :: i

do i=1,n
  a_high(i)=shiftr(a(i),var_prec_half_eff_size)
enddo
do i=1,n
  a_low(i)=a(i)-a_high(i)*var_prec_half_base
enddo

end subroutine

subroutine split(a,a_high,a_low,n)

integer :: n
integer (kind=int_kind), intent(in) :: a(1:n)
integer (kind=int_kind), intent(out) :: a_high(1:n)
integer (kind=int_kind), intent(out) :: a_low(1:n)
integer :: i

do i=1,n
  a_high(i)=shiftr(a(i),var_prec_eff_size)
enddo
do i=1,n
  a_low(i)=a(i)-a_high(i)*var_prec_base
enddo

end subroutine

end module

