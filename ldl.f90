module ldl_mod

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
use hyper_dual_mod

implicit none

contains

subroutine ldl(a)

type (hyper_dual), intent(inout) :: a(:,:)

integer :: i,j,n
type (hyper_dual) :: v(size(a,1))
type (hyper_dual) :: tmp

n=size(a,1)

v=hyper_dual(vp_zero,vp_zero,vp_zero,vp_zero)
do j=2,n
  do i=1,j-1
    v(i)=a(j,i)*a(i,i)
  enddo
  v(j)=a(j,j)-dot_product_hd(a(j,1:j-1),v(1:j-1))
  a(j,j)=v(j)
  tmp=hd_one/v(j)
  do i=j+1,n
    a(i,j)=(a(i,j)-dot_product_hd(a(i,1:j-1),v(1:j-1)))*tmp
  enddo
enddo


end subroutine

end module
