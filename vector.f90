module vector_mod

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

integer, parameter :: kfold_working=2

contains

function twosum(a,b)
real (kind=real_kind) :: twosum(2)
real (kind=real_kind), intent(in) :: a,b

real (kind=real_kind) :: z

twosum(2)=a+b
z=twosum(2)-a
twosum(1)=((a-(twosum(2)-z))+(b-z))

end function

function split(a)
real (kind=real_kind) :: split(2)
real (kind=real_kind), intent(in) :: a

real (kind=real_kind) :: c

c=(2**27+1)*a
split(1)=c-(c-a)
split(2)=a-split(1)

end function

function twoproduct(a,b)
real (kind=real_kind) :: twoproduct(2)
real (kind=real_kind), intent(in) :: a,b

real (kind=real_kind) :: as(2),bs(2)

twoproduct(1)=a*b
as=split(a)
bs=split(b)
twoproduct(2)=as(2)*bs(2)-(((twoproduct(1)-as(1)*bs(1))-as(2)*bs(1))-as(1)*bs(2))

end function

function vectorsum(p)
real (kind=real_kind), intent(in) :: p(:)
real (kind=real_kind) :: vectorsum(1:size(p,1))

real (kind=real_kind) :: q(1:size(p,1))

integer :: i,n

n=size(p,1)
q=p

do i=2,n
  q(i-1:i)=twosum(p(i),q(i-1))
enddo

vectorsum=q

end function

function sumkk(p)
real (kind=real_kind), intent(in) :: p(:)
real (kind=real_kind) :: sumkk(1:kfold_working)

integer :: i,n
real (kind=real_kind) :: q(size(p,1))

n=size(p,1)

sumkk(1:kfold_working)=0d0
q=p
do i=0,kfold_working-2
  q(1:n-i)=vectorsum(q(1:n-i))
  sumkk(i+1)=q(n-i)
enddo
sumkk(kfold_working)=sum(q(1:n-kfold_working+1))

end function

function kfold_dot_product(x,y)

real (kind=real_kind) :: kfold_dot_product
real (kind=real_kind), intent(in) :: x(:),y(:)
real (kind=real_kind) :: prod(size(x,1)),sm(2)

integer :: i,n

n=size(x,1)

do i=1,n
  prod(i)=x(i)*y(i)
enddo

sm=sumkk(prod)

kfold_dot_product=sm(1)+sm(2)

end function

function kfold_sum(x)

real (kind=real_kind) :: kfold_sum
real (kind=real_kind), intent(in) :: x(:)
real (kind=real_kind) :: sm(2)


sm=sumkk(x)

kfold_sum=sm(1)+sm(2)

end function

end module
