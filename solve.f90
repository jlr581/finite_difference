module solve_mod

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
use var_prec_mod
use hyper_dual_mod
use ldl_mod

implicit none

contains

subroutine ldl_solve(a,b,x,h)

type (hyper_dual), intent(inout) :: a(:,:)
type (hyper_dual), intent(in) :: b(:,:)
real (kind=real_kind), intent(out) :: x(:,:,0:)
type (hyper_dual), intent(in) :: h

integer :: i,j,n,m,k
type (hyper_dual) :: y(size(a,1),size(b,2))
type (var_prec) :: x_vp
type (hyper_dual) :: inv_a(size(a,1)), inv_single
type (hyper_dual) :: best_y(size(a,1),size(b,2))
type (hyper_dual) :: orig_a(size(a,1),size(a,2)),orig_b(size(b,1),size(b,2))
type (hyper_dual) :: part_a(size(a,1),size(a,2))
type (hyper_dual) :: resid(size(b,1),size(b,2))
real (kind=real_kind) :: norm(1:4,size(b,2)) 
real (kind=real_kind) :: t_norm,last_tnorm,tmp,best_tnorm
real (kind=real_kind) :: r_r,r_e1,r_e1e2,r_e1e2e3
real (kind=real_kind) :: last_x(size(x,1),size(x,2),0:3),diff

n=size(a,1)
m=size(b,2)

do j=1,n
  do i=1,n
    orig_a(i,j)=a(i,j)
  enddo
enddo
do j=1,m
  do i=1,n
    orig_b(i,j)=b(i,j)
  enddo
enddo

call ldl(a) ! Cholesky decomposition of a
do j=1,n
  do i=1,n
    part_a(i,j)=a(i,j)
  enddo
enddo

! forward substution solution of Ly=b
y=b
inv_single=hd_one/a(1,1)
do j=1,m
  y(1,j)=y(1,j)*inv_single
  do i=2,n
    y(i,j)=y(i,j)-dot_product_hd(a(i,1:i-1),y(1:i-1,j))
  enddo
enddo

! back substution solution of DL(T)x=y
do i=1,n-1
  do j=i+1,n
    a(i,j)=a(i,i)*a(j,i)
  enddo
enddo

do i=1,n
  inv_a(i)=hd_one/a(i,i)
enddo

do j=1,m
  y(n,j)=y(n,j)*inv_a(n)
  do i=n-1,1,-1
    y(i,j)=(y(i,j)-dot_product_hd(a(i,i+1:n),y(i+1:n,j)))*inv_a(i)
  enddo
enddo

! iterative improvement
last_tnorm=huge(last_tnorm)
best_tnorm=huge(best_tnorm)
best_y=y
last_x=0d0
do k=1,n

  var_prec_working=var_prec_working+1
  norm(:,:)=0d0
  do j=1,m
    do i=1,n
      resid(i,j)=orig_b(i,j)-dot_product_hd(orig_a(:,i),y(:,j))
      r_r=resid(i,j)%r
      r_e1=resid(i,j)%e1
      r_e1e2=resid(i,j)%e1e2
      r_e1e2e3=resid(i,j)%e1e2e3
      norm(1,j)=norm(1,j)+abs(r_r)
      norm(2,j)=norm(2,j)+abs(r_e1)
      norm(3,j)=norm(3,j)+abs(r_e1e2)
      norm(4,j)=norm(4,j)+abs(r_e1e2e3)
    enddo
  enddo
  var_prec_working=var_prec_working-1
  t_norm=0d0
  do j=1,m
    t_norm=t_norm+sum(norm(:,j))
  enddo
  if ((t_norm.gt.2d0*last_tnorm).and.(k.gt.2)) then
    do j=1,m
      do i=1,n
        x_vp=best_y(i,j)%r
        x(i,j,0)=x_vp
        x_vp=x_vp-h%r*best_y(i,j)%e1
        x(i,j,1)=x_vp
        x_vp=2*x_vp+h%r*h%r*best_y(i,j)%e1e2
        x(i,j,2)=x_vp
        x(i,j,2)=x(i,j,2)/2d0
        x_vp=6*x_vp-h%r*h%r*h%r*best_y(i,j)%e1e2e3
        x(i,j,3)=x_vp
        x(i,j,3)=x(i,j,3)/12d0
      enddo
    enddo
    exit
  endif
  if (t_norm.lt.best_tnorm) then
    best_y=y
    best_tnorm=t_norm
  endif
  last_tnorm=t_norm

! forward substution solution of Ly=b
  do j=1,m
    resid(1,j)=resid(1,j)*inv_single
    do i=2,n
      resid(i,j)=resid(i,j)-dot_product_hd(part_a(i,1:i-1),resid(1:i-1,j))
    enddo
  enddo

  do j=1,m
    resid(n,j)=resid(n,j)*inv_a(n)
    do i=n-1,1,-1
      resid(i,j)=(resid(i,j)-dot_product_hd(a(i,i+1:n),resid(i+1:n,j)))*inv_a(i)
    enddo
  enddo

  do j=1,m
    do i=1,n
      y(i,j)=y(i,j)+resid(i,j)
    enddo
  enddo

  do j=1,m
    do i=1,n
      x_vp=y(i,j)%r
      x(i,j,0)=x_vp
      x_vp=x_vp-h%r*y(i,j)%e1
      x(i,j,1)=x_vp
      x_vp=2*x_vp+h%r*h%r*y(i,j)%e1e2
      x(i,j,2)=x_vp
      x(i,j,2)=x(i,j,2)/2d0
      x_vp=3*x_vp-h%r*h%r*h%r*y(i,j)%e1e2e3
      x(i,j,3)=x_vp
      x(i,j,3)=x(i,j,3)/6d0
    enddo
  enddo

  diff=maxval(abs(x-last_x))
  last_x=x
  if (diff.eq.0) exit
enddo

end subroutine

end module
