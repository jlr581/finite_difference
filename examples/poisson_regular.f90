program poisson

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

use dq_coeff_mod
use vector_mod

implicit none

integer, parameter :: nx=10,nit=81,nei=16,ntot=nx*4+nit

real*8 :: points(ntot,3),coeff(nei,2),phi(ntot),rnd(2)
real*8 :: dist(nei),d,pi,n_points(nei,3),err,dis
real*8 :: local_err(3),local_rms(nx*4:ntot),dx,local_phi(nei)
integer :: stat,i,j,neighbour(nei),k,seed(33)

pi=abs(atan2(0d0,-1d0))

! points on boundary
points(:,:)=0d0
do i=1,nx
  points(i,1)=2d0*real(i-1,real_kind)/real(nx,real_kind)-1d0
  points(i,2)=1d0
  points(i+nx,1)=2d0*real(i,real_kind)/real(nx,real_kind)-1d0
  points(i+nx,2)=-1d0
  points(i+2*nx,1)=1d0
  points(i+2*nx,2)=2*real(i,real_kind)/real(nx,real_kind)-1d0
  points(i+3*nx,1)=-1d0
  points(i+3*nx,2)=2*real(i-1,real_kind)/real(nx,real_kind)-1d0
enddo

dx=abs(points(2,1)-points(1,1))

! points in interior
do j=1,nx-1
  do i=1,nx-1
    points(nx*4+(j-1)*(nx-1)+i,1)=2d0*real(i,real_kind)/real(nx,real_kind)-1d0
    points(nx*4+(j-1)*(nx-1)+i,2)=2d0*real(j,real_kind)/real(nx,real_kind)-1d0
  enddo
enddo

! calculate correct solution
do i=1,ntot
  phi(i)=cos(points(i,1)*pi/2d0)*cos(points(i,2)*pi/2d0)
enddo

! calculate rms error
print *,'x                         y                       local error'

!$OMP parallel do default(shared), private(i,dist,d,j,k,neighbour,coeff,n_points,err,local_err,stat,local_phi)
do i=4*nx+1,ntot
  ! find nearest neighbours
  dist(:)=huge(dist(1))
  do k=1,ntot
    d=(points(i,1)-points(k,1))**2+(points(i,2)-points(k,2))**2
    if (d.lt.dist(nei)) then
      do j=1,nei
        if (d.lt.dist(j)) then
          dist(j+1:nei)=dist(j:nei-1)
          neighbour(j+1:nei)=neighbour(j:nei-1)
          dist(j)=d
          neighbour(j)=k
          exit
        endif
      enddo
    endif
  enddo

  do j=1,nei
    n_points(j,:)=points(neighbour(j),:)
  enddo

  ! calculate weighting coefficients
  call calc_dq_coeff(points(i,:),n_points(:,:),(/"d2_dx2","d2_dy2"/), &
    coeff,stat)

  ! calculate local error
  local_err(1)=pi*pi/2d0*cos(points(i,1)*pi/2)*cos(points(i,2)*pi/2)
  do j=1,nei
    local_phi(j)=phi(neighbour(j))
  enddo
  local_err(2)=kfold_dot_product(coeff(:,1),local_phi)
  local_err(3)=kfold_dot_product(coeff(:,2),local_phi)
  err=kfold_sum(local_err)
  ! add to global error
  local_rms(i)=err**2
  print *,points(i,1:2),err
enddo
!$OMP end parallel do

print *,"global rms error",sqrt(kfold_sum(local_rms)/(ntot-4*nx))

end program
