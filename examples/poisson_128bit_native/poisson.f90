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

integer, parameter :: nx=30,nit=841,nei=36,ntot=nx*4+nit

real (kind=real_kind):: points(ntot,3),coeff(nei,2),phi(ntot),rnd(2)
real (kind=real_kind):: dist(nei),d,pi,n_points(nei,3),rms,err,dis,e2
real (kind=real_kind):: local_err(3),local_rms(nx*4:ntot),dx,local_phi(nei)
integer :: stat,i,j,neighbour(nei),k,seed(33)

pi=abs(atan2(0d0,-1d0))

print *,'!Enter epsilon2'
read *,e2

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

seed=1
call random_seed(put=seed)

dx=abs(points(2,1)-points(1,1))

do i=1,nit
  do
    call random_number(rnd)
    points(i+nx*4,1:2)=(rnd-0.5d0)*2d0
    dis=huge(dis)
    do j=1,i-1+nx*4
      dis=min(dis,(points(i+nx*4,1)-points(j,1))**2+(points(i+nx*4,2)-points(j,2))**2)
    enddo
    if (dis<4d0/nit/2) cycle
    dx=min(dx,sqrt(dis))
    exit
  enddo
enddo

do i=1,ntot
  phi(i)=cos(points(i,1)*pi/2d0)*cos(points(i,2)*pi/2d0)
enddo

rms=0d0
!$OMP parallel do default(shared), private(i,dist,d,j,k,neighbour,coeff,n_points,err,local_err,stat,local_phi)
do i=4*nx+1,ntot
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

 call calc_dq_coeff(points(i,:),n_points(:,:),(/"d2_dx2","d2_dy2"/),coeff,stat,set_epsilon2=e2,taylor_series=0)
  
  local_err(1)=pi*pi/2d0*cos(points(i,1)*pi/2)*cos(points(i,2)*pi/2)
  do j=1,nei
    local_phi(j)=phi(neighbour(j))
  enddo
  local_err(2)=kfold_dot_product(coeff(:,1),local_phi)
  local_err(3)=kfold_dot_product(coeff(:,2),local_phi)
  err=kfold_sum(local_err)
  local_rms(i)=err**2
  rms=rms+err**2
  if (rms.ne.rms) stop
enddo
!$OMP end parallel do

rms=kfold_sum(local_rms)
print *,"rms",sqrt(rms/(ntot-4*nx)),dx

end program
