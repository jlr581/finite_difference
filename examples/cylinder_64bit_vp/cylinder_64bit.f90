program main

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

integer :: nc, nb, ninter, np
integer :: n,info,lwork
integer :: i,j,k,m,np1,np2,stat,seed(33)
integer, allocatable :: nei(:)
real (kind=real_kind) :: pi,rad,dist,nu
real (kind=real_kind), allocatable :: min_dist(:)
real (kind=real_kind), allocatable :: u(:),v(:),p(:)
real (kind=real_kind), allocatable :: neigb(:,:),local_var(:)
real (kind=real_kind), allocatable :: a(:,:),b(:),d(:,:),e(:),x(:),work(:)
real (kind=real_kind) :: dr,dt,theta,txx,txy,theta1
real (kind=real_kind), allocatable :: wsum(:),csum(:,:)
integer :: i1,i2

character(7) :: derivative(5)
integer, allocatable :: coeff_neigh(:,:)
integer, allocatable :: coeff_num_neigh(:)
real (kind=real_kind), allocatable :: points(:,:)
integer, allocatable :: coeff_prec(:)
real (kind=real_kind), allocatable :: coeff(:,:,:)

print *,'Enter the number of points on the cylinder'
read *,nc
print *,'Enter the number of points on the boundary (divisible by 6)'
read *,nb
nb=6*int(nb/6)
print *,'Enter the number of points in the interior'
read *,ninter
print *,'Enter the number of neighbours for each point'
read *,np

n=nc+nb+ninter  ! total number of points
np1=sqrt(real(np))*2-1  ! number of neighbours that can be on boundary
np2=sqrt(real(np))      ! number of neighbours that can be on cylinder
pi=abs(atan2(0d0,-1d0))
nu=1d0

lwork=3*n+3*n+2*(nb+nc)+1  ! work space for constrained linear algebra solver


allocate (points(n,3))          ! points
allocate (coeff(n,np,5))    ! weighting coefficients
allocate (coeff_prec(n))      ! precision of weighting coefficients
allocate (coeff_num_neigh(n)) ! number of neighbours for each point - could
                              ! be different for each point, constant in this
                              ! example
allocate (coeff_neigh(n,np))  ! index of each neighbour for each point
allocate (u(n),v(n),p(n))     ! velocity components and pressure
allocate (a(3*n,3*n),b(3*n),x(3*n)) ! matrices/vectors for coupled solution
                                    ! to Stokes equation
allocate (d(2*(nc+nb)+1,3*n),e(2*(nc+nb)+1)) ! constrains for boundary conditions
allocate (work(lwork))              ! work space for linear algebra routine
allocate (nei(np),min_dist(np),neigb(np,3)) ! work space used when finding
                                          ! neighbours
allocate (wsum(np),csum(nc,3),local_var(np))

points(:,:)=0.0

! define cylinder
do i=1,nc/4
  theta=(dble(i)/(nc/4+0.5))**2*pi/2d0
  points(i,1)=0.5d0*cos(theta)
  points(i,2)=0.5d0*sin(theta)
  points(nc-i+1,1)=points(i,1)
  points(nc-i+1,2)=-points(i,2)
  points(nc/2-i+1,1)=-points(i,1)
  points(nc/2-i+1,2)=points(i,2)
  points(nc/2+i,1)=-points(i,1)
  points(nc/2+i,2)=-points(i,2)
enddo

! define boundary
do i=1,nb/3
  ! bottom
  points(nc+i,1)=real(i-1,real_kind)/(nb/3-1)*10d0-5d0
  points(nc+i,2)=-2.5d0
  ! top
  points(nc+i+nb/3,1)=real(i-1,real_kind)/(nb/3-1)*10d0-5d0
  points(nc+i+nb/3,2)=2.5d0
enddo

do i=1,nb/6
  ! inflow
  points(nc+i+2*nb/3,1)=-5d0
  points(nc+i+2*nb/3,2)=real(i,real_kind)/(nb/6+1)*5d0-2.5d0
  ! outflow
  points(nc+i+2*nb/3+nb/6,1)=5d0
  points(nc+i+2*nb/3+nb/6,2)=real(i,real_kind)/(nb/6+1)*5d0-2.5d0
enddo

seed=1
call random_seed(put=seed)

! interior points
do i=nc+nb+1,n
  do
    call random_number(points(i,1))
    call random_number(points(i,2))
    points(i,1)=(points(i,1)-0.5d0)*10d0
    points(i,2)=(points(i,2)-0.5d0)*5d0
    rad=points(i,1)**2+points(i,2)**2
    if (rad<=0.25d0) cycle ! exclude points inside cylinder
    dist=huge(dist)
    ! ensure minimum distance between points
    do j=1,i-1
      dist=min(dist,(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2)
    enddo
    if (dist<50d0/ninter/4.1*sqrt(rad)) cycle 
    exit
  enddo
enddo

! set boundary conditions
u(:)=1.0
v(:)=0.0
u(1:nc)=0.0

! set which derivatives to calculate
derivative(1)="d_dx   "
derivative(2)="d_dy   "
derivative(3)="d2_dx2 "
derivative(4)="d2_dy2 "
derivative(5)="d2_dxdy"


! find neigbours and weighting coefficients
!$OMP parallel do schedule(dynamic) default(shared), private(i,dist,j,min_dist,k,m,nei,stat,neigb,required_epsilon2)
do i=1,n
  coeff_num_neigh(i)=np
  if (i<=nc) then ! points on cylinder
    min_dist(:)=huge(min_dist(:))
    do j=1,nc 
      dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
      if (dist<min_dist(np)) then
        do k=1,np
          if (dist<min_dist(k)) then
            do m=np,k+1,-1
              min_dist(m)=min_dist(m-1)
              nei(m)=nei(m-1)
            enddo
            min_dist(k)=dist
            nei(k)=j
            exit
          endif
        enddo
      endif
    enddo
    min_dist(np2+1:np)=huge(min_dist(np2+1:np)) ! only allow np2 points on cylinder
    do j=nc+1,n ! rest of points must be from boundary or interior
      dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
      if (dist<min_dist(np)) then
        do k=1,np
          if (dist<min_dist(k)) then
            do m=np,k+1,-1
              min_dist(m)=min_dist(m-1)
              nei(m)=nei(m-1)
            enddo
            min_dist(k)=dist
            nei(k)=j
            exit
          endif
        enddo
      endif
    enddo
  elseif (i<=nc+nb) then ! points on boundary
    min_dist(:)=huge(min_dist(:))
    do j=1,nc+nb
      dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
      if (dist<min_dist(np)) then
        do k=1,np
          if (dist<min_dist(k)) then
            do m=np,k+1,-1
              min_dist(m)=min_dist(m-1)
              nei(m)=nei(m-1)
            enddo
            min_dist(k)=dist
            nei(k)=j
            exit
          endif
        enddo
      endif
    enddo
    min_dist(np1+1:np)=huge(min_dist(np1+1:np)) ! only allow np1 points on boundary
    do j=nc+nb+1,n ! rest from interior
      dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
      if (dist<min_dist(np)) then
        do k=1,np
          if (dist<min_dist(k)) then
            do m=np,k+1,-1
              min_dist(m)=min_dist(m-1)
              nei(m)=nei(m-1)
            enddo
            min_dist(k)=dist
            nei(k)=j
            exit
          endif
        enddo
      endif
    enddo
  else ! interior points
    min_dist(:)=huge(min_dist(:))
    do j=1,nc 
      dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
      if (dist<min_dist(np)) then
        do k=1,np
          if (dist<min_dist(k)) then
            do m=np,k+1,-1
              min_dist(m)=min_dist(m-1)
              nei(m)=nei(m-1)
            enddo
            min_dist(k)=dist
            nei(k)=j
            exit
          endif
        enddo
      endif
    enddo
    min_dist(np2+1:np)=huge(min_dist(np2+1:np)) ! only allow np2 points from cylinder
    do j=nc+1,n ! rest from boundary or interior
      dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
      if (dist<min_dist(np)) then
        do k=1,np
          if (dist<min_dist(k)) then
            do m=np,k+1,-1
              min_dist(m)=min_dist(m-1)
              nei(m)=nei(m-1)
            enddo
            min_dist(k)=dist
            nei(k)=j
            exit
          endif
        enddo
      endif
    enddo
  endif
  coeff_neigh(i,:)=nei(:)
  do j=1,np
    neigb(j,:)=points(nei(j),:)
  enddo
  ! calcualate weighting coefficients
  call calc_dq_coeff(points(i,:),neigb(:,:),derivative,coeff(i,:,:),stat) 
  if (stat/=0) then
    print *,'error calculating weighting coefficients'
    stop
  endif
enddo
!$OMP end parallel do
print *,'coefficients calculated'

! solve fluid flow
a(:,:)=0
b(:)=0
do i=1,n
  do j=1,coeff_num_neigh(i)
    ! x momentum
    a(i,coeff_neigh(i,j))=-nu*(coeff(i,j,3)+coeff(i,j,4)) 
    a(i,2*n+coeff_neigh(i,j))=coeff(i,j,1)
    ! y momentum
    a(n+i,n+coeff_neigh(i,j))=-nu*(coeff(i,j,3)+coeff(i,j,4)) 
    a(n+i,2*n+coeff_neigh(i,j))=coeff(i,j,2)
  enddo
  do j=1,coeff_num_neigh(i)
    ! continuity
    a(2*n+i,coeff_neigh(i,j))=coeff(i,j,1)
    a(2*n+i,n+coeff_neigh(i,j))=coeff(i,j,2)
  enddo
enddo
! boundary condition constraints
! velocity boundary conditions
do i=1,nb+nc
  d(i,i)=1
  e(i)=u(i)
  d(nb+nc+i,n+i)=1
  e(nb+nc+i)=v(i)
enddo
! pressure boundary condition
d(2*(nb+nc)+1,2*n+1)=1
e(2*(nb+nc)+1)=0
! LAPACK constained least squares solution (in this case not over determined)
call dgglse(3*n,3*n,2*(nc+nb)+1,a,3*n,d,2*(nc+nb)+1,b,e,x,work,lwork,info)
if (info/=0) then
  print *,'error in soln'
  stop
endif
! extract velocity and pressure from solution vector
do i=1,n
  u(i)=x(i)
  v(i)=x(n+i)
  p(i)=x(2*n+i)
enddo

write(*,'(5a15)')'x','y','u','v','p'
do i=1,n
  write(*,'(5e15.7)')points(i,1),points(i,2),u(i),v(i),p(i)
enddo

print *,'on surface of cylinder'
write(*,'(a5,5a20)')'point','x','y','p','tau_xx','tau_xy'
! calculate drag on cylinder
dr=0
csum=0d0
do i=1,nc
  do j=1,np
    local_var(j)=u(coeff_neigh(i,j))
  enddo
  csum(i,1)=kfold_dot_product(coeff(i,1:np,1),local_var(1:np))
  csum(i,2)=kfold_dot_product(coeff(i,1:np,2),local_var(1:np))
  do j=1,np
    local_var(j)=v(coeff_neigh(i,j))
  enddo
  csum(i,3)=kfold_dot_product(coeff(i,1:np,1),local_var(1:np))
   write(*,'(i5,5e20.8)')i,points(i,1),points(i,2),p(i),nu*csum(i,1),nu*(csum(i,2)+csum(i,3))
  i1=i+1
  if (i1.gt.nc) i1=1
  theta=atan2(points(i1,2),points(i1,1))
  if (theta.lt.0) theta=theta+2d0*pi
  i2=i-1
  if (i2.lt.1) i2=nc
  theta1=atan2(points(i2,2),points(i2,1))
  if (theta1.lt.0) theta1=theta1+2d0*pi
  if ((theta-theta1).lt.0d0) theta1=theta1-2d0*pi
  dt=(theta-theta1)/2d0
   theta=atan2(points(i,2),points(i,1))
   dr=dr+(-p(i)*cos(theta)+2d0*nu*csum(i,1)*cos(theta)+nu*(csum(i,2)+csum(i,3))*sin(theta))*dt*0.5d0
enddo
write(*,'(a,f7.3,a)')'total drag per unit length of cylinder',dr,' N'
  
deallocate (points,coeff,coeff_prec,coeff_num_neigh,coeff_neigh)
deallocate (u,v,p)
deallocate (a,b,x,d,e,work,local_var)
deallocate (nei,min_dist,neigb,wsum,csum)


end program
