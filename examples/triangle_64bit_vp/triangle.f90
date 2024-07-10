program main

use dq_coeff_mod
use vector_mod

implicit none

integer :: no, ni, ninter, np
integer :: n,lwork,error
integer :: i,j,k,m,np1,stat,seed(33)
integer, allocatable :: nei(:)
real (kind=real_kind) :: pi,dist,annular_area
real (kind=real_kind), allocatable :: min_dist(:)
real (kind=real_kind), allocatable :: T(:)
real (kind=real_kind), allocatable :: neigb(:,:)
real (kind=real_kind) :: xi,yi,li
real (kind=real_kind) :: d1,d2,ext_length
real (kind=real_kind), allocatable :: A(:,:),b(:),x(:),bv(:,:)
real (kind=real_kind) :: required_epsilon2
real (kind=real_kind) :: Qn, Qn_x,Qn_y,l1,n1,bp(3)
integer, allocatable :: work(:)
integer :: required_precision
logical :: interior

character(7) :: derivative(2)
integer, allocatable :: coeff_neigh(:,:)
integer, allocatable :: coeff_num_neigh(:)
real (kind=real_kind), allocatable :: points(:,:)
integer, allocatable :: coeff_prec(:)
real (kind=real_kind), allocatable :: coeff(:,:,:)

print *,'!Enter the number of points in the interior'
read *,ninter
print *,'!Enter the number of neighbours for each point'
read *,np

pi=abs(atan2(0d0,-1d0))
li=4d0*sqrt(1d0-0.5d0**2)
annular_area=li*(li*sin(60d0*pi/180d0))/2d0*(1d0-0.5d0**2)
no=ceiling(li/sqrt(annular_area/ninter)*2d0)
ni=no/2
ni=3*ceiling(ni/3d0)
no=3*ceiling(no/3d0)

n=ni+no+ninter  ! total number of points
np1=sqrt(real(np))*2-1  ! number of neighbours that can be on boundary

lwork=3*n+3*n+2*(ni+no)+1  ! work space for constrained linear algebra solver


allocate (points(n,3))          ! points
allocate (coeff(n,np,2))    ! weighting coefficients
allocate (coeff_prec(n))      ! precision of weighting coefficients
allocate (coeff_num_neigh(n)) ! number of neighbours for each point - could
                              ! be different for each point, constant in this
                              ! example
allocate (coeff_neigh(n,np))  ! index of each neighbour for each point
allocate (T(n))               ! Temperature
allocate (nei(np),min_dist(np),neigb(np,3)) ! work space used when finding
                                          ! neighbours
allocate (A(ninter,ninter),b(ninter),x(ninter))
allocate (bv(ninter,1))

points(:,:)=0.0

! exscribed circle
!xc=2d0*cos(pi*30d0/180d0)
!yc=2d0*sin(pi*30d0/180d0)
!
!do i=0,359,5
!  print *,xc+cos(pi*dble(i)/180d0),yc+sin(pi*dble(i)/180d0),"0.0"
!  print *,xc+2d0*cos(pi*dble(i)/180d0),yc+2d0*sin(pi*dble(i)/180d0),"0.0"
!enddo


! define inner triangle
li=2d0*sqrt(1d0-0.5d0**2)
xi=1d0*cos(pi*30d0/180d0)
yi=1d0*sin(pi*30d0/180d0)
do i=1,ni/3
  points(ni/3-i+1,1)=xi+li*dble(i-1)/dble(ni/3)
  points(ni/3-i+1,2)=yi
enddo
do i=ni/3+1,2*ni/3
  points(i,1)=xi+li*dble(i-ni/3)/dble(ni/3)*cos(pi*60d0/180d0)
  points(i,2)=yi+li*dble(i-ni/3)/dble(ni/3)*sin(pi*60d0/180d0)
enddo
xi=xi+li*cos(pi*60d0/180d0)
yi=yi+li*sin(pi*60d0/180d0)
do i=2*ni/3+1,ni
  points(i,1)=xi+li*dble(i-2*ni/3)/dble(ni/3)*cos(pi*60d0/180d0)
  points(i,2)=yi-li*dble(i-2*ni/3)/dble(ni/3)*sin(pi*60d0/180d0)
enddo

! define outer triangle
xi=0d0
yi=0d0
li=4d0*sqrt(1d0-0.5d0**2)
annular_area=li*(li*sin(60d0*pi/180d0))/2d0*(1d0-0.5d0**2)
do i=no/3,1,-1
  points(ni+no/3-i+1,1)=xi+li*dble(i-1)/dble(no/3)
  points(ni+no/3-i+1,2)=yi
enddo
do i=no/3+1,2*no/3
  points(ni+i,1)=xi+li*dble(i-no/3)/dble(no/3)*cos(pi*60d0/180d0)
  points(ni+i,2)=yi+li*dble(i-no/3)/dble(no/3)*sin(pi*60d0/180d0)
enddo
xi=xi+li*cos(pi*60d0/180d0)
yi=yi+li*sin(pi*60d0/180d0)
do i=2*no/3+1,no
  points(ni+i,1)=xi+li*dble(i-2*no/3)/dble(no/3)*cos(pi*60d0/180d0)
  points(ni+i,2)=yi-li*dble(i-2*no/3)/dble(no/3)*sin(pi*60d0/180d0)
enddo

seed=1
call random_seed(put=seed)

! interior points
do i=ni+no+1,n
  do
    call random_number(points(i,1))
    call random_number(points(i,2))
    points(i,1)=points(i,1)*li
    points(i,2)=points(i,2)*li
! check if inside inner triangle
    interior=.true.
    d1=(points(i,1)-points(1,1))*(points(2,2)-points(1,2))- &
         (points(i,2)-points(1,2))*(points(2,1)-points(1,1))
    do j=2,ni
      d2=(points(i,1)-points(j,1))*(points(mod(j,ni)+1,2)-points(j,2))- &
        (points(i,2)-points(j,2))*(points(mod(j,ni)+1,1)-points(j,1))
      if (d1*d2.le.0d0) interior=.false. 
      d1=d2
    enddo
    if (interior) cycle
! check if outside outter triangle
    interior=.true.
    d1=(points(i,1)-points(ni+1,1))*(points(ni+2,2)-points(ni+1,2))- &
         (points(i,2)-points(ni+1,2))*(points(ni+2,1)-points(ni+1,1))
    do j=ni+2,ni+no
      d2=(points(i,1)-points(j,1))*(points(mod(j-ni,no)+ni+1,2)-points(j,2))- &
        (points(i,2)-points(j,2))*(points(mod(j-ni,no)+ni+1,1)-points(j,1))
      if (d1*d2.le.0d0) interior=.false. 
      d1=d2
    enddo
    if (.not.interior) cycle
    dist=huge(dist)
! ensure minimum distance between points
    do j=1,i-1
      dist=min(dist,(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2)
    enddo
    ! update min dist to be related to area of annulus and number of points
    if (dist<annular_area/2d0/(ninter+(ni+no)/2d0)) cycle 
    exit
  enddo
enddo

print *,'!points generated'

! Assign boundary conditions
T(1:ni)=0d0
T(ni+1:ni+no)=1d0

! Assign initial guess at temperature distribution in interior
do i=ni+no+1,n
  d1=huge(d1)
  do j=1,ni
    dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
    if (dist.lt.d1) then
      d1=dist
      k=j
    endif
  enddo
  d2=huge(d2)
  do j=ni+1,ni+no
    dist=(points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2
    if (dist.lt.d2) then
      d2=dist
      m=j
    endif
  enddo
  d1=1d0/d1
  d2=1d0/d2
  T(i)=(T(k)*d1+T(m)*d2)/(d1+d2)
enddo


derivative(1)="d2_dx2 "
derivative(2)="d2_dy2 "


! find neigbours and coefficients
!$OMP parallel do schedule(dynamic) default(shared), private(i,dist,j,min_dist,k,m,nei,stat,neigb,required_precision,required_epsilon2)
do i=ni+no+1,n
  coeff_num_neigh(i)=np
  min_dist(:)=huge(min_dist(:))
  do j=1,ni+no 
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
  min_dist(np1+1:np)=huge(min_dist(np1+1:np)) ! only allow np2 points from boundary
  do j=ni+no+1,n ! rest of points must be from interior
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
  coeff_neigh(i,:)=nei(:)
  do j=1,np
    neigb(j,:)=points(nei(j),:)
  enddo
! calcualate weighting coefficients
  call calc_dq_coeff(points(i,:),neigb(:,:),derivative(1:2),coeff(i,:,1:2), &
    stat,required_precision=required_precision, &
    required_epsilon2=required_epsilon2)
!  print *,i,required_precision,required_epsilon2
  if (stat.ne.0) then
    print *,'error calculating weighting coefficients'
    stop
  endif
enddo
!$OMP end parallel do
print *,'coefficients calculated'

! solve heat flow
A(:,:)=0d0
b(:)=0d0
do i=ni+no+1,n
  do j=1,coeff_num_neigh(i)
    if (coeff_neigh(i,j).gt.ni+no) then
      A(i-no-ni,coeff_neigh(i,j)-ni-no)=(coeff(i,j,1)+coeff(i,j,2))
    else
      B(i-no-ni)=B(i-no-ni)-T(coeff_neigh(i,j))*(coeff(i,j,1)+coeff(i,j,2))
    endif
  enddo
enddo

!x=Precon_Richardson(A,b)

bv(:,1)=b(:)
lwork=10*ninter
allocate(work(ninter))

call dgesv(ninter,1,A,ninter,work,Bv,ninter,error)
!call dgels('N',ninter,ninter,1,A,ninter,bv,ninter,work,lwork,error)
x=bv(:,1)

T(ni+no+1:n)=x(:)

write(*,'(3a15)')'x','y','T'
do i=1,n
  write(*,'(5e15.7)')points(i,1),points(i,2),T(i)
enddo

print *,'on outer surface'
write(*,'(a5,5a20)')'point','x','y','Qn'
! calculate heat flux on outer surface
Qn=0d0
ext_length=0d0
derivative(1)="d_dx   "
derivative(2)="d_dy   "
!$OMP parallel do schedule(dynamic) default(shared), private(i,Qn_x,Qn_y,bp,l1,n1,coeff_num_neigh,min_dist,j,dist,k,m,nei,coeff_neigh,neigb,coeff)
do i=ni+1,ni+no
  Qn_x=0d0
  Qn_y=0d0
  ! define points midway between boundary nodes
  bp(:)=(points(i,:)+points(mod(i-ni,no)+ni+1,:))/2d0
  l1=sqrt((points(mod(i-ni,no)+ni+1,1)-points(i,1))**2+ &
    (points(mod(i-ni,no)+ni+1,2)-points(i,2))**2)
  n1=atan2(points(mod(i-ni,no)+ni+1,2)-points(i,2), &
    points(mod(i-ni,no)+ni+1,1)-points(i,1))+pi/2d0
  ext_length=ext_length+l1
  ! find neighbours 
  coeff_num_neigh(i)=np
  min_dist(:)=huge(min_dist(:))
  do j=1,ni+no 
    dist=(bp(1)-points(j,1))**2+(bp(2)-points(j,2))**2
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
  min_dist(3:np)=huge(min_dist(3:np)) ! only allow 2 points from boundary
  do j=ni+no+1,n ! rest of points must be from interior
    dist=(bp(1)-points(j,1))**2+(bp(2)-points(j,2))**2
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
  coeff_neigh(i,:)=nei(:)
  do j=1,np
    neigb(j,:)=points(nei(j),:)
  enddo
! calcualate weighting coefficients
  call calc_dq_coeff(points(i,:),neigb(:,:),derivative(1:2),coeff(i,:,1:2), &
    stat) 
  if (stat.ne.0) then
    print *,'error calculating weighting coefficients'
    stop
  endif
 do j=1,coeff_num_neigh(i)
      Qn_x=Qn_x+T(coeff_neigh(i,j))*coeff(i,j,1)
      Qn_y=Qn_y+T(coeff_neigh(i,j))*coeff(i,j,2)
  enddo
  Qn=Qn+l1*(Qn_x*cos(n1)+Qn_y*sin(n1))
  write(*,'(i5,9e20.8)')i,bp(1),bp(2),Qn_x,Qn_y,Qn_x*cos(n1),Qn_y*sin(n1),l1,l1*(Qn_x*cos(n1)+Qn_y*sin(n1)),Qn
enddo
!$OMP end parallel do
print *,'S',Qn
  
deallocate (points,coeff,coeff_prec,coeff_num_neigh,coeff_neigh)
deallocate (T)
deallocate (A,b,x)
deallocate (nei,min_dist,neigb,work,bv)

contains

function iterative_inverse(A,tol,v_init)

real*8, intent(in) :: A(:,:)
real*8 :: iterative_inverse(size(A,1),size(A,1))
real*8, intent(in), optional :: tol
real*8, intent(in), optional :: v_init(:,:)
real*8 :: norm_1,norm_inf,temp,diag,off_diag,eps
real*8 :: phi(size(A,1),size(A,1)),eta(size(A,1),size(A,1))
real*8 :: kn(size(A,1),size(A,1))
real*8 :: work0(size(A,1),size(A,1)),work1(size(A,1),size(A,1))
integer :: i,j,n,iter

n=size(A,1)
if (present(tol)) then
  eps=tol
else
  eps=epsilon(eps)*n
endif

if (present(v_init)) then
  iterative_inverse=v_init
else
  norm_1=0d0
  do i=1,n
    temp=sum(abs(A(:,i)))
    norm_1=max(norm_1,temp)
  enddo
  norm_inf=0d0
  do i=1,n
    temp=sum(abs(A(i,:)))
    norm_inf=max(norm_inf,temp)
  enddo
  iterative_inverse=transpose(A)/(norm_1*norm_inf)
endif

iter=0
do
  iter=iter+1
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      phi(i,j)=dot_product(A(i,:),iterative_inverse(:,j))
    enddo
  enddo
!$OMP end parallel do
  diag=0d0
  do i=1,n
    diag=max(diag,abs(1d0-phi(i,i)))
  enddo
  if (diag.lt.eps) then
    off_diag=0d0
    do i=1,n
      temp=sum(abs(phi(i,:)))-abs(phi(i,i))
      off_diag=max(off_diag,temp)
    enddo
    if (off_diag.lt.eps) exit
  endif
  print *,'!iter inverse',iter,diag,off_diag,eps
  work0=phi
  do i=1,n
    work0(i,i)=work0(i,i)-8d0
  enddo
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      work1(i,j)=dot_product(phi(i,:),work0(:,j))
    enddo
  enddo
!$OMP end parallel do
  do i=1,n
    work1(i,i)=work1(i,i)+22d0
  enddo
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      work0(i,j)=dot_product(phi(i,:),work1(:,j))
    enddo
  enddo
!$OMP end parallel do
  do i=1,n
    work0(i,i)=work0(i,i)-28d0
  enddo
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      eta(i,j)=dot_product(phi(i,:),work0(:,j))
    enddo
  enddo
!$OMP end parallel do
  do i=1,n
    eta(i,i)=eta(i,i)+17d0
  enddo
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      kn(i,j)=dot_product(phi(i,:),eta(:,j))
    enddo
  enddo
!$OMP end parallel do
  work0=kn
  do i=1,n
    work0(i,i)=work0(i,i)-12d0
  enddo
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      work1(i,j)=dot_product(kn(i,:),work0(:,j))
    enddo
  enddo
!$OMP end parallel do
  do i=1,n
    work1(i,i)=work1(i,i)+48d0
  enddo
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      work0(i,j)=dot_product(eta(i,:),work1(:,j))
    enddo
  enddo
!$OMP end parallel do
!$OMP parallel do default(shared), private(i,j)
  do i=1,n
    do j=1,n
      work1(i,j)=dot_product(iterative_inverse(i,:),work0(:,j))/64d0
    enddo
  enddo
!$OMP end parallel do
  iterative_inverse=work1
enddo

end function

function Precon_Richardson(A,b,tol)

real*8, intent(in) :: A(:,:),b(:)
real*8 :: Precon_Richardson(size(b,1))
real*8,intent(in), optional :: tol
real*8 :: ainv(size(A,1),size(A,1)),r(size(b,1)),eps,inv_eps
real*8 :: x(size(b,1)),rnorm,last_rnorm,last_x(size(b,1))
real*8 :: r_last(size(b,1))
integer :: i,n,iter

if (present(tol)) then
  eps=tol**2
else
  eps=(epsilon(eps)*n)**2
endif

n=size(b,1)

inv_eps=1d-2

ainv=iterative_inverse(A,inv_eps)
print *,'!approx inverse calculated'
last_x=0d0
r(:)=b(:)
last_rnorm=dot_product(r,r)
iter=0
do
  iter=iter+1
  print *,'!Precon_Richardson iteration ',iter
!$OMP parallel do default(shared), private(i)
  do i=1,n
    x(i)=last_x(i)+dot_product(ainv(i,:),r)
  enddo
!$OMP end parallel do
!$OMP parallel do default(shared), private(i)
  do i=1,n
    r(i)=b(i)-dot_product(A(i,:),x)
  enddo
!$OMP end parallel do
  rnorm=dot_product(r,r)
  if (rnorm.le.eps) exit
  if (rnorm.gt.last_rnorm/2d0) then
    inv_eps=inv_eps/1d2
    if (inv_eps.lt.1d-10) exit
    ainv=iterative_inverse(A,inv_eps,ainv)
    r=r_last
  else
    last_rnorm=rnorm
    last_x=x
    r_last=r
  endif
enddo

Precon_Richardson(:)=x(:)

end function

end program
