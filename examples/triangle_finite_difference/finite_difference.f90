program main

!use vector_mod

implicit none

integer, parameter :: real_kind=selected_real_kind(10,100)

integer :: no, ni, nx
integer :: iter,n_int
integer :: i,j,k,l,i_offset
real (kind=real_kind) :: pi,cf
real (kind=real_kind) :: Qn, Qn_x,Qn_y,l1,rms,error,sten
real (kind=real_kind) :: li,xi,yi,lambda,r,bp(6,2),ext_length,n1,r2,n2
real (kind=real_kind) :: dx,dy,x0,y0,l2,thres
real (kind=real_kind) :: d1,d2,dist,h,last_rms
real (kind=real_kind), allocatable :: t(:,:),x(:,:),y(:,:),Ti(:,:)
logical,allocatable :: calc(:,:)
logical :: interior

print *,'!number of points on x axis'
read *,nx

allocate (T(nx,nx))
allocate (Ti(nx,nx))
allocate (calc(nx,nx))
allocate (x(nx,nx))
allocate (y(nx,nx))

ni=3
no=3

pi=abs(atan2(0d0,-1d0))

! define inner triangle
li=2d0*sqrt(1d0-0.5d0**2)
xi=1d0*cos(pi*30d0/180d0)
yi=1d0*sin(pi*30d0/180d0)
  bp(1,1)=sqrt(3d0)/2d0
  bp(1,2)=0.5d0
  bp(2,1)=sqrt(3d0)
  bp(2,2)=2d0
xi=xi+li*cos(pi*60d0/180d0)
yi=yi+li*sin(pi*60d0/180d0)
  bp(3,1)=sqrt(3d0)*1.5d0
  bp(3,2)=0.5d0

! define outer triangle
xi=0d0
yi=0d0
  bp(4,1)=xi
  bp(4,2)=yi
  bp(5,1)=xi+sqrt(3d0)
  bp(5,2)=yi+3d0
  bp(6,1)=xi+sqrt(3d0)*2d0
  bp(6,2)=yi

!  print *,bp(1,:) 
!  print *,bp(2,:) 
!  print *,bp(3,:) 
!  print *,bp(1,:) 
!  print *,'* *'
!  print *,bp(4,:) 
!  print *,bp(5,:) 
!  print *,bp(6,:) 
!  print *,bp(4,:) 
!  stop

dx=(bp(6,1)-bp(4,1))/dble(nx-1)
dy=dx*sqrt(3d0)/2d0
calc(:,:)=.true.
T(:,:)=0.5d0
n_int=0
cf=huge(cf)
do i=1,nx
  do j=1,nx
    x0=0d0+mod(j+1,2)*dx/2d0
    x(i,j)=x0+(i-1)*dx
    y(i,j)=(j-1)*dy
    interior=.true.
    d1=(x(i,j)-bp(1,1))*(bp(2,2)-bp(1,2))- &
         (y(i,j)-bp(1,2))*(bp(2,1)-bp(1,1))
    do k=2,ni
      d2=(x(i,j)-bp(k,1))*(bp(mod(k,ni)+1,2)-bp(k,2))- &
        (y(i,j)-bp(k,2))*(bp(mod(k,ni)+1,1)-bp(k,1))
      if (d1*d2.le.0d0) interior=.false. 
      d1=d2
    enddo
    if (interior) then
      calc(i,j)=.false.
      T(i,j)=0d0
      cf=min(cf,y(i,j))
    endif
    interior=.true.
    d1=(x(i,j)-bp(ni+1,1))*(bp(ni+2,2)-bp(ni+1,2))- &
         (y(i,j)-bp(ni+1,2))*(bp(ni+2,1)-bp(ni+1,1))
    do k=ni+2,ni+no
      d2=(x(i,j)-bp(k,1))*(bp(mod(k-ni,no)+ni+1,2)-bp(k,2))- &
        (y(i,j)-bp(k,2))*(bp(mod(k-ni,no)+ni+1,1)-bp(k,1))
      if (d1*d2.le.1d-8) interior=.false.
      d1=d2
    enddo
    if (.not.interior) then
      calc(i,j)=.false.
      T(i,j)=1d0
    endif
    if (T(i,j).eq.0.5d0) n_int=n_int+1
  enddo
enddo

cf=cf/0.5d0

! initial guess
do i=1,nx
  do j=1,nx
    d1=huge(d1)
    d2=huge(d2)
    if (calc(i,j)) then
      do k=1,nx
        do l=1,nx 
          if (.not.calc(k,l)) then
            dist=(x(i,j)-x(k,l))**2+(y(i,j)-y(k,l))**2
            if (T(k,l).lt.0.5) then
              d1=min(d1,dist)
            else
              d2=min(d2,dist)
            endif
          endif
        enddo
      enddo
      T(i,j)=sqrt(d1)/(sqrt(d1)+sqrt(d2))
    endif
  enddo
enddo

! iterative solve
h=3d0*dx**2/2d0
last_rms=huge(last_rms)
thres=epsilon(thres)*dble(n_int)**2
do iter=1,n_int**2
  rms=0d0
  do i=1,nx
    do j=1,nx
      if (calc(i,j)) then
        i_offset=-mod(j,2)
!        print *,x(i,j)
!        print *,x(i_offset+i,j-1),x(i_offset+i+1,j-1),x(i-1,j)
        sten=(T(i_offset+i,j-1)+T(i_offset+i+1,j-1)+ &
          T(i-1,j)+T(i+1,j)+ &
          T(i_offset+i,j+1)+T(i_offset+i+1,j+1))
        error=sten-6d0*T(i,j)
        rms=rms+error**2
        T(i,j)=T(i,j)+error*h/6d0*1.5d0
!        print *,Ti(i,j),T(i,j),error
      endif
    enddo
  enddo
  rms=sqrt(rms/n_int)
  if (last_rms.lt.rms) exit
  last_rms=rms
  if (rms.lt.thres) exit
  if (mod(iter,100).eq.0) print *,rms,iter,thres
enddo

Qn=0d0
ext_length=0d0
i_offset=-mod(1,2)
h=sqrt(3d0)*dx
do i=1,nx-1
  l1=x(i+1,1)-x(i,1)
  if ((i.eq.1).or.(i.eq.nx-1)) l1=l1/2d0
  ext_length=ext_length+l1
  Qn_y=(2d0*T(i+i_offset,2)-T(i,1)-T(i+1,1))/h
  Qn=Qn-Qn_y*l1
  print *,Qn_y,l1
enddo
  
print *,'S',3d0*Qn*cf,n_int


deallocate (T,calc,x,y,Ti)     

end program
