module dq_coeff_mod

use data_types_mod
use hyper_dual_mod
use solve_mod

implicit none

contains

subroutine calc_dq_coeff(target_point,points,derivative,coeff,stat,epsilon2_step,set_epsilon2,required_epsilon2,required_precision,initial_epsilon2,initial_precision,taylor_series)

real (kind=real_kind), intent(in) :: target_point(3)
real (kind=real_kind), intent(in) :: points(:,:)
character(*), intent(in) :: derivative(:)
real (kind=real_kind), intent(out) :: coeff(:,:)
integer, intent(out) :: stat
real (kind=real_kind), optional, intent(in) :: epsilon2_step
real (kind=real_kind), optional, intent(in) :: set_epsilon2
real (kind=real_kind), optional, intent(out) :: required_epsilon2
integer, optional, intent(out) :: required_precision
real (kind=real_kind), optional, intent(in) :: initial_epsilon2
integer, optional, intent(in) :: initial_precision
integer, optional, intent(in) :: taylor_series

type (hyper_dual) :: a(size(points,1),size(points,1)),b(size(points,1),size(derivative,1))
type (hyper_dual) :: last_a(size(points,1),size(points,1)),last_b(size(points,1),size(derivative,1))
real (kind=real_kind) :: x(size(points,1),size(derivative,1),0:3)
real (kind=real_kind) :: last_x(size(points,1),size(derivative,1),0:3)
real (kind=real_kind) :: scaled_points(size(points,1),size(points,2))
real (kind=real_kind) :: scaled_target_point(3),scl,e2_step
integer :: int_scl,n,num_der
real (kind=real_kind) :: d2_vp,d_vp,wrk,epsilon2_vp
type (hyper_dual) :: h,d_hd
real (kind=real_kind) :: p1,p2
real (kind=real_kind) :: eps,epsilon2,max_diff(size(b,2)),best_diff

integer :: i,j
logical :: converged,test_both,test_prec,iterate_on_epsilon2
integer :: taylor
logical :: var_taylor

real (kind=real_kind) :: t
type (hyper_dual) :: tmp

stat=0
n=size(points,1)
num_der=size(derivative,1)

int_scl=0
do i=1,n-1
  do j=i+1,n
    int_scl=max(int_scl,floor(log(sqrt((points(i,1)-points(j,1))**2+(points(i,2)-points(j,2))**2+(points(i,3)-points(j,3))**2))/log2r+0.5d0))
  enddo
enddo

scl=2d0**int_scl

scaled_points=points/scl
scaled_target_point=target_point/scl

i=min(n,100)
var_prec_working=var_prec_default(i)
epsilon2=2d0**epsilon2_default(i)

if (present(initial_epsilon2)) then
  epsilon2=1d0/initial_epsilon2
endif

if (present(initial_precision)) then
  var_prec_working=initial_precision
endif

e2_step=4d0
if (var_prec_working.gt.11) e2_step=2d0

if (present(epsilon2_step)) then
  e2_step=epsilon2_step
endif

iterate_on_epsilon2=.true.
if (present(set_epsilon2)) then
  epsilon2=1d0/set_epsilon2
  iterate_on_epsilon2=.false.
endif

var_taylor=.true.
taylor=3
if (present(taylor_series)) then
  taylor=taylor_series
  var_taylor=.false.
endif

last_x=huge(last_x)

test_both=.false.
test_prec=.true.

eps=epsilon(eps)

if (iterate_on_epsilon2) epsilon2=epsilon2/e2_step
var_prec_working=max(var_prec_working-1,2)


do
  test_both=(.not.test_both).and.iterate_on_epsilon2
  if (iterate_on_epsilon2) test_prec=.not.test_prec
  
  if (test_prec.and.iterate_on_epsilon2) var_prec_working=var_prec_working-1
  if (test_both) then
    epsilon2=epsilon2*e2_step
    var_prec_working=var_prec_working+1
  endif

  if (var_prec_working.gt.var_prec_max) then
    stat=1
    exit
  endif

  epsilon2_vp=1d0/epsilon2
  h=hyper_dual(epsilon2_vp,1d0,0d0,0d0)

!  print *,'var_prec',var_prec_working,epsilon2

  if (test_prec.and.(.not.test_both).and.iterate_on_epsilon2) then
     a=last_a
     b=last_b
  else
    p1=0d0
    p2=0d0
    do i=1,n
     do j=1,i
      d2_vp=0d0
      p1=scaled_points(i,1)
      p2=scaled_points(j,1)
      d_vp=p1-p2
      d2_vp=d_vp*d_vp
      p1=scaled_points(i,2)
      p2=scaled_points(j,2)
      d_vp=p1-p2
      d2_vp=d2_vp+d_vp*d_vp
      p1=scaled_points(i,3)
      p2=scaled_points(j,3)
      d_vp=p1-p2
      d2_vp=d2_vp+d_vp*d_vp
      a(i,j)=hyper_dual_exp(-d2_vp*h)
      a(j,i)=a(i,j)
     enddo
    enddo
    last_a=a

    do j=1,num_der
      do i=1,n
        d2_vp=0d0
        p1=scaled_points(i,1)
        p2=scaled_target_point(1)
        wrk=p2-p1
        d2_vp=wrk*wrk
        if (trim(derivative(j)).eq."interp") d_hd=hyper_dual(1d0,0d0,0d0,0d0) 
        if (trim(derivative(j)).eq."d_dx") d_hd=-2*wrk*h
        if (trim(derivative(j)).eq."d2_dx2") d_hd=-2*h+4*wrk*wrk*h*h
        if (derivative(j).eq."d2_dxdy") d_hd=4*wrk*h*h
        if (derivative(j).eq."d2_dxdz") d_hd=4*wrk*h*h
        p1=scaled_points(i,2)
        p2=scaled_target_point(2)
        wrk=p2-p1
        d2_vp=d2_vp+wrk*wrk
        if (trim(derivative(j)).eq."d_dy") d_hd=-2*wrk*h
        if (trim(derivative(j)).eq."d2_dy2") d_hd=-2*h+4*wrk*wrk*h*h
        if (derivative(j).eq."d2_dxdy") d_hd=wrk*d_hd
        if (derivative(j).eq."d2_dydz") d_hd=4*wrk*h*h
        p1=scaled_points(i,3)
        p2=scaled_target_point(3)
        wrk=p2-p1
        d2_vp=d2_vp+wrk*wrk
        if (trim(derivative(j)).eq."d_dz") d_hd=-2*wrk*h
        if (trim(derivative(j)).eq."d2_dz2") d_hd=-2*h+4*wrk*wrk*h*h
        if (derivative(j).eq."d2_dxdz") d_hd=wrk*d_hd
        if (derivative(j).eq."d2_dydz") d_hd=wrk*d_hd
        b(i,j)=d_hd*hyper_dual_exp(-d2_vp*h)
      enddo
    enddo
    last_b=b
  endif

  call ldl_solve(a,b,x,h)

  converged=.false.
  if (var_taylor) then
    do taylor=3,0,-1
      do j=1,num_der
        max_diff(j)=maxval(abs(x(:,j,taylor)-last_x(:,j,taylor)))/(maxval(abs(x(:,j,taylor)))*eps)
      enddo
      best_diff=maxval(max_diff)
      if (best_diff.le.0.5d0) then 
        converged=.true.
        exit
      endif
    enddo
  else
    do j=1,num_der
      max_diff(j)=maxval(abs(x(:,j,taylor)-last_x(:,j,taylor)))/(maxval(abs(x(:,j,taylor)))*eps)
    enddo
    best_diff=maxval(max_diff)
    if (best_diff.le.0.5d0) then 
      converged=.true.
    endif
  endif
  if (taylor.lt.0) taylor=3

  if (test_both.and.converged) exit
  if ((.not.iterate_on_epsilon2).and.converged) exit
  if (test_prec) then
    if (.not.converged) var_prec_working=var_prec_working+1
    if (.not.iterate_on_epsilon2) last_x=x
    cycle
  endif 

  last_x=x
  
enddo

x=x/scl

do j=1,num_der
  if (trim(derivative(j)).eq."interp") x(:,j,taylor)=x(:,j,taylor)*scl
  if ((trim(derivative(j)).eq."d2_dx2").or.(trim(derivative(j)).eq."d2_dy2").or.(trim(derivative(j)).eq."d2_dz2").or.(derivative(j).eq."d2_dxdy").or.(derivative(j).eq."d2_dxdz").or.(derivative(j).eq."d2_dydz")) x(:,j,taylor)=x(:,j,taylor)/scl
enddo

coeff=x(:,:,taylor)

if (present(required_epsilon2)) then
  if (iterate_on_epsilon2) then
    required_epsilon2=e2_step/epsilon2
    if (test_prec.and.iterate_on_epsilon2) required_epsilon2=required_epsilon2/e2_step
  else
    required_epsilon2=1d0/epsilon2
  endif
endif

if (present(required_precision)) then
  required_precision=var_prec_working-1
  if (test_prec.and.iterate_on_epsilon2) required_precision=required_precision+1
endif

end subroutine

end module
