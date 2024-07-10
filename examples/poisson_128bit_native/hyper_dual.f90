module hyper_dual_mod

use data_types_mod

implicit none

type hyper_dual
  real (kind=real_kind) :: r,e1,e1e2,e1e2e3
end type hyper_dual

type (hyper_dual), parameter :: hd_one=hyper_dual(1d0,0d0,0d0,0d0)

interface operator (+)
  module procedure hyper_dual_add
end interface

interface operator (-)
  module procedure hyper_dual_minus, hyper_dual_negate
end interface

interface operator (*)
  module procedure hyper_dual_times, real_times_hyper_dual, int_times_hyper_dual
end interface

interface operator (/)
  module procedure hyper_dual_div
end interface

contains

function hyper_dual_add(a,b)

type (hyper_dual) :: hyper_dual_add
type (hyper_dual), intent(in) :: a,b

hyper_dual_add=hyper_dual(a%r+b%r,a%e1+b%e1,a%e1e2+b%e1e2,a%e1e2e3+b%e1e2e3)

end function

function hyper_dual_minus(a,b)

type (hyper_dual) :: hyper_dual_minus
type (hyper_dual), intent(in) :: a,b

hyper_dual_minus=hyper_dual(a%r-b%r,a%e1-b%e1,a%e1e2-b%e1e2,a%e1e2e3-b%e1e2e3)

end function

function hyper_dual_negate(a)

type (hyper_dual) :: hyper_dual_negate
type (hyper_dual), intent(in) :: a

hyper_dual_negate=hyper_dual(-a%r,-a%e1,-a%e1e2,-a%e1e2e3)

end function

function hyper_dual_times(a,b)

type (hyper_dual) :: hyper_dual_times
type (hyper_dual), intent(in) :: a,b

hyper_dual_times=hyper_dual(a%r*b%r,a%r*b%e1+a%e1*b%r,a%r*b%e1e2+2*a%e1*b%e1+a%e1e2*b%r,a%r*b%e1e2e3+a%e1e2e3*b%r+3*(a%e1*b%e1e2+a%e1e2*b%e1))

end function

function real_times_hyper_dual(r,a)

type (hyper_dual) :: real_times_hyper_dual
type (hyper_dual), intent(in) :: a
real (kind=real_kind), intent(in) :: r
real (kind=real_kind) :: vp

vp=r

real_times_hyper_dual=hyper_dual(vp*a%r,vp*a%e1,vp*a%e1e2,vp*a%e1e2e3)

end function

function vp_times_hyper_dual(r,a)

type (hyper_dual) :: vp_times_hyper_dual
type (hyper_dual), intent(in) :: a
real (kind=real_kind), intent(in) :: r

vp_times_hyper_dual=hyper_dual(r*a%r,r*a%e1,r*a%e1e2,r*a%e1e2e3)

end function

function int_times_hyper_dual(i,a)

type (hyper_dual) :: int_times_hyper_dual
type (hyper_dual), intent(in) :: a
integer, intent(in) :: i

int_times_hyper_dual=hyper_dual(i*a%r,i*a%e1,i*a%e1e2,i*a%e1e2e3)

end function

function hyper_dual_div(a,b)

type (hyper_dual) :: hyper_dual_div
type (hyper_dual), intent(in) :: a,b

real (kind=real_kind) :: r,e1,e1e2,e1e2e3
real (kind=real_kind) :: inv_br,inv_br2,inv_br3,arbe1,be1be1

real (kind=real_kind) :: tmp

inv_br=1d0/b%r

inv_br2=inv_br*inv_br
inv_br3=inv_br2*inv_br
arbe1=a%r*b%e1
be1be1=b%e1*b%e1

r=a%r*inv_br
e1=a%e1*inv_br-arbe1*inv_br2
e1e2=2*arbe1*b%e1*inv_br3-a%r*b%e1e2*inv_br2+a%e1e2*inv_br-2*a%e1*b%e1*inv_br2
e1e2e3=6*arbe1*b%e1e2*inv_br3-a%r*b%e1e2e3*inv_br2-6*arbe1*be1be1*inv_br2*inv_br2-3*a%e1*b%e1e2*inv_br2+6*a%e1*be1be1*inv_br3-3*a%e1e2*b%e1*inv_br2+a%e1e2e3*inv_br

hyper_dual_div=hyper_dual(r,e1,e1e2,e1e2e3)

end function

function hyper_dual_exp(a)

type (hyper_dual) :: hyper_dual_exp
type (hyper_dual), intent(in) :: a

real (kind=real_kind) :: r,ae1ae1

r=exp(a%r)
ae1ae1=a%e1*a%e1

hyper_dual_exp=hyper_dual(r,r*a%e1,r*(a%e1e2+ae1ae1),r*(a%e1e2e3+3*a%e1e2*a%e1+a%e1*ae1ae1))

end function

function dot_product_hd(x,y)

type (hyper_dual) :: dot_product_hd
type (hyper_dual), intent(in) :: x(:),y(:)

integer :: i,n
real (kind=real_kind) :: r,e1,e1e2,e1e2e3

n=size(x,1)

r=dot_product(x%r,y%r)
e1=dot_product(x%r,y%e1)+dot_product(x%e1,y%r)
e1e2=dot_product(x%r,y%e1e2)+2*dot_product(x%e1,y%e1)+dot_product(x%e1e2,y%r)
e1e2e3=dot_product(x%r,y%e1e2e3)+dot_product(x%e1e2e3,y%r)+3*(dot_product(x%e1,y%e1e2)+dot_product(x%e1e2,y%e1))

dot_product_hd=hyper_dual(r,e1,e1e2,e1e2e3)

end function

end module
