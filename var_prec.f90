module var_prec_mod

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
use var_prec_extra_mod

implicit none

interface assignment (=)
  module procedure real_to_vp, vp_to_real
end interface

interface operator (+)
  module procedure vp_add_vp
end interface

interface operator (-)
  module procedure vp_sub_vp, vp_negate
end interface

interface operator (*)
  module procedure int_mult_vp, vp_mult_vp
end interface

interface operator (/)
  module procedure vp_div_vp
end interface

contains

subroutine real_to_vp(a,b)

type (var_prec), intent(out) :: a
real (kind=real_kind), intent(in) :: b
real (kind=real_kind) :: c

a%sgn=int(sign(1d0,b))

if (b.eq.0d0) then
  a%expn=-1
  a%words(1:var_prec_working)=0
else
  a%expn=floor(real(exponent(b))/var_prec_eff_size)
  do
    c=abs(b)*var_prec_inv_base_real**a%expn
    a%words(1)=floor(c,kind=int_kind)
    if (a%words(1).ne.0) exit
    a%expn=a%expn-1
  enddo
  c=(c-a%words(1))*var_prec_base
  a%words(2)=floor(c+0.5d0,kind=int_kind)
  a%words(3:var_prec_working)=0
endif

a%ops=1

end subroutine

subroutine vp_to_real(a,b)

type (var_prec), intent(in) :: b
real (kind=real_kind), intent(out) :: a

a=real(b%words(2),kind=real_kind)/var_prec_base
a=a+b%words(1)
a=b%sgn*a*var_prec_base_real**b%expn

end subroutine

function vp_add_vp(a,b)

type (var_prec), intent(in) :: a,b
type (var_prec) :: vp_add_vp
integer (kind=int_kind) :: work(1:2*var_prec_working+1)
integer (kind=int_kind), parameter :: one=1
integer :: carry(1:var_prec_working+1),carry_sign
integer :: renorm
integer :: max_ops
integer :: i,j,offset

vp_add_vp%sgn=1
vp_add_vp%ops=1
max_ops=max(a%ops,b%ops)

offset=a%expn-b%expn
if (offset.ge.0) then
  do i=1,var_prec_working
    work(i)=a%sgn*a%words(i)
  enddo
  work(var_prec_working+1:)=0
  j=min(var_prec_working,2*var_prec_working-offset) 
  do i=1,j
    work(i+offset)=work(i+offset)+b%sgn*b%words(i)
  enddo
  vp_add_vp%expn=a%expn
else
  do i=1,var_prec_working
    work(i)=b%sgn*b%words(i)
  enddo
  work(var_prec_working+1:)=0
  vp_add_vp%expn=b%expn
  j=min(var_prec_working,2*var_prec_working+offset) 
  do i=1,j
    work(i-offset)=work(i-offset)+a%sgn*a%words(i)
  enddo
endif

do i=1,var_prec_working+1
  vp_add_vp%words(i)=work(i)
enddo
do
  if (vp_add_vp%words(1).lt.0) then
    vp_add_vp%sgn=-1*vp_add_vp%sgn
    do i=1,var_prec_working+1
      vp_add_vp%words(i)=-vp_add_vp%words(i)
    enddo
  endif
  do i=var_prec_working+1,2,-1
    carry(i)=(1-sign(one,vp_add_vp%words(i)))/2
    carry(i)=carry(i)*max_ops
    vp_add_vp%words(i)=vp_add_vp%words(i)+carry(i)*var_prec_base
    renorm=shiftr(vp_add_vp%words(i),var_prec_eff_size)
    vp_add_vp%words(i)=vp_add_vp%words(i)-renorm*var_prec_base
    carry(i)=carry(i)-renorm
    vp_add_vp%words(i-1)=vp_add_vp%words(i-1)-carry(i)
  enddo
  if (vp_add_vp%words(1).ge.0) exit
enddo
carry_sign=int((sign(one,vp_add_vp%words(1))-1)/2)
carry(1)=int((vp_add_vp%words(1)-carry_sign)/var_prec_base+carry_sign)
vp_add_vp%words(1)=vp_add_vp%words(1)-carry(1)*var_prec_base
if (carry(1).ne.0) then
  vp_add_vp%sgn=sign(1,carry(1))*vp_add_vp%sgn
  vp_add_vp%expn=vp_add_vp%expn+min(1,abs(carry(1)))
  do i=var_prec_working,1,-1
    vp_add_vp%words(i+1)=vp_add_vp%words(i)
  enddo
  vp_add_vp%words(1)=carry(1)
endif
do i=1,var_prec_working+1
  if (vp_add_vp%words(i).ne.0) exit
enddo
if (i.eq.var_prec_working+2) i=0
if (i.gt.1) then
  do j=1,var_prec_working-i+2
    vp_add_vp%words(j)=vp_add_vp%words(j+i-1)
  enddo
  vp_add_vp%words(var_prec_working-i+3:var_prec_working+1)=0
  vp_add_vp%expn=vp_add_vp%expn-i+1
endif
carry(var_prec_working+1)=shiftr(vp_add_vp%words(var_prec_working+1)+var_prec_base_div2,var_prec_eff_size)
vp_add_vp%words(var_prec_working)=vp_add_vp%words(var_prec_working)+carry(var_prec_working+1)

end function

function vp_sub_vp(a,b)

type (var_prec), intent(in) :: a,b
type (var_prec) :: vp_sub_vp
type (var_prec) :: tmp
integer :: i

tmp%sgn=-b%sgn
tmp%ops=b%ops
tmp%expn=b%expn
do i=1,var_prec_working
  tmp%words(i)=b%words(i)
enddo

vp_sub_vp=vp_add_vp(a,tmp)

end function

function vp_negate(a)

type (var_prec), intent(in) :: a
type (var_prec) :: vp_negate
integer :: i

vp_negate%ops=a%ops
vp_negate%expn=a%expn
vp_negate%sgn=-a%sgn

do i=1,var_prec_working
  vp_negate%words(i)=a%words(i)
enddo

end function

function vp_mult_vp(a,b)

type (var_prec), intent(in) :: a,b
type (var_prec) :: vp_mult_vp
integer :: i,k,current_prec
integer (kind=int_kind) :: a_high(1:var_prec_working+2)
integer (kind=int_kind) :: a_low(1:var_prec_working+2)
integer (kind=int_kind) :: b_high(1:var_prec_working+2)
integer (kind=int_kind) :: b_low(1:var_prec_working+2)
integer (kind=int_kind) :: high(0:var_prec_working+2)
integer (kind=int_kind) :: mid(0:var_prec_working+2)
integer (kind=int_kind) :: low(0:var_prec_working+2)
integer (kind=int_kind) :: carry(1:var_prec_working+2)
integer (kind=int_kind) :: remaind(1:var_prec_working+2)
integer (kind=int_kind) :: a_renorm(1:var_prec_working+2)
integer (kind=int_kind) :: b_renorm(1:var_prec_working+2)

vp_mult_vp%expn=a%expn+b%expn
vp_mult_vp%ops=a%ops*b%ops*var_prec_working+1
vp_mult_vp%sgn=a%sgn*b%sgn

if (vp_mult_vp%ops.ge.var_prec_renorm) then
  call split(a%words(1:var_prec_working),carry(1:var_prec_working),a_renorm(1:var_prec_working),var_prec_working)
  a_renorm(var_prec_working+1)=0
  do i=1,var_prec_working-1
    a_renorm(i)=a_renorm(i)+carry(i+1)
  enddo
  if (carry(1).ne.0) then
    a_renorm(2:var_prec_working+1)=a_renorm(1:var_prec_working)
    a_renorm(1)=carry(1)
    vp_mult_vp%expn=vp_mult_vp%expn+1
  endif
  call split(b%words(1:var_prec_working),carry(1:var_prec_working),b_renorm(1:var_prec_working),var_prec_working)
  b_renorm(var_prec_working+1)=0
  do i=1,var_prec_working-1
    b_renorm(i)=b_renorm(i)+carry(i+1)
  enddo
  if (carry(1).ne.0) then
    b_renorm(2:var_prec_working+1)=b_renorm(1:var_prec_working)
    b_renorm(1)=carry(1)
    vp_mult_vp%expn=vp_mult_vp%expn+1
  endif
  call split_half(a_renorm(1:var_prec_working+1),a_high(1:var_prec_working+1),a_low(1:var_prec_working+1),var_prec_working+1)
  call split_half(b_renorm(1:var_prec_working+1),b_high(1:var_prec_working+1),b_low(1:var_prec_working+1),var_prec_working+1)
  vp_mult_vp%ops=var_prec_working+1
else
  call split_half(a%words(1:var_prec_working),a_high(1:var_prec_working),a_low(1:var_prec_working),var_prec_working)
  call split_half(b%words(1:var_prec_working),b_high(1:var_prec_working),b_low(1:var_prec_working),var_prec_working)
  a_high(var_prec_working+1)=0
  a_low(var_prec_working+1)=0
  b_high(var_prec_working+1)=0
  b_low(var_prec_working+1)=0
endif

! short product

current_prec=var_prec_working
var_prec_working=var_prec_working+1

high(1:var_prec_working+1)=0
mid(1:var_prec_working+1)=0
low(1:var_prec_working+1)=0

do i=1,var_prec_working
  high(i)=high(i)+a_high(1)*b_high(i)
  mid(i)=mid(i)+a_high(1)*b_low(i)
  mid(i)=mid(i)+a_low(1)*b_high(i)
  low(i)=low(i)+a_low(1)*b_low(i)
enddo
do i=1,var_prec_working
  high(i+1)=high(i+1)+a_high(2)*b_high(i)
  mid(i+1)=mid(i+1)+a_high(2)*b_low(i)
  mid(i+1)=mid(i+1)+a_low(2)*b_high(i)
  low(i+1)=low(i+1)+a_low(2)*b_low(i)
enddo

do k=3,var_prec_working
  do i=1,var_prec_working-k+2
    high(k+i-1)=high(k+i-1)+a_high(k)*b_high(i)
    mid(k+i-1)=mid(k+i-1)+a_high(k)*b_low(i)
    mid(k+i-1)=mid(k+i-1)+a_low(k)*b_high(i)
    low(k+i-1)=low(k+i-1)+a_low(k)*b_low(i)
  enddo
enddo

vp_mult_vp%words(var_prec_working+1)=low(var_prec_working+1)
do i=1,var_prec_working
  vp_mult_vp%words(i)=high(i+1)+low(i)
enddo
call split_half(mid(1:var_prec_working+1),carry(1:var_prec_working+1),remaind(1:var_prec_working+1),var_prec_working+1)
do i=1,var_prec_working
  vp_mult_vp%words(i)=vp_mult_vp%words(i)+carry(i+1)
enddo
do i=1,var_prec_working+1
  vp_mult_vp%words(i)=vp_mult_vp%words(i)+remaind(i)*var_prec_half_base
enddo
carry(1)=carry(1)+high(1)
if (carry(1).ne.0) then
  vp_mult_vp%words(2:var_prec_working+1)=vp_mult_vp%words(1:var_prec_working)
  vp_mult_vp%words(1)=carry(1)
  vp_mult_vp%expn=vp_mult_vp%expn+1
endif
var_prec_working=current_prec
carry(var_prec_working+1)=shiftr(vp_mult_vp%words(var_prec_working+1)+var_prec_base_div2,var_prec_eff_size)
vp_mult_vp%words(var_prec_working)=vp_mult_vp%words(var_prec_working)+carry(var_prec_working+1)

end function

function int_mult_vp(a,b)

integer, intent(in) :: a
type (var_prec), intent(in) :: b
type (var_prec) :: int_mult_vp
integer :: i
integer (kind=int_kind) :: carry(1:var_prec_working+1)
integer (kind=int_kind) :: remaind(1:var_prec_working+1)
integer (kind=int_kind) :: b_renorm(1:var_prec_working+1)

int_mult_vp%expn=b%expn
int_mult_vp%sgn=b%sgn*sign(1,a)
int_mult_vp%ops=b%ops*abs(a)

if (int_mult_vp%ops.ge.var_prec_renorm) then
  call split(b%words(1:var_prec_working),carry(1:var_prec_working),b_renorm(1:var_prec_working),var_prec_working)
  do i=1,var_prec_working-1
    b_renorm(i)=b_renorm(i)+carry(i+1)
  enddo
  if (carry(1).ne.0) then
    b_renorm(2:var_prec_working)=b_renorm(1:var_prec_working-1)
    b_renorm(1)=carry(1)
    int_mult_vp%expn=int_mult_vp%expn+1
  endif
  do i=1,var_prec_working+1
    int_mult_vp%words(i)=a*b_renorm(i)
  enddo
  carry(var_prec_working+1)=shiftr(int_mult_vp%words(var_prec_working+1)+var_prec_base_div2,var_prec_eff_size)
  int_mult_vp%words(var_prec_working)=int_mult_vp%words(var_prec_working)+carry(var_prec_working+1)
  int_mult_vp%ops=abs(a)
else
  do i=1,var_prec_working
    int_mult_vp%words(i)=a*b%words(i)
  enddo
endif

end function 

function vp_div_vp(a,b)

type (var_prec), intent(in) :: a,b
type (var_prec) :: vp_div_vp
real (kind=real_kind) :: x0
type (var_prec) :: c,b_abs
integer :: i,current_prec
integer (kind=int_kind) :: carry

x0=b
if (x0.eq.0) then
  vp_div_vp=vp_zero
else
  x0=1d0/abs(x0)
  c=x0

  b_abs%sgn=1
  b_abs%ops=b%ops
  b_abs%expn=b%expn
  do i=1,var_prec_working
    b_abs%words(i)=b%words(i)
  enddo
  current_prec=var_prec_working
  var_prec_working=var_prec_working+1
  b_abs%words(var_prec_working)=0
  c%words(var_prec_working)=0
  do i=1,floor(sqrt(real(var_prec_working)))+1
    c=c+c*(vp_one-c*b_abs)
  enddo

  c%sgn=b%sgn
  vp_div_vp=c*a
  carry=shiftr(vp_div_vp%words(var_prec_working)+var_prec_base_div2,var_prec_eff_size)
  var_prec_working=current_prec
  vp_div_vp%words(var_prec_working)=vp_div_vp%words(var_prec_working)+carry
endif

end function

function vp_dot_product(a,b)

type (var_prec), intent(in) :: a(:),b(:)
type (var_prec) :: vp_dot_product
integer :: i,j,k,n,offset,max_expn,max_ops,current_prec
integer (kind=int_kind) :: a_high(1:var_prec_working+2)
integer (kind=int_kind) :: a_low(1:var_prec_working+2)
integer (kind=int_kind) :: b_high(1:var_prec_working+2)
integer (kind=int_kind) :: b_low(1:var_prec_working+2)
integer (kind=int_kind) :: high(0:var_prec_working+2)
integer (kind=int_kind) :: mid(0:var_prec_working+2)
integer (kind=int_kind) :: low(0:var_prec_working+2)
integer (kind=int_kind) :: carry(1:var_prec_working+2)
integer (kind=int_kind) :: remaind(1:var_prec_working+2)
integer (kind=int_kind) :: a_renorm(1:var_prec_working+2)
integer (kind=int_kind) :: b_renorm(1:var_prec_working+2)
integer (kind=int_kind) :: work(1:2*var_prec_working+3)
integer (kind=int_kind), parameter :: one=1
integer :: carry_sign,renorm,crry(1:var_prec_working+3)
integer :: products_sgn(1:size(a,1))
integer :: products_expn(1:size(a,1))
integer (kind=int_kind) :: products_words(1:var_prec_max+3,1:size(a,1))

n=size(a,1)

do j=1,n

  products_expn(j)=a(j)%expn+b(j)%expn
  products_sgn(j)=a(j)%sgn*b(j)%sgn

  call split(a(j)%words(1:var_prec_working),carry(1:var_prec_working),a_renorm(1:var_prec_working),var_prec_working)
  a_renorm(var_prec_working+1)=0
  do i=1,var_prec_working-1
    a_renorm(i)=a_renorm(i)+carry(i+1)
  enddo
  if (carry(1).ne.0) then
    a_renorm(2:var_prec_working+1)=a_renorm(1:var_prec_working)
    a_renorm(1)=carry(1)
    products_expn(j)=products_expn(j)+1
  endif
  call split(b(j)%words(1:var_prec_working),carry(1:var_prec_working),b_renorm(1:var_prec_working),var_prec_working)
  b_renorm(var_prec_working+1)=0
  do i=1,var_prec_working-1
    b_renorm(i)=b_renorm(i)+carry(i+1)
  enddo
  if (carry(1).ne.0) then
    b_renorm(2:var_prec_working+1)=b_renorm(1:var_prec_working)
    b_renorm(1)=carry(1)
    products_expn(j)=products_expn(j)+1
  endif
  call split_half(a_renorm(1:var_prec_working+1),a_high(1:var_prec_working+1),a_low(1:var_prec_working+1),var_prec_working+1)
  call split_half(b_renorm(1:var_prec_working+1),b_high(1:var_prec_working+1),b_low(1:var_prec_working+1),var_prec_working+1)

! short product

  current_prec=var_prec_working
  var_prec_working=var_prec_working+1

  high(1:var_prec_working+1)=0
  mid(1:var_prec_working+1)=0
  low(1:var_prec_working+1)=0

  do i=1,var_prec_working
    high(i)=high(i)+a_high(1)*b_high(i)
    mid(i)=mid(i)+a_high(1)*b_low(i)
    mid(i)=mid(i)+a_low(1)*b_high(i)
    low(i)=low(i)+a_low(1)*b_low(i)
  enddo
  do i=1,var_prec_working
    high(i+1)=high(i+1)+a_high(2)*b_high(i)
    mid(i+1)=mid(i+1)+a_high(2)*b_low(i)
    mid(i+1)=mid(i+1)+a_low(2)*b_high(i)
    low(i+1)=low(i+1)+a_low(2)*b_low(i)
  enddo

  do k=3,var_prec_working
    do i=1,var_prec_working-k+2
      high(k+i-1)=high(k+i-1)+a_high(k)*b_high(i)
      mid(k+i-1)=mid(k+i-1)+a_high(k)*b_low(i)
      mid(k+i-1)=mid(k+i-1)+a_low(k)*b_high(i)
      low(k+i-1)=low(k+i-1)+a_low(k)*b_low(i)
    enddo
  enddo

  products_words(var_prec_working+1,j)=low(var_prec_working+1)
  do i=1,var_prec_working
    products_words(i,j)=high(i+1)+low(i)
  enddo
  call split_half(mid(1:var_prec_working+1),carry(1:var_prec_working+1),remaind(1:var_prec_working+1),var_prec_working+1)
  do i=1,var_prec_working
    products_words(i,j)=products_words(i,j)+carry(i+1)
  enddo
  do i=1,var_prec_working+1
    products_words(i,j)=products_words(i,j)+remaind(i)*var_prec_half_base
  enddo
  carry(1)=carry(1)+high(1)
  if (carry(1).ne.0) then
    products_words(2:var_prec_working+1,j)=products_words(1:var_prec_working,j)
    products_words(1,j)=carry(1)
    products_expn(j)=products_expn(j)+1
  endif
  var_prec_working=current_prec
enddo
  
current_prec=var_prec_working
var_prec_working=var_prec_working+1

vp_dot_product%sgn=1
vp_dot_product%ops=1
max_ops=(var_prec_working+1)*n

max_expn=maxval(products_expn(:))
vp_dot_product%expn=max_expn

work(:)=0
offset=max_expn-products_expn(1)
do k=1,var_prec_working+1
  j=min(var_prec_working+1,2*var_prec_working+1-offset) 
  do i=1,j
    work(i+offset)=products_sgn(1)*products_words(i,1)
  enddo
enddo

do k=2,n
  offset=max_expn-products_expn(k)
  j=min(var_prec_working+1,2*var_prec_working+1-offset) 
  do i=1,j
    work(i+offset)=work(i+offset)+products_sgn(k)*products_words(i,k)
  enddo
enddo

do i=1,var_prec_working+1
  vp_dot_product%words(i)=work(i)
enddo
do
  if (vp_dot_product%words(1).lt.0) then
    vp_dot_product%sgn=-1*vp_dot_product%sgn
    do i=1,var_prec_working+1
      vp_dot_product%words(i)=-vp_dot_product%words(i)
    enddo
  endif
  do i=var_prec_working+1,2,-1
    crry(i)=(1-sign(one,vp_dot_product%words(i)))/2
    crry(i)=crry(i)*max_ops
    vp_dot_product%words(i)=vp_dot_product%words(i)+crry(i)*var_prec_base
    renorm=shiftr(vp_dot_product%words(i),var_prec_eff_size)
    vp_dot_product%words(i)=vp_dot_product%words(i)-renorm*var_prec_base
    crry(i)=crry(i)-renorm
    vp_dot_product%words(i-1)=vp_dot_product%words(i-1)-crry(i)
  enddo
  if (vp_dot_product%words(1).ge.0) exit
enddo
carry_sign=int((sign(one,vp_dot_product%words(1))-1)/2)
crry(1)=int((vp_dot_product%words(1)-carry_sign)/var_prec_base+carry_sign)
vp_dot_product%words(1)=vp_dot_product%words(1)-crry(1)*var_prec_base
if (crry(1).ne.0) then
  vp_dot_product%sgn=sign(1,crry(1))*vp_dot_product%sgn
  vp_dot_product%expn=vp_dot_product%expn+min(1,abs(crry(1)))
  do i=var_prec_working,1,-1
    vp_dot_product%words(i+1)=vp_dot_product%words(i)
  enddo
  vp_dot_product%words(1)=crry(1)
endif
do i=1,var_prec_working+1
  if (vp_dot_product%words(i).ne.0) exit
enddo
if (i.eq.var_prec_working+2) i=0
if (i.gt.1) then
  do j=1,var_prec_working-i+2
    vp_dot_product%words(j)=vp_dot_product%words(j+i-1)
  enddo
  vp_dot_product%words(var_prec_working-i+3:var_prec_working+1)=0
  vp_dot_product%expn=vp_dot_product%expn-i+1
endif
crry(var_prec_working+1)=shiftr(vp_dot_product%words(var_prec_working+1)+var_prec_base_div2,var_prec_eff_size)
vp_dot_product%words(var_prec_working)=vp_dot_product%words(var_prec_working)+crry(var_prec_working+1)
var_prec_working=current_prec
crry(var_prec_working+1)=shiftr(vp_dot_product%words(var_prec_working+1)+var_prec_base_div2,var_prec_eff_size)
vp_dot_product%words(var_prec_working)=vp_dot_product%words(var_prec_working)+crry(var_prec_working+1)

end function

function vp_exp(a)

type (var_prec), intent(in) :: a
type (var_prec) :: vp_exp
type (var_prec) :: r,s,rp
real (kind=real_kind) :: ar
integer :: i,m,current_prec,terms_indx,carry
type (var_prec) :: tmp,fact

ar=a
m=int(ar/log2r)
current_prec=var_prec_working
var_prec_working=var_prec_working+1
tmp=a
tmp%words(var_prec_working)=0
r=tmp+m*log2
ar=r
terms_indx=max(exponent(ar),-65)

s%ops=1
s%sgn=1
s%expn=-1
s%words(2:var_prec_working)=0
s%words(1)=var_prec_base
rp=r

fact=vp_one
do i=1,vp_n_terms(terms_indx,current_prec)
  s=s+rp*vp_inv_fact(i)
  rp=rp*r
enddo

ar=2d0**m
tmp=ar
vp_exp=tmp*s
var_prec_working=current_prec
carry=shiftr(vp_exp%words(var_prec_working+1)+var_prec_base_div2,var_prec_eff_size)
vp_exp%words(var_prec_working)=vp_exp%words(var_prec_working)+carry

end function

end module

