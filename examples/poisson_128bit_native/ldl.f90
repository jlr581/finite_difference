module ldl_mod


use data_types_mod
use hyper_dual_mod

implicit none

contains

subroutine ldl(a)

type (hyper_dual), intent(inout) :: a(:,:)

integer :: i,j,n
type (hyper_dual) :: v(size(a,1))
type (hyper_dual) :: tmp

n=size(a,1)

v=hyper_dual(0d0,0d0,0d0,0d0)
do j=2,n
  do i=1,j-1
    v(i)=a(j,i)*a(i,i)
  enddo
  v(j)=a(j,j)-dot_product_hd(a(j,1:j-1),v(1:j-1))
  a(j,j)=v(j)
  tmp=hd_one/v(j)
  do i=j+1,n
    a(i,j)=(a(i,j)-dot_product_hd(a(i,1:j-1),v(1:j-1)))*tmp
  enddo
enddo


end subroutine

end module
