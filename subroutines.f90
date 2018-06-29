! collection of subroutines used in the main program

subroutine sort(n,m,e,order)
  ! linear sorting

  implicit none

  integer,intent(in)  :: n,m
  real(8),intent(in)  :: e(m)
  integer,intent(out) :: order(m)

  integer :: k,i,j

  do i=1,m
     order=i
  end do

  do i=1,m-1
     k=i

     do j=i+1,m
        if (e(j)<e(k)) then
           k=j
        end if
     end do

     if (k/=i) then
        order((/i,k/))=order((/k,i/))
     end if
  end do

end subroutine sort


subroutine symtridmatvec(n,alp,bet,m,x,y)
  ! symmetric tridiagonal mat-vec y=A*x

  implicit none

  integer :: n,m
  real(8) :: alp(n),bet(n-1),x(n,m),y(n,m)

  integer :: i

  y=0

  do i=1,n
     y(i,:)=y(i,:)+alp(i)*x(i,:)

     if (i/=1) then
        y(i,:)=y(i,:)+bet(i-1)*x(i-1,:)
     end if

     if (i/=n) then
        y(i,:)=y(i,:)+bet(i)*x(i+1,:)
     end if
  end do

end subroutine symtridmatvec


subroutine inv_iter(n,alp,bet,lam,m,v,tau)
  ! inverse iteration

  use omp_lib
  implicit none
  include 'mkl.fi'

  integer,intent(in)  :: n,m
  real(8),intent(in)  :: alp(n),bet(n-1),lam,tau
  real(8),intent(out) :: v(n,m)

  integer :: info,lwork
  integer,allocatable :: ipiv(:)
  real(8) :: a
  real(8),allocatable :: dl(:),d(:),du(:),du2(:),tau2(:),work(:)

  call random_number(v)

  allocate(dl(n-1),d(n),du(n-1),du2(n-2),ipiv(n))
  d=alp-lam
  dl=bet
  du=bet

  call dgttrf(n,dl,d,du,du2,ipiv,info)

  do while (1)
     call dgttrs('N',n,m,dl,d,du,du2,ipiv,v,n,info)

     a=maxval(abs(v))

     if (a*tau>=1) then
        exit
     end if

  end do

  deallocate(dl,d,du,du2,ipiv)

  if (m==1) then
     a=dnrm2(n,v,1)
     v=v/a
  else
     !--- orthonormalize v ---!
     lwork=m   
     allocate(tau2(m),work(lwork))

     call dgeqrf(n,m,v,n,tau2,work,lwork,info)
     call dorgqr(n,m,m,v,n,tau2,work,lwork,info)

     deallocate(tau2,work)
  end if

end subroutine inv_iter


subroutine sturm(n,alp,bet,lam,c)
!!!--- Sturm Sequence ---

!!!--- input ---
!!!      n: matrix size of tridiagonal matrix A
!!!    alp: diagonal entries of A
!!!    bet: subdiagonal entries of A
!!!    lam:
!!!
!!!--- output ---
!!!      c: number of eigenvalues of A that are smaller than lam

  implicit none

  integer,intent(in)  :: n
  real(8),intent(in)  :: alp(n),bet(n-1),lam
  integer,intent(out) :: c

  real(8) :: q
  integer :: i

  q=alp(1)-lam

  if (q<0) then
     c=1
  else
     c=0
  end if

  do i=2,n
     q=alp(i)-lam-bet(i-1)**2/q

     if (q<0) then
        c=c+1
     end if
  end do

end subroutine sturm


subroutine bisection(n,alp,bet,p1,p2,c1,c2,tau,lam)
  ! find simple eigenvalue in interval (c1,c2) via bisection

  implicit none

  integer,intent(in)  :: n,c1,c2
  real(8),intent(in)  :: alp(n),bet(n-1),p1,p2,tau
  real(8),intent(out) :: lam

  real(8) :: pp1,pp2,pp3
  integer :: cc1,cc2,cc3

  pp1=p1
  pp2=p2

  cc1=c1
  cc2=c2

  do while ((pp2-pp1)>tau)
     pp3=(pp1+pp2)/2

     call sturm(n,alp,bet,pp3,cc3)

     if (cc3==cc1) then
        cc1=cc3
        pp1=pp3
     else
        cc2=cc3
        pp2=pp3
     end if
  end do

  lam=(pp1+pp2)/2

end subroutine bisection


subroutine disk(n,alp,bet,a,b)
  ! locate eigenvalue range by Gershgorin disk

  integer,intent(in)  :: n
  real(8),intent(in)  :: alp(n),bet(n-1)
  real(8),intent(out) :: a,b

  integer :: i

  a=min(alp(1)-abs(bet(1)),alp(n)-abs(bet(n-1)))
  b=max(alp(1)+abs(bet(1)),alp(n)+abs(bet(n-1)))

  do i=2,n-1
     a=min(a,alp(i)-abs(bet(i-1))-abs(bet(i)))
     b=min(a,alp(i)+abs(bet(i-1))+abs(bet(i)))
  end do

end subroutine disk