subroutine sym_poly(degree,np)
! calculate the symmetric polynomials
! in degree N, save the corresponding
! monomials to mons in module
! output: mons(in module), np
use mod_sympoly, only : npoly
implicit none
integer,intent(in) :: degree
integer,intent(out) :: np
integer :: a(degree)
npoly=0
call intsplit(degree,1,a)
np=npoly
end subroutine sym_poly

recursive subroutine intsplit(n,index,a)
use mod_reynolds, only : nr
implicit none
integer,intent(in) :: n,index
integer,intent(inout) :: a(n)
integer :: i
integer :: x(nr) ! 
if(n.le.0) then
  call combinations(1,nr-index+1+1,1,index-1,nr,a(1:index-1),x)
!  write(*,'(100(1x,i3))') a(1:index-1)
  return
endif
do i=n,1,-1
  a(index)=i
  call intsplit(n-i,index+1,a)
enddo
end subroutine intsplit
recursive subroutine combinations(ia,ib,id,m,n,coef,x)
use mod_reynolds, only : nr
use mod_sympoly, only : rtest,npoly
implicit none
integer,intent(in) :: ia,ib,id,m,n,coef(m)
integer,intent(inout) :: x(nr)
integer :: i,flag
real(kind=8) :: v
if(id.gt.m) then
 call reduce(id,x,coef,flag)
 !call poly_value(id,x,coef,rtest,v)
 !call compare_value(v,flag)
 if(flag.eq.0) then
   call save_monomials(id,x,coef)
 endif
 !if(flag.eq.1) then
 !  npoly=npoly+1
 !  call save_monomials(id,x,coef)
 !endif
 return
endif
do i=ia,ib
  x(id)=i
  call combinations(i+1,ib+1,id+1,m,n,coef,x)
enddo
end subroutine combinations

subroutine reduce(id,xid,coef,flag)
use mod_reynolds, only : nr,ngroup,mat_rey
use mod_sympoly, only : npoly,vrecx,maxv
implicit none
integer,intent(in) :: id,xid(id-1),coef(id-1)
integer,intent(out) :: flag
integer :: i,j,k
integer :: vtmp(nr),tmp
flag=1
if(npoly.eq.0) then
  flag=0
  npoly=1
  vrecx(:,npoly)=0
  do j=1,id-1
    vrecx(mat_rey(xid(j),1),npoly)=coef(j)
  enddo
elseif(npoly.ge.1) then
  flag=0
  do j=1,ngroup
    vtmp=0
    do k=1,id-1
      vtmp(mat_rey(xid(k),j))=coef(k)
    enddo
    do i=1,npoly
      tmp=sum(abs(vtmp-vrecx(:,i)))
      if(tmp.eq.0) then
        flag=1
        exit
      endif
    enddo
    if(flag.eq.1) then
      exit
    endif
  enddo
  if(flag.eq.0) then
    npoly=npoly+1
    vrecx(:,npoly)=0
    do j=1,id-1
      vrecx(mat_rey(xid(j),1),npoly)=coef(j)
    enddo
    if(npoly.gt.maxv) then
      write(*,*) " !!! Error: number of polynomials exceeds..."
    endif
  endif
endif
end subroutine reduce

subroutine poly_value(id,xid,coef,r0,v)
use mod_reynolds, only : nr,ngroup,mat_rey
implicit none
integer,intent(in) :: id,xid(id-1),coef(id-1)
real(kind=8),intent(in) :: r0(nr)
real(kind=8),intent(out) :: v
real(kind=8) :: tmp
integer :: j,k
v=0.d0
do k=1,ngroup
  tmp=1.d0
  do j=1,id-1
    tmp=tmp*r0(mat_rey(xid(j),k))**(real(coef(j)))
  enddo
  v=v+tmp
enddo
end subroutine poly_value
subroutine compare_value(v,flag)
use mod_sympoly, only : npoly,vrec,maxv
implicit none
real(kind=8),intent(in) :: v
integer,intent(out) :: flag
real(kind=8) :: tmp
integer :: i
!------------------------
! compare the value of the new polynomial
! with the rest existing polynomials
flag=1
if(npoly.eq.0) then
  flag=0
  npoly=1
  vrec(1)=v
elseif(npoly.ge.1) then
  flag=0
  do i=1,npoly
    tmp=abs(vrec(i)-v)
    if(tmp.lt.1.d-9) then
      flag=1
    endif
  enddo
  if(flag.eq.0) then
    npoly=npoly+1
    vrec(npoly)=v
    if(npoly.gt.maxv) then
      write(*,*) " !!! Error: number of polynomials exceeds..."
    endif
  endif
endif
end subroutine compare_value
subroutine save_monomials(id,x,coef)
use mod_sympoly, only : mons,npoly,monid
implicit none
integer,intent(in) :: id,x(id-1),coef(id-1)
integer :: j
monid(npoly)=id
do j=1,id-1
  mons(1,j,npoly)=x(j)
  mons(2,j,npoly)=coef(j)
enddo
end subroutine save_monomials
