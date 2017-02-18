subroutine gen_poly(degree,np)
use mod_genpoly, only : ntotb
implicit none
integer,intent(in) :: degree
integer,intent(out) :: np
integer :: a(degree)
ntotb=0
call intsplitx(degree,1,a,degree)
np=ntotb
end subroutine gen_poly

recursive subroutine intsplitx(n,index,ax,degree)
! split integer with decrement number
! remove the first one
use mod_genpoly, only : b,norder,btotal,ntotb,maxb
implicit none
integer,intent(in) :: n,index,degree
integer,intent(inout) :: ax(n)
integer :: i,j,k,id,id1,id2,ntmp(degree)
integer :: indexb,ntotal,nx
integer :: atmp(degree) ! the highest degree of ntmp(i) is n-1
if(n.le.0.and.index-1.ge.2) then ! index-1.ge.2 to remove the first one
  ntmp=0
  do i=1,index-1
    ntmp(ax(i))=ntmp(ax(i))+1 ! the number of 
  enddo
  ! ntmp(1) = 3  degree 1 occurs 3 times in ax(1:n)
  ! ntmp(2) = 1  degree 2 occurs 1 time  in ax(1:n)
  ! ntmp(3) = 2  degree 3 occurs 2 times in ax(1:n)
  ! ...
  ntotal=0
  nx=0
  do i=1,maxval(ax(1:index-1))
    if(ntmp(i).gt.0) then
      indexb=0
      call comb_order(1,1,indexb,ntmp(i),norder(i),atmp)
      if(ntotal.eq.0) then
        ntotal=indexb; nx=ntmp(i)
        do j=1,indexb
          do k=1,ntmp(i)
            btotal(k,j+ntotb)=b(k,j)+sum(norder(1:i-1))
          enddo
        enddo
      elseif(ntotal.gt.0) then
        if(ntotb+ntotal*indexb.gt.maxb) write(*,*) " Error: ntotb"
        do j=1,indexb
          id1=ntotal*(j-1)+1+ntotb; id2=ntotal*j+ntotb
          btotal(1:nx,id1:id2)=btotal(1:nx,1+ntotb:ntotal+ntotb)
          do k=1,ntmp(i)
            btotal(nx+k,id1:id2)=b(k,j)+sum(norder(1:i-1))
          enddo
        enddo
        ntotal=ntotal*indexb
        nx=nx+ntmp(i)
      endif
    endif
  enddo
  btotal(0,ntotb+1:ntotb+ntotal)=nx
  ntotb=ntotb+ntotal
  return
endif
if(index.eq.1) then
  id=n
elseif(index.gt.1) then
  if(n.ge.ax(index-1)) then
    id=ax(index-1)
  else
    id=n
  endif
endif
do i=id,1,-1
  ax(index)=i
  call intsplitx(n-i,index+1,ax,degree)
enddo
end subroutine intsplitx

recursive subroutine comb_order(t,index,indexb,m,n,a)
use mod_genpoly, only : b
! select M elements in a set
! contains N elements
implicit none
integer,intent(in) :: t,m,n,index
integer,intent(inout) :: a(m),indexb
integer :: i
if(index-1.eq.m) then
  indexb=indexb+1
  b(1:index-1,indexb)=a(1:index-1)
  return
endif
do i=t,n
  a(index)=i
  call comb_order(i,index+1,indexb,m,n,a)
enddo
end subroutine comb_order
