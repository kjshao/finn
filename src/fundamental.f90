subroutine fundamental
use mod_reynolds, only : nr,ngroup
use mod_sympoly, only : rtest,mons,maxv,monid,polys,vrecx
use mod_genpoly, only : norder,b,btotal,maxb
use mod_fi, only : nfi,fmons,fmonid,nexit
implicit none
integer :: degree,np,nb
integer :: i,j,nx0,ny0
real(kind=8),allocatable :: frtest(:,:)

! degree of the fundamental invariants to be calculated
read(*,*); read(*,*) degree,i
write(601,'("degree and number of exit point:")'); write(601,'((i2,1x,i2))') degree,i
allocate(nexit(0:i))
nexit(0)=i
read(*,*); read(*,*) nexit(1:i)
write(601,'("exit index of FI:")')
do i=1,nexit(0)
  write(601,'(i4\)') nexit(i)
enddo
allocate(norder(degree))
allocate(rtest(nr),mons(2,nr,maxv),polys(0:nr,ngroup))
allocate(fmons(2,nr,maxv),fmonid(maxv),frtest(nr,maxv))
allocate(b(degree,maxb),btotal(0:degree,maxb))
allocate(vrecx(nr,maxv))
call random_number(rtest)
call random_number(frtest)
frtest=frtest*5.d0
nfi=0
norder=0
open(701,file='dat.log',status='replace')
open(702,file='dat.err',status='replace')
open(703,file='dat.polys',status='replace')
open(705,file='dat.fpolys',status='replace')
open(707,file='dat.nfi',status='replace')
open(711,file='dat.cpolys',status='replace')
open(717,file='dat.latex',status='replace')
!-----------------
do i=1,10
  call sym_poly(i,np)
  write(*,*) "*",np
  do j=1,np
  call write_poly(monid(j),mons(1,1:monid(j)-1,j),mons(2,1:monid(j)-1,j),j)
  enddo
enddo
stop
!-----------------
do i=1,degree
  call sym_poly(i,np)
  if(i.eq.1) then
    fmons(:,:,1:np)=mons(:,:,1:np)
    fmonid(1:np)=monid(1:np)
    norder(1)=np
    nfi=np
    do j=1,nfi
      call write_poly(fmonid(j),fmons(1,1:fmonid(j)-1,j),fmons(2,1:fmonid(j)-1,j),j)
      write(702,'(i3,"*  ",i5,",",i5," in degree ",(i2),",")') j,j,np,i
    enddo
  else if(i.gt.1) then
    call gen_poly(i,nb)
    nx0=nb+np+100 ! nx >= ny+1 in the dgels, here add 100
    ny0=nfi+nb+np
    call cal_error(nx0,ny0,np,nb,frtest,i)
  endif
  write(701,'(i3," fundamental invariants, ",i4," symmetric polynomials in degree: ",i3)') norder(i),np,i
  if(norder(i).eq.0) exit;
enddo
call write_rsingle
close(701); close(702); close(703); close(705); close(707); close(711); close(717)
deallocate(norder,rtest,mons,polys)
deallocate(fmons,fmonid,frtest)
deallocate(b,btotal)
deallocate(vrecx)
end subroutine fundamental
