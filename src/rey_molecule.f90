subroutine rey_molecule
implicit none
integer :: natom,nkind,ir
character(len=1),allocatable :: label(:)
integer,allocatable :: id(:),idx(:),nperm(:)
integer,allocatable :: idr(:,:)
integer :: i
!---------------------------------------------------------
! calculate the reynolds operator for the molecular system
! AxBy...Mp, the corresponding group is the direct product
! of symmetric group SxSy...Sp
read(*,*); read(*,*) natom
write(601,'("number of atoms:")'); write(601,'(i2)') natom
allocate(label(natom),id(natom),idx(natom),nperm(natom))
allocate(idr(natom*(natom-1)/2,2))
read(*,*); read(*,*) label(1:natom) ! label of atoms, as A A A B B C ...
write(601,'("labels of atoms:")')
do i=1,natom
  write(601,'(a2\)') adjustl(label(i))
enddo
write(601,*) ""
call kinds(natom,label,nkind,id,idx,nperm)
call index_r(natom,nperm,idx,ir,idr)
call permutation(natom,nkind,ir,id,nperm,idr)
deallocate(label,id,idx,nperm,idr)
!---------------------------------------------------------
end subroutine rey_molecule
subroutine kinds(natom,label,nkind,id,idx,nperm)
use mod_reynolds, only : ngroup
! Count the atom number of different kinds of atoms
! Calculate the order of the group SxSy...Sp
implicit none
integer,intent(in) :: natom
character(len=1),intent(in) :: label(natom)
integer,intent(out) :: nkind,id(natom),idx(natom),nperm(natom)
character(len=1) :: atmp
integer :: i,j
nkind=1
atmp=label(1)
id(nkind)=1
do i=1,natom
  if(trim(label(i)).ne.trim(atmp)) then
    nkind=nkind+1 ! kinds of atoms
    atmp=label(i)
    id(nkind)=i ! index of atoms kinds
  endif
idx(i)=nkind
enddo
id(nkind+1)=natom+1
ngroup=1
do i=1,nkind
  nperm(i)=id(i+1)-id(i) ! number of N-th kind atom
  do j=1,nperm(i)
    ngroup=ngroup*j
  enddo 
enddo
write(*,*) " Group order:",ngroup
end subroutine kinds
subroutine index_r(natom,nperm,idx,ir,idr)
use mod_reynolds, only : nr
use mod_molecule, only : nrs,idrs,idrf
implicit none
integer,intent(in) :: natom,nperm(natom),idx(natom)
integer,intent(out) :: ir,idr(natom*(natom-1)/2,2)
integer :: i,j,id
ir=0
nrs=0
id=0
if(natom.gt.20) then
  write(*,*) "Error: idrf exceed"
endif
open(1,file='datx.rij',status='replace')
do i=1,natom-1
  do j=i+1,natom
      id=id+1
      if(nperm(idx(i)).eq.1.and.nperm(idx(j)).eq.1) then
        ! only have one atom in the different 2 kinds
        nrs=nrs+1; idrs(nrs)=id 
        cycle
      elseif(idx(i).eq.idx(j)) then
        if(nperm(idx(i)).eq.2) then
        ! only 2 atoms in a same kind, remove it
          nrs=nrs+1; idrs(nrs)=id
          cycle
        else
          ir=ir+1
          idrf(ir)=id
          idr(ir,1)=i ! index of the atoms of r_ij
          idr(ir,2)=j
          write(1,*) idr(ir,1),idr(ir,2)
        endif
      else
        ir=ir+1
        idrf(ir)=id
        idr(ir,1)=i
        idr(ir,2)=j
        write(1,*) idr(ir,1),idr(ir,2)
      endif
  enddo
enddo
close(1)
nr=ir
end subroutine index_r
subroutine permutation(natom,nkind,ir,id,nperm,idr)
implicit none
integer,intent(in) :: natom,nkind,ir,id(natom),nperm(natom)
integer,intent(in) :: idr(natom*(natom-1)/2,2)
integer :: ntotal(-1:nkind+1),ipall(natom)
integer :: i,j
integer,allocatable :: ip(:,:,:)
i=maxval(nperm(1:nkind))
ntotal=1
ntotal(-1)=i
ntotal(0)=1
do j=i,1,-1
  ntotal(0)=ntotal(0)*j ! max number of permutation of one kind atom
enddo
allocate(ip(i,ntotal(0),nkind)) ! all permutations of atom A, atom B, atom C...
do i=1,nkind
    ntotal(i)=1
    do j=nperm(i),1,-1
      ntotal(i)=ntotal(i)*j
    enddo
    ntotal(nkind+1)=ntotal(nkind+1)*ntotal(i)
    call cperm(nperm(i),ntotal(i),ip(1:nperm(i),1:ntotal(i),i))
    do j=1,ntotal(i)
      ip(1:nperm(i),j,i)=ip(1:nperm(i),j,i)+id(i)-1
    enddo
enddo
! direct product of ip(:,i)
open(200,file='datx.perm',status='replace')
call directpro(1,nkind,ntotal,natom,id,ip,ipall)
close(200)
call get_reynolds(natom,nkind,ntotal,ir,idr)
deallocate(ip)
end subroutine permutation
subroutine cperm(nperm,ntotal,px)
implicit none
integer,intent(in) :: nperm,ntotal
integer,intent(inout) :: px(nperm,ntotal)
integer :: i,itmp,itotal
integer :: n,flag,id,idmax(2),imax
integer,allocatable :: p(:,:)
! p(:,1) the elements
! p(:,2) the directions, -1 is left, 1 is right
itotal=0
n=nperm
!open(100,file='cperm.dat',status='replace')
allocate(p(n,2))
do i=1,n
  p(i,1)=i
  p(i,2)=-1
enddo
!do iloop=1,n
!  write(100,'((1x,i5\))') p(iloop,1)
!enddo
!write(100,*) ""
itotal=itotal+1
px(1:n,itotal)=p(1:n,1)
flag=1
do while(flag.eq.1)
  imax=-1
  do i=1,n
    id=i+p(i,2) ! index after moving
    if(id.ge.1.and.id.le.n) then
      if(p(i,1).gt.p(id,1)) then
        if(p(i,1).gt.imax) then
          idmax(1)=i
          idmax(2)=id
          imax=p(i,1)
        endif
      endif
    endif
  enddo
  if(imax.eq.-1) then
    flag=0
  else ! change elements and directions
    do i=1,n
      if(p(i,1).gt.p(idmax(1),1)) then
        p(i,2)=-p(i,2)
      endif
    enddo
    itmp=p(idmax(1),1)
    p(idmax(1),1)=p(idmax(2),1)
    p(idmax(2),1)=itmp
    itmp=p(idmax(1),2)
    p(idmax(1),2)=p(idmax(2),2)
    p(idmax(2),2)=itmp
    !do iloop=1,n
    !  write(100,'((1x,i5\))') p(iloop,1)
    !enddo
    !write(100,*) ""
    itotal=itotal+1
    px(1:n,itotal)=p(1:n,1)
  endif
enddo
!close(100)
deallocate(p)
end subroutine cperm
recursive subroutine directpro(k,nkind,ntotal,natom,id,ip,ipall)
implicit none
integer,intent(in) :: k,nkind,natom,ntotal(-1:nkind),id(nkind+1)
integer,intent(in) :: ip(ntotal(-1),ntotal(0),nkind)
integer,intent(inout) :: ipall(natom)
integer :: i,nperm,iloop
if(k.gt.nkind) then
  do iloop=1,natom
    write(200,'((i3\))') ipall(iloop)
  enddo
  write(200,*) ""
  return
endif
do i=1,ntotal(k)
  nperm=id(k+1)-id(k)
  ipall(id(k):id(k+1)-1)=ip(1:nperm,i,k)
  call directpro(k+1,nkind,ntotal,natom,id,ip,ipall)
enddo
end subroutine directpro
subroutine get_reynolds(natom,nkind,ntotal,ir,idr)
use mod_reynolds, only : nr,ngroup,mat_rey
implicit none
integer,intent(in) :: natom,nkind,ntotal(-1:nkind+1),ir
integer,intent(in) :: idr(natom*(natom-1)/2,2)
integer :: i,j,k,iloop,ida,idb,i1,i2,ipall(natom)
integer,allocatable :: idout(:),idmat(:,:)
allocate(idout(natom*(natom-1)/2))
allocate(idmat(ir,ir),mat_rey(nr,ngroup))
open(1,file='datx.perm',status='old')
open(100,file='datx.mata',status='replace')
open(200,file='datx.matb',status='replace')
do k=1,ntotal(nkind+1)
  read(1,*) ipall(1:natom)
  i2=0
  do i=1,natom-1
    do j=i+1,natom
      ida=ipall(i)
      idb=ipall(j)
      if(ida.gt.idb) then
        i1=ida
        ida=idb
        idb=i1
      endif  
      do i1=1,ir
        if(ida.eq.idr(i1,1).and.idb.eq.idr(i1,2)) then
          i2=i2+1
          idout(i2)=i1
        endif
      enddo
    enddo
  enddo
  do iloop=1,ir
    write(100,'((i3\))') idout(iloop)
    mat_rey(iloop,k)=idout(iloop)
  enddo
    write(100,*) ""
  do i=1,ir
    do j=1,ir
      if(idout(i).eq.j) then
        idmat(j,i)=1
      else
        idmat(j,i)=0
      endif
    enddo
    do iloop=1,ir
      write(200,'((1x,i1\))') idmat(iloop,i)
    enddo
      write(200,*) ""
  enddo
  write(200,*) ""
enddo
close(1)
close(100)
close(200)
deallocate(idout,idmat)
end subroutine get_reynolds
