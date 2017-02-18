subroutine write_poly(id,x,coef,np)
use mod_reynolds, only : nr,ngroup,mat_rey
use mod_sympoly, only : polys
use mod_molecule, only : idrf
implicit none
integer,intent(in) :: id,coef(id-1),x(id-1),np
character(len=100) :: charc,charx,charp,charf
integer :: i,j,k,npoly,i0,idr,idc,maxp
integer :: nchar
write(707,'(1(i4)\)') id
do j=1,id-1
  write(707,'(2(i4)\)') x(j),coef(j)
enddo
write(707,*) ""
polys=0
npoly=0
do k=1,ngroup
    do j=1,id-1
      i=mat_rey(x(j),k)
      polys(i,npoly+1)=coef(j)
    enddo
    if(npoly.eq.0) then
      npoly=npoly+1
      polys(0,npoly)=1
    else
      i0=0
      do j=1,npoly
        i=sum(abs(polys(1:nr,j)-polys(1:nr,npoly+1)))
        if(i.eq.0) then
          polys(0,j)=polys(0,j)+1
          polys(0:nr,npoly+1)=0
          exit
        else
          i0=i0+1
        endif
      enddo
      if(i0.eq.npoly) then
        npoly=npoly+1
        polys(0,npoly)=1
      endif
    endif
enddo
write(charc,*) np-1
write(charf,*) np
nchar=0
if(np.le.9) then
  write(703,'(("p"a1"="\))') trim(adjustl(charf))
  write(705,'(("p("a1")="\))') trim(adjustl(charf))
  write(711,'(("p["a1"]="\))') trim(adjustl(charc))
  write(717,'(("p_"a1"="\))') trim(adjustl(charf))
  nchar=5
elseif(np.ge.10.and.np.le.99) then
  write(703,'(("p"a2"="\))') trim(adjustl(charf))
  write(705,'(("p("a2")="\))') trim(adjustl(charf))
  write(711,'(("p["a2"]="\))') trim(adjustl(charc))
  write(717,'(("p_{"a2"}="\))') trim(adjustl(charf))
  nchar=6
elseif(np.ge.100.and.np.le.999) then
  write(703,'(("p"a3"="\))') trim(adjustl(charf))
  write(705,'(("p("a3")="\))') trim(adjustl(charf))
  write(711,'(("p["a3"]="\))') trim(adjustl(charc))
  write(717,'(("p_{"a3"}="\))') trim(adjustl(charf))
  nchar=7
elseif(np.ge.1000.and.np.lt.9999) then
  write(703,'(("p"a4"="\))') trim(adjustl(charf))
  write(705,'(("p("a4")="\))') trim(adjustl(charf))
  write(711,'(("p["a4"]="\))') trim(adjustl(charc))
  write(717,'(("p_{"a4"}="\))') trim(adjustl(charf))
  nchar=8
endif
do k=1,npoly
  ! The coefs of monomials should be same in a same polynomial,
  ! whose sum equals the order of the group
  ! write(charc,*) polys(0,k)
  ! write(100,'((a3\))') trim(adjustl(charc))
  maxp=0
  do j=1,nr
    if(polys(j,k).gt.0) then
      maxp=j
    endif
  enddo
  do j=1,nr
  idr=idrf(j);idc=polys(j,k)
  if(polys(j,k).gt.0) then
    write(charx,*) idr
    write(charc,*) idr-1
    write(charp,*) idc
    if(idr.lt.10.and.idc.lt.10) then
      if(idc.eq.1) then
        write(703,'(("r"a1\))') trim(adjustl(charx))
        write(717,'(("r_"a1\))') trim(adjustl(charx))
        if(j.lt.maxp) then
          write(705,'(("r("a1")*"\))') trim(adjustl(charx))
          nchar=nchar+5
        else
          write(705,'(("r("a1")"\))') trim(adjustl(charx))
          nchar=nchar+4
        endif
      else
        write(703,'(("r"a1"^"a1\))') trim(adjustl(charx)),trim(adjustl(charp))
        write(717,'(("r_"a1"^"a1\))') trim(adjustl(charx)),trim(adjustl(charp))
        if(j.lt.maxp) then
          write(705,'(("r("a1")**"a1"*"\))') trim(adjustl(charx)),trim(adjustl(charp))
          nchar=nchar+8
        else
          write(705,'(("r("a1")**"a1\))') trim(adjustl(charx)),trim(adjustl(charp))
          nchar=nchar+7
        endif
      endif
    endif
    if(idr-1.lt.10.and.idc.lt.10) then
      if(idc.eq.1) then
        if(j.lt.maxp) then
          write(711,'(("r["a1"]*"\))') trim(adjustl(charc))
        else
          write(711,'(("r["a1"]"\))') trim(adjustl(charc))
        endif
      else
        if(j.lt.maxp) then
          write(711,'(("r["a1"]^"a1"*"\))') trim(adjustl(charc)),trim(adjustl(charp))
        else
          write(711,'(("r["a1"]^"a1\))') trim(adjustl(charc)),trim(adjustl(charp))
        endif
      endif
    endif
    if(idr.lt.10.and.idc.ge.10) then
      write(703,'(("r"a1"^"a2\))') trim(adjustl(charx)),trim(adjustl(charp))
      write(717,'(("r_"a1"^{"a2"}"\))') trim(adjustl(charx)),trim(adjustl(charp))
      if(j.lt.maxp) then
        write(705,'(("r("a1")**"a2"*"\))') trim(adjustl(charx)),trim(adjustl(charp))
        nchar=nchar+8
      else
        write(705,'(("r("a1")**"a2\))') trim(adjustl(charx)),trim(adjustl(charp))
        nchar=nchar+7
      endif
    endif
    if(idr-1.lt.10.and.idc.ge.10) then
      if(j.lt.maxp) then
        write(711,'(("r["a1"]^"a2"*"\))') trim(adjustl(charc)),trim(adjustl(charp))
      else
        write(711,'(("r["a1"]^"a2\))') trim(adjustl(charc)),trim(adjustl(charp))
      endif
    endif
    if(idr.ge.10.and.idc.lt.10) then
      if(idc.eq.1) then
        write(703,'(("r"a2\))') trim(adjustl(charx))
        write(717,'(("r_{"a2"}"\))') trim(adjustl(charx))
        if(j.lt.maxp) then
          write(705,'(("r("a2")*"\))') trim(adjustl(charx))
          nchar=nchar+6
        else
          write(705,'(("r("a2")"\))') trim(adjustl(charx))
          nchar=nchar+5
        endif
      else
        write(703,'(("r"a2"^"a1\))') trim(adjustl(charx)),trim(adjustl(charp))
        write(717,'(("r_{"a2"}^"a1\))') trim(adjustl(charx)),trim(adjustl(charp))
        if(j.lt.maxp) then
          write(705,'(("r("a2")**"a1"*"\))') trim(adjustl(charx)),trim(adjustl(charp))
          nchar=nchar+9
        else
          write(705,'(("r("a2")**"a1\))') trim(adjustl(charx)),trim(adjustl(charp))
          nchar=nchar+8
        endif
      endif
    endif
    if(idr-1.ge.10.and.idc.lt.10) then
      if(idc.eq.1) then
        if(j.lt.maxp) then
          write(711,'(("r["a2"]*"\))') trim(adjustl(charc))
        else
          write(711,'(("r["a2"]"\))') trim(adjustl(charc))
        endif
      else
        if(j.lt.maxp) then
          write(711,'(("r["a2"]^"a1"*"\))') trim(adjustl(charc)),trim(adjustl(charp))
        else
          write(711,'(("r["a2"]^"a1\))') trim(adjustl(charc)),trim(adjustl(charp))
        endif
      endif
    endif
    if(idr.ge.10.and.idc.ge.10) then
        write(703,'(("r"a2"^"a2\))') trim(adjustl(charx)),trim(adjustl(charp))
        write(717,'(("r_{"a2"}^{"a2"}"\))') trim(adjustl(charx)),trim(adjustl(charp))
        if(j.lt.maxp) then
          write(705,'(("r("a2")**"a2"*"\))') trim(adjustl(charx)),trim(adjustl(charp))
          nchar=nchar+9
        else
          write(705,'(("r("a2")**"a2\))') trim(adjustl(charx)),trim(adjustl(charp))
          nchar=nchar+8
        endif
    endif
    if(idr-1.ge.10.and.idc.ge.10) then
        if(j.lt.maxp) then
          write(711,'(("r["a2"]^"a2"*"\))') trim(adjustl(charc)),trim(adjustl(charp))
        else
          write(711,'(("r["a2"]^"a2\))') trim(adjustl(charc)),trim(adjustl(charp))
        endif
    endif
  endif
  enddo
  if(k.lt.npoly) then
    write(703,'("+"\)')
    write(717,'("+"\)')
    if(nchar.gt.80) then
      nchar=0
      write(705,*) "&"
    endif
    write(705,'("+"\)')
    write(711,'("+"\)')
  endif
enddo
write(703,*) ""
write(717,*) ""
write(705,*) ""
write(711,*) ""
end subroutine write_poly
subroutine write_rsingle
use mod_molecule, only : nrs,idrs
use mod_fi, only : nfi
implicit none
character(len=100) :: charc,charr
character(len=100) :: charcx,charrx
integer :: i,id
do i=1,nrs
id=nfi+i
write(charc,*) id
write(charr,*) idrs(i)
write(charcx,*) id-1
write(charrx,*) idrs(i)-1
if(id.le.9) then
  write(705,'(("p("a1")="\))') trim(adjustl(charc))
  write(717,'(("p_"a1"="\))') trim(adjustl(charc))
elseif(id.ge.10.and.id.le.99) then
  write(705,'(("p("a2")="\))') trim(adjustl(charc))
  write(717,'(("p_{"a2"}="\))') trim(adjustl(charc))
elseif(id.ge.100.and.id.le.999) then
  write(705,'(("p("a3")="\))') trim(adjustl(charc))
  write(717,'(("p_{"a3"}="\))') trim(adjustl(charc))
elseif(id.ge.1000.and.id.le.9999) then
  write(705,'(("p("a4")="\))') trim(adjustl(charc))
  write(717,'(("p_{"a4"}="\))') trim(adjustl(charc))
endif
if(idrs(i).lt.10) then
  write(705,'(("r("a1")"))') trim(adjustl(charr))
  write(717,'(("r_"a1""))') trim(adjustl(charr))
else
  write(705,'(("r("a2")"))') trim(adjustl(charr))
  write(717,'(("r_{"a2"}"))') trim(adjustl(charr))
endif
if(id.le.10) then
  write(711,'(("p["a1"]="\))') trim(adjustl(charcx))
elseif(id.gt.10.and.id.le.100) then
  write(711,'(("p["a2"]="\))') trim(adjustl(charcx))
elseif(id.gt.100.and.id.le.1000) then
  write(711,'(("p["a3"]="\))') trim(adjustl(charcx))
elseif(id.gt.1000.and.id.le.10000) then
  write(711,'(("p["a4"]="\))') trim(adjustl(charcx))
endif
if(idrs(i).le.10) then
  write(711,'(("r["a1"]"))') trim(adjustl(charrx))
else
  write(711,'(("r["a2"]"))') trim(adjustl(charrx))
endif
enddo
end subroutine write_rsingle
