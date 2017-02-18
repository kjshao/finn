module mod_reynolds
implicit none
integer :: nr,ngroup
integer,allocatable :: mat_rey(:,:)
end module mod_reynolds

module mod_sympoly
implicit none
integer :: npoly
integer,parameter :: maxv=100000
real(kind=8) :: vrec(maxv)
real(kind=8),allocatable :: rtest(:)
integer,allocatable :: mons(:,:,:),polys(:,:)
integer,allocatable :: vrecx(:,:)
integer :: monid(maxv)
end module mod_sympoly

module mod_genpoly
implicit none
integer,parameter :: maxb=100000
integer,allocatable :: b(:,:),norder(:),btotal(:,:)
integer :: ntotb
end module mod_genpoly

module mod_fi
implicit none
integer :: nfi
integer, allocatable :: fmons(:,:,:),fmonid(:)
integer, allocatable :: nexit(:)
end module mod_fi

module mod_molecule
implicit none
integer :: nrs
integer :: idrf(190),idrs(190)
end module mod_molecule
