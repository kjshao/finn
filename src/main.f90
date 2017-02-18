program generators
use mod_reynolds, only : mat_rey
implicit none
!single precesion
!dgels
call random_seed()
open(601,file="dat.input",status='replace')
! calculate the Reynolds Operator
call rey_molecule
! calculate the Fundamental Invariants (Minimal Generating Sets)
! of degree N
call fundamental
close(601)
deallocate(mat_rey)
end program generators
