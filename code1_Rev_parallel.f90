PROGRAM LINEAR1
use MPI
!include 'mpif.h'

IMPLICIT NONE

!!!!DECLARATION OF VARIABLES!!!!
REAL :: LB     ! left boundary LOCATION
REAL :: RB     ! RIGHT BOUNDARY LOCATION
REAL :: LENGTH ! DOMAIN LENGTH
REAL :: DX     ! GRID SPACING
REAL :: DT     ! TIME STEP SIZE
REAL :: U	     ! VELOCITY
REAL :: T_C    ! CURRENT TIME LEVEL
REAL :: T_END  ! FINAL TIME LEVEL
REAL :: CFL    ! CFL NUMBER
INTEGER :: POINTS	! NUMBER OF GRID POINTS
INTEGER :: IT	! ITERATIONS NUMBER
REAL, ALLOCATABLE, DIMENSION(:) :: X ! COORDINATES FOR EACH GRID POINT
REAL, ALLOCATABLE, DIMENSION(:,:) :: FI ! SOLUTION SCALAR
INTEGER :: I    ! LOOP COUNTER
REAL :: pi      ! pi
INTEGER :: parts	! number of elements per cpu
INTEGER :: istart, iend	! beginning and end of array
INTEGER :: myid, ierr		! error integer 
INTEGER, ALLOCATABLE, DIMENSION(:) :: reqs	! requests for non-blocking
INTEGER :: my_rank	! my rank (id)
INTEGER :: n_cpus		! number of processes
INTEGER :: method
CHARACTER(LEN=20) :: PROC3, outfile
INTEGER :: STATUS(MPI_STATUS_SIZE)
REAL :: sentinfo0, sentinfo1, sentinfo2  ! Declarations for MPI communication sentinfo
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status14  ! Declaration of MPI_Status
character(len=20) :: filename
    integer :: iounit
INTEGER :: t_start_exchange, t_end_exchange, t_exchange
INTEGER :: t_start_computation, t_end_computation, t_computation

	CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpus, ierr)


open(30,file="input.dat", form="formatted",action="read")
read(30,*)
read(30,*)POINTS
read(30,*)t_end
read(30,*)cfl
read(30,*)method
 close(30)
call mpi_barrier(MPI_COMM_WORLD, ierr)


pi=4.0*atan(1.0d0)
LB=-10.0
RB=10.0
u=1
LENGTH=ABS(RB)+ABS(LB)


DX=LENGTH/(POINTS-1)
dt=cfl*dx/u
parts=points/n_cpus



istart=(my_rank*parts)+1
iend=(my_rank+1)*parts

if (my_rank.eq.n_cpus-1)then
	iend=POINTS
end if
print*,"here",istart,iend,my_rank
call mpi_barrier(MPI_COMM_WORLD, ierr)


!!!MEMORY ALLOCATION

ALLOCATE(X(istart:iend)); X=0.0D0
ALLOCATE(FI(istart-2:iend+2,1:2)); FI=0.0D0



DO I=istart,iend

X(I)=lb+(DX*I)

END DO

!COORDINATES ALREADY SETUP



DO I=istart,iend
FI(I, 1) = SIN(2.0 * pi * X(I))
END DO


WRITE(PROC3,FMT='(I10)') my_rank
OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".dat"

open(30,file=outfile, form="formatted",status="replace")
DO I=istart,iend
WRITE(30,*)X(I), FI(I,1)
END DO
 close(30)
call mpi_barrier(MPI_COMM_WORLD, ierr)
stop


T_C=0.0		!current time
it=0

do 	!start time loop

	dt=min(dt,t_end-t_c)
	! Start timer for boundary exchange
    	t_start_exchange = MPI_Wtime()
	
	!Multiprocessing start
	!Rank 0
	if(my_rank == 0) then
		x(0) = -10.0d0
		FI(0,1) = 0
		sentinfo0 = FI(iend, 1)
		call MPI_SEND(sentinfo0, 1, MPI_COMPLEX, 1, 1, MPI_COMM_WORLD, ierr) !Send to 1
		call MPI_RECV(sentinfo0, 1, MPI_COMPLEX, 1, 2, MPI_COMM_WORLD, status14, ierr) !Rec from 1
		FI(iend+1,1)=sentinfo0
	end if

	!Rank between 0 and ncp-1
	if(my_rank > 0 .and. my_rank < n_cpus-1) then
		call MPI_RECV(sentinfo1, 1, MPI_COMPLEX, my_rank-1, 1, MPI_COMM_WORLD, status14, ierr) !Request to left
		FI(istart-1,1) = sentinfo1
		sentinfo1 = FI(iend+1,1)
		call MPI_SEND(sentinfo1, 1, MPI_COMPLEX, my_rank+1, 1, MPI_COMM_WORLD, ierr) !Send to right
		sentinfo1 = FI(istart,1)
		call MPI_SEND(sentinfo1, 1, MPI_COMPLEX, my_rank-1, 2, MPI_COMM_WORLD, ierr) !Send to left
		call MPI_RECV(sentinfo1, 1, MPI_COMPLEX, my_rank+1, 2, MPI_COMM_WORLD, status14, ierr) !Rec from right
		FI(iend+1,1) = sentinfo1
	end if

	!If rank=ncp-1
	if(my_rank == n_cpus-1) then
		x(iend+1) = 10.0d0
		FI(iend+1,1) = 0
		FI(istart-1,1) = sentinfo2
		call MPI_RECV(sentinfo2, 1, MPI_COMPLEX, my_rank-1, 1, MPI_COMM_WORLD, status14, ierr) !Rec from left
 		sentinfo2 = FI(istart,1)
		call MPI_SEND(sentinfo2, 1, MPI_COMPLEX, my_rank-1, 2, MPI_COMM_WORLD, ierr) !Send to left
	end if
	! End timer for boundary exchange
    	t_end_exchange = MPI_Wtime()	
	t_exchange = t_end_exchange - t_start_exchange

	call mpi_barrier(MPI_COMM_WORLD, ierr)
	! Start timer for time-step computation
    	t_start_computation = MPI_Wtime()

	!FTCS scheme
	do i=istart, iend-1
		fi(i,2)=fi(i,1)-(0.5*U*dt/dx)*(fi(i+1,1)-fi(i-1,1))
	end do
	fi(istart,2)=fi(iend-1,2)

  	! Current time update
  	t_c = t_c + dt

  	! Update solution (next time level copied to current)
  	Fi(istart:iend,1) = fi(istart:iend,2) 

	

	if (t_c.ge.t_end)then
	! End timer for time-step computation
    	t_end_computation = MPI_Wtime()
    	t_computation = t_end_computation - t_start_computation

	!Print
	print*, 'Time Spent in Boundary Exchange:', t_exchange
	print*, 'Time Spent in time-step computation:', t_computation
	print*, X, FI
		exit		!exit loop
	end if

end do

open(unit=31, file='output.dat', status="replace")
write(31, *) 'Time spent in boundary exchange: ', t_exchange
write(31, *) 'Time spent in time-step computation: ', t_computation
 close(31)


!write final solution
WRITE(PROC3,FMT='(I10)') my_rank
OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".dat"
open(30,file=outfile, form="formatted",status="replace")

DO I=istart-1,iend+1
WRITE(30,*)X(I), FI(I,2)
END DO
 close(30)
call mpi_barrier(MPI_COMM_WORLD, ierr)
call mpi_finalize(ierr)






END PROGRAM LINEAR1
