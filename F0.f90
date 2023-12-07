PROGRAM LINEAR1

IMPLICIT NONE


!!!!DECLARATION OF VARIABLES!!!!
!NUMERICAL METHOD
REAL::LB     !LEFT BOUNDARY LOCATION
REAL::RB     !RIGHT BOUNDARY LOCATION
REAL::LENGTH !DOMAIN LENGTH
REAL::DX     !GRID SPACING
REAL::DT     !TIME STEP SIZE
REAL::U	     !VELOCITY
REAL::T_C    !CURRENT TIME LEVEL
REAL::T_END  !FINAL TIME LEVEL
REAL::CFL    !CFL NUMBER
REAL:: PI    !PI (3.14...)
INTEGER::POINTS	!NUMBER OF GRID POINTS
INTEGER::IT	!ITERATIONS NUMBER
INTEGER :: ierr
REAL,ALLOCATABLE, DIMENSION(:)::X !COORDINATES FOR EACH GRID POINT
REAL, ALLOCATABLE, DIMENSION (:,:)::FI !SOLUTION SCALAR
INTEGER::I       !LOOP COUNTER
INTEGER::SCHEME  !SCHEME (CHOSEN BY USER)
REAL:: EPS = 1e-6, alpha(3), beta(3), w(3), sum_w, omega(3) !WENO SCHEME
!ANALYTICAL METHOD
REAL,ALLOCATABLE, DIMENSION(:)::X_A !COORDINATES FOR EACH GRID POINT
REAL, ALLOCATABLE, DIMENSION (:,:)::FI_A !SOLUTION SCALAR
REAL :: T_CA
!ERROR CALCULATION
real :: l1_norm, l2_norm, linf_norm
real :: L1, L2, Linfmax
REAL, ALLOCATABLE, DIMENSION (:,:)::Linf !L INF
!EXECUTION TIME
real:: start_time, end_time, elapsed_time

!START TIME
call CPU_TIME(start_time)

!CALL INPUT FILE
open(30,file="input_noparallel.dat", form="formatted",action="read")
read(30,*)
read(30,*)POINTS
read(30,*)t_end
read(30,*)cfl
read(30,*)scheme
 close(30)

pi = 4.0 * ATAN(1.0)
LB = -10.0
RB = 10.0
U = 1
LENGTH = ABS(RB) + ABS(LB)
DX = LENGTH / REAL(POINTS-1)
DT = CFL * DX / ABS(U)


! Memory allocation
ALLOCATE(X(0:POINTS-1))
ALLOCATE(FI(0:POINTS-1, 1:2))

! Initialize X using a loop
DO I = 0, POINTS-1
  X(I) = LB + REAL(I) * (RB - LB) / REAL(POINTS-1)
END DO

! Initialize FI using a loop
DO I = 0, POINTS-1
  FI(I, 1) = SIN(2.0 * pi * X(I))
END DO



!ANALYTICAL SOLUTION
ALLOCATE(X_A(0:points-1)); X_A=0.0D0
ALLOCATE(FI_A(0:points-1,1:2)); FI_A=0.0D0
ALLOCATE(LINF(0:points,1:2)); LINF=0.0D0

X_A(0)=LB
X_A(POINTS-1)=RB

DO I=1,points-1
	FI_A(I,1)=sin(2*pi*(x(I)-1*t_end))
END DO

open(32,file="f0_analytical.dat", form="formatted",status="replace")
DO I=0,POINTS-1
	WRITE(32,*)X(I), FI_A(I,1)
END DO
 close(32)




!NUMERICAL CALCULATION
if (scheme==1) then
T_C=0.0			!current time
  DO
  	! Ensure the time step does not exceed the final time level
  	DT = MIN(DT, T_END - T_C)
	OPEN(31, file="L2_ftcs_f0.dat", form="formatted")

  	! FTCS Scheme
  	DO I = 1, POINTS-1
   		fi(i,2)=fi(i,1)-(0.5*U*dt/dx)*(fi(i+1,1)-fi(i-1,1))
		l1_norm = sum(abs(fi(:,2)-fi_a(:,1)))
		l2_norm = sqrt(sum((fi(:, 2) - fi_a(:, 1))**2))
    		linf_norm = maxval(abs(fi(:, 2) - fi_a(:, 1)))
	END DO
	fi(0,2)=fi(points-1,2)
	
	! Write the output to the file
    	WRITE(31, *) t_c, l2_norm
	

  	! Current time update
  	t_c = t_c + dt

  	! Update solution (next time level copied to current)
  	fi(:, 1) = fi(:, 2)

	print*,t_c,dt,l1_norm,l2_norm,linf_norm

  	! If the final time level is reached, exit the loop
  	if (t_c.ge.t_end)then
	CLOSE(31)
		exit		!exit loop
	end if
   END DO
print*,"FTCS was used"


! Write the final solution to the output file
OPEN(31, file="f0_ftcs_400.dat", form="formatted", status="replace")
DO I = 0, POINTS-1
  WRITE(31, *) X(I), FI(I, 1)
END DO
 CLOSE(31)

else if (scheme==2) then
T_C=0.0			!current time
  DO
  	! Ensure the time step does not exceed the final time level
  	DT = MIN(DT, T_END - T_C)
	OPEN(31, file="Linf_uw_f0.dat", form="formatted")

  	! Upwind Scheme
  	DO I = 1, POINTS-1
   		fi(i,2)=fi(i,1)-((U*dt/dx)*(fi(i,1)-fi(i-1,1)))
		l1_norm = sum(abs(fi(:,2)-fi_a(:,1)))
		l2_norm = sqrt(sum((fi(:, 2) - fi_a(:, 1))**2))
    		linf_norm = maxval(abs(fi(:, 2) - fi_a(:, 1)))
	END DO
	fi(0,2)=fi(points-1,2)
	
	! Write the output to the file
    	WRITE(31, *) t_c, linf_norm
	

  	! Current time update
  	t_c = t_c + dt

  	! Update solution (next time level copied to current)
  	fi(:, 1) = fi(:, 2)

	print*,t_c,dt,l1_norm,l2_norm,linf_norm

  	! If the final time level is reached, exit the loop
  	if (t_c.ge.t_end)then
	CLOSE(31)
		exit		!exit loop
	end if
   END DO
print*,"Upwind was used"


! Write the final solution to the output file
OPEN(31, file="f0_uw_400.dat", form="formatted", status="replace")
DO I = 0, POINTS-1
  WRITE(31, *) X(I), FI(I, 1)
END DO
 CLOSE(31)


else if(scheme == 3) then
  T_C=0.0			!current time
  DO
  	! Ensure the time step does not exceed the final time level
  	DT = MIN(DT, T_END - T_C)
	OPEN(31, file="Linf_lf_f0.dat", form="formatted")

  	! Lax-Friedrichs Scheme
  	DO I = 1, POINTS-1
   		fi(i,2)=0.5*((u*fi(i+1,1))+(u*fi(i-1,1)))-(u*dt/(2*dx))*(fi(i+1,1)-fi(i-1,1))
		l1_norm = sum(abs(fi(:,2)-fi_a(:,1)))
		l2_norm = sqrt(sum((fi(:, 2) - fi_a(:, 1))**2))
    		linf_norm = maxval(abs(fi(:, 2) - fi_a(:, 1)))S
	END DO
	fi(0,2)=fi(points-1,2)
	
	! Write the output to the file
    	WRITE(31, *) t_c, linf_norm
	

  	! Current time update
  	t_c = t_c + dt

  	! Update solution (next time level copied to current)
  	fi(:, 1) = fi(:, 2)

	print*,t_c,dt,l1_norm,l2_norm,linf_norm

  	! If the final time level is reached, exit the loop
  	if (t_c.ge.t_end)then
	CLOSE(31)
		exit		!exit loop
	end if
   END DO
print*,"Lax-Friedrichs was used"
  

! Write the final solution to the output file
OPEN(31, file="f0_lf_400.dat", form="formatted", status="replace")
DO I = 0, POINTS-1
  WRITE(31, *) X(I), FI(I, 1)
END DO
 CLOSE(31)



else if(scheme == 4) then
  T_C=0.0			!current time
  DO
  	! Ensure the time step does not exceed the final time level
  	DT = MIN(DT, T_END - T_C)
	OPEN(31, file="L2_lw_f0.dat", form="formatted")

  	! Lax-Wendroff Scheme
  	DO I = 1, POINTS-1
   		FI(I, 2) = FI(I, 1) - (U * DT / (2.0 * DX)) * (FI(I+1, 1) - FI(I-1, 1)) + &
                             (U**2 * DT**2 / (2.0 * DX**2)) * (FI(I+1, 1) - 2.0 * FI(I, 1) + FI(I-1, 1))
		l1_norm = sum(abs(fi(:,2)-fi_a(:,1)))
		l2_norm = sqrt(sum((fi(:, 2) - fi_a(:, 1))**2))
    		linf_norm = maxval(abs(fi(:, 2) - fi_a(:, 1)))
	END DO
	fi(0,2)=fi(points-1,2)
	
	! Write the output to the file
    	WRITE(31, *) t_c, l2_norm
	

  	! Current time update
  	t_c = t_c + dt

  	! Update solution (next time level copied to current)
  	fi(:, 1) = fi(:, 2)

	print*,t_c,dt,l1_norm,l2_norm,linf_norm

  	! If the final time level is reached, exit the loop
  	if (t_c.ge.t_end)then
	CLOSE(31)
		exit		!exit loop
	end if
   END DO
print*,"Lax-Wendroff was used"


! Write the final solution to the output file
OPEN(31, file="f0_lw_400.dat", form="formatted", status="replace")
DO I = 0, POINTS-1
  WRITE(31, *) X(I), FI(I, 1)
END DO
 CLOSE(31)


else if(scheme == 5) then
print*,"Using WENO Scheme"


print*,"WENO was used"
else
   print*,"Invalid Scheme"
   stop
end if











!EXECUTION TIME CALCULATION
call CPU_TIME(end_time)
elapsed_time = end_time - start_time

write(*,*)"Elapsed time: ", elapsed_time, "seconds"



END PROGRAM LINEAR1

