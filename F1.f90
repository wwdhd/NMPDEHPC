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
REAL :: L1, L2, Linfmax
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

!VARIABLES DEFINTION
pi = 4.0 * ATAN(1.0)
LB=-10
RB=10.0D0
u=1
LENGTH=ABS(RB)+ABS(LB)

DX=LENGTH/(POINTS-1)

!TIMESTEP INITIALIZATION
dt=cfl*dx/u

!!!MEMORY ALLOCATION
ALLOCATE(X(0:points-1)); X=0.0D0
ALLOCATE(FI(0:points-1,1:2)); FI=0.0D0

X(0)=LB
X(POINTS-1)=RB

DO I=1,POINTS-2
X(I)=X(0)+(DX*I)
END DO

!BOUNDARY CONDITION

FI(0,:)=-1.0D0
FI(points-1,:)=1.0D0

!INITIAL VALUE DISTRIBUTION
DO I=1,POINTS-2
FI(I,1)=(SIGN(1.0,X(I)))
END DO

open(30,file="initial.dat", form="formatted",status="replace")
DO I=0,POINTS-1
WRITE(30,*)X(I), FI(I,1)
END DO
 close(30)

!Analytical Solution
ALLOCATE(X_A(0:points-1)); X_A=0.0D0
ALLOCATE(FI_A(0:points-1,1:2)); FI_A=0.0D0
ALLOCATE(LINF(0:points,1:2)); LINF=0.0D0

X_A(0)=LB
X_A(POINTS-1)=RB

DO I=0,POINTS-1
FI_A(I,1)=SIGN(1.0,X_A(I)-1*t_end)
END DO
open(32,file="analytical.dat", form="formatted",status="replace")
DO I=0,POINTS-1
WRITE(32,*)X(I), FI_A(I,1)
END DO
 close(32)


!NUMERICAL CALCULATION
if (scheme == 1) then
   T_C=0.0			!current time
   open(31,file="Linf_FTCS_F1.dat", form="formatted",status="replace")

   do 	!start time loop

	dt=min(dt,t_end-t_c)
	print*,t_c,dt,L1,L2

	!FTCS Scheme
	do i=1,points-2
		fi(i,2)=fi(i,1)-(0.5*U*dt/dx)*(fi(i+1,1)-fi(i-1,1))
		if (fi(i,1)< fi_a(i,1)) then
	   	   Linf(i,1) = fi_a(i,1) - fi(i,1)
		else
	   	   Linf(i,1) = fi(i,1) - fi_a (i,1)
		end if
		!Compute the errors
		L1 = dx*sum(abs(fi-fi_a))
		L2 = sqrt(dx*sum((fi-fi_a)**2))
		Linfmax = maxval(Linf(:,1))
		WRITE(31,*)t_c,Linfmax
	end do


	!current time update
	t_c=t_c+dt

	!update solution (next time level copied to current)
	fi(:,1)=fi(:,2)

	print*,t_c,dt,L1,L2,Linfmax
	if (t_c.ge.t_end)then
		exit		!exit loop
	end if

   end do

   print*,"Euler FTCS is used"

   !write final solution
   open(31,file="final_FTCS_F1.dat", form="formatted",status="replace")
   DO I=0,POINTS-1
   WRITE(31,*)X(I), FI(I,1)
   END DO
   close(31)

else if (scheme == 2) then
   T_C=0.0			!current time
   open(31,file="L2_UW_F1.dat", form="formatted",status="replace")


   do 	!start time loop
	
	dt=min(dt,t_end-t_c)

	!Upwind Scheme
	do i=1,points-2
		fi(i,2)=fi(i,1)-((U*dt/dx)*(fi(i,1)-fi(i-1,1)))
		if (fi(i,1)< fi_a(i,1)) then
	   	   Linf(i,1) = fi_a(i,1) - fi(i,1)
		else
	   	   Linf(i,1) = fi(i,1) - fi_a (i,1)
		end if
		!Compute the errors
		L1 = dx*sum(abs(fi-fi_a))
		L2 = sqrt(dx*sum((fi-fi_a)**2))
		Linfmax = maxval(Linf(:,1))
		WRITE(31,*)t_c,L2
	end do

	!current time update
	t_c=t_c+dt

	!update solution (next time level copied to current)
	fi(:,1)=fi(:,2)

	!Compute the errors
	L1 = dx*sum(abs(fi-fi_a))
	L2 = sqrt(dx*sum((fi-fi_a)**2))
	Linfmax = maxval(Linf(:,1))

	print*,t_c,dt,L1,L2,Linfmax
	if (t_c.ge.t_end)then
		exit		!exit loop
	end if

   end do

   print*,"Upwind Scheme is used"

   !write final solution
   open(31,file="final_UW_F1.dat", form="formatted",status="replace")
   DO I=0,POINTS-1
   WRITE(31,*)X(I), FI(I,1)
   END DO
   close(31)

else if (scheme== 3) then
   T_C=0.0			!current time
   open(31,file="Linf_LF_F1.dat", form="formatted",status="replace")

   do 	!start time loop
	
	dt=min(dt,t_end-t_c)
	print*,t_c,dt,L1,L2

	!Lax-Friedrichs Scheme
	do i=1,points-2
		fi(i,2)=0.5*((u*fi(i+1,1))+(u*fi(i-1,1)))-(u*dt/(2*dx))*(fi(i+1,1)-fi(i-1,1))
		if (fi(i,1)< fi_a(i,1)) then
	   	   Linf(i,1) = fi_a(i,1) - fi(i,1)
		else
	   	   Linf(i,1) = fi(i,1) - fi_a (i,1)
		end if
		!Compute the errors
		L1 = dx*sum(abs(fi-fi_a))
		L2 = sqrt(dx*sum((fi-fi_a)**2))
		Linfmax = maxval(Linf(:,1))
		WRITE(31,*)t_c,Linfmax
	end do

	!current time update
	t_c=t_c+dt

	!update solution (next time level copied to current)
	fi(:,1)=fi(:,2)

	!Compute the errors
	L1 = dx*sum(abs(fi-fi_a))
	L2 = sqrt(dx*sum((fi-fi_a)**2))
	Linfmax = maxval(Linf(:,1))

	print*,t_c,dt,L1,L2,Linfmax
	if (t_c.ge.t_end)then
		exit		!exit loop
	end if

   end do

   print*,"Lax-Friedrichs Scheme is used"

   !write final solution
   open(31,file="final_LF_F1.dat", form="formatted",status="replace")
   DO I=0,POINTS-1
   WRITE(31,*)X(I), FI(I,1)
   END DO
   close(31)

else if (scheme == 4) then
   T_C=0.0			!current time
   open(31,file="Linf_LW_F1.dat", form="formatted",status="replace")

   do 	!start time loop

	dt=min(dt,t_end-t_c)
	print*,t_c,dt,L1,L2

	!Lax-Wendroff Scheme
	do i=1,points-2
		fi(i,2)=fi(i,1)-(u*dt/(2*dx))*(fi(i+1,1)-fi(i-1,1))+(u*u*dt*dt/(dx*dx*2))*(fi(i+1,1)-2*fi(i,1)+fi(i-1,1))
		if (fi(i,1)< fi_a(i,1)) then
	   	   Linf(i,1) = fi_a(i,1) - fi(i,1)
		else
	   	   Linf(i,1) = fi(i,1) - fi_a (i,1)
		end if
		!Compute the errors
		L1 = dx*sum(abs(fi-fi_a))
		L2 = sqrt(dx*sum((fi-fi_a)**2))
		Linfmax = maxval(Linf(:,1))
		WRITE(31,*)t_c,Linfmax
	end do

	!current time update
	t_c=t_c+dt

	!update solution (next time level copied to current)
	fi(:,1)=fi(:,2)
	
	!Compute the errors

	print*,t_c,dt,L1,L2,Linfmax
	if (t_c.ge.t_end)then
		exit		!exit loop
	end if

   end do

   print*,"Lax-Wendroff Scheme is used"

   !write final solution
   open(31,file="final_LW_F1.dat", form="formatted",status="replace")
   DO I=0,POINTS-1
   WRITE(31,*)X(I), FI(I,1)
   END DO
   close(31)

else if(scheme == 5) then
print*,"Using WENO Scheme"
   T_C=0.0			!current time
   open(31,file="L1_WENO_F1.dat", form="formatted",status="replace")

   do 	!start time loop
	dt=min(dt,t_end-t_c)
	print*,t_c,dt

   	do i=2,points-2
   	!compute the smoothness indicators
  	alpha(1) = FI(i+1,1)-FI(i,1)**2
  	alpha(2) = FI(i,1)-FI(i-1,1)**2
   	alpha(3) = FI(i+2,1)-FI(i-1,1)**2

   	!compute the smoothness indicators
   	beta(1) = (eps + alpha(1))**(-2)
   	beta(2) = (eps + alpha(2))**(-2)
   	beta(3) = (eps + alpha(3))**(-2)

   	!compute the linear weights
   	w(1) = beta(1) / (beta(1) + beta(2) + beta(3))
   	w(2) = beta(2) / (beta(1) + beta(2) + beta(3))
   	w(3) = beta(3) / (beta(1) + beta(2) + beta(3))

   	!Smoothness Measures
   	sum_w = w(1) + w(2) + w(3)
   	omega (1) = w(1) / sum_w
   	omega (2) = w(2) / sum_w
   	omega (3) = w(3) / sum_w

	!Compute the numerical flux using WENO scheme
	fi(i+0.5,2)=omega(1)*(fi(i+1,1)-fi(i,1))+omega(2)*(fi(i,1)-fi(i-1,1))+omega(3)*(fi(i-1,1)-fi(i-2,1))

	!Update the function
	fi(i+0.5,2)=0.5*(fi(i+0.5,2)+sign(1.0,fi(i+0.5,2))*sqrt(abs(fi(i+0.5,2)**2+eps)))

	!Compute the numerical flux using WENO Scheme with limited derivatives
	fi(i+1,2)=fi(i,1)+dx/dt*(fi(i+0.5,2)-fi(i-0.5,2))
        !Compute the errors
	L1 = dx*sum(abs(fi-fi_a))
	L2 = sqrt(dx*sum((fi-fi_a)**2))
	Linfmax = maxval(Linf(:,1))
	WRITE(31,*)t_c,L1
   	end do

   !current time update
   t_c=t_c+dt

   !update solution (next time level copied to current)
   fi(:,1)=fi(:,2)

   if (t_c.ge.t_end)then
	exit		!exit loop
	end if
	
end do
print*,"Using WENO Scheme"

!write final solution
   open(31,file="final_WENO_F1.dat", form="formatted",status="replace")
   DO I=0,POINTS-1
   WRITE(31,*)X(I), FI(I,1)
   END DO
   close(31)



else
   print*,"Invalid Scheme"
   stop
end if


!EXECUTION TIME CALCULATION
call CPU_TIME(end_time)
elapsed_time = end_time - start_time

write(*,*)"Elapsed time: ", elapsed_time, "seconds"




END PROGRAM LINEAR1
