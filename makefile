
# ***************************** MAKEFILE ***************************** #

# $ make main_program

obj_linalg = array_memory.o 	  \
			 array_operations.o   \
			 matrix_operations.o  \
			 rmatrix_operations.o \
		   	 iterative_solver.o   \
			 tridiagonal_solver.o

obj_gp = $(obj_linalg)   \
		 time_routine.o  \
		 itime_routine.o \
		 calculus.o      \
		 NewtonCG.o      \
		 rk4.o

obj_mctdhb = $(obj_linalg)           \
		 	 MCTDHB_configurations.o \
			 calculus.o	             \
			 MCTDHB_datatype.o       \
			 MCTDHB_observables.o    \
		 	 MCTDHB_integrator.o

linalg_header = include/array.h 			 \
				include/array_memory.h		 \
				include/array_operations.h   \
				include/matrix_operations.h  \
				include/rmatrix_operations.h \
	  		    include/tridiagonal_solver.h \
				include/iterative_solver.h

gp_header = $(linalg_header) 		\
			include/time_routine.h 	\
			include/itime_routine.h \
			include/calculus.h      \
			include/NewtonCG.h      \
			include/rk4.h

mctdhb_header = $(linalg_header) 	    	    \
				include/MCTDHB_configurations.h \
				include/calculus.h		        \
				include/MCTDHB_datatype.h       \
				include/MCTDHB_observables.h    \
				include/MCTDHB_integrator.h




# Main program
# ************

time_evolution : libgp.a exe/time_evolution.c $(gp_header)
	gcc -o time_evolution exe/time_evolution.c          \
		-L${MKLROOT}/lib/intel64                        \
		-Wl,--no-as-needed                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core    \
		-lm -fopenmp                                    \
		-L./lib -I./include -lgp -O3

itime_propagate : libgp.a exe/itime_propagate.c $(gp_header)
	gcc -o itime_propagate exe/itime_propagate.c        \
		-L${MKLROOT}/lib/intel64                        \
		-Wl,--no-as-needed                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core    \
		-lm -fopenmp                                    \
		-L./lib -I./include -lgp -O3

mu_steady : libgp.a exe/mu_steady.c $(gp_header)
	gcc -o mu_steady exe/mu_steady.c                    \
		-L${MKLROOT}/lib/intel64                        \
		-Wl,--no-as-needed                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core    \
		-lm -fopenmp                                    \
		-L./lib -I./include -lgp -O3

MCTDHB_time : libmctdhb.a exe/MCTDHB_time.c include/MCTDHB_integrator.h
	gcc -o MCTDHB_time exe/MCTDHB_time.c				\
		-L${MKLROOT}/lib/intel64                        \
		-Wl,--no-as-needed                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core    \
		-lm -fopenmp                                    \
		-L./lib -lmctdhb -O3




# Libraries to be linked
# **********************

libgp.a : $(obj_gp)
	ar rcs libgp.a $(obj_gp)
	mv libgp.a lib
	mv $(obj_gp) build

libmctdhb.a : $(obj_mctdhb)
	ar rcs libmctdhb.a $(obj_mctdhb)
	mv libmctdhb.a lib
	mv $(obj_mctdhb) build





# Object files to the library
# ***************************

array_memory.o : src/array_memory.c include/array_memory.h
	gcc -c -O3 src/array_memory.c





array_operations.o : src/array_operations.c include/array_operations.h
	gcc -c -O3 -fopenmp src/array_operations.c





matrix_operations.o : src/matrix_operations.c include/matrix_operations.h
	gcc -c -O3 -fopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp \
	src/matrix_operations.c





rmatrix_operations.o : src/rmatrix_operations.c include/rmatrix_operations.h
	gcc -c -O3 -fopenmp src/rmatrix_operations.c





tridiagonal_solver.o : src/tridiagonal_solver.c include/tridiagonal_solver.h
	gcc -c -O3 -fopenmp src/tridiagonal_solver.c





iterative_solver.o : src/iterative_solver.c	include/iterative_solver.h
	gcc -c -O3 -fopenmp src/iterative_solver.c





time_routine.o : src/time_routine.c include/time_routine.h
	gcc -c -O3 -fopenmp                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		src/time_routine.c





itime_routine.o : src/itime_routine.c include/itime_routine.h
	gcc -c -O3 -fopenmp                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		src/itime_routine.c


rk4.o : src/rk4.c include/rk4.h
	gcc -c -O3 -fopenmp src/rk4.c


calculus.o : src/calculus.c include/calculus.h
	gcc -c -O3 -fopenmp                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		src/calculus.c

NewtonCG.o : src/NewtonCG.c include/NewtonCG.h
	gcc -c -O3 -fopenmp src/NewtonCG.c

MCTDHB_integrator.o : src/MCTDHB_integrator.c include/MCTDHB_integrator.h
	gcc -c -O3 -fopenmp src/MCTDHB_integrator.c

MCTDHB_configurations.o : src/MCTDHB_configurations.c include/MCTDHB_configurations.h
	gcc -c -O3 -fopenmp src/MCTDHB_configurations.c

MCTDHB_observables.o : src/MCTDHB_observables.c include/MCTDHB_observables.h
	gcc -c -O3 -fopenmp src/MCTDHB_observables.c

MCTDHB_datatype.o : src/MCTDHB_datatype.c include/MCTDHB_datatype.h
	gcc -c -O3 src/MCTDHB_datatype.c

clean :
	-rm build/*.o
	-rm lib/lib*
