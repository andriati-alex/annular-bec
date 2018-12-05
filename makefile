
# ***************************** MAKEFILE ***************************** #

# $ make main_program



obj_linalg = array_memory.o 	  \
			 array_operations.o   \
			 matrix_operations.o  \
			 rmatrix_operations.o \
		   	 iterative_solver.o   \
			 tridiagonal_solver.o



obj_gp = $(obj_linalg) \
		 calculus.o    \
		 NewtonCG.o    \
		 rk4.o         \
		 linear_potential.o       \
		 GP_realtime_integrator.o \
		 GP_imagtime_integrator.o



obj_mctdhb = $(obj_linalg)           \
			 calculus.o              \
			 linear_potential.o      \
		 	 MCTDHB_configurations.o \
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



gp_header = $(linalg_header) 	\
			include/calculus.h  \
			include/NewtonCG.h  \
			include/rk4.h		\
		 	include/linear_potential.h       \
			include/GP_realtime_integrator.h \
			include/GP_imagtime_integrator.h



mctdhb_header = $(linalg_header) 	    	    \
				include/calculus.h		        \
		 		include/linear_potential.h      \
				include/MCTDHB_configurations.h \
				include/MCTDHB_datatype.h       \
				include/MCTDHB_observables.h    \
				include/MCTDHB_integrator.h




# Main program
# ************



GP_timepropagate : libgp.a exe/GP_timepropagate.c $(gp_header)
	icc -o GP_timepropagate exe/GP_timepropagate.c   \
		-L${MKLROOT}/lib/intel64                     \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		-lm -qopenmp                                 \
		-L./lib -I./include -lgp -O3





mu_steady : libgp.a exe/mu_steady.c $(gp_header)
	icc -o mu_steady exe/mu_steady.c                    \
		-L${MKLROOT}/lib/intel64                        \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core    \
		-lm -qopenmp                                    \
		-L./lib -I./include -lgp -O3





MCTDHB_time : libmctdhb.a exe/MCTDHB_time.c include/MCTDHB_integrator.h
	icc -o MCTDHB_time exe/MCTDHB_time.c			 \
		-L${MKLROOT}/lib/intel64                     \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		-lm -qopenmp                                 \
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

array_memory.o : src/array_memory.c
	icc -c -O3 src/array_memory.c





array_operations.o : src/array_operations.c
	icc -c -O3 -qopenmp src/array_operations.c





matrix_operations.o : src/matrix_operations.c
	icc -c -O3 -qopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp \
	src/matrix_operations.c





rmatrix_operations.o : src/rmatrix_operations.c
	icc -c -O3 -qopenmp src/rmatrix_operations.c





tridiagonal_solver.o : src/tridiagonal_solver.c
	icc -c -O3 -qopenmp src/tridiagonal_solver.c





iterative_solver.o : src/iterative_solver.c
	icc -c -O3 -qopenmp src/iterative_solver.c



linear_potential.o : src/linear_potential.c
	icc -c -O3 src/linear_potential.c





GP_realtime_integrator.o : src/GP_realtime_integrator.c
	icc -c -O3 -qopenmp                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		src/GP_realtime_integrator.c





GP_imagtime_integrator.o : src/GP_imagtime_integrator.c
	icc -c -O3 -qopenmp                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		src/GP_imagtime_integrator.c





rk4.o : src/rk4.c
	icc -c -O3 -qopenmp src/rk4.c





calculus.o : src/calculus.c
	icc -c -O3 -qopenmp                              \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		src/calculus.c





NewtonCG.o : src/NewtonCG.c
	icc -c -O3 -qopenmp src/NewtonCG.c





MCTDHB_integrator.o : src/MCTDHB_integrator.c
	icc -c -O3 -qopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp \
		src/MCTDHB_integrator.c





MCTDHB_configurations.o : src/MCTDHB_configurations.c
	icc -c -O3 -qopenmp src/MCTDHB_configurations.c





MCTDHB_observables.o : src/MCTDHB_observables.c
	icc -c -O3 -qopenmp src/MCTDHB_observables.c





MCTDHB_datatype.o : src/MCTDHB_datatype.c
	icc -c -O3 src/MCTDHB_datatype.c





clean :
	-rm build/*.o
	-rm lib/lib*
