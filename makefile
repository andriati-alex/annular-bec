
# ***************************** MAKEFILE ***************************** #

# $ make main_program



obj_linalg = inout.o              \
			 array_memory.o 	  \
			 array_operations.o   \
			 matrix_operations.o  \
		   	 iterative_solver.o   \
			 tridiagonal_solver.o



obj_gp = $(obj_linalg)         \
		 data_structure.o      \
		 calculus.o            \
		 rk4.o                 \
		 linear_potential.o    \
		 observables.o         \
		 realtime_integrator.o \
		 imagtime_integrator.o



linalg_header = include/inout.h              \
				include/array.h 			 \
				include/array_memory.h		 \
				include/array_operations.h   \
				include/matrix_operations.h  \
	  		    include/tridiagonal_solver.h \
				include/iterative_solver.h



gp_header = $(linalg_header) 	          \
			include/data_structure.h      \
			include/calculus.h            \
			include/rk4.h		          \
		 	include/linear_potential.h    \
			include/realtime_integrator.h \
			include/imagtime_integrator.h



   # ------------------------------------------------------------------ #

                         ###     EXECUTABLES     ###

   # ------------------------------------------------------------------ #



time_evolution : libgp.a exe/time_evolution.c $(gp_header)
	icc -o time_evolution exe/time_evolution.c -L${MKLROOT}/lib/intel64 \
		-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lm -qopenmp \
		-L./lib -I./include -lgp -O3





# Libraries to be linked
# ----------------------

libgp.a : $(obj_gp)
	ar rcs libgp.a $(obj_gp)
	mv libgp.a lib
	mv $(obj_gp) build





# Object files to the library
# ---------------------------

inout.o : src/inout.c
	icc -c -O3 -I./include src/inout.c



data_structure.o : src/data_structure.c
	icc -c -O3 -I./include src/data_structure.c



array_memory.o : src/array_memory.c
	icc -c -O3 -I./include src/array_memory.c



array_operations.o : src/array_operations.c
	icc -c -O3 -qopenmp -I./include src/array_operations.c



matrix_operations.o : src/matrix_operations.c
	icc -c -O3 -qopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp \
		-I./include src/matrix_operations.c



tridiagonal_solver.o : src/tridiagonal_solver.c
	icc -c -O3 -qopenmp -I./include src/tridiagonal_solver.c



iterative_solver.o : src/iterative_solver.c
	icc -c -O3 -qopenmp -I./include src/iterative_solver.c



linear_potential.o : src/linear_potential.c
	icc -c -O3 -I./include src/linear_potential.c



rk4.o : src/rk4.c
	icc -c -O3 -qopenmp -I./include src/rk4.c



calculus.o : src/calculus.c
	icc -c -O3 -qopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		-I./include src/calculus.c



observables.o : src/observables.c
	icc -c -O3 -qopenmp -I./include src/observables.c



realtime_integrator.o : src/realtime_integrator.c
	icc -c -O3 -qopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		-I./include src/realtime_integrator.c



imagtime_integrator.o : src/imagtime_integrator.c
	icc -c -O3 -qopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
		-I./include src/imagtime_integrator.c



clean :
	-rm build/*.o
	-rm lib/lib*
	-rm time_evolution
