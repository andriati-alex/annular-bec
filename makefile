# run in command line: $ make

objects =  array_memory.o system_solvers.o matrix_operations.o \
		   array_operations.o

#client : libnavalbattle.so client.c datastructure.h navallib.h wait.h
#	gcc -o client client.c -L$(HOME)/NavalBattle/lib -lnavalbattle

# Prepara a biblioteca para compilar o client
libAAlinalg.so : $(objects)
	gcc -o libAAlinalg.so -shared $(objects)
	export LD_LIBRARY_PATH=$(HOME)/AndriatiLibrary/linear-algebra:$$LD_LIBRARY_PATH
#	mkdir ~/NavalBattle
#	mkdir ~/NavalBattle/lib
#	mv libnavalbattle.so ~/NavalBattle/lib

array_memory.o : array_memory.c array_memory.h
	gcc -c -fPIC array_memory.c

array_operations.o : array_operations.c array_operations.h
	gcc -c -fopenmp -fPIC array_operations.c

matrix_operations.o : matrix_operations.c matrix_operations.h
	gcc -c -fPIC matrix_operations.c

system_solvers.o : system_solvers.c system_solvers.h \
	array_operations.o matrix_operations.o
	gcc -c -fPIC system_solvers.c

# O traco antes do rm (-rm) indica que o comando sera executado mesmo
# na presenca de erros.
clean :
	-rm $(objects)
	-rm libAAlinalg.so
#	-rm ~/NavalBattle/lib/libnavalbattle.so
#	-rmdir ~/NavalBattle/lib
#	-rmdir ~/NavalBattle
