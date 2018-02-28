# run in command line: $ make

objects =  array_memory.o system_solvers.o matrix_operations.o \
		   array_operations.o

TestLarge : libAAlinalg.a TestLarge.c array_operations.h system_solvers.h \
			matrix_operations.h array_memory.h
	icc -o TestLarge -static -O3 TestLarge.c \
		-L$(HOME)/AndriatiLibrary/linear-algebra -lAAlinalg

# Prepara a biblioteca para compilar o client
libAAlinalg.a : $(objects)
	ar rcs libAAlinalg.a $(objects)
#	mkdir ~/NavalBattle
#	mkdir ~/NavalBattle/lib
#	mv libnavalbattle.so ~/NavalBattle/lib

array_memory.o : array_memory.c array_memory.h
	icc -c -O3 -fPIC array_memory.c

array_operations.o : array_operations.c array_operations.h
	icc -c -O3 -fPIC array_operations.c

matrix_operations.o : matrix_operations.c matrix_operations.h
	icc -c -O3 -fPIC matrix_operations.c

system_solvers.o : system_solvers.c system_solvers.h \
	array_operations.o matrix_operations.o
	icc -c -O3 -fPIC system_solvers.c

# O traco antes do rm (-rm) indica que o comando sera executado mesmo
# na presenca de erros.
clean :
	-rm $(objects)
	-rm libAAlinalg.a
	-rm TestLarge
#	-rm ~/NavalBattle/lib/libnavalbattle.so
#	-rmdir ~/NavalBattle/lib
#	-rmdir ~/NavalBattle
