all: test1

test1: test1.o
	g++ test1.o GePUP.o GePUP_IMEX.o amgSolver.o matrix.o sparseMatrix.o avl.o norm.o TimeFunction2D.o -o test1 -O2 -O3 -Ofast

test1.o: test1.cpp GePUP_IMEX.o
	g++ -c test1.cpp -Iinclude -O2 -O3 -Ofast

GePUP.o: include/GePUP.h src/GePUP.cpp include/idpair.h include/Field.h amgSolver.o matrix.o sparseMatrix.o avl.o norm.o TimeFunction2D.o
	g++ -c src/GePUP.cpp -Iinclude -O2 -O3 -Ofast

GePUP_IMEX.o: include/GePUP_IMEX.h src/GePUP_IMEX.cpp GePUP.o
	g++ -c src/GePUP_IMEX.cpp -Iinclude -O2 -O3 -Ofast

amgSolver.o: include/amgSolver.h src/amgSolver.cpp sparseMatrix.o avl.o
	g++ -c src/amgSolver.cpp -Iinclude -O2 -O3 -Ofast

matrix.o: include/matrix.h src/matrix.cpp
	g++ -c src/matrix.cpp -Iinclude -O2 -O3 -Ofast

sparseMatrix.o: include/sparseMatrix.h src/sparseMatrix.cpp matrix.o
	g++ -c src/sparseMatrix.cpp -Iinclude -O2 -O3 -Ofast

avl.o: include/avl.h src/avl.cpp
	g++ -c src/avl.cpp -Iinclude -O2 -O3 -Ofast

norm.o: include/norm.h src/norm.cpp
	g++ -c src/norm.cpp -Iinclude -O2 -O3 -Ofast

TimeFunction2D.o: include/TimeFunction2D.h src/TimeFunction2D.cpp
	g++ -c src/TimeFunction2D.cpp -Iinclude -O2 -O3 -Ofast

clean:
	rm *.o
