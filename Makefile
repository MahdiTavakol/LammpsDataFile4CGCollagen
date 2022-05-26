ld4cc:
	g++ -c main.cpp Atoms.cpp Bonds.cpp Angles.cpp
	g++ -o ld4cc Atoms.o Bonds.o Angles.o main.o
	rm -rf *.o
clean:
	rm -rf *.o
