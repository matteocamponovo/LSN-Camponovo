CC = g++
CFLAGS = -Wall -O3 --std=c++11

chiquadro.exe : chiquadro.o random.o
	$(CC) random.o chiquadro.o -o chiquadro.exe
chiquadro.o : chiquadro.cpp
	$(CC) -c chiquadro.cpp -o chiquadro.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o chiquadro.exe seed.out
