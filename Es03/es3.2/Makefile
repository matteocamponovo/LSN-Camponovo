CC = g++
CFLAGS = -Wall -O3 --std=c++11

GBM2.exe : GBM2.o random.o
	$(CC) random.o GBM2.o -o GBM2.exe
GBM2.o : GBM2.cpp
	$(CC) -c GBM2.cpp -o GBM2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o GBM2.exe seed.out
