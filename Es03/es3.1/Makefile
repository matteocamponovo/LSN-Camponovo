CC = g++
CFLAGS = -Wall -O3 --std=c++11

GBM1.exe : GBM1.o random.o
	$(CC) random.o GBM1.o -o GBM1.exe
GBM1.o : GBM1.cpp
	$(CC) -c GBM1.cpp -o GBM1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o GBM1.exe seed.out
