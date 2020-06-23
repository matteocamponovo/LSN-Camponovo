CC = g++
CFLAGS = -Wall -O3 --std=c++11

RWdiscreto.exe : RWdiscreto.o random.o
	$(CC) random.o RWdiscreto.o -o RWdiscreto.exe
RWdiscreto.o : RWdiscreto.cpp
	$(CC) -c RWdiscreto.cpp -o RWdiscreto.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o RWdiscreto.exe seed.out
