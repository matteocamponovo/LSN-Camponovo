CC = g++
CFLAGS = -Wall -O3 --std=c++11

RWcontinuo.exe : RWcontinuo.o random.o
	$(CC) random.o RWcontinuo.o -o RWcontinuo.exe
RWcontinuo.o : RWcontinuo.cpp
	$(CC) -c RWcontinuo.cpp -o RWcontinuo.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o RWcontinuo.exe seed.out
