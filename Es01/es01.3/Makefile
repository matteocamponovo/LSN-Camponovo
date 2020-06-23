CC = g++
CFLAGS = -Wall -O3 --std=c++11

buffon.exe : buffon.o random.o
	$(CC) random.o buffon.o -o buffon.exe
buffon.o : buffon.cpp
	$(CC) -c buffon.cpp -o buffon.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o buffon.exe seed.out
