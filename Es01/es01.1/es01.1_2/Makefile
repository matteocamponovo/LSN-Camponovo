CC = g++
CFLAGS = -Wall -O3 --std=c++11

sigmamedio.exe : sigmamedio.o random.o
	$(CC) random.o sigmamedio.o -o sigmamedio.exe
sigmamedio.o : sigmamedio.cpp
	$(CC) -c sigmamedio.cpp -o sigmamedio.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o sigmamedio.exe seed.out
