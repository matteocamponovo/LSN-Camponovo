CC = g++
CFLAGS = -Wall -O3 --std=c++11

rmedio.exe : rmedio.o random.o
	$(CC) random.o rmedio.o -o rmedio.exe
rmedio.o : rmedio.cpp
	$(CC) -c rmedio.cpp -o rmedio.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o rmedio.exe seed.out
