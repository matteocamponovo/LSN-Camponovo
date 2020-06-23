CC = g++
CFLAGS = -Wall -O3 --std=c++11

es9.exe : es9.o random.o Maps.o Chrom.o
	$(CC) random.o Maps.o Chrom.o es9.o -o es9.exe

es9.o : es9.cpp es9.h
	$(CC) -c es9.cpp -o es9.o $(CFLAGS)

Chrom.o : Chrom.cpp Chrom.h 
	$(CC) -c Chrom.cpp -o Chrom.o $(CFLAGS)

Maps.o : Maps.cpp Maps.h 
	$(CC) -c Maps.cpp -o Maps.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es9.exe
