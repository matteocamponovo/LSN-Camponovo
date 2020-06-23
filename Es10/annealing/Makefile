CC = g++
CFLAGS = -Wall -O3 --std=c++11

es10.1.exe : es10.1.o random.o Maps.o Chrom.o
	$(CC) random.o Maps.o Chrom.o es10.1.o -o es10.1.exe

es10.1.o : es10.1.cpp es10.1.h
	$(CC) -c es10.1.cpp -o es10.1.o $(CFLAGS)

Chrom.o : Chrom.cpp Chrom.h 
	$(CC) -c Chrom.cpp -o Chrom.o $(CFLAGS)

Maps.o : Maps.cpp Maps.h 
	$(CC) -c Maps.cpp -o Maps.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es10.1.exe
