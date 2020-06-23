CC = g++
CFLAGS = -Wall -O3 --std=c++11

integrale1.exe : integrale1.o random.o
	$(CC) random.o integrale1.o -o integrale1.exe
integrale1.o : integrale1.cpp
	$(CC) -c integrale1.cpp -o integrale1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o integrale1.exe seed.out
