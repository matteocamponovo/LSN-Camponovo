CC = g++
CFLAGS = -Wall -O3 --std=c++11

integrale2.exe : integrale2.o random.o
	$(CC) random.o integrale2.o -o integrale2.exe
integrale2.o : integrale2.cpp
	$(CC) -c integrale2.cpp -o integrale2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o integrale2.exe seed.out
