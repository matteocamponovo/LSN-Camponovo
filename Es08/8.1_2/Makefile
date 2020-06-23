CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo.exe : Monte_Carlo.o random.o
	$(CC) random.o Monte_Carlo.o -o Monte_Carlo.exe
Monte_Carlo.o : Monte_Carlo.cpp Monte_Carlo.h
	$(CC) -c Monte_Carlo.cpp -o Monte_Carlo.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo.exe
