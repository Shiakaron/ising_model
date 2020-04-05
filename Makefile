CC = g++

output: main.o bootstrap.o general.o metropolis.o wolff.o autocorrelation.o
	$(CC) -o output main.o bootstrap.o general.o metropolis.o wolff.o autocorrelation.o

main.o: main.cpp
	$(CC) -c main.cpp

bootstrap.o: bootstrap.cpp
	$(CC) -c bootstrap.cpp

autocorrelation.o: autocorrelation.cpp
	$(CC) -c autocorrelation.cpp

metropolis.o: metropolis.cpp
	$(CC) -c metropolis.cpp

wolff.o: wolff.cpp
	$(CC) -c wolff.cpp

general.o: general.cpp
	$(CC) -c general.cpp

clean:
	del *.o output
