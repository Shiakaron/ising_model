CC = g++

output: main.o bootstrap.o general.o metropolis.o wolff.o autocorrelation.o
	$(CC) -o output main.o bootstrap.o general.o metropolis.o wolff.o autocorrelation.o

main.o: main.cpp header.h
	$(CC) -c main.cpp

bootstrap.o: bootstrap.cpp header.h
	$(CC) -c bootstrap.cpp

autocorrelation.o: autocorrelation.cpp header.h
	$(CC) -c autocorrelation.cpp

metropolis.o: metropolis.cpp header.h
	$(CC) -c metropolis.cpp

wolff.o: wolff.cpp header.h
	$(CC) -c wolff.cpp

general.o: general.cpp header.h
	$(CC) -c general.cpp

clean:
	del *.o output
