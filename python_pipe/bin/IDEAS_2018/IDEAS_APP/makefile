ideas:	main.cpp tensorHMMbase.cpp tensorHMMbase.h genomicTensor.cpp genomicTensor.h MixGauss.cpp MixGauss.h datastructure.h
	g++  -I. -I/usr/local/include/gsl -c main.cpp tensorHMMbase.cpp genomicTensor.cpp MixGauss.cpp
	g++ -L/usr/local/lib main.o tensorHMMbase.o genomicTensor.o MixGauss.o -lgsl -lgslcblas -lm -o ideas
