CC=g++

output: MainGradient.o ELTensors.o
	$(CC) MainGradient.o ELTensors.o -o output
	./output

MainGradient.o: MainGradient.cpp EmbedInitial.h CSVinto4Darray.h CSVintomatrix.h CSVintoVector.h ELInitialization.h EmbedInitial.h Volumefraction.h Volfnminmax.h 
	$(CC) -c MainGradient.cpp EmbedInitial.h CSVinto4Darray.h CSVintomatrix.h CSVintoVector.h ELInitialization.h EmbedInitial.h Volumefraction.h Volfnminmax.h 

ELTensors.o: ELTensors.h ELTensors.cpp ELGradient.cpp ELDerivTensors.cpp
	$(CC) -c ELTensors.h ELTensors.cpp ELGradient.cpp ELDerivTensors.cpp
