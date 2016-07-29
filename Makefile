CC=g++

output: MainGradient.o ELTensors.o
	$(CC) MainGradient.o ELTensors.o -o output
	./output

MainGradient.o: MainGradient.cpp EmbedInitial.h CSVinto4Darray.h CSVintomatrix.h CSVintoVector.h ELInitialization.h EmbedInitial.h Volumefraction.h Volfnminmax.h ELTensors.h AhatInitializing.h
	$(CC) -c MainGradient.cpp 

ELTensors.o: ELTensors.h ELTensors.cpp ELInitialization.h
	$(CC) -c ELTensors.cpp


clean:
	rm -rf *.o
