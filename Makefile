CC=g++

output: main.o GameMaster.o Player.o Field.o Deck.o RandomPlayer.o LogWriter.o
	$(CC) main.o GameMaster.o Player.o Field.o Deck.o RandomPlayer.o LogWriter.o -o output
	./output
