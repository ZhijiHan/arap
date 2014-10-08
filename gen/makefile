all: model_gen

model_gen: main.o modelgenerator.o
	g++ main.o modelgenerator.o -o model_gen

main.o: main.cc
	g++ -c main.cc

modelgenerator.o: modelgenerator.cc
	g++ -c modelgenerator.cc

clean:
	rm -rf *o model_gen
