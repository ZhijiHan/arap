all: rotate_bar

rotate_bar: rotate_bar.o
	g++ rotate_bar.o -o rotate_bar

rotate_bar.o: rotate_bar.cc
	g++ -c rotate_bar.cc

clean:
	rm -rf *o rotate_bar
