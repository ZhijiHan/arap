all: rotate_bar stretch_bar

rotate_bar: rotate_bar.o
	g++ rotate_bar.o -o rotate_bar

rotate_bar.o: rotate_bar.cc
	g++ -c rotate_bar.cc

stretch_bar: stretch_bar.o
	g++ stretch_bar.o -o stretch_bar

stretch_bar.o: stretch_bar.cc
	g++ -c stretch_bar.cc

clean:
	rm -rf *o rotate_bar stretch_bar
