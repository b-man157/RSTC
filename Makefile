all: a.out

a.out: TI.cpp
	g++ -g TI.cpp

clean:
	rm a.out
