all: a.out

a.out: RSTC.cpp
	g++ -g RSTC.cpp -o RSTC

clean:
	rm RSTC
