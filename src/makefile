LDIR=/usr/lib/

example: example.o hyp.o mpreal.o dlmalloc.o
	g++ -ggdb example.o mpreal.o dlmalloc.o hyp.o $(LDIR)libmpfr.a $(LDIR)libgmp.a -o example

mpreal.o: mpreal.cpp mpreal.h
	g++ -ggdb -c mpreal.cpp

example.o: example.cpp mpreal.h
	g++ -ggdb -c example.cpp

hyp.o: hyp.cpp hyp.h
	g++ -ggdb -c hyp.cpp

dlmalloc.o: dlmalloc.c dlmalloc.h
	g++ -c dlmalloc.c

clean:
	rm -f *.o 
