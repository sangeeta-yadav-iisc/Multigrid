CC = g++ -Wall
LIBS = -lm -lblas
problem1.exe:problem1.o SOR.o csr_mv.o norm.o
	$(CC) $^ -o $@ $(LIBS)
problem1.o:problem1.cpp  *.h
	$(CC) -c $^ $(LIBS)
SOR.o:SOR.cpp *.h
	$(CC) -c $< $(LIBS) 
csr_mv.o: csr_mv.cpp *.h
	$(CC) -c $< $(LIBS)
norm.o: norm.cpp *.h
	$(CC) -c $< $(LIBS)
clean:
	rm -f *.o a.out *.gch *.exe
scrub: clean
	rm -f *.dat


