
CFLAGS= -fPIC -O3 

all: librevolve.a

python: revolve.cpp wrap_revolve.cpp revolve_python.c
	python setup.py build
	@echo "Run sudo python setup.py install to install the revolve module."

checkpoints: librevolve.a checkpoints.o
	g++ $(CFLAGS) -o checkpoints checkpoints.o librevolve.a 

checkpoints.o: checkpoints.c
	gcc -c $(CFLAGS) checkpoints.c

librevolve.a: revolve.cpp revolve.h wrap_revolve.cpp
	g++ -c $(CFLAGS) revolve.cpp
	g++ -c $(CFLAGS) wrap_revolve.cpp
	ar -rs librevolve.a  revolve.o wrap_revolve.o
clean: 
	rm -rf checkpoints librevolve.a *.o build

.PHONY: clean
