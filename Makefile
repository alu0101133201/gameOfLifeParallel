FC := mpifort
# FC_FLAGS := 

EXECUTABLE := main
TESTDIR := tests

all: $(EXECUTABLE)

run: clean all
		clear
		mpirun -np 12 --oversubscribe $(EXECUTABLE) < $(TESTDIR)/diehard.in

$(EXECUTABLE): ./game_of_life_parallel.f90
	$(FC) $^ -o $@

clean:
	rm -f $(EXECUTABLE)
