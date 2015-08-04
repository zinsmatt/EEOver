CC      = g++
CFLAGS  = -Wall -g 
LDFLAGS = -lgsl -lgslcblas

OBJ = ellipse_ellipse_overlap.o call_ee.o gsl_poly.o zsolve_quartic.o

prog: $(OBJ)
	$(CC) $(CFLAGS) -o overlap $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ)
