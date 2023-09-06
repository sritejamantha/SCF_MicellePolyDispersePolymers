CFLAGS=-O3 -g "-I/fs01/ma03/homes/smantha/Softwares/GSL_Install/include" -Wall 
LDLIBS =-O3 -g "-L/fs01/ma03/homes/smantha/Softwares/GSL_Install/lib" -lgsl -lgslcblas -lm 
CC=gcc

SRC=$(wildcard *.c)

OBJ=$(patsubst %c,%o,$(SRC) )

out: $(OBJ)
	gcc $^ -o $@ $(LDLIBS)

include $(SRC:.c=.d)

%.d: %.c
	set -e; rm -f $@; \
	gcc -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@:,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

clean:
	rm -f *.o *.d out
