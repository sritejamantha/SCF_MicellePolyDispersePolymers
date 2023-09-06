CFLAGS= #-ansi -pedantic -Wshadow -Wundef -Wcast-align -Wpointer-arith -O3 -D_XOPEN_SOURCE=600
LDLIBS= -lm -Wall
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
