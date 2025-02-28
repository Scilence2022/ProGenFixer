CFLAGS=-g -Wall -Ofast 
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz -lm
PROG=QuickVar kcov

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)

QuickVar:QuickVar.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ QuickVar.c kthread.c $(LIBS) -lpthread

kcov:kcov.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ kcov.c kthread.c $(LIBS) -lpthread

clean:
	rm -fr *.dSYM $(PROG)

