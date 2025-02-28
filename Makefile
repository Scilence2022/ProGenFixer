CFLAGS=-g -Wall -Ofast 
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz -lm
PROG=kcov ProGenFixer

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)

kcov:kcov.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ kcov.c kthread.c $(LIBS) -lpthread

ProGenFixer:ProGenFixer.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ ProGenFixer.c kthread.c $(LIBS) -lpthread

clean:
	rm -fr *.dSYM $(PROG)

