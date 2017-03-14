CFLAGS=-std=c99 -w
OPT=

NW_PATH=seq-align
NW_LIBS_PATH=$(NW_PATH)/libs
NW_INCS=-I $(NW_LIBS_PATH) -I $(NW_PATH)/src
NW_LIBS=-L $(NW_LIBS_PATH)/string_buffer -L $(NW_PATH)/src
NW_LINK=-lalign -lstrbuf -lpthread -lz
NW_DEPS=$(NW_PATH)/src/libalign.a

PQ_PATH=libpqueue
PQ_INCS=-I $(PQ_PATH)/src
PQ_COMP=$(PQ_PATH)/src/pqueue.c

INCS=$(NW_INCS) $(PQ_INCS) -I ./src
LIBS=$(NW_LIBS)
LINK=$(NW_LINK)
DEPS=$(NW_DEPS)
COMP=$(PQ_COMP) src/msa.c

all: nw_build msa

nw_build:
	make -C $(NW_PATH)

msa: main.c src/msa.c $(DEPS)
	$(CC) $< $(COMP) -o $@ $(CFLAGS) $(INCS) $(LIBS) $(LINK) $(OPT)

nw_clean:
	make clean -C $(NW_PATH)

clean: nw_clean
	$(RM) *.o msa

.PHONY: all clean
