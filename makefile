CC = gcc
LINK = gcc
AR = ar

ARFLAGS = -rv
CFLAGS_DBG = -DFP_TYPE=float -Og -g -std=c99 -Wall -fPIC -fopenmp
CFLAGS_REL = -DFP_TYPE=float -Ofast -std=c99 -Wall -fPIC -fopenmp

OUT_DIR = ./build
OBJS = $(OUT_DIR)/container.o \
  $(OUT_DIR)/frame.o \
  $(OUT_DIR)/dsputils.o \
  $(OUT_DIR)/layer0.o \
  $(OUT_DIR)/layer1.o
TARGET_A = $(OUT_DIR)/libllsm.a

CIGLET_DIR = ./external/ciglet
CIGLET_SRC = $(CIGLET_DIR)/ciglet.c
CIGLET_O   = $(OUT_DIR)/ciglet.o
LIBGVPS_DIR = ./external/libgvps
LIBGVPS_A = ./external/libgvps/build/libgvps.a
LIBPYIN_DIR = ./external/libpyin
LIBPYIN_A = ./external/libpyin/build/libpyin.a

default: $(TARGET_A)

test: $(OUT_DIR)/test-structs \
	  $(OUT_DIR)/test-dsputils \
	  $(OUT_DIR)/test-harmonic \
	  $(OUT_DIR)/test-layer0-anasynth \
	  $(OUT_DIR)/test-layer0-edgecase \
	  $(OUT_DIR)/test-layer1-anasynth
	$(OUT_DIR)/test-structs
	$(OUT_DIR)/test-dsputils
	$(OUT_DIR)/test-harmonic
	$(OUT_DIR)/test-layer0-anasynth pp
	$(OUT_DIR)/test-layer0-anasynth czt
	$(OUT_DIR)/test-layer0-edgecase
	$(OUT_DIR)/test-layer1-anasynth

test-layer0: $(OUT_DIR)/test-layer0-anasynth \
	  $(OUT_DIR)/test-layer0-edgecase
	$(OUT_DIR)/test-layer0-anasynth pp
	$(OUT_DIR)/test-layer0-anasynth czt
	$(OUT_DIR)/test-layer0-edgecase

test-layer1: $(OUT_DIR)/test-layer1-anasynth
	$(OUT_DIR)/test-layer1-anasynth

test-dsputils: $(OUT_DIR)/test-structs \
	  $(OUT_DIR)/test-dsputils \
	  $(OUT_DIR)/test-harmonic
	$(OUT_DIR)/test-structs
	$(OUT_DIR)/test-dsputils
	$(OUT_DIR)/test-harmonic

$(OUT_DIR)/test-%: test/test-%.c $(TARGET_A) \
	  $(LIBPYIN_A) $(LIBGVPS_A) $(CIGLET_O) test/verify-utils.h
	$(CC) $(CFLAGS_DBG) -o $(OUT_DIR)/test-$* test/test-$*.c \
	  $(TARGET_A) $(LIBPYIN_A) $(LIBGVPS_A) $(CIGLET_O) -lm

$(TARGET_A): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET_A) $(OBJS)

$(OUT_DIR)/frame.o: llsm.h
$(OUT_DIR)/container.o: llsm.h
$(OUT_DIR)/dsputils.o: dsputils.h llsm.h
$(OUT_DIR)/layer0.o: dsputils.h llsm.h
$(OUT_DIR)/layer1.o: dsputils.h llsm.h

$(CIGLET_O): $(CIGLET_SRC)
	mkdir -p $(OUT_DIR)
	$(CC) $(CFLAGS_REL) -o $(CIGLET_O) -c $(CIGLET_SRC)

$(OUT_DIR)/%.o: %.c
	mkdir -p $(OUT_DIR)
	$(CC) $(CFLAGS_DBG) -o $(OUT_DIR)/$*.o -c $*.c

clean:
	rm -f build/*
