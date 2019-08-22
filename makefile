PREFIX = /usr
CC = gcc
LINK = gcc
AR = ar

OUT_DIR = ./build
OBJS = $(OUT_DIR)/container.o \
  $(OUT_DIR)/frame.o \
  $(OUT_DIR)/dsputils.o \
  $(OUT_DIR)/llsmutils.o \
  $(OUT_DIR)/layer0.o \
  $(OUT_DIR)/layer1.o \
  $(OUT_DIR)/coder.o \
  $(OUT_DIR)/llsmrt.o
TARGET_A = $(OUT_DIR)/libllsm.a

CIGLET_DIR = ./external/ciglet
CIGLET_SRC = $(CIGLET_DIR)/ciglet.c
CIGLET_O   = $(OUT_DIR)/ciglet.o
LIBCIGLET_A = # use it if building entire libciglet is more appropriate
LIBGVPS_DIR = ./external/libgvps
LIBGVPS_A = ./external/libgvps/build/libgvps.a
LIBGVPS_INCLUDE = ./external/libgvps/
LIBPYIN_DIR = ./external/libpyin
LIBPYIN_A = ./external/libpyin/build/libpyin.a
LIBPYIN_INCLUDE = ./external/libpyin/
LIBNEBULA_DIR = ./external/libnebula
LIBNEBULA_A = ./external/libnebula/build/libnebula.a
LIBIHNM_DIR = ./external/libihnm
LIBIHNM_A = ./external/libihnm/build/libihnm.a

ARFLAGS = -rv
CFLAGS_COMMON = -DFP_TYPE=float -std=c99 -Wall -fPIC -pthread -DUSE_PTHREAD \
	-I$(LIBGVPS_INCLUDE) -I$(LIBPYIN_INCLUDE)
CFLAGS_DBG = $(CFLAGS_COMMON) -Og -g
CFLAGS_REL = $(CFLAGS_COMMON) -Ofast
CFLAGS = $(CFLAGS_DBG)

default: $(TARGET_A)

test: $(OUT_DIR)/test-structs \
	  $(OUT_DIR)/test-dsputils \
	  $(OUT_DIR)/test-harmonic \
	  $(OUT_DIR)/test-layer0-anasynth \
	  $(OUT_DIR)/test-layer0-edgecase \
	  $(OUT_DIR)/test-layer1-anasynth \
	  $(OUT_DIR)/test-llsmrt \
	  $(OUT_DIR)/test-coder
	$(OUT_DIR)/test-structs
	$(OUT_DIR)/test-dsputils
	$(OUT_DIR)/test-harmonic
	$(OUT_DIR)/test-layer0-anasynth pp
	$(OUT_DIR)/test-layer0-anasynth czt
	$(OUT_DIR)/test-layer0-edgecase
	$(OUT_DIR)/test-layer1-anasynth
	$(OUT_DIR)/test-llsmrt
	$(OUT_DIR)/test-coder

test-layer0: $(OUT_DIR)/test-layer0-anasynth \
	  $(OUT_DIR)/test-layer0-edgecase
	$(OUT_DIR)/test-layer0-anasynth pp
	$(OUT_DIR)/test-layer0-anasynth czt
	$(OUT_DIR)/test-layer0-edgecase

test-layer1: $(OUT_DIR)/test-layer1-anasynth
	$(OUT_DIR)/test-layer1-anasynth

test-coder: $(OUT_DIR)/test-coder
	$(OUT_DIR)/test-coder

test-dsputils: $(OUT_DIR)/test-structs \
	  $(OUT_DIR)/test-dsputils \
	  $(OUT_DIR)/test-harmonic
	$(OUT_DIR)/test-structs
	$(OUT_DIR)/test-dsputils
	$(OUT_DIR)/test-harmonic

test-llsmrt: $(OUT_DIR)/test-llsmrt
	$(OUT_DIR)/test-llsmrt

test-pbpeffects: $(OUT_DIR)/test-pbpeffects
	$(OUT_DIR)/test-pbpeffects

$(OUT_DIR)/test-structs: buffer.h

$(OUT_DIR)/test-%: test/test-%.c $(TARGET_A) \
	  $(LIBPYIN_A) $(LIBGVPS_A) $(CIGLET_O) test/verify-utils.h
	$(CC) $(CFLAGS) -o $(OUT_DIR)/test-$* test/test-$*.c \
	  $(TARGET_A) $(LIBPYIN_A) $(LIBGVPS_A) $(LIBIHNM_A) $(CIGLET_O) $(LIBCIGLET_A) -lm

$(OUT_DIR)/demo-%: test/demo-%.c $(TARGET_A) \
	  $(LIBPYIN_A) $(LIBNEBULA_A) $(LIBGVPS_A) $(LIBIHNM_A) $(CIGLET_O)
	$(CC) $(CFLAGS) -o $(OUT_DIR)/demo-$* test/demo-$*.c -fopenmp \
	  $(TARGET_A) $(LIBPYIN_A) $(LIBNEBULA_A) $(LIBGVPS_A) $(LIBIHNM_A) $(CIGLET_O) $(LIBCIGLET_A) -lm

$(TARGET_A): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET_A) $(OBJS)

$(OUT_DIR)/frame.o: llsm.h dsputils.h
$(OUT_DIR)/container.o: llsm.h
$(OUT_DIR)/dsputils.o: dsputils.h llsm.h
$(OUT_DIR)/llsmutils.o: llsmutils.h dsputils.h llsm.h
$(OUT_DIR)/layer0.o: llsmutils.h dsputils.h llsm.h
$(OUT_DIR)/layer1.o: llsmutils.h dsputils.h llsm.h
$(OUT_DIR)/coder.o: dsputils.h llsm.h
$(OUT_DIR)/llsmrt.o: buffer.h llsmutils.h dsputils.h llsm.h llsmrt.h

$(CIGLET_O): $(CIGLET_SRC)
	mkdir -p $(OUT_DIR)
	$(CC) $(CFLAGS_REL) -o $(CIGLET_O) -c $(CIGLET_SRC)

$(OUT_DIR)/%.o: %.c
	mkdir -p $(OUT_DIR)
	$(CC) $(CFLAGS) -o $(OUT_DIR)/$*.o -c $*.c

install: $(OUT_DIR)/libllsm.a 
	mkdir -p $(PREFIX)/lib $(PREFIX)/include
	cp $(OUT_DIR)/libllsm.a $(PREFIX)/lib
	cp llsm.h llsmrt.h llsmutils.h $(PREFIX)/include

clean:
	rm -f build/*
