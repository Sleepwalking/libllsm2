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
TARGET_A = $(OUT_DIR)/libllsm2.a

CIGLET_PREFIX = /usr
GVPS_PREFIX = /usr
PYIN_PREFIX = /usr
NEBULA_PREFIX = /usr
IHNM_PREFIX = /usr

FP_TYPE = float
CONFIG  = Debug

CIGLET_A = $(CIGLET_PREFIX)/lib/libciglet.a
CIGLET_INCLUDE = $(CIGLET_PREFIX)/include/
GVPS_A = $(GVPS_PREFIX)/lib/libgvps.a
GVPS_INCLUDE = $(GVPS_PREFIX)/include/
PYIN_A = $(PYIN_PREFIX)/lib/libpyin.a
PYIN_INCLUDE = $(PYIN_PREFIX)/include/
NEBULA_A = $(NEBULA_PREFIX)/lib/libnebula.a
NEBULA_INCLUDE = $(NEBULA_PREFIX)/include
IHNM_A = $(IHNM_PREFIX)/lib/libihnm.a
IHNM_INCLUDE = $(IHNM_PREFIX)/include

ifeq 'Darwin' '$(shell uname)'
  CFLAGS_PLAT =
else
  CFLAGS_PLAT = -fopenmp
endif

ARFLAGS = -rv
CFLAGS_COMMON = -DFP_TYPE=$(FP_TYPE) -std=c99 -Wall -fPIC -pthread -DUSE_PTHREAD $(CFLAGS_PLAT) \
	-I$(CIGLET_INCLUDE) -I$(GVPS_INCLUDE) -I$(PYIN_INCLUDE) -I$(NEBULA_INCLUDE) -I$(IHNM_INCLUDE)
CFLAGS_DBG = $(CFLAGS_COMMON) -Og -g
CFLAGS_REL = $(CFLAGS_COMMON) -Ofast
ifeq ($(CONFIG), Debug)
  CFLAGS = $(CFLAGS_DBG)
else
  CFLAGS = $(CFLAGS_REL)
endif

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
	  $(CIGLET_A) $(PYIN_A) $(GVPS_A) $(IHNM_A) test/verify-utils.h
	$(CC) $(CFLAGS) -o $(OUT_DIR)/test-$* test/test-$*.c \
	  $(TARGET_A) $(CIGLET_A) $(PYIN_A) $(NEBULA_A) $(GVPS_A) $(IHNM_A) -lm

$(OUT_DIR)/demo-%: test/demo-%.c $(TARGET_A) \
	  $(CIGLET_A) $(PYIN_A) $(NEBULA_A) $(GVPS_A) $(IHNM_A)
	$(CC) $(CFLAGS) -o $(OUT_DIR)/demo-$* test/demo-$*.c $(CFLAGS_PLAT) \
	  $(TARGET_A) $(CIGLET_A) $(PYIN_A) $(NEBULA_A) $(GVPS_A) $(IHNM_A) -lm

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

$(OUT_DIR)/%.o: %.c constants.h
	mkdir -p $(OUT_DIR)
	$(CC) $(CFLAGS) -o $(OUT_DIR)/$*.o -c $*.c

install: $(OUT_DIR)/libllsm2.a 
	mkdir -p $(PREFIX)/lib $(PREFIX)/include/libllsm2
	cp $(OUT_DIR)/libllsm2.a $(PREFIX)/lib
	cp llsm.h llsmrt.h llsmutils.h dsputils.h buffer.h $(PREFIX)/include/libllsm2

clean:
	rm -f build/*
