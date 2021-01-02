ifndef COMPILER
COMPILER=x64
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean

all: libvecjac_$(TEST)$(DEBUG)_$(ABI).a

help:
	@echo "make [COMPILER=x64|x200] [NDEBUG=0|1|2|3|4|5] [TEST=0|1|2|3|4|5|6|7] [all|clean|help]"

libvecjac_$(TEST)$(DEBUG)_$(ABI).a: common_$(TEST)$(DEBUG)_$(ABI).o $(MKFS)
	$(AR) $(ARFLAGS) $@ common_$(TEST)$(DEBUG)_$(ABI).o

common_$(TEST)$(DEBUG)_$(ABI).o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-$(RM) common_$(TEST)$(DEBUG)_$(ABI).o
	-$(RM) libvecjac_$(TEST)$(DEBUG)_$(ABI).a
