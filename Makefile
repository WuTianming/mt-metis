#######################################################################
# Makefile generated by './configure' at Thu 16 Feb 2023 06:07:22 PM CST
# Using flags:
#	-DCMAKE_VERBOSE_MAKEFILE=1
#	-DBIGVERTICES=1
#	-DBIGEDGES=1
#	-DBIGWEIGHTS=1
#	-DCMAKE_C_COMPILER=gcc-9
#	-DDEBUG=1
#	-DASSERT=1
#	-DDOMLIB_PATH=domlib
#	-DWILDRIVER_PATH=wildriver
#	-DMETIS_PATH=metis
#######################################################################

all test clean install:
	make -C build/Linux-x86_64 $@ $(MAKEFLAGS)

distclean:
	rm -rvf build/Linux-x86_64 Makefile

