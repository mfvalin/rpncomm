include $(VPATH)/RPN_COMM_version.inc
$(info this is RPN_COMM version $(RPN_COMM_version_s))

# general building rules
include $(VPATH)/Makefile.common
# sources specific (mechanically generated) dependencies
include $(VPATH)/dependencies.mk

LIB      = rpn_comm
CLEAN    = rpn_comm_fortran_stubs.f rpn_comm_c_stubs.c \
           $(STUB_LIBRARY) $(LIBRARY) $(VPATH)/rpn-comm_$(RPN_COMM_version_s)_multi.ssm
CLEANDIRS= $(VPATH)/rpn-comm_$(RPN_COMM_version_s)_multi $(LIBDIR)
TESTS    = TEST_000.Abs TEST_001.Abs TEST_002.Abs TEST_004.Abs TEST_005.Abs TEST_006.Abs TEST_007.Abs TEST_008.Abs
FMODULES = RPN_COMM_mod.o
LIBNAME  = $(LIB)_$(RPN_COMM_version)
LIBRARY  = $(LIBDIR)/lib$(LIBNAME).a
STUB_LIBRARY = $(LIBDIR)/lib$(LIB)stubs_$(RPN_COMM_version).a
SOURCES  = $(INCDECKS) $(CDECKS) $(FDECKS) $(HDECKS) $(F90DECKS)

ALL:  lib stublib tests $(VPATH)/includes $(VPATH)/RPN_COMM_interfaces.inc

tests:	$(TESTS)

stublib: $(STUB_LIBRARY)

lib: $(LIBRARY) $(VPATH)/includes

$(VPATH)/RPN_COMM_interfaces.inc: $(wildcard $(VPATH)/*.f) $(wildcard $(VPATH)/*.f90) $(wildcard $(VPATH)/*.c)
	(cd $(VPATH) ; cat RPN_COMM_*.f RPN_COMM_*.f90 RPN_COMM_*.c | ../tools/extract_interface.sh >RPN_COMM_interfaces.inc )

$(VPATH)/dependencies.mk: $(wildcard $(VPATH)/*.f) $(wildcard $(VPATH)/*.f90) $(wildcard $(VPATH)/*.c) $(wildcard $(VPATH)/*.h)
	-which gnu_find 2>/dev/null 1>/dev/null || (cd $(VPATH) ; find . -maxdepth 1 -type f | grep -v TEST_0 | ../tools/mk.dependencies.pl >dependencies.mk )
	-which gnu_find 2>/dev/null 1>/dev/null && (cd $(VPATH) ; gnu_find . -maxdepth 1 -type f | grep -v TEST_0 ../tools/mk.dependencies.pl >dependencies.mk )

ssm-package:
	rm -rf $(VPATH)/rpn-comm_${RPN_COMM_version_s}_multi
	(cd $(VPATH) ; tar zxf ssmtemplate_1.0_all.ssm ; mv ssmtemplate_1.0_all rpn-comm_$(RPN_COMM_version_s)_multi )
	(cd $(VPATH) ; cp rpn_comm_stubs.sh $(SOURCES) rpn-comm_$(RPN_COMM_version_s)_multi/src/.)
	(cd $(VPATH) ; cp Makefile Makefile.common Makefile.default Makefile.ECsetup *.mk \
	    rpn-comm_$(RPN_COMM_version_s)_multi/src/.)
	(cd $(VPATH) ; tar zcf rpn-comm_$(RPN_COMM_version_s)_multi.ssm rpn-comm_$(RPN_COMM_version_s)_multi)

rpn_comm_fortran_stubs.f: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh fortran

rpn_comm_c_stubs.c: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh c

$(STUB_LIBRARY): rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	mkdir -p $(LIBDIR)
	ar rcv $(STUB_LIBRARY) rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	(cd $(LIBDIR) ; ln -sf lib$(LIB)stubs_$(RPN_COMM_version).a lib$(LIB)stubs.a)

$(LIBRARY): $(OBJECTS)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY)_ $(OBJECTS)
	ar d $(LIBRARY)_ TEST_stubs.o rpn_comm_c_stubs.o rpn_comm_fortran_stubs.o
	mv $(LIBRARY)_ $(LIBRARY)
	(cd $(LIBDIR) ; ln -sf lib$(LIB)_$(RPN_COMM_version).a  lib$(LIB).a)
	mkdir -p $(INCDIR)
	cp *.mod $(INCDIR)

.PHONY:	includes
$(VPATH)/includes:
	cp $(VPATH)/*.inc $(INCDIR)
	touch $(VPATH)/includes

TEST_000.Abs: $(LIBRARY) TEST_000.o
TEST_001.Abs: $(LIBRARY) TEST_001.o
TEST_002.Abs: $(LIBRARY) TEST_002.o
TEST_003.Abs: $(LIBRARY) TEST_003.o
TEST_004.Abs: $(LIBRARY) TEST_004.o
TEST_005.Abs: $(LIBRARY) TEST_005.o
TEST_006.Abs: $(LIBRARY) TEST_006.o
TEST_007.Abs: $(LIBRARY) TEST_007.o
TEST_008.Abs: $(LIBRARY) TEST_008.o

