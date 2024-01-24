# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

#include Makefile.arch

#pour une architecture standard avec command use root
#include $(ROOTSYS)/share/doc/root/test/Makefile.arch
include $(ROOTSYS)/etc/Makefile.arch

#------------------------------------------------------------------------------

CONVERTO        = Convert_Narval.$(ObjSuf)
CONVERTS        = Convert_Narval.$(SrcSuf)

FONCTIONO	= Fonctions_convert.$(ObjSuf)
FONCTIONS	= Fonctions_convert.$(SrcSuf)	

ifeq ($(PLATFORM),win32)
CLASSLIB      = libClasslib.lib
else
CLASSLIB      = $(CLASSSO)
endif

MAINO    = Convert_Narval_Main.$(ObjSuf) 
MAINS    = Convert_Narval_Main.$(SrcSuf)
MAIN     = Convert$(ExeSuf)

OBJS          = $(CONVERTO) $(FONCTIONO) $(MAINO)

PROGRAMS      = $(MAIN)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    Aclock Hello Tetris

all:            $(PROGRAMS)

$(CLASSSO):     $(CLASSO) $(CONVERTO) $(FONCTIONO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "$@ done"

$(MAIN):        $(MAINO) $(CONVERTO) $(FONCTIONO) $(CLASSSO)
		$(LD) $(LDFLAGS) $(MAINO)$(CONVERTO) $(FONCTIONO) $(CLASSLIB) $(LIBS) \
		   $(OutPutOpt)$(MAIN)
		@echo "$@ done"

clean:
		@rm -f $(OBJS)  *Dict* core

.SUFFIXES: .$(SrcSuf)

###
.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
###
