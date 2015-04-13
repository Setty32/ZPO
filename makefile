################################################################################
#
# Project name: zpo_projekt
# File name:    Makefile
# Author:       Pavol Eldes
# Login:        xeldes00
#
################################################################################

.PHONY = all clean debug release prepare all-clean pack

################################################################################
# compilation variables
################################################################################

CPP = g++ -std=c++11 -Wall -Wextra -Werror -pedantic -pedantic-errors

DIR = .

NAMES = main
OBJS = $(addsuffix .o,$(NAMES))
BIN = zpo_projekt

OBJDIR = obj
BINDIR = bin

# put library directories here
ifeq ($(OS), Windows_NT)
  INCDIR =
  LIBDIR =
  LINKDIR =
else
  INCDIR = -I/usr/local/include
  LIBDIR = -L
  LINKDIR = -lopencv_core -lopencv_imgproc -lopencv_highgui
endif

# directory structure written here is expected to be present in the project dir
DBGDIR = debug
DBGOBJS = $(addprefix $(OBJDIR)/$(DBGDIR)/,$(OBJS))
DBGBIN = $(addprefix $(BINDIR)/$(DBGDIR)/,$(BIN))
CPPDBG = $(CPP) -g -pg

RLSDIR = release
RLSOBJS = $(addprefix $(OBJDIR)/$(RLSDIR)/,$(OBJS))
RLSBIN = $(addprefix $(BINDIR)/$(RLSDIR)/,$(BIN))
CPPRLS = $(CPP) -O3

################################################################################
# default compile target
################################################################################

all: release

################################################################################
# debug compilation
################################################################################

debug: prepare $(DBGBIN)

#$(BINDIR)/$(DBGDIR)/%: $(OBJDIR)/$(DBGDIR)/%.o
#	$(CPPDBG) $< $(LIBDIR) $(LINKDIR) -o $@

$(DBGBIN): $(DBGOBJS)
	$(CPPDBG) $^ $(LIBDIR) $(LINKDIR) -o $@

$(OBJDIR)/$(DBGDIR)/%.o: $(DIR)/%.cpp $(DIR)/%.h
	$(CPPDBG) -c $< -o $@ $(INCDIR)

$(OBJDIR)/$(DBGDIR)/main.o: $(DIR)/main.cpp
	$(CPPDBG) -c $< -o $@ $(INCDIR)

################################################################################
# release compilation
################################################################################

release: prepare $(RLSBIN)

#$(BINDIR)/$(RLSDIR)/%: $(OBJDIR)/$(RLSDIR)/%.o
#	$(CPPRLS) $< $(LIBDIR) $(LINKDIR) -o $@

$(RLSBIN): $(RLSOBJS)
	$(CPPRLS) $^ $(LIBDIR) $(LINKDIR) -o $@

$(OBJDIR)/$(RLSDIR)/%.o: $(DIR)/%.cpp $(DIR)/%.h
	$(CPPRLS) -c $< -o $@ $(INCDIR)

$(OBJDIR)/$(RLSDIR)/main.o: $(DIR)/main.cpp
	$(CPPRLS) -c $< -o $@ $(INCDIR)

################################################################################
# creates needed directory structure
################################################################################

prepare:
	@mkdir -p $(OBJDIR)/$(DBGDIR) $(OBJDIR)/$(RLSDIR) $(BINDIR)/$(DBGDIR) \
	$(BINDIR)/$(RLSDIR)

################################################################################
# removes all make-made object and binary files together with their directory 
# structure
################################################################################

clean:
	@rm -rf $(OBJDIR)/$(DBGDIR)/ $(OBJDIR)/$(RLSDIR)/ $(BINDIR)/$(DBGDIR)/ \
	$(BINDIR)/$(RLSDIR)/

################################################################################
# removes all object and binary files together with their directory structure
################################################################################

all-clean:
	@rm -rf $(OBJDIR) $(BINDIR)

################################################################################
# creates an archive containing all project's source files
################################################################################

pack:
	@tar -czf $(DIR).tar.gz *.cpp *.h makefile

################################################################################
