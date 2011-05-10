CP               = @cp@
RM               = @rm@

#-----------------------------------------------------------------------------#
LEX              = @lex@

SRCCOMMONDIR     = ../src/common

.SUFFIXES:
.SUFFIXES: .c .l

MSGSRC = \
	hecmw_msg_table.c \
	hecmw_msgno.h \
	hecmw_msgno_f.f90

MSGSETUP = ./msg_setup.pl

MSGDEFMASTER = msg_master.xml

LEXSRCCOMMON = \
	hecmw_ctrllex.c \
	hecmw_heclex.c \
	hecmw_ablex.c \
	hecmw_gflex.c

LEXSRC = $(LEXSRCCOMMON)

all: setup-msg setup-lex

setup-msg:
	$(MSGSETUP) $(MSGDEFMASTER)
	for i in $(MSGSRC); do $(CP) $$i $(SRCCOMMONDIR)/.; done

setup-lex: $(LEXSRC)
	for i in $(LEXSRCCOMMON); do $(CP) $$i $(SRCCOMMONDIR)/.; done

hecmw_ctrllex.c: $(SRCCOMMONDIR)/hecmw_ctrllex.l
	$(LEX) -t -PHECMW_ctrl $(SRCCOMMONDIR)/hecmw_ctrllex.l > $@

hecmw_gflex.c: $(SRCCOMMONDIR)/hecmw_gflex.l
	$(LEX) -t -PHECMW_gf $(SRCCOMMONDIR)/hecmw_gflex.l > $@

hecmw_heclex.c: $(SRCCOMMONDIR)/hecmw_heclex.l
	$(LEX) -t -PHECMW_hec $(SRCCOMMONDIR)/hecmw_heclex.l > $@

hecmw_ablex.c: $(SRCCOMMONDIR)/hecmw_ablex.l
	$(LEX) -t -PHECMW_ab -i $(SRCCOMMONDIR)/hecmw_ablex.l > $@

install:

clean:
	$(RM) $(MSGSRC) $(LEXSRC)

distclean: clean
