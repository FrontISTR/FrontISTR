CP               = @cp@
MV               = @mv@
RM               = @rm@

#-----------------------------------------------------------------------------#
LEX              = @lex@

SRCCOMMONDIR     = ../src/common

SRCPARTITIONERDIR = ../tools/partitioner

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

LEXSRCPARTITIONER = \
	hecmw_partlex.c

LEXSRC = $(LEXSRCCOMMON) $(LEXSRCPARTITIONER)

all: setup-msg setup-lex

setup-msg:
	$(MSGSETUP) $(MSGDEFMASTER)
	for i in $(MSGSRC); do $(MV) $$i $(SRCCOMMONDIR)/.; done

setup-lex: $(LEXSRC)
	for i in $(LEXSRCCOMMON); do $(MV) $$i $(SRCCOMMONDIR)/.; done
	for i in $(LEXSRCPARTITIONER); do $(MV) $$i $(SRCPARTITIONERDIR)/.; done

hecmw_ctrllex.c: $(SRCCOMMONDIR)/hecmw_ctrllex.l
	$(LEX) -t -PHECMW_ctrl $(SRCCOMMONDIR)/hecmw_ctrllex.l > $@

hecmw_gflex.c: $(SRCCOMMONDIR)/hecmw_gflex.l
	$(LEX) -t -PHECMW_gf $(SRCCOMMONDIR)/hecmw_gflex.l > $@

hecmw_heclex.c: $(SRCCOMMONDIR)/hecmw_heclex.l
	$(LEX) -t -PHECMW_hec $(SRCCOMMONDIR)/hecmw_heclex.l > $@

hecmw_ablex.c: $(SRCCOMMONDIR)/hecmw_ablex.l
	$(LEX) -t -PHECMW_ab -i $(SRCCOMMONDIR)/hecmw_ablex.l > $@

hecmw_partlex.c: $(SRCPARTITIONERDIR)/hecmw_partlex.l
	$(LEX) -t -PHECMW_part $(SRCPARTITIONERDIR)/hecmw_partlex.l > $@

install:

clean:
	$(RM) $(MSGSRC) $(LEXSRC)

distclean: clean

depend:
