PREFIX     = @prefix@
BINDIR     = @bindir@
LIBDIR     = @libdir@
INCLUDEDIR = @includedir@

MAKE       = @make@
AR         = @ar@
CP         = @cp@
RM         = @rm@
RANLIB     = @ranlib@
MKDIR      = @mkdir@

all: @all_build_target@

@build_target@:
	if [ ! -d include ]; then \
		$(MKDIR) include; \
	fi
	if [ ! -d lib ]; then \
		$(MKDIR) lib; \
	fi
	if [ ! -d bin ]; then \
		$(MKDIR) bin; \
	fi
	$(RM) ./include/*
	$(RM) ./lib/*
	@cd src && $(MAKE) @build_target@

build-tools:
#	@cd tools && $(MAKE)

clean:
	@cd src && $(MAKE) clean
#	@cd tools && $(MAKE) clean

distclean:
	@cd src && $(MAKE) distclean
#	@cd tools && $(MAKE) distclean
	$(RM) ./include/*
	$(RM) ./bin/*
	$(RM) ./lib/*

install:
	if [ ! -d $(PREFIX) ]; then \
		$(MKDIR) $(PREFIX); \
	fi
	if [ ! -d $(BINDIR) ]; then \
		$(MKDIR) $(BINDIR); \
	fi
	$(CP) ./bin/* $(BINDIR)
