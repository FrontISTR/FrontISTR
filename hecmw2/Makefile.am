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

create-include-dir:
	if [ ! -d include ]; then \
		$(MKDIR) include; \
	fi

create-lib-dir:
	if [ ! -d lib ]; then \
		$(MKDIR) lib; \
	fi

create-bin-dir:
	if [ ! -d bin ]; then \
		$(MKDIR) bin; \
	fi

build-libs: clean-includes clean-libs
	@cd src && $(MAKE) @build_target@

build-tools:
	@cd tools && $(MAKE) @build_target@

clean-includes:
	$(RM) ./include/*

clean-libs:
	$(RM) ./lib/*

clean-bins:
	$(RM) ./bin/*

clean:
	@cd src && $(MAKE) clean
	@cd tools && $(MAKE) clean

distclean: clean-includes clean-libs clean-bins
	@cd src && $(MAKE) distclean
	@cd tools && $(MAKE) distclean

install:
	if [ ! -d $(PREFIX) ]; then \
		$(MKDIR) $(PREFIX); \
	fi
	if [ ! -d $(BINDIR) ]; then \
		$(MKDIR) $(BINDIR); \
	fi
	$(CP) ./bin/* $(BINDIR)
