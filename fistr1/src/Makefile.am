MAKE = @make@

all: @build_target@

@build_target@:
	@cd common   && $(MAKE) install-header
	@cd lib      && $(MAKE) @build_target@
	@cd common   && $(MAKE) @build_target@
	@cd analysis && $(MAKE) @build_target@
	@cd main     && $(MAKE) @build_target@

clean:
	@cd lib      && $(MAKE) clean
	@cd common   && $(MAKE) clean
	@cd analysis && $(MAKE) clean
	@cd main     && $(MAKE) clean

distclean:
	@cd lib      && $(MAKE) distclean
	@cd common   && $(MAKE) distclean
	@cd analysis && $(MAKE) distclean
	@cd main     && $(MAKE) distclean

install:

depend:
	@cd lib      && $(MAKE) depend
	@cd common   && $(MAKE) depend
	@cd analysis && $(MAKE) depend
	@cd main     && $(MAKE) depend
