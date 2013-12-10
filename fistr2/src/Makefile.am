MAKE = @make@

all: @build_target@

@build_target@:
	@cd common   && $(MAKE) install-header
	@cd hecmw2   && $(MAKE) @build_target@
	@cd lib      && $(MAKE) @build_target@
	@cd common   && $(MAKE) @build_target@
	@cd analysis && $(MAKE) @build_target@
	@cd main     && $(MAKE) @build_target@

clean:
	@cd hecmw2   && $(MAKE) clean
	@cd lib      && $(MAKE) clean
	@cd common   && $(MAKE) clean
	@cd analysis && $(MAKE) clean
	@cd main     && $(MAKE) clean

distclean:
	@cd hecmw2   && $(MAKE) distclean
	@cd lib      && $(MAKE) distclean
	@cd common   && $(MAKE) distclean
	@cd analysis && $(MAKE) distclean
	@cd main     && $(MAKE) distclean

install:
