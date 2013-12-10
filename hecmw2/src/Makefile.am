MAKE = @make@

all: @build_target@

@build_target@:
	@cd hecmw                  && $(MAKE) @build_target@
	@cd visualizer             && $(MAKE) @build_target@

clean:
	@cd hecmw                  && $(MAKE) clean
	@cd visualizer             && $(MAKE) clean

distclean:
	@cd hecmw                  && $(MAKE) distclean
	@cd visualizer             && $(MAKE) distclean

install:
	@cd hecmw                  && $(MAKE) install
	@cd visualizer             && $(MAKE) install
