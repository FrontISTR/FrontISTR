MAKE = @make@

all: @build_target@

@build_target@:
	@cd hecmw                  && $(MAKE) @build_target@
	@cd solver_direct_parallel && $(MAKE) @build_target@
#	@cd partitioner            && $(MAKE) @build_target@
	@cd visualizer             && $(MAKE) @build_target@

clean:
	@cd hecmw                  && $(MAKE) clean
	@cd solver_direct_parallel && $(MAKE) clean
	@cd partitioner            && $(MAKE) clean
	@cd visualizer             && $(MAKE) clean

distclean:
	@cd hecmw                  && $(MAKE) distclean
	@cd solver_direct_parallel && $(MAKE) distclean
	@cd partitioner            && $(MAKE) distclean
	@cd visualizer             && $(MAKE) distclean

install:
	@cd hecmw                  && $(MAKE) install
	@cd solver_direct_parallel && $(MAKE) install
	@cd partitioner            && $(MAKE) install
	@cd visualizer             && $(MAKE) install
