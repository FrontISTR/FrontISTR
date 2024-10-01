MAKE = @make@

all: @build_target@

install-header:
	@cd dynamic_load_balancing && $(MAKE) install_header

@build_target@:
	@cd adaptation             && $(MAKE) @build_target@
	@cd dynamic_load_balancing && $(MAKE) @build_target@
	@cd element_smoothing      && $(MAKE) @build_target@
	@cd jacobian               && $(MAKE) @build_target@

clean:
	@cd adaptation             && $(MAKE) clean
	@cd dynamic_load_balancing && $(MAKE) clean
	@cd element_smoothing      && $(MAKE) clean
	@cd jacobian               && $(MAKE) clean

distclean:
	@cd adaptation             && $(MAKE) distclean
	@cd dynamic_load_balancing && $(MAKE) distclean
	@cd element_smoothing      && $(MAKE) distclean
	@cd jacobian               && $(MAKE) distclean

install:
	@cd adaptation             && $(MAKE) install
	@cd dynamic_load_balancing && $(MAKE) install
	@cd element_smoothing      && $(MAKE) install
	@cd jacobian               && $(MAKE) install

depend:
	@cd adaptation             && $(MAKE) depend
	@cd dynamic_load_balancing && $(MAKE) depend
	@cd element_smoothing      && $(MAKE) depend
	@cd jacobian               && $(MAKE) depend
