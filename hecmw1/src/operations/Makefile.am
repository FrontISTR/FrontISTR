MAKE = @make@

all: @build_target@

install-header:
	@cd dynamic_load_balancing && $(MAKE) install_header

@build_target@:
	@cd adaptation             && $(MAKE) @build_target@
	@cd dynamic_load_balancing && $(MAKE) @build_target@
	@cd jacobian               && $(MAKE) @build_target@

clean:
	@cd adaptation             && $(MAKE) clean
	@cd dynamic_load_balancing && $(MAKE) clean
	@cd jacobian               && $(MAKE) clean

distclean:
	@cd adaptation             && $(MAKE) distclean
	@cd dynamic_load_balancing && $(MAKE) distclean
	@cd jacobian               && $(MAKE) distclean

install:
	@cd adaptaion              && $(MAKE) install
	@cd dynamic_load_balancing && $(MAKE) install
	@cd jacobian               && $(MAKE) install
