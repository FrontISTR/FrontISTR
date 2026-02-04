MAKE = @make@

all: @build_target@

@build_target@:
	@cd interaction && $(MAKE) @build_target@
	@cd static  && $(MAKE) @build_target@
	@cd dynamic && $(MAKE) @build_target@
	@cd heat    && $(MAKE) @build_target@

clean:
	@cd interaction && $(MAKE) clean
	@cd static  && $(MAKE) clean
	@cd dynamic && $(MAKE) clean
	@cd heat    && $(MAKE) clean

distclean:
	@cd interaction && $(MAKE) distclean
	@cd static  && $(MAKE) distclean
	@cd dynamic && $(MAKE) distclean
	@cd heat    && $(MAKE) distclean

install:

depend:
	@cd interaction && $(MAKE) depend
	@cd static  && $(MAKE) depend
	@cd dynamic && $(MAKE) depend
	@cd heat    && $(MAKE) depend
