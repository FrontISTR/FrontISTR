MAKE = @make@

all: @build_target@

@build_target@:
	@cd static  && $(MAKE) @build_target_mkl@
	@cd dynamic && $(MAKE) @build_target@
	@cd heat    && $(MAKE) @build_target@

clean:
	@cd static  && $(MAKE) clean
	@cd dynamic && $(MAKE) clean
	@cd heat    && $(MAKE) clean

distclean:
	@cd static  && $(MAKE) distclean
	@cd dynamic && $(MAKE) distclean
	@cd heat    && $(MAKE) distclean

install:
