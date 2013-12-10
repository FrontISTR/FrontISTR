MAKE = @make@

all: @build_target@

@build_target@:
	@cd conv2mw3    && $(MAKE) @build_target@
	@cd partitioner && $(MAKE) @build_target@
	@cd visualizer  && $(MAKE) @build_target@

clean:
	@cd conv2mw3    && $(MAKE) clean
	@cd partitioner && $(MAKE) clean
	@cd visualizer  && $(MAKE) clean

distclean:
	@cd conv2mw3    && $(MAKE) distclean
	@cd partitioner && $(MAKE) distclean
	@cd visualizer  && $(MAKE) distclean

install:
	@cd conv2mw3    && $(MAKE) install
	@cd partitioner && $(MAKE) install
	@cd visualizer  && $(MAKE) install
