MAKE = @make@

all: @build_target@

@build_target@:
	@cd partitioner           && $(MAKE) @build_target@
	@cd visualizer            && $(MAKE) @build_target@
	@cd hec2rcap              && $(MAKE) @build_target@
	@cd result_type_converter && $(MAKE) @build_target@
	@cd result_file_merger    && $(MAKE) @build_target@

clean:
	@cd partitioner           && $(MAKE) clean
	@cd visualizer            && $(MAKE) clean
	@cd hec2rcap              && $(MAKE) clean
	@cd result_type_converter && $(MAKE) clean
	@cd result_file_merger    && $(MAKE) clean

distclean:
	@cd partitioner           && $(MAKE) distclean
	@cd visualizer            && $(MAKE) distclean
	@cd hec2rcap              && $(MAKE) distclean
	@cd result_type_converter && $(MAKE) distclean
	@cd result_file_merger    && $(MAKE) distclean

install:
	@cd partitioner           && $(MAKE) install
	@cd visualizer            && $(MAKE) install
	@cd hec2rcap              && $(MAKE) install
	@cd result_type_converter && $(MAKE) install
	@cd result_file_merger    && $(MAKE) install
