MAKE = @make@

all: @build_target@

install-header:

@build_target@:
	@cd communication && $(MAKE) @build_target@
	@cd matrix        && $(MAKE) @build_target@
	@cd init          && $(MAKE) @build_target@
	@cd sparse_matrix && $(MAKE) @build_target@
	@cd solver_direct_mumps && $(MAKE) @build_target_mumps@
	@cd solver_11     && $(MAKE) @build_target@
	@cd solver_22     && $(MAKE) @build_target@
	@cd solver_33     && $(MAKE) @build_target@
	@cd solver_direct && $(MAKE) @build_target@
	@cd solver_direct_parallel && $(MAKE) @build_target@
	@cd solver_direct_lag && $(MAKE) @build_target@

clean:
	@cd communication && $(MAKE) clean
	@cd matrix        && $(MAKE) clean
	@cd init          && $(MAKE) clean
	@cd sparse_matrix && $(MAKE) clean
	@cd solver_direct_mumps && $(MAKE) clean
	@cd solver_11     && $(MAKE) clean
	@cd solver_22     && $(MAKE) clean
	@cd solver_33     && $(MAKE) clean
	@cd solver_direct && $(MAKE) clean
	@cd solver_direct_parallel && $(MAKE) clean
	@cd solver_direct_lag && $(MAKE) clean

distclean:
	@cd communication && $(MAKE) distclean
	@cd matrix        && $(MAKE) distclean
	@cd init          && $(MAKE) distclean
	@cd sparse_matrix && $(MAKE) distclean
	@cd solver_direct_mumps && $(MAKE) distclean
	@cd solver_11     && $(MAKE) distclean
	@cd solver_22     && $(MAKE) distclean
	@cd solver_33     && $(MAKE) distclean
	@cd solver_direct && $(MAKE) distclean
	@cd solver_direct_parallel && $(MAKE) distclean
	@cd solver_direct_lag && $(MAKE) distclean

maintainer-clean:
	@cd communication && $(MAKE) maintainer-clean
	@cd matrix        && $(MAKE) maintainer-clean
	@cd init          && $(MAKE) maintainer-clean
	@cd sparse_matrix && $(MAKE) maintainer_clean
	@cd solver_direct_mumps && $(MAKE) maintainer_clean
	@cd solver_11     && $(MAKE) maintainer-clean
	@cd solver_22     && $(MAKE) maintainer-clean
	@cd solver_33     && $(MAKE) maintainer-clean
	@cd solver_direct && $(MAKE) maintainer-clean
	@cd solver_direct_parallel && $(MAKE) maintainer-clean
	@cd solver_direct_lag && $(MAKE) maintainer-clean

install:
	@cd communication && $(MAKE) install
	@cd matrix        && $(MAKE) install
	@cd init          && $(MAKE) install
	@cd sparse_matrix && $(MAKE) install
	@cd solver_direct_mumps && $(MAKE) install
	@cd solver_11     && $(MAKE) install
	@cd solver_22     && $(MAKE) install
	@cd solver_33     && $(MAKE) install
	@cd solver_direct && $(MAKE) install
	@cd solver_direct_parallel && $(MAKE) install
	@cd solver_direct_lag && $(MAKE) install
