MAKE = @make@

all: @build_target@

@build_target@:
	@cd communication && $(MAKE) @build_target@
	@cd matrix        && $(MAKE) @build_target@
	@cd init          && $(MAKE) @build_target@
	@cd sparse_matrix && $(MAKE) @build_target@
	@cd las           && $(MAKE) @build_target@
	@cd mumps         && $(MAKE) @build_target@
	@cd mkl           && $(MAKE) @build_target@
	@cd clustermkl    && $(MAKE) @build_target@
	@cd precond       && $(MAKE) @build_target@
	@cd iterative     && $(MAKE) @build_target@
	@cd solver_direct && $(MAKE) @build_target@
	@cd solver_direct_parallel && $(MAKE) @build_target@
	@cd solver_direct_lag && $(MAKE) @build_target@
	@cd mpc           && $(MAKE) @build_target@
	@cd main          && $(MAKE) @build_target@
	@cd contact       && $(MAKE) @build_target@

install-header:
	@cd precond       && $(MAKE) install-header

install:
	@cd communication && $(MAKE) install
	@cd matrix        && $(MAKE) install
	@cd init          && $(MAKE) install
	@cd sparse_matrix && $(MAKE) install
	@cd las           && $(MAKE) install
	@cd mumps         && $(MAKE) install
	@cd mkl           && $(MAKE) install
	@cd clustermkl    && $(MAKE) install
	@cd precond       && $(MAKE) install
	@cd iterative     && $(MAKE) install
	@cd solver_direct && $(MAKE) install
	@cd solver_direct_parallel && $(MAKE) install
	@cd solver_direct_lag && $(MAKE) install
	@cd mpc           && $(MAKE) install
	@cd main          && $(MAKE) install
	@cd contact       && $(MAKE) install

clean:
	@cd communication && $(MAKE) clean
	@cd matrix        && $(MAKE) clean
	@cd init          && $(MAKE) clean
	@cd sparse_matrix && $(MAKE) clean
	@cd las           && $(MAKE) clean
	@cd mumps         && $(MAKE) clean
	@cd mkl           && $(MAKE) clean
	@cd clustermkl    && $(MAKE) clean
	@cd precond       && $(MAKE) clean
	@cd iterative     && $(MAKE) clean
	@cd solver_direct && $(MAKE) clean
	@cd solver_direct_parallel && $(MAKE) clean
	@cd solver_direct_lag && $(MAKE) clean
	@cd mpc           && $(MAKE) clean
	@cd main          && $(MAKE) clean
	@cd contact       && $(MAKE) clean

distclean: clean
	@cd communication && $(MAKE) distclean
	@cd matrix        && $(MAKE) distclean
	@cd init          && $(MAKE) distclean
	@cd sparse_matrix && $(MAKE) distclean
	@cd las           && $(MAKE) distclean
	@cd mumps         && $(MAKE) distclean
	@cd mkl           && $(MAKE) distclean
	@cd clustermkl    && $(MAKE) distclean
	@cd precond       && $(MAKE) distclean
	@cd iterative     && $(MAKE) distclean
	@cd solver_direct && $(MAKE) distclean
	@cd solver_direct_parallel && $(MAKE) distclean
	@cd solver_direct_lag && $(MAKE) distclean
	@cd mpc           && $(MAKE) distclean
	@cd main          && $(MAKE) distclean
	@cd contact       && $(MAKE) distclean

maintainer-clean:
	@cd communication && $(MAKE) maintainer-clean
	@cd matrix        && $(MAKE) maintainer-clean
	@cd init          && $(MAKE) maintainer-clean
	@cd sparse_matrix && $(MAKE) maintainer_clean
	@cd las           && $(MAKE) maintainer-clean
	@cd mumps         && $(MAKE) maintainer_clean
	@cd mkl           && $(MAKE) maintainer_clean
	@cd clustermkl    && $(MAKE) maintainer_clean
	@cd precond       && $(MAKE) maintainer-clean
	@cd iterative     && $(MAKE) maintainer-clean
	@cd solver_direct && $(MAKE) maintainer-clean
	@cd solver_direct_parallel && $(MAKE) maintainer-clean
	@cd solver_direct_lag && $(MAKE) maintainer-clean
	@cd mpc           && $(MAKE) maintainer-clean
	@cd main          && $(MAKE) maintainer-clean
	@cd contact       && $(MAKE) maintainer-clean

depend:
	@cd communication && $(MAKE) depend
	@cd matrix        && $(MAKE) depend
	@cd init          && $(MAKE) depend
	@cd sparse_matrix && $(MAKE) depend
	@cd las           && $(MAKE) depend
	@cd mumps         && $(MAKE) depend
	@cd mkl           && $(MAKE) depend
	@cd clustermkl    && $(MAKE) depend
	@cd precond       && $(MAKE) depend
	@cd iterative     && $(MAKE) depend
	@cd solver_direct && $(MAKE) depend
	@cd solver_direct_parallel && $(MAKE) depend
	@cd solver_direct_lag && $(MAKE) depend
	@cd mpc           && $(MAKE) depend
	@cd main          && $(MAKE) depend
	@cd contact       && $(MAKE) depend
