CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32/openblas

build: serial mpi openmp hybrid
push: login push_base push_serial push_mpi push_openmp push_hybrid

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

base: base_serial base_openmp base_mpi base_hybrid

base_serial:
	docker build -t $(CI_REGISTRY_IMAGE):base_serial       -f Dockerfile.base . --target serial
base_openmp:
	docker build -t $(CI_REGISTRY_IMAGE):base_openmp       -f Dockerfile.base . --target openmp
base_mpi:    base_serial
	docker build -t $(CI_REGISTRY_IMAGE):base_mpi          -f Dockerfile.base . --target mpi
base_hybrid: base_openmp
	docker build -t $(CI_REGISTRY_IMAGE):base_hybrid       -f Dockerfile.base . --target hybrid

serial: base_serial
	docker build -t $(CI_REGISTRY_IMAGE):serial            -f Dockerfile.serial . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):serial_trilinos12 -f Dockerfile.serial . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):serial_metis4     -f Dockerfile.serial . --target lib-metis4
mpi:    base_mpi
	docker build -t $(CI_REGISTRY_IMAGE):mpi               -f Dockerfile.mpi    . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):mpi_trilinos12    -f Dockerfile.mpi    . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):mpi_metis4        -f Dockerfile.mpi    . --target lib-metis4
openmp: base_openmp
	docker build -t $(CI_REGISTRY_IMAGE):openmp            -f Dockerfile.openmp . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):openmp_trilinos12 -f Dockerfile.openmp . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):openmp_metis4     -f Dockerfile.openmp . --target lib-metis4
hybrid: base_hybrid
	docker build -t $(CI_REGISTRY_IMAGE):hybrid            -f Dockerfile.hybrid . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):hybrid_trilinos12 -f Dockerfile.hybrid . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):hybrid_metis4     -f Dockerfile.hybrid . --target lib-metis4

push_base: push_base_serial push_base_openmp push_base_mpi push_base_hybrid
push_base_serial: base_serial
	docker push $(CI_REGISTRY_IMAGE):base_serial
push_base_openmp: base_openmp
	docker push $(CI_REGISTRY_IMAGE):base_openmp
push_base_mpi:    base_mpi    push_base_serial
	docker push $(CI_REGISTRY_IMAGE):base_mpi
push_base_hybrid: base_hybrid push_base_openmp
	docker push $(CI_REGISTRY_IMAGE):base_hybrid
push_serial: serial push_base_serial
	docker push $(CI_REGISTRY_IMAGE):serial
	docker push $(CI_REGISTRY_IMAGE):serial_trilinos12
	docker push $(CI_REGISTRY_IMAGE):serial_metis4
push_openmp: openmp push_base_openmp
	docker push $(CI_REGISTRY_IMAGE):openmp
	docker push $(CI_REGISTRY_IMAGE):openmp_trilinos12
	docker push $(CI_REGISTRY_IMAGE):openmp_metis4
push_mpi:    mpi push_base_mpi
	docker push $(CI_REGISTRY_IMAGE):mpi
	docker push $(CI_REGISTRY_IMAGE):mpi_trilinos12
	docker push $(CI_REGISTRY_IMAGE):mpi_metis4
push_hybrid: hyb push_base_hybrid
	docker push $(CI_REGISTRY_IMAGE):hybrid
	docker push $(CI_REGISTRY_IMAGE):hybrid_trilinos12
	docker push $(CI_REGISTRY_IMAGE):hybrid_metis4

