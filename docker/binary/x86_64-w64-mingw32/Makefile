CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32

all: ser mpi omp hyb

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

base:
	docker build -t $(CI_REGISTRY_IMAGE)/base           -f Dockerfile.ser . --target base
	docker build -t $(CI_REGISTRY_IMAGE)/common         -f Dockerfile.ser . --target common
ser: base
	docker build -t $(CI_REGISTRY_IMAGE)/ser            -f Dockerfile.ser . --target lib
	docker build -t $(CI_REGISTRY_IMAGE)/ser:trilinos12 -f Dockerfile.ser . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE)/ser:metis4     -f Dockerfile.ser . --target lib-metis4
mpi: base
	docker build -t $(CI_REGISTRY_IMAGE)/mpi            -f Dockerfile.mpi . --target lib
	docker build -t $(CI_REGISTRY_IMAGE)/mpi:trilinos12 -f Dockerfile.mpi . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE)/mpi:metis4     -f Dockerfile.mpi . --target lib-metis4
omp: base
	docker build -t $(CI_REGISTRY_IMAGE)/omp            -f Dockerfile.omp . --target lib
	docker build -t $(CI_REGISTRY_IMAGE)/omp:trilinos12 -f Dockerfile.omp . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE)/omp:metis4     -f Dockerfile.omp . --target lib-metis4
hyb: base
	docker build -t $(CI_REGISTRY_IMAGE)/hyb            -f Dockerfile.hyb . --target lib
	docker build -t $(CI_REGISTRY_IMAGE)/hyb:trilinos12 -f Dockerfile.hyb . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE)/hyb:metis4     -f Dockerfile.hyb . --target lib-metis4

push_base: base login
	docker push $(CI_REGISTRY_IMAGE)/base
	docker push $(CI_REGISTRY_IMAGE)/common
push_ser: ser login push_base
	docker push $(CI_REGISTRY_IMAGE)/ser
	docker push $(CI_REGISTRY_IMAGE)/ser:trilinos12
	docker push $(CI_REGISTRY_IMAGE)/ser:metis4
push_omp: omp login push_base
	docker push $(CI_REGISTRY_IMAGE)/omp
	docker push $(CI_REGISTRY_IMAGE)/omp:trilinos12
	docker push $(CI_REGISTRY_IMAGE)/omp:metis4
push_mpi: mpi login push_base
	docker push $(CI_REGISTRY_IMAGE)/mpi
	docker push $(CI_REGISTRY_IMAGE)/mpi:trilinos12
	docker push $(CI_REGISTRY_IMAGE)/mpi:metis4
push_hyb: hyb login push_base
	docker push $(CI_REGISTRY_IMAGE)/hyb
	docker push $(CI_REGISTRY_IMAGE)/hyb:trilinos12
	docker push $(CI_REGISTRY_IMAGE)/hyb:metis4

push: push_ser push_mpi push_omp push_hyb
