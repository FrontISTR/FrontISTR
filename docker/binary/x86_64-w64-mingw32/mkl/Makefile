CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32/mkl

build: ser mpi omp hyb

push: push_ser push_mpi push_omp push_hyb

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

base:
	docker build -t $(CI_REGISTRY_IMAGE):base_ser       -f base.Dockerfile . --target ser
	docker build -t $(CI_REGISTRY_IMAGE):base_omp       -f base.Dockerfile . --target omp
	docker build -t $(CI_REGISTRY_IMAGE):base_mpi       -f base.Dockerfile . --target mpi
	docker build -t $(CI_REGISTRY_IMAGE):base_hyb       -f base.Dockerfile . --target hyb

ser: base
	docker build -t $(CI_REGISTRY_IMAGE):ser            -f ser.Dockerfile . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):ser_trilinos12 -f ser.Dockerfile . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):ser_metis4     -f ser.Dockerfile . --target lib-metis4
mpi: base
	docker build -t $(CI_REGISTRY_IMAGE):mpi            -f mpi.Dockerfile . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):mpi_trilinos12 -f mpi.Dockerfile . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):mpi_metis4     -f mpi.Dockerfile . --target lib-metis4
omp: base
	docker build -t $(CI_REGISTRY_IMAGE):omp            -f omp.Dockerfile . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):omp_trilinos12 -f omp.Dockerfile . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):omp_metis4     -f omp.Dockerfile . --target lib-metis4
hyb: base
	docker build -t $(CI_REGISTRY_IMAGE):hyb            -f hyb.Dockerfile . --target lib
	docker build -t $(CI_REGISTRY_IMAGE):hyb_trilinos12 -f hyb.Dockerfile . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE):hyb_metis4     -f hyb.Dockerfile . --target lib-metis4

push_ser: ser login
	docker push $(CI_REGISTRY_IMAGE):ser
	docker push $(CI_REGISTRY_IMAGE):ser_trilinos12
	docker push $(CI_REGISTRY_IMAGE):ser_metis4
push_omp: omp login
	docker push $(CI_REGISTRY_IMAGE):omp
	docker push $(CI_REGISTRY_IMAGE):omp_trilinos12
	docker push $(CI_REGISTRY_IMAGE):omp_metis4
push_mpi: mpi login
	docker push $(CI_REGISTRY_IMAGE):mpi
	docker push $(CI_REGISTRY_IMAGE):mpi_trilinos12
	docker push $(CI_REGISTRY_IMAGE):mpi_metis4
push_hyb: hyb login
	docker push $(CI_REGISTRY_IMAGE):hyb
	docker push $(CI_REGISTRY_IMAGE):hyb_trilinos12
	docker push $(CI_REGISTRY_IMAGE):hyb_metis4

