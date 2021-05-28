CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
REGISTRY_IMAGE     ?= $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/hybrid

build: build_msmpi build_impi
build_msmpi: build_msmpi_openblas build_msmpi_mkl
build_msmpi_openblas:
	docker build -t $(REGISTRY_IMAGE):msmpi_openblas            -f Dockerfile.msmpi . --target openblas
	docker build -t $(REGISTRY_IMAGE):msmpi_openblas_trilinos12 -f Dockerfile.msmpi . --target openblas-trilinos12
	docker build -t $(REGISTRY_IMAGE):msmpi_openblas_metis4     -f Dockerfile.msmpi . --target openblas-metis4
build_msmpi_mkl:
	docker build -t $(REGISTRY_IMAGE):msmpi_mkl                 -f Dockerfile.msmpi . --target mkl
	docker build -t $(REGISTRY_IMAGE):msmpi_mkl_trilinos12      -f Dockerfile.msmpi . --target mkl-trilinos12
	docker build -t $(REGISTRY_IMAGE):msmpi_mkl_metis4          -f Dockerfile.msmpi . --target mkl-metis4
build_impi: build_impi_openblas build_impi_mkl
build_impi_openblas:
	docker build -t $(REGISTRY_IMAGE):impi_openblas             -f Dockerfile.impi  . --target openblas
	docker build -t $(REGISTRY_IMAGE):impi_openblas_trilinos12  -f Dockerfile.impi  . --target openblas-trilinos12
	docker build -t $(REGISTRY_IMAGE):impi_openblas_metis4      -f Dockerfile.impi  . --target openblas-metis4
build_impi_mkl:
	docker build -t $(REGISTRY_IMAGE):impi_mkl                  -f Dockerfile.impi  . --target mkl
	docker build -t $(REGISTRY_IMAGE):impi_mkl_trilinos12       -f Dockerfile.impi  . --target mkl-trilinos12
	docker build -t $(REGISTRY_IMAGE):impi_mkl_metis4           -f Dockerfile.impi  . --target mkl-metis4

push: push_msmpi push_impi
push_msmpi:  push_msmpi_openblas push_msmpi_mkl
push_msmpi_openblas: build_msmpi_openblas
	docker push $(REGISTRY_IMAGE):msmpi_openblas
	docker push $(REGISTRY_IMAGE):msmpi_openblas_trilinos12
	docker push $(REGISTRY_IMAGE):msmpi_openblas_metis4
push_msmpi_mkl: build_msmpi_mkl
	docker push $(REGISTRY_IMAGE):msmpi_mkl
	docker push $(REGISTRY_IMAGE):msmpi_mkl_trilinos12
	docker push $(REGISTRY_IMAGE):msmpi_mkl_metis4
push_impi:  push_impi_openblas push_impi_mkl
push_impi_openblas: build_impi_openblas
	docker push $(REGISTRY_IMAGE):impi_openblas
	docker push $(REGISTRY_IMAGE):impi_openblas_trilinos12
	docker push $(REGISTRY_IMAGE):impi_openblas_metis4
push_impi_mkl: build_impi_mkl
	docker push $(REGISTRY_IMAGE):impi_mkl
	docker push $(REGISTRY_IMAGE):impi_mkl_trilinos12
	docker push $(REGISTRY_IMAGE):impi_mkl_metis4
