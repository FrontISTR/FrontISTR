CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
REGISTRY_IMAGE     ?= $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/serial

build: build_openblas build_mkl
build_openblas:
	docker build -t $(REGISTRY_IMAGE):openblas            -f Dockerfile . --target openblas
	docker build -t $(REGISTRY_IMAGE):openblas_trilinos12 -f Dockerfile . --target openblas-trilinos12
	docker build -t $(REGISTRY_IMAGE):openblas_metis4     -f Dockerfile . --target openblas-metis4
build_mkl:
	docker build -t $(REGISTRY_IMAGE):mkl                 -f Dockerfile . --target mkl
	docker build -t $(REGISTRY_IMAGE):mkl_trilinos12      -f Dockerfile . --target mkl-trilinos12
	docker build -t $(REGISTRY_IMAGE):mkl_metis4          -f Dockerfile . --target mkl-metis4

push: push_openblas push_mkl
push_openblas: build_openblas
	docker push $(REGISTRY_IMAGE):openblas
	docker push $(REGISTRY_IMAGE):openblas_trilinos12
	docker push $(REGISTRY_IMAGE):openblas_metis4
push_mkl: build_mkl
	docker push $(REGISTRY_IMAGE):mkl
	docker push $(REGISTRY_IMAGE):mkl_trilinos12
	docker push $(REGISTRY_IMAGE):mkl_metis4
