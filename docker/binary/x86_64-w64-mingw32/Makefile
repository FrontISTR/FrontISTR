CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32

all: build lib mpilib

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

build:
	docker build -t $(CI_REGISTRY_IMAGE)/build  . --target build
lib:
	docker build -t $(CI_REGISTRY_IMAGE)/lib               . --target lib
	docker build -t $(CI_REGISTRY_IMAGE)/lib:trilinos12    . --target lib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE)/lib:metis4        . --target lib-metis4
mpilib:
	docker build -t $(CI_REGISTRY_IMAGE)/mpilib            . --target mpilib
	docker build -t $(CI_REGISTRY_IMAGE)/mpilib:trilinos12 . --target mpilib-trilinos12
	docker build -t $(CI_REGISTRY_IMAGE)/mpilib:metis4     . --target mpilib-metis4

push: login build lib mpilib
	docker push $(CI_REGISTRY_IMAGE)/build
	docker push $(CI_REGISTRY_IMAGE)/lib
	docker push $(CI_REGISTRY_IMAGE)/mpilib
	docker push $(CI_REGISTRY_IMAGE)/lib:trilinos12
	docker push $(CI_REGISTRY_IMAGE)/mpilib:trilinos12
	docker push $(CI_REGISTRY_IMAGE)/lib:metis4
	docker push $(CI_REGISTRY_IMAGE)/mpilib:metis4

