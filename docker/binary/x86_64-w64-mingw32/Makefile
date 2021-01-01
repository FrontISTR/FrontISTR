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
	docker build -t $(CI_REGISTRY_IMAGE)/lib    . --target lib
mpilib:
	docker build -t $(CI_REGISTRY_IMAGE)/mpilib . --target mpilib

push: login build lib mpilib
	docker push $(CI_REGISTRY_IMAGE)/build
	docker push $(CI_REGISTRY_IMAGE)/lib
	docker push $(CI_REGISTRY_IMAGE)/mpilib

