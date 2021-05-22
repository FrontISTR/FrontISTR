CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32

SUBDIRS := openblas mkl

all: build

build: $(SUBDIRS)

push: $(SUBDIRS)


login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

$(SUBDIRS): base
	make -C $@ $(MAKECMDGOALS)

base:
	docker build -t $(CI_REGISTRY_IMAGE)/core  -f Dockerfile . --target core
	docker build -t $(CI_REGISTRY_IMAGE)/base  -f Dockerfile . --target base

push_base: base login
	docker push $(CI_REGISTRY_IMAGE)/base




