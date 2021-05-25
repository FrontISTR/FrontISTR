CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
SUBDIRS := serial thread process hybrid

.PHONY: clean all build push $(SUBDIRS)

build: build_serial build_thread
	$(MAKE) subdir MAKECMDGOALS=$@
push: login push_serial push_thread
	$(MAKE) subdir MAKECMDGOALS=$@

subdir: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

build_core:
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/env   -f Dockerfile . --target env
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/core  -f Dockerfile . --target core
build_serial: build_core
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/base:serial -f Dockerfile . --target serial
build_thread: build_core
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/base:thread -f Dockerfile . --target thread

push_core: build_core login
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/env
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/core
push_serial: build_serial
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/base:serial
push_thread: build_thread
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32/base:thread
