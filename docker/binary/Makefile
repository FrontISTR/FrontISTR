CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
CI_COMMIT_REF_NAME ?= _no_branch

all: x86_64-w64-mingw32-build x86_64-w64-mingw32-lib

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(REGISTRY)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

x86_64-w64-mingw32-build:
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-build . --target x86_64-w64-mingw32-build

x86_64-w64-mingw32-lib:
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-lib . --target x86_64-w64-mingw32-lib

push: x86_64-w64-mingw32-build x86_64-w64-mingw32-lib
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-build
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-lib

