CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
CI_COMMIT_REF_NAME ?= _no_branch

all: cross mingw64lib fistr1

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(REGISTRY)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

fistr1:
	docker build \
		-t $(CI_REGISTRY_IMAGE)/fistr1:$(CI_COMMIT_REF_NAME) \
		-f fistr1.Dockerfile \
		--target fistr1 \
		../..

push_fistr1: login fistr1
	docker push $(CI_REGISTRY_IMAGE)/fistr1:$(CI_COMMIT_REF_NAME)

cross:
	docker build \
		-t $(CI_REGISTRY_IMAGE)/cross:$(CI_COMMIT_REF_NAME) \
		-f mingw64.Dockerfile \
		--target cross \
		.

mingw64lib:
	docker build \
		-t $(CI_REGISTRY_IMAGE)/mingw64lib:$(CI_COMMIT_REF_NAME) \
		-f mingw64.Dockerfile \
		--target mingw64lib \
		.

push_mingw64: cross mingw64lib login
	docker push $(CI_REGISTRY_IMAGE)/cross:$(CI_COMMIT_REF_NAME)
	docker push $(CI_REGISTRY_IMAGE)/mingw64lib:$(CI_COMMIT_REF_NAME)
