CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
CI_COMMIT_REF_NAME ?= _no_branch

all: build

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

build:
	docker build \
		-t $(CI_REGISTRY_IMAGE)/fistr1:$(CI_COMMIT_REF_NAME) \
		-f Dockerfile \
		--target fistr1 \
		../..

push: login build
	docker push $(CI_REGISTRY_IMAGE)/fistr1:$(CI_COMMIT_REF_NAME)

