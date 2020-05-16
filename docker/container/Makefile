CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
CI_COMMIT_REF_SLUG ?= _no_branch
CI_COMMIT_REF_SLUG := $(shell echo $(CI_COMMIT_REF_SLUG) | sed -e 's/\//__/')

all: build

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

build:
	docker build \
		-t $(CI_REGISTRY_IMAGE)/fistr1:$(CI_COMMIT_REF_SLUG) \
		-f Dockerfile \
		--target fistr1 \
		../..

push: login build
	docker push $(CI_REGISTRY_IMAGE)/fistr1:$(CI_COMMIT_REF_SLUG)

