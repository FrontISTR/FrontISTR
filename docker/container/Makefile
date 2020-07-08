CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr

ifdef CI_COMMIT_TAG
  TAG ?= $(CI_COMMIT_TAG)
else
  ifdef CI_COMMIT_REF_SLUG
    TAG ?= $(CI_COMMIT_REF_SLUG)
    TAG := $(shell echo $(TAG) | sed -e 's/\//__/')
  else
    TAG := _no_branch
  endif
endif

all: build

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

build:
	docker pull registry.gitlab.com/frontistr-commons/frontistr/build:latest
	docker pull registry.gitlab.com/frontistr-commons/frontistr/runtime:latest
	docker build \
		-t $(CI_REGISTRY_IMAGE)/fistr1:$(TAG) \
		-f Dockerfile \
		--target fistr1 \
		../..

push: login build
	docker push $(CI_REGISTRY_IMAGE)/fistr1:$(TAG)

