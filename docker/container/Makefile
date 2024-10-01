CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
ARCH  ?= x86_64

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

all: push

login:
ifeq ($(CI_JOB_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_JOB_TOKEN) $(CI_REGISTRY_IMAGE)
endif

push: login
	-docker buildx create --name multiarch --driver docker-container
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/fistr1:$(TAG) -f Dockerfile     --target fistr1 ../..
ifeq ($(TAG),master)
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/fistr1:latest -f Dockerfile     --target fistr1 ../..
endif
