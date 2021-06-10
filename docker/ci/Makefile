CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
ARCH  ?= x86_64,aarch64,ppc64le,s390x
all: runtime build document

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

runtime: builder
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/runtime -f Dockerfile . --target runtime
build: builder
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/build -f Dockerfile . --target build
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/build:metis4 -f Dockerfile . --target build-metis4
document: builder
	docker buildx build --push --builder multiarch --platform x86_64 -t $(CI_REGISTRY_IMAGE)/document      -f Dockerfile      . --target document
builder:
	-docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
	-docker buildx create --name multiarch --driver docker-container
