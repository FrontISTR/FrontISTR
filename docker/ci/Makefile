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
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/runtime:ubuntu2004 -f Dockerfile.ubuntu2004 . --target runtime
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/runtime:ubuntu2204 -f Dockerfile.ubuntu2204 . --target runtime
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/runtime:ubuntu2004-metis4 -f Dockerfile.ubuntu2004 . --target runtime
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/runtime:ubuntu2204-metis4 -f Dockerfile.ubuntu2204 . --target runtime
build: builder
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/build:ubuntu2004 -f Dockerfile.ubuntu2004 . --target build
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/build:ubuntu2004-metis4 -f Dockerfile.ubuntu2004 . --target build-metis4
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/build:ubuntu2204 -f Dockerfile.ubuntu2204 . --target build
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/build:ubuntu2204-metis4 -f Dockerfile.ubuntu2204 . --target build-metis4
document: builder
	docker buildx build --push --builder multiarch --platform x86_64 -t $(CI_REGISTRY_IMAGE)/document      -f Dockerfile.ubuntu2004      . --target document
builder:
	-docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
	-docker buildx create --name multiarch --driver docker-container
