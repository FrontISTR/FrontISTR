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
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/focal/runtime -f Dockerfile.focal . --target runtime
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/jammy/runtime -f Dockerfile.jammy . --target runtime
build: builder
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/focal/build -f Dockerfile.focal . --target build
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/focal/build:metis4 -f Dockerfile.focal . --target build-metis4
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/jammy/build -f Dockerfile.jammy . --target build
	docker buildx build --push --builder multiarch --platform $(ARCH) -t $(CI_REGISTRY_IMAGE)/jammy/build:metis4 -f Dockerfile.jammy . --target build-metis4
document: builder
	docker buildx build --push --builder multiarch --platform x86_64 -t $(CI_REGISTRY_IMAGE)/document      -f Dockerfile      . --target document
builder:
	-docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
	-docker buildx create --name multiarch --driver docker-container
