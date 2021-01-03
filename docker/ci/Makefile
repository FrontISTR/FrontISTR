CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr
all: runtime build document

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

runtime:
	docker build -t $(CI_REGISTRY_IMAGE)/runtime . --target runtime
build:
	docker build -t $(CI_REGISTRY_IMAGE)/build . --target build
	docker build -t $(CI_REGISTRY_IMAGE)/build:metis4 . --target build-metis4
document:
	docker build -t $(CI_REGISTRY_IMAGE)/document . --target document
push: login build document
	docker push $(CI_REGISTRY_IMAGE)/document
	docker push $(CI_REGISTRY_IMAGE)/build
	docker push $(CI_REGISTRY_IMAGE)/runtime
