CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr

all: build document

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

build:
	docker build -t $(CI_REGISTRY_IMAGE)/build . --target build

document:
	docker build -t $(CI_REGISTRY_IMAGE)/document . --target document

push: login build document
	docker push $(CI_REGISTRY_IMAGE)/build
	docker push $(CI_REGISTRY_IMAGE)/document
