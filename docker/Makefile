REGISTRY := registry.ritc.jp/ricos/frontistr
CI_COMMIT_REF_NAME ?= _manual

all: build document

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(REGISTRY)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(REGISTRY)
endif

# Image for building FrontISTR binary
build:
	docker build -t $(REGISTRY)/build . --target build

# Image for document generation
document:
	docker build -t $(REGISTRY)/document . --target document

push: build document
	docker push $(REGISTRY)/build
	docker push $(REGISTRY)/document

# Image which includes fistr1 binary
fistr1:
	docker build -t $(REGISTRY)/fistr1:$(CI_COMMIT_REF_NAME) .. --target fistr1 -f $(PWD)/Dockerfile

push_fistr1: login fistr1
	docker push $(REGISTRY)/fistr1:$(CI_COMMIT_REF_NAME)
