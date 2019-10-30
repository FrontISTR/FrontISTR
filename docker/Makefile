REGISTORY := registry.ritc.jp/ricos/frontistr
CI_COMMIT_REF_NAME ?= manual_deploy

all: build document

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(REGISTRY)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(REGISTRY)
endif

# Image for building FrontISTR binary
build:
	docker build -t $(REGISTORY)/build . --target build

# Image for document generation
document:
	docker build -t $(REGISTORY)/document . --target document

push: build document
	docker push $(REGISTORY)/build
	docker push $(REGISTORY)/document

# Image which includes fistr1 binary
fistr1:
	docker build -t $(REGISTORY)/fistr1:$(CI_COMMIT_REF_NAME) .. --target fistr1 -f $(PWD)/Dockerfile

push_fistr1: login fistr1
	docker push $(REGISTORY)/fistr1:$(CI_COMMIT_REF_NAME)
