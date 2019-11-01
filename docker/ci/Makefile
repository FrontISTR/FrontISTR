REGISTRY := registry.gitlab.com/frontistr-commons/frontistr

all: build document

login:
	docker login $(REGISTRY)

build:
	docker build -t $(REGISTRY)/build . --target build

document:
	docker build -t $(REGISTRY)/document . --target document

push: build document
	docker push $(REGISTRY)/build
	docker push $(REGISTRY)/document
