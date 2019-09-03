REGISTORY := registry.ritc.jp/ricos/frontistr

all: build document

login:
	docker login $(REGISTORY)

build:
	docker build -t $(REGISTORY)/build . --target build

document:
	docker build -t $(REGISTORY)/document . --target document

push: build document
	docker push $(REGISTORY)/build
	docker push $(REGISTORY)/document
