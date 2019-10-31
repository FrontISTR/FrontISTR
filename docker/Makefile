REGISTORY := registry.gitlab.com/frontistr-commons/frontistr

all: build document crossmake

login:
	docker login $(REGISTORY)

build:
	docker build -t $(REGISTORY)/build . --target build

document:
	docker build -t $(REGISTORY)/document . --target document

cross:
	docker build -t $(REGISTORY)/cross . --target cross

push: build document
	docker push $(REGISTORY)/build
	docker push $(REGISTORY)/document
	docker push $(REGISTORY)/cross
