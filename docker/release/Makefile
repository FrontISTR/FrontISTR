REGISTORY := registry.gitlab.com/frontistr-commons/frontistr

all: cross mingw64lib

login:
	docker login $(REGISTORY)

cross:
	docker build -t $(REGISTORY)/cross . --target cross

mingw64lib:
	docker build -t $(REGISTORY)/mingw64lib . --target mingw64lib

push: cross mingw64lib
	docker push $(REGISTORY)/cross
	docker push $(REGISTORY)/mingw64lib
