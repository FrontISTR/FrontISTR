REGISTORY := registry.gitlab.com/frontistr-commons/frontistr

all: cross mingw64lib

login:
	docker login $(REGISTORY)

cross:
	docker build -t $(REGISTORY)/cross -f mingw64.Dockerfile . --target cross

mingw64lib:
	docker build -t $(REGISTORY)/mingw64lib mingw64.Dockerfile . --target mingw64lib

push: cross mingw64lib login
	docker push $(REGISTORY)/cross
	docker push $(REGISTORY)/mingw64lib
