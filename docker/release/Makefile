REGISTORY := registry.gitlab.com/frontistr-commons/frontistr

all: x86_64-w64-mingw32-build x86_64-w64-mingw32-lib

login:
	docker login $(REGISTORY)

x86_64-w64-mingw32-build:
	docker build -t $(REGISTORY)/x86_64-w64-mingw32-build . --target x86_64-w64-mingw32-build

x86_64-w64-mingw32-lib:
	docker build -t $(REGISTORY)/x86_64-w64-mingw32-lib . --target x86_64-w64-mingw32-lib

push: x86_64-w64-mingw32-build x86_64-w64-mingw32-lib
	docker push $(REGISTORY)/x86_64-w64-mingw32-build
	docker push $(REGISTORY)/x86_64-w64-mingw32-lib

