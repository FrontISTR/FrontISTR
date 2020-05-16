CI_REGISTRY_IMAGE  ?= registry.gitlab.com/frontistr-commons/frontistr

all: x86_64-w64-mingw32-build x86_64-w64-mingw32-lib x86_64-w64-mingw32-mpilib

login:
ifeq ($(CI_BUILD_TOKEN),)
	docker login $(CI_REGISTRY_IMAGE)
else
	docker login -u gitlab-ci-token -p $(CI_BUILD_TOKEN) $(CI_REGISTRY_IMAGE)
endif

x86_64-w64-mingw32-build:
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-build -f x86_64-w64-mingw32.Dockerfile . --target x86_64-w64-mingw32-build
x86_64-w64-mingw32-lib:
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-lib   -f x86_64-w64-mingw32.Dockerfile . --target x86_64-w64-mingw32-lib
x86_64-w64-mingw32-mpilib:
#	curl -O https://download.microsoft.com/download/a/5/2/a5207ca5-1203-491a-8fb8-906fd68ae623/msmpisetup.exe
#	curl -O https://download.microsoft.com/download/a/5/2/a5207ca5-1203-491a-8fb8-906fd68ae623/msmpisdk.msi
	docker build -t $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-mpilib   -f x86_64-w64-mingw32.Dockerfile . --target x86_64-w64-mingw32-mpilib

push: login x86_64-w64-mingw32-build x86_64-w64-mingw32-lib x86_64-w64-mingw32-mpilib
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-build
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-lib
	docker push $(CI_REGISTRY_IMAGE)/x86_64-w64-mingw32-mpilib

