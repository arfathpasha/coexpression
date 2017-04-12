TOP_DIR ?= /kb/dev_container
TOP_DIR_NAME = $(shell basename $(TOP_DIR))
DEPLOY_RUNTIME?=/kb/runtime
TARGET ?= /kb/deployment
SERVICE_SPEC = CoExpression.spec
SERVICE_NAME = CoExpression
SERVICE_DIR = $(TARGET)/services/$(SERVICE_NAME)
#SERVICE = CoExpressionService
SERVICE_URL = https://kbase.us/services/CoExpression
SERVICE_PORT = 7063

TPAGE = $(DEPLOY_RUNTIME)/bin/tpage
TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE_NAME) \
        --define kb_service_port=$(SERVICE_PORT)

DIR = $(shell pwd)
CONTAINER_DIR = $(shell readlink -f ../.. )
CONTAINER_DIR_NAME = $(shell basename $(CONTAINER_DIR))

LIB_DIR = lib
LBIN_DIR = bin
EXECUTABLE_SCRIPT_NAME = run_$(SERVICE_NAME).sh
TEST_DIR = ltest


default: compile build-executable-script-python

include $(TOP_DIR)/tools/Makefile.common
include $(TOP_DIR)/tools/Makefile.common.rules
####
# Warning: Inside Docker, KB_TOP is KB_TARGET, /kb/deployment instead of /kb/dev_container
#          But, the following assumes KB_TOP contains dev_container/bin for succussful make behavior
#          If everything is transitioned to kb_sdk, then it is fine...
build-executable-script-python: setup-local-dev-kb-py-libs
	mkdir -p $(LBIN_DIR)
	echo '#!/bin/bash' > $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_DEPLOYMENT_CONFIG="$(DIR)/deploy.cfg"' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_SERVICE_NAME="$(SERVICE_NAME)"' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export PYTHONPATH="$(DIR)/$(LIB_DIR)"' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'python $(DIR)/lib/biokbase/$(SERVICE_NAME)/$(SERVICE_NAME).py $$1 $$2 $$3' \
		>> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	chmod +x $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
ifeq ($(CONTAINER_DIR_NAME), dev_container)
	cp $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME) $(TOP_DIR)/bin/.
else
	cp $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME) $(KB_TOP)/bin/.
endif

setup-local-dev-kb-py-libs:
	touch lib/biokbase/__init__.py
	touch lib/biokbase/$(SERVICE_NAME)/__init__.py
ifeq ($(CONTAINER_DIR_NAME), dev_container)
	-rsync -vrh $(TOP_DIR)/modules/kbapi_common/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(TOP_DIR)/modules/auth/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(TOP_DIR)/modules/handle_service/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(TOP_DIR)/modules/workspace_deluxe/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(TOP_DIR)/modules/transform/lib/biokbase/* lib/biokbase/.
else
# Assume dev_container is on the same level of $KB_TOP, which works for KB_TOP=/mnt/kb/runtime and TOP_DIR=/mnt/kb/dev_container
# NEXT version will be using submodule, then no assumption
	-rsync -vrh $(KB_TOP)/../dev_container/modules/kbapi_common/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(KB_TOP)/../dev_container/modules/auth/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(KB_TOP)/../dev_container/modules/handle_service/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(KB_TOP)/../dev_container/modules/workspace_deluxe/lib/biokbase/* lib/biokbase/.
	-rsync -vrh $(KB_TOP)/../dev_container/modules/transform/lib/biokbase/* lib/biokbase/.
endif

update-R:
	export TARGET=$(TARGET) && export R_LIBS=$(TARGET)/lib && R --vanilla --slave -e "library('WGCNA')" &&  \
        if [ $$? -ne 0 ] ; then  $(DIR)/deps/WGCNA/install-r-packages.sh ; else echo "WGCNA is installed on $(TARGET)/lib"; fi
	
dk-build:
	docker build -t kbase/coex:test .

dk-bash:
	docker run -it --rm --entrypoint bash kbase/coex:test

deploy: deploy-scripts

deploy-scripts: deploy-libs deploy-executable-script 

deploy-service: deploy-libs deploy-executable-script deploy-service-scripts deploy-cfg

#deploy-libs:
#	@echo "Deploying libs to target: $(TARGET)"
#	mkdir -p $(TARGET)/lib/biokbase
#	rsync -vrh lib/biokbase/$(MODULE) $(TARGET)/lib/biokbase/.

deploy-executable-script:
	@echo "Installing executable scripts to target: $(TARGET)/bin"
	echo '#!/bin/bash' > $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_RUNTIME=$(DEPLOY_RUNTIME)' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export PYTHONPATH="$(TARGET)/lib"' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'python $(TARGET)/lib/biokbase/$(SERVICE_NAME)/$(SERVICE_NAME).py $$1 $$2 $$3' \
		>> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	chmod +x $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)

# Test Section


test: test-impl create-test-wrapper


test-impl: create-test-wrapper
	./$(TEST_DIR)/script_test/run_tests.sh
	coverage report -m
	cp .coverage work/
	mkdir -p work/kb/deployment/lib
	cp -R /kb/deployment/lib/biokbase work/kb/deployment/lib

create-test-wrapper:
	@echo "Creating test script wrapper"
ifeq ($(CONTAINER_DIR_NAME), dev_container)
	echo '#!/bin/bash' > $(TEST_DIR)/script_test/run_tests.sh
	echo 'export PATH=$(PATH):$(TARGET)/bin' >> $(TEST_DIR)/script_test/run_tests.sh
	echo 'export KB_RUNTIME=$(DEPLOY_RUNTIME)' >> $(TEST_DIR)/script_test/run_tests.sh
	echo 'export PYTHONPATH="$(DIR)/$(LIB_DIR)"' >> $(TEST_DIR)/script_test/run_tests.sh
else
	echo '#!/bin/bash' > $(TEST_DIR)/script_test/run_tests.sh
	echo 'export PATH=$(PATH):$(TARGET)/bin' >> $(TEST_DIR)/script_test/run_tests.sh
	echo 'export KB_RUNTIME=$(DEPLOY_RUNTIME)' >> $(TEST_DIR)/script_test/run_tests.sh
	echo 'export PYTHONPATH="$(TARGET)/lib"' >> $(TEST_DIR)/script_test/run_tests.sh
	echo 'export KB_SERVICE_NAME="$(SERVICE_NAME)"' >> $(TEST_DIR)/script_test/run_tests.sh
	echo 'export KB_DEPLOYMENT_CONFIG="$(DIR)/deploy.cfg"' >> $(TEST_DIR)/script_test/run_tests.sh # TODO: not sure about this line
endif
	echo 'coverage run $(DIR)/$(TEST_DIR)/script_test/basic_test.py $$1 $$2 $$3' \
		>> $(TEST_DIR)/script_test/run_tests.sh
#	chmod +x $(TEST_DIR)/script_test/run_tests.sh


compile:
	mkdir -p scripts; kb-sdk compile $(SERVICE_SPEC)\
		--out lib\
		--pyclname biokbase.$(SERVICE_NAME).$(SERVICE_NAME)Client \
		--pysrvname biokbase.$(SERVICE_NAME).$(SERVICE_NAME) \
		--pyimplname biokbase.$(SERVICE_NAME).$(SERVICE_NAME)Impl;
