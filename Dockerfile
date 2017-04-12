FROM kbase/kbase:sdkbase.latest 
MAINTAINER KBase Developer 

RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf transform && \
  git clone https://github.com/kbase/transform && \
  rm -rf workspace_deluxe && \
  git clone https://github.com/kbase/workspace_deluxe && \
  cd /kb/dev_container/modules/workspace_deluxe && \
  cp -rv lib/* /kb/deployment/lib/

RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf coex_helper && \
  git clone https://github.com/sjyoo/coex_helper && \
  cd /kb/dev_container/modules/coex_helper && \
  make update-R

####END OF KBASE #############################
#apt-get update && apt-get install -y ant && \
# -----------------------------------------
# Insert apt-get instructions here to install
# any required dependencies for your module.
# -----------------------------------------
#RUN \
#  . /kb/dev_container/user-env.sh && \
#  cd /kb/dev_container/modules/coex_helper && \
#  make update-R 

RUN \
  apt-get update && \
  apt-get install -y bzip2 \
                     gcc \
		     ncurses-dev \
		     tofrodos \
		     unzip 
RUN pip install mpipe
RUN pip install pandas numpy
RUN pip install coverage
WORKDIR /kb/module
COPY ./deps /kb/deps
COPY ./ /kb/module
RUN chmod -R a+rw /kb/module
RUN chmod +x /kb/module/ltest/script_test/run_tests.sh
# Windows compatibility line
#RUN bash -c "for i in `find . -name '*.sh'`; do dos2unix -v $i; done"
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/module && \
  make && make deploy
#  rm -rf narrative_method_store \
ENV PATH=$PATH:/kb/dev_container/modules/kb_sdk/bin
WORKDIR /kb/module
RUN mkdir -p /kb/module/work
ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]
