FROM ubuntu:16.04

##S3 awsmfa --identity-profile bridge.prod.jgockley --target-profile jgockley@sandbox --token $658059
#_From Dockerfile containing folder run: CMD# docker build -t test ~/Software/Docker/twas/ CMD# docker run -it 
#<IMAGE_ID> /bin/bash

RUN apt-get update 
RUN apt-get install -y --no-install-recommends apt-utils libblas-dev liblapack-dev gfortran r-base aptitude build-essential libcurl4-openssl-dev libxml2-dev libxml2-dev libssl-dev nano python-dev plink git curl
#R Packages
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile 

RUN Rscript -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && Rscript -e "install.packages('RColorBrewer', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && Rscript -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && Rscript -e "install.packages('glmnet', dependencies=TRUE)" \
 && Rscript -e "install.packages('methods', dependencies=TRUE)" \
 && Rscript -e "devtools::install_github('RcppCore/RcppEigen', dependencies=TRUE)" \
 && Rscript -e "devtools::install_github('gabraham/plink2R/plink2R', dependencies=TRUE)" \
 && Rscript -e "install.packages('data.table', dependencies=TRUE)" \
 && Rscript -e "install.packages('foreach', dependencies=TRUE)" \
 && Rscript -e "install.packages('doMC', dependencies=TRUE)"\
 && Rscript -e "source('http://depot.sagebase.org/CRAN.R')" \
 && Rscript -e "install.packages('glmnet', repos = 'http://cran.us.r-project.org')" \
 && Rscript -e "install.packages('synapseClient')" 

RUN apt-get install -y wget nano cdda2wav cdrecord mkisofs python-subprocess32 

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \ 
 && python get-pip.py 
#RUN apt-get install -y multiprocessing

RUN pip install synapseclient \
 && pip install awscli \
 && pip install aws-mfa \
 && pip install awsmfa \
 && pip install pyliftover \
 && pip install pprint \
 && pip install Biopython \
 && pip install bitarray \
 && pip install scipy \
 && pip install pandas
 #&& pip install pybedtools

#Command Line plink2 capability - Works! (need to try binary on AWS for compatability test)
RUN wget http://s3.amazonaws.com/plink2-assets/alpha1/plink2_linux_x86_64.zip \
  && unzip plink2_linux_x86_64.zip -d /usr/local/bin/ \
  && rm -f plink2_linux_x86_64.zip
RUN wget https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz \
  && gunzip GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz \
  && wget https://cnsgenomics.com/software/gcta/bin/gcta_1.92.4beta.zip \
  && unzip gcta_1.92.4beta.zip

# Get the fusion software - Works!
RUN wget https://github.com/gusevlab/fusion_twas/archive/master.zip \
 && unzip master.zip \
 && cd fusion_twas-master \
 && wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2 \
 && tar xjvf LDREF.tar.bz2 \
 && wget https://github.com/gabraham/plink2R/archive/master.zip \
 && unzip master.zip

#Get predefined weights - Works!
RUN wget https://data.broadinstitute.org/alkesgroup/FUSION/SUM/PGC2.SCZ.sumstats \
 && mkdir WEIGHTS \
 && cd WEIGHTS \
 && wget https://data.broadinstitute.org/alkesgroup/FUSION/WGT/GTEx.Whole_Blood.tar.bz2 \
 && tar -xjf GTEx.Whole_Blood.tar.bz2

#Get Gemma - Works! ENV GEMMA_URL 
#https://github.com/genetics-statistics/GEMMA/releases/download/v0.97/gemma-git-gn2-0.97-c760aa0-gener$ ENV 
#GEMMA_NAME_ARCHIVE gemma-git-gn2-0.97-c760aa0-generic-x86_64.tar.bz2 ENV GEMMA_PATH /opt/gemma/ RUN mkdir /opt RUN 
#mkdir /opt/gemma WORKDIR ${GEMMA_PATH} RUN wget ${GEMMA_URL} && \
# tar xvjf ${GEMMA_NAME_ARCHIVE} && \ sh install.sh ${GEMMA_PATH}/install && \ cp ${GEMMA_PATH}/install/bin/* 
# /usr/local/bin

RUN git clone https://github.com/gusevlab/fusion_twas.git \
 && wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2 \
 && tar -xvf LDREF.tar.bz2 \
 && git clone https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline.git \
 && wget https://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz \
 && tar -xvf gtool_v0.7.5_x86_64.tgz \
 && wget http://sourceforge.net/projects/fcgene/files/fcgene-1.0.7.tar.gz \
 && tar -zxvf fcgene-1.0.7.tar.gz \
 && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip \ 
 && unzip plink_linux_x86_64_20190304.zip 

RUN wget https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz 

RUN wget -O /fusion_twas-master/gemma-0.98.1-linux-static.gz https://github.com/genetics-statistics/GEMMA/releases/download/0.98.1/gemma-0.98.1-linux-static.gz \
 && gunzip /fusion_twas-master/gemma-0.98.1-linux-static.gz \
 && chmod 777 /fusion_twas-master/gemma-0.98.1-linux-static

RUN wget -O /fusion_twas-master/plink_linux_x86_64_20190304.zip http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip \
 && unzip -o /fusion_twas-master/plink_linux_x86_64_20190304.zip 
# && rm /fusion_twas-master/toy.*


