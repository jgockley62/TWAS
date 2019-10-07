# TWAS
## Preprocessing of TWAS Data:

### Reports on all preprocessing steps can be seen in the respective subfolders in the synapse project here:
[Synapse Data Project](https://www.synapse.org/#!Synapse:syn18936948/files/)

### All Expression data is re-processed to correct for diagnosis with the following scripts
* code/MSBB_RegressDiagnosis_4_TWAS.Rmd
* code/Mayo_RegressDiagnosis_4_TWAS.Rmd
* code/ROSMAP_RegressDiagnosis_4_TWAS.Rmd

### All Variant data processing and Ancestry Culstering is preformed with:
* code/Ancestry_PCA.Rmd

### Harmonization of Sample names and RNA-Seq Dataframe Finalization can be found here:
* code/ProcessData.Rmd

## Running Code

Log into EC2 instance
```{bash}
 #Pull Repo
 git clone https://github.com/jgockley62/TWAS.git
 
 #Build Image from Dockerfile
 docker image build --build-arg USER_ID=$(id -u ${USER}) --build-arg GROUP_ID=$(id -g ${USER}) -t twas /home/${USER}/TWAS/
 
 #Open Docker container and RInstance
 docker run -v "/home/${USER}/TWAS/:/home/${USER}/TWAS/" -e USER=${USER} -e PASSWORD=<PASSWORD> -d -p 8787:8787 <IMAGE_ID>
```

Now loginto RStudio Instance
https://<AWS_Instance_IP>:8787

```{r}
setwd("~/TWAS/code/")

#Sample Code Run Command:

synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

library(data.table)

source("~/TWAS/utilityFunctions/knitfile2synapseClient.R")
source("~/TWAS/utilityFunctions/hook_synapseMdSyntax_plot.R")

createAndKnitToFolderEntityClient(file = "ProcessData.Rmd",
                                          parentId ="syn18936948",
                                          folderName = 'Processed_TWAS_Training_Data')


```

## Train Weights on AWS Instance

docker build -t train /home/${USER}/TWAS_Training/

screen

docker run -it <Image_ID> /bin/bash

pip install awscli
pip install aws-mfa
pip install awsmfa
apt-get install -y python-subprocess32
apt-get install -y multiprocessing
pip install scipy
pip install pandas
apt-get install python-dev
pip install bitarray
pip install pyliftover
pip install pprint
pip install Biopython


wget https://github.com/genetics-statistics/GEMMA/releases/download/0.98.1/gemma-0.98.1-linux-static.gz
gunzip gemma-0.98.1-linux-static.gz
chmod 777 gemma-0.98.1-linux-static


wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip
unzip plink_linux_x86_64_20190304.zip
rm toy.*

