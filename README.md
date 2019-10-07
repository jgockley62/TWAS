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
 docker run -v "/home/${USER}/TWAS/:/home/${USER}/TWAS/" -e USER=${USER} -e PASSWORD=HornetHockey -d -p 8787:8787 <IMAGE_ID>
```
