---
title: "CMC_DataCleaning"
author: 'JKG'
output: html_notebook
editor_options:
  chunk_output_type: console
---

Date of analysis update: "Tue Dec 17 16:12:51 2019"

| *Synapse ID* | *File Name* |
|  -----------------------------------------  |   ---------                      |
| syn21385083 | | 
| syn20835005 | |
| syn20835007 | |
| syn20834996 | |
| syn10887100 | |
| syn10887068 | |
| syn10887117 | |
| syn10887048 | |
| syn10887106 | |
| syn10887111 | |
| syn10887027 | |
| syn10887033 | |
| syn10887101 | |
| syn10887029 | |
| syn10887044 | |
| syn11273049 | |
| syn10887084 | |
| syn10887091 | |
| syn10887078 | |
| syn10887115 | |
| syn10887041 | |
| syn10887071 | |
| syn10887037 | |
| syn10887103 | |
| syn10887102 | |
| syn10887030 | |
| syn11273048 | |
| syn10887059 | |
| syn10887040 | |
| syn10887057 | |
| syn10887065 | |
| syn10887055 | |
| syn10887119 | |
| syn10887073 | |
| syn10887050 | |
| syn10887014 | |
| syn10887021 | |
| syn10887062 | |
| syn10887107 | |
| syn10887096 | |
| syn11273047 | |
| syn10887049 | |
| syn10887058 | |
| syn10887020 | |
| syn10887031 | |
| syn10887108 | |
| syn10887056 | |
| syn11273044 | |
| syn10887032 | |
| syn10887094 | |
| syn10887069 | |
| syn10887077 | |
| syn10887019 | |
| syn10887104 | |


##Before running 
Run the following in docker exec -it <CONTAINER> /bin/bash 
chmod 777 /root/
mv TWAS/ /home/<USR>/TWAS/
chown -hR <GID>:<USR> /home/<USR>/TWAS/
Create a ~/.synapseConfig file in the container, but do not push a container containing your credentials file to Docker Hub!!!

## Login to Synapse Will need to replace with your credentials

```bash
#source /root/.bashrc
synapse login -u <USR> -p <PSWD> --rememberMe

#Alternativly you can setup a credentials file as such:
touch ~/.synapseConfig
echo "[authentication]" >> ~/.synapseConfig
echo "username = <USR>" >> ~/.synapseConfig
echo "password = <PASWD>" >> ~/.synapseConfig
```


```bash
mkdir ~/TWAS/bin
wget -q -P ~/TWAS/bin http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip 
unzip -d ~/TWAS/bin/ ~/TWAS/bin/plink_linux_x86_64_20190617.zip
wget -q -P ~/TWAS/bin http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20191015.zip
unzip -d ~/TWAS/bin/ ~/TWAS/bin/plink2_linux_x86_64_20191015.zip

rm ~/TWAS/bin/*.zip
rm ~/TWAS/bin/toy.*
rm ~/TWAS/bin/LICENSE
```

## Pull CMC Data From Synapse and Prep for Imputation - 
This data is protected, you wil be required to apply for access through PhsycEncode

```bash
CORES=`expr $(nproc) - 2`
#Pull Data with needed SNPS only and fiter GTEx for that SNP List
synapse get syn20835005
synapse get syn20835007
synapse get syn20834996
~/TWAS/bin/plink --threads $CORES --bfile All_CEU_ToTrainForTWAS --write-snplist

mkdir CMC_Genos
synapse get syn10537112 --recursive --downloadLocation CMC_Genos/

mkdir CMC_Filt_Genos
touch merge.txt
#Pull out Snps
for i in {1..22}
  do
    #--positions plink.snplist 
    bcftools norm --threads $CORES --rm-dup snps -Oz CMC_Genos/CMC_chr$i.dose.vcf.gz > CMC_Filt_Genos/Filt_CMC_chr$i.dose.vcf.gz
    vcftools --gzvcf CMC_Filt_Genos/Filt_CMC_chr$i.dose.vcf.gz --plink --out CMC_Filt_Genos/Init_Filtered_chr$i
    echo 'CMC_Filt_Genos/Init_Filtered_chr'$i'.ped CMC_Filt_Genos/Init_Filtered_chr'$i'.map' >> merge.txt
  done
grep -v 'chr1.ped' merge.txt > foo
mv foo merge.txt

#Merge Genotype
~/TWAS/bin/plink --threads $CORES --file CMC_Filt_Genos/Init_Filtered_chr1 --merge-list merge.txt --make-bed --out Merged_Initial_CMC_Geno

#Remove old Map/Ped Files
rm CMC_Filt_Genos/Init_Filtered_chr*.map
rm CMC_Filt_Genos/Init_Filtered_chr*.ped
```

### ReName the MSBB SNPS so the can be combined


## Pull Imputation SNPS - 

```bash
CORES=`expr $(nproc) - 2`

mv Merged_Initial_CMC_GenoRenamed.bim Merged_Initial_CMC_Geno.bim
~/TWAS/bin/plink --threads $CORES --bfile Merged_Initial_CMC_Geno --extract plink.snplist --make-bed --out CMC_Genos_For_Impute

```





## Clean


```bash
rm *.fam
rm *.bed
rm *.bim
rm *.fam
rm *.nosex
rm *.ped
rm *.map
rm *.missnp
rm *.log
rm *.snplist
rm merge.txt
rm -r CMC_Genos/
rm -r CMC_Filt_Genos/
```

### R Source Code
[Github](https://github.com/jgockley62/TWAS/blob/4c7b0da35fd2d06639b59f360717c8e9bd5b6cc1/code/CMC_DataProcessing.Rmd)


