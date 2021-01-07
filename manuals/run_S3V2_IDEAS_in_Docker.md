# Running S3V2-IDEAS in operating system set in Docker

### Outline
**[(1) Set up Docker package in a laptop or desktop.](#Set-up-Docker-package-in-a-laptop-or-desktop.)**<br>
#####
**[(2) Running S3V2-IDEAS in Docker](#Running-S3V2-IDEAS-in-Docker)**<br>
#####
**[(3) Extract the outputs in Docker](#Extract-the-outputs-in-Docker)**<br>
#####
**[(4) Download whole genome data and run the pipeline in Docker](#Download-whole-genome-data-and-run-the-pipeline-in-Docker)**<br>
#####


####

## Set up Docker package in a laptop or desktop.

###### Download and Install Docker in the (https://www.docker.com/get-started)
###### Set the DResources that can be assigned to Docker
###### After open Docker package in MAC OS system, you should be able to see the following page.
###### Got to Setting resource page (highlighted by the red dash box).
<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/set_docker_resource1.png" width="500"/>

###### Change the default resource to the following settings (Or the resource your laptop or desktop can support).
#### For memory, at least 5GB is required to run the pipeline on the testing data.
###### Click the "Apply & Restart" button at the bottom right corner (highlighted by the red dash box).

<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/set_docker_resource2.png" width="500"/>

###### Open a Terminal in MAC and try the following command to see if the Docker has been successfully installed.

```
Docker system
```
###### If the Docker has been successfully installed, you should be able to see something similar to the followings:

<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/set_docker_resource3.png" width="300"/>


## 
## Running S3V2-IDEAS in Docker

#### Step1: Open a Terminal in MAC. --> Clone the S3V2_IDEAS_ESMP pipeline from GitHub. --> Then Enter the "test_run_S3V2_in_Docker_container" folder include the Docker image file in the S3V2_IDEAS_ESMP pipeline folder.
```
### Enter a local folder
cd /Path_to_S3V2_folder/

### clone S3V2_IDEAS_ESMP pipeline
git clone https://github.com/guanjue/S3V2_IDEAS_ESMP.git

### Enter the "test_run_S3V2_in_Docker_container" folder with Docker image file
cd /Path_to_S3V2_folder/S3V2_IDEAS_ESMP/test_run_S3V2_in_Docker_container/
```

#### Step2: Set up Docker container. 
##### This step may take a few minutes to install all of the Dependencies in the Linux-64 operating system in Docker container
```
docker image build -t test_s3v2 .
```

#### Step3: Start the Docker container
##### This command will start interaction mode within the Linux-64 operating system in Docker container
```
docker container run --rm -it -m 15G test_s3v2
```
##### Noted! For request 5GB is enough for just running the testing data. So user can simply change the "15G" to "5G" in this command.


#### Step4: run S3V2-IDEAS pipeline on the testing datasets in Docker container.
##### A Docker container is similar to a Terminal. 
##### The user can type the following command to run the S3V2-IDEAS pipeline for the testing data in the Docker container
```
### set the GSL PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/gsl/gsl/lib/
### enter the 'test_data' folder
cd S3V2_IDEAS_ESMP/test_data
### run the S3V2-IDEAS pipeline
time bash run_S3V2_IDEAS_ESMP.sh
```
##### For the testing dataset, it will take around 2 hours to finish one run of the S3V2-IDEAS pipeline.


## Extract the outputs in Docker 
##### After the S3V2-IDEAS run is finished, you can extract the outputs generated in the Docker container.
##### Open a new tab in the Terminal and Get Docker container ID by the following command:
```
docker container ls      
```
##### You should see something like the following containing the "CONTAINER ID":
```
CONTAINER ID   IMAGE          COMMAND       CREATED          STATUS          PORTS     NAMES
350c8c6d806c   4d3a38fbe43f   "/bin/bash"   27 minutes ago   Up 27 minutes             wizardly_lamport
```

##### Extract output from Docker container 
###### Here, the container ID of my test run is: 350c8c6d806c
###### Then, you can use the "docker cp 350c8c6d806c:..." command to extract the outputs to a local folder in your laptop or desktop. 
```
### mkdir ~/Documents/test_S3V2_docker_outputs/
### Get the bigWig files for the S3V2 normalized signals
docker cp 350c8c6d806c:/app/S3V2_IDEAS_ESMP/test_data/outputs/test_S3V2_IDEAS_pipeline_bws_RC ~/Documents/test_S3V2_docker_outputs/

### Get the IDEAS genome segmentation results (Heatmap for the Epigenetic states and the corresponding state genome browser Tracks)
docker cp 350c8c6d806c:/app/S3V2_IDEAS_ESMP/test_data/outputs/test_S3V2_IDEAS_pipeline_IDEAS_output ~/Documents/test_S3V2_docker_outputs/
```

## Download whole genome data and run the pipeline in Docker 
##### For testing the whole genome data analysis, user can also use the script in "get_hg38wgbw.sh" in the following link to download some whole genome data into the Docker container. (https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/test_data/)
##### Then, use the script in "run_S3V2_IDEAS_ESMP.hg38wg.sh" to run the S3V2-IDEAS pipeline in the whole genome data
```
### set the GSL PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/gsl/gsl/lib/
### enter the 'test_data' folder
cd S3V2_IDEAS_ESMP/test_data
### Download the bigWig files into a "input_bw_files_wg/" folder in the Docker container
time bash get_hg38wgbw.sh
### run the S3V2-IDEAS pipeline
time bash run_S3V2_IDEAS_ESMP.hg38wg.sh
```
##### For the user's own data analysis, user can also modified the "get_hg38wgbw.sh" script to download the data into Docker container.
##### Or directly cp into the Docker container from the local folder:
```
### cp local bigWig files into the Docker container
docker cp ~/local_folder/some_ct_mk.bigWig 350c8c6d806c:/app/S3V2_IDEAS_ESMP/test_data/input_bw_files_wg/ 
```
##### Then, user need to modified the "metadata.forEScall.hg38wg.txt" and "run_S3V2_IDEAS_ESMP.hg38wg.sh" files to run the pipeline on their own datasets.
##### The instruction for setting up the "metadata.forEScall.hg38wg.txt" and "run_S3V2_IDEAS_ESMP.hg38wg.sh" files can be found in the following page:
**[Inputs-ParameterSettings-RunningSteps-Outputs for S3V2_IDEAS_ESMP](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/manuals/inoutput_for_S3V2_IDEAS_pipeline.md)**<br>




