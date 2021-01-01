# Running S3V2-IDEAS in operating system set in Docker

### Outline
**[(1) Set up Docker package in a laptop or desktop.](#Set-up-Docker-package-in-a-laptop-or-desktop.)**<br>
#####
**[(2) Running S3V2-IDEAS in Docker](#Running-S3V2-IDEAS-in-Docker)**<br>
#####
**[(3) Extract the outputs in Docker](#Extract-the-outputs-in-Docker)**<br>
#####


####

## Set up Docker package in a laptop or desktop.

###### Download and Install Docker in the (https://www.docker.com/get-started)
###### Set the DResources that can be assigned to Docker
###### After open Docker package in MAC OS system, you should be able to see the following page.
###### Got to Setting resource page (highlighted by the red dash box).
<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/set_docker_resource1.png" width="800"/>

###### Change the default resource to the following settings (Or the resource your laptop or desktop can support).
###### Click the "Apply & Restart" button at the bottom right corner (highlighted by the red dash box).

<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/set_docker_resource2.png" width="800"/>

###### Open a Terminal in MAC and try the following command to see if the Docker has been successfully installed.

```
Docker system
```
###### If the Docker has been successfully installed, you should be able to see something similar to the followings:

<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/set_docker_resource3.png" width="800"/>


## 
## Running S3V2-IDEAS in Docker

#### Step1: Open a Terminal in MAC and Enter the Dockerfile folder
```
cd /Path_to_S3V2_folder/S3V2_IDEAS_ESMP/test_docker/
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

#### Step4: run S3V2-IDEAS pipeline on the testing datasets in Docker container
##### For the testing dataset, it will take around 2 hours to finish one run of the S3V2-IDEAS pipeline.
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/gsl/gsl/lib/
cd S3V2_IDEAS_ESMP/test_data
time bash run_S3V2_IDEAS_ESMP.sh
```

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
##### Here, the container ID of my test run is: 350c8c6d806c
##### Then, you can use the "docker cp 350c8c6d806c:..." command to extract the outputs to a local folder in your laptop or desktop. 
```
### mkdir ~/Documents/test_S3V2_docker_outputs/
### Get the bigWig files for the S3V2 normalized signals
docker cp 350c8c6d806c:/app/S3V2_IDEAS_ESMP/test_data/outputs/test_S3V2_IDEAS_pipeline_bws_RC ~/Documents/test_S3V2_docker_outputs/

### Get the IDEAS genome segmentation results (Heatmap for the Epigenetic states and the corresponding state genome browser Tracks)
docker cp 350c8c6d806c:/app/S3V2_IDEAS_ESMP/test_data/outputs/test_S3V2_IDEAS_pipeline_IDEAS_output ~/Documents/test_S3V2_docker_outputs/
```

