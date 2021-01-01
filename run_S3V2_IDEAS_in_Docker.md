# Running S3V2-IDEAS in Docker

### Outline
**[(1) Step0: Set up Docker package in a laptop or desktop.](#Step0-Set-up-Docker-package-in-a-laptop-or-desktop.)**<br>
#####
**[(2) Running S3V2-IDEAS in Docker](#Running-S3V2-IDEAS-in-Docker)**<br>
#####


####

## Step0 Set up Docker package in a laptop or desktop.

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

#### Step2: Set up docker container
```
docker image build -t test_s3v2 .
```

#### Step3: Start docker container
```
docker container run --rm -it -m 15G test_s3v2
```

#### Step4: run S3V2 in docker container
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/gsl/gsl/lib/
cd S3V2_IDEAS_ESMP/test_data
time bash run_S3V2_IDEAS_ESMP.sh
```

### get Docker container ID
docker container ls           
CONTAINER ID   IMAGE          COMMAND       CREATED          STATUS          PORTS     NAMES
350c8c6d806c   4d3a38fbe43f   "/bin/bash"   27 minutes ago   Up 27 minutes             wizardly_lamport

### Extract output from Docker container (container ID: 350c8c6d806c)
docker cp 350c8c6d806c:/app/S3V2_IDEAS_ESMP/test_data/outputs/test_S3V2_IDEAS_pipeline_IDEAS_output ~/Documents/2020_BG/test_S3V2_docker/
docker cp 350c8c6d806c:/app/S3V2_IDEAS_ESMP/test_data/outputs/test_S3V2_IDEAS_pipeline_bws ~/Documents/2020_BG/test_S3V2_docker/

#real	60m22.976s
#user	272m45.208s
#sys	2m27.606s

docker cp 725076b71d65:/app/S3V2_IDEAS_ESMP/test_data/input_bw_files_wg ~/Documents/2020_BG/test_S3V2_docker/
docker cp ~/Documents/2020_BG/test_S3V2_docker/input_bw_files_wg 725076b71d65:/app/S3V2_IDEAS_ESMP/test_data/