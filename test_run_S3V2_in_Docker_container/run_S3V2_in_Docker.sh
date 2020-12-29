### enter the Dockerfile folder
cd /Path_to_S3V2_folder/S3V2_IDEAS_ESMP/test_docker/

### set up docker container
docker image build -t test_s3v2 .

### start docker container
docker container run --rm -it -m 15G test_s3v2

### run S3V2 in docker container
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/gsl/gsl/lib/
cd S3V2_IDEAS_ESMP/test_data
time bash run_S3V2_IDEAS_ESMP.sh

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