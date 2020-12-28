### enter the Dockerfile folder
cd /Path_to_S3V2_folder/S3V2_IDEAS_ESMP/test_docker/

### set up docker container
docker image build -t test_s3v2 .

### start docker container
docker container run --rm -it -m 10G test_s3v2

### run S3V2 in docker container
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/gsl/gsl/lib/
cd S3V2_IDEAS_ESMP/test_data
time bash run_S3V2_IDEAS_ESMP.sh

### get Docker container ID
docker container ls

### Extract output from Docker container (container ID: 5f5a44199278)
docker cp 5f5a44199278:/app/S3V2_IDEAS_ESMP/test_data/outputs/test_S3V2_IDEAS_pipeline_IDEAS_output ~/Documents/2020_BG/test_S3V2_docker/

#real	60m22.976s
#user	272m45.208s
#sys	2m27.606s
