# JointAnalysis

## docker
```shell
# OpenSwath
docker pull ghcr.io/openms/openms-executables:3.1.0

# Pyprophet
docker pull pyprophet/pyprophet:2.2.5

# MRGD
docker pull meiyulab/mrgd:v1
```
## JointAnalysis acquisition
JointAnalysis is containerized by Docker into an image, the installation tutorial of Docker is described in the [Docker documentation](https://docs.docker.com/engine) (both for Linux and Windows). On your machine, please start a Terminal (PowerShell) session and then execute the following command within the console:
```shell
docker pull meiyulab/mrgd:v1
```
This will take a few minutes to pull the Diamond image from [Docker Hub](https://hub.docker.com/r/zeroli/diamond/) to your machine. You can check whether the image `meiyulab/mrgd:v1` is successfully pulled by executing `docker images`, and if successfully, it will appear in the images list.  

## Container creation and startup
Create a container (named JointAnalysis_test) based on the image `meiyulab/mrgd:v1` and simultaneously mount the local folder `/path/to/JointAnalysis` to the folder `/path/to/JointAnalysis` (in the container) by running the following command in your terminal:
```shell
docker run -it --name JointAnalysis_test -v /path/to/JointAnalysis:/path/to/JointAnalysis meiyulab/mrgd:v1 bash
```
Please change /path/to/JointAnalysis to your own path. The path on your machine should be exactly the same as the path in the container.

After the above command is executed, you will enter the container. Please switch to the folder `/path/to/JointAnalysis` by executing `cd /path/to/JointAnalysis` in your terminal.

**Note:** Type in `exit` and press `Enter`, or hit `Ctrl+D` to exit the container. To re-enter the container after exiting, please follow the commands below :
```shell
docker start JointAnalysis_test
docker exec -it JointAnalysis_test bash
```
## Data analysis
The Nextflow script is saved as a `pipeline.nf` file in the `JointAnalysis` folder. JointAnalysis's execution commands are as follows.

Execute the following command in your terminal to start the analysis of MS data by providing an assay library:
```shell
nextflow run /path/to/JointAnalysis/pipeline.nf  --rawData "/path/to/JointAnalysis/data/rawdata/*.xzML" --library "/path/to/JointAnalysis/data/library/*.pqp" --tr_irt "/path/to/JointAnalysis/data/irt/*.tsv" --swath_windows_file "/path/to/JointAnalysis/data/win/*.tsv" --outputDir "/path/to/JointAnalysis/results"
```
Please change /path/to/JointAnalysis to your own path. Also the filename.

**Note:**  The `--outputDir` parameter specifies the storage location of the data processing intermediate results. The final peptide identification results are saved in the outdir folder, named `jointAnalysis_results.tsv`. Please refer to the **Help Message** section or execute `nextflow run /path/to/JointAnalysis/pipeline.nf --help` in the container to view the detailed information of parameter passing.

# Help Message
```shell
nextflow run /path/to/JointAnalysis/pipeline.nf --help
```
## Command: 
```
nextflow run /path/to/JointAnalysis/pipeline.nf --rawData "" --library "" --tr_irt "" --swath_windows_file "" --outputDir "" <Options> <Functions>
```


