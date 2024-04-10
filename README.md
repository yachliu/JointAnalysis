# JointAnalysis

## docker
```shell
# JointAnalysis
docker pull meiyulab/jointanalysis:v1

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
docker pull meiyulab/jointanalysis:v1
```
This will take a few minutes to pull the Diamond image from [Docker Hub](https://hub.docker.com/r/zeroli/diamond/) to your machine. You can check whether the image `meiyulab/jointanalysis:v1` is successfully pulled by executing `docker images`, and if successfully, it will appear in the images list.  

## Container creation and startup
Create a container (named JointAnalysis_test) based on the image `meiyulab/analysis:v1` and simultaneously mount the local folder `/path/to/JointAnalysis` to the folder `/path/to/JointAnalysis` (in the container) by running the following command in your terminal:
```shell
docker run -it --name JointAnalysis_test -v /path/to/JointAnalysis:/path/to/JointAnalysis meiyulab/analysis:v1 bash
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

## Help Message
```shell
nextflow run /path/to/JointAnalysis/pipeline.nf --help
```
## Command: 
```
nextflow run /path/to/JointAnalysis/pipeline.nf --rawData "" --library "" --tr_irt "" --swath_windows_file "" --outputDir "" <Options> <Functions>
```
## Parameters descriptions

### Mandatory arguments
|parameters|descriptions|
|---|---|
|--rawData|Path to raw data. For example: --rawData "/path/to/JointAnalysis/data/rawdata/*.xzML".|
|--library|Library file path. For example: --library "/path/to/JointAnalysis/data/library/*.pqp".|
|--tr_irt|IRT file path. For example: --tr_irt "/path/to/JointAnalysis/data/irt/*.tsv".|
|--swath_windows_file|Swath windows file path. For example: --swath_windows_file "/path/to/JointAnalysis/data/win/*.tsv".|
|--outputDir|Output directory for results. For example: --outputDir "/path/to/JointAnalysis/results".|
|--threads|Number of threads. (default: 48).|

### OpenSwath arguments
|parameters|descriptions|
|---|---|
|--openSWATH_paraNumber|Specify the maximum number of parallel data processing for openSWATH (Default: "10").|

**Note:** We process the MS data on a machine with a 64-core CPU and 256G memory. The greater the number of parallel data processing, the higher the memory and CPU resources consumed. If the memory is insufficient, you can appropriately reduce the number of parallel data processing.

### Pyprophet arguments
|parameters|descriptions|
|---|---|
|--classifier|Either a "LDA" or "XGBoost" classifier is used for semi-supervised learning.(default: XGBoost).|
|--lever|Either "ms1", "ms2", "ms1ms2" or "transition". the data level selected for scoring. "ms1ms2 integrates both MS1- andMS2-level scores and can be used instead of "ms2"-level results." (default: "ms1ms2").|

### MRGDiscirm arguments
|parameters|descriptions|
|---|---|
|--seed|Random seed for decoy generation (default: 123).|
|--map_size|The size of the temporary database (default: 32).|
|--fdr_precursor|FDR of precursor level (default: 0.01).|
|--n_mrg|The number of candidate MRGroup (default: 3).|
|--min_nuf|The minimum value of unique features in MRGroup (default: 2).| 
|--nrt_interval_percent|Percentage of the smallest interval in normalized retention time (default: 0.0005).| 
|--nrt_width_percent|Percentage of the search range in normalized retention time (default: 0.02).| 

## Citation
Please cite this.


