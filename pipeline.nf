#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 定义帮助信息
help = """
Usage:

The basic usage of the workflow:

nextflow run pipeline.nf --rawData <path> --library <path> --tr_irt <path> --swath_windows_file <path> [Options]

Options:
--rawData                Path to raw data
--library                Library file path
--tr_irt                 IRT file path
--swath_windows_file     Swath windows file path
--outputDir              Output directory for results (default: /mnt/data_nas/mlf/openswath-pyprophet/results)
--openSWATH_paraNumber   Number of parallel processes for OpenSWATH (default: 4)
--threads                Number of threads (default: 48)
--batchsize              Batch size for processing
--classifier             Classifier for PyProphet(XGBoost / LDA) (default: XGBoost)
--lever                  Leveraging level: (default: ms1ms2)
--seed                   Random seed for decoy generation (default: 123)
--map_size               The size of the temporary database (default: 32)
--fdr_precursor          FDR of precursor level (default: 0.01)
--n_mrg                  The number of candidate MRGroup (default: 3)
--min_nuf                The minimum value of unique features in MRGroup (default: 2)
--nrt_interval_percent   Percentage of the smallest interval in normalized retention time (default: 0.0005)
--nrt_width_percent      Percentage of the search range in normalized retention time (default: 0.02)

For more information, visit [URL to your documentation or GitHub repository].
"""

// 检测是否使用了 --help 参数
if (params.help) {
    println(help)
    exit 0
}

params.help = " "
params.rawData = " "
params.library = " "
params.tr_irt = " "
params.swath_windows_file = " "
params.outputDir = " "
params.openSWATH_paraNumber = 10
params.threads = 48

// openswath
// params.batchsize = 1000

// pyprophet
params.classifier = "XGBoost"
params.lever = "ms1ms2"

// mrgd
params.seed =  123
params.map_size = 32
params.fdr_precursor = 0.01
params.n_mrg = 3
params.min_nuf = 2
params.nrt_interval_percent = 0.0005
params.nrt_width_percent = 0.02


// 定义通道
Channel.fromPath(params.rawData).set { profile_mzMLs }
Channel.fromPath(params.swath_windows_file).set { windows_file }
Channel.fromPath(params.tr_irt).set { input_irt }
Channel.fromPath(params.library).set { input_lib }


process OpenSWATH {
    publishDir path : "${params.outputDir}/openswath", mode: 'copy'
    container 'ghcr.io/openms/openms-executables'

    input:
    path rawData
    path irt
    path library
    path windows_file

    output:

    path "*.osw"
    path "*.chrom.sqMass"

    script:
    
    """
    find ${params.rawData} -name '*.mzML' | while read i ; do
        j=`echo \${i} | awk -F '/' '{print \$NF;}'`
        echo \${j}
    done | xargs -P ${params.openSWATH_paraNumber} -I{} \\
    OpenSwathWorkflow -in ${params.rawData}/{} \\
    -tr $library -tr_irt $irt -swath_windows_file $windows_file -out_osw {}.osw \\
    -out_chrom {}.chrom.sqMass -threads ${params.threads}

    """
}

process PyProphet {
    publishDir "${params.outputDir}/pyprophet", mode: 'copy'
    container 'pyprophet/pyprophet:2.2.5'

    input:
    path mzML_osw
    path chrom
    path library

    output:
    path "*.osw"
    path "*.tsv"

    script:
    """

    pyprophet merge --template $library --out merged.osw $mzML_osw
    pyprophet score --in merged.osw --classifier ${params.classifier} --level ${params.lever} --threads ${params.threads}
    pyprophet peptide --in merged.osw --context experiment-wide
    pyprophet peptide --in merged.osw --context global
    pyprophet protein --in merged.osw --context experiment-wide
    pyprophet protein --in merged.osw --context global
    pyprophet export --in merged.osw --format legacy_split --max_global_peptide_qvalue 1 --max_rs_peakgroup_qvalue 1 --max_global_protein_qvalue 1 --out results.tsv
    
    """
}

process JointAnalyse {
    publishDir "${params.outputDir}/mrgd", mode: 'copy'
    container 'meiyulab/mrgd:v1'

    input:
    path os_osw
    path os_chrom_sqMass
    path pyp_osw
    path pyp_tsv

    output:

    path "*.tsv"

    script:
    
    """
    source /home/user/miniconda3/bin/activate MRGD
    mrgd --db_fpath $pyp_osw --chrom_dpath ./ --work_dpath . --seed ${params.seed} --map_size ${params.map_size}  \\
    --fdr_precursor ${params.fdr_precursor} --n_mrg ${params.n_mrg} --min_nuf ${params.min_nuf}  \\
    --nrt_interval_percent ${params.nrt_interval_percent} --nrt_width_percent ${params.nrt_width_percent}

    """
}


workflow {
    // 将通道和过程连接起来
    OS_output = OpenSWATH(params.rawData, input_irt, input_lib, windows_file)
    PyP_output = PyProphet(OS_output , input_lib)
    JointAnalyse(OS_output, PyP_output)
}


