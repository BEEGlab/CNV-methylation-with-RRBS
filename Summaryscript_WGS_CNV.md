# WGS scripts used for analysis related to the manuscript 'Epigenetic diversity of genes with copy number variations among natural populations of the three-spined stickleback'
# by Frédéric J. J. Chain,*, Britta S. Meyer,*, Melanie J. Heckwolf, Sören Franzenburg, Christophe Eizaguirre and Thorsten B.H. Reusch.

# by Britta S. Meyer 
# 26th October 2023

#################################
#1 Preprocessing of Fastq files #
#################################

# 1.0 Fastqc on raw reads
```bash
for file in *.fastq.gz
do
fastqc ${file} -t 16 -o /01FASTQCresults
done
```

# 1.1 FastqtoSAM
for i in *R1_001.fastq.gz; do
    [[ $i =~ (.*)_(.*)_(.*)_R1_001.fastq.gz ]]
    java -Xmx16G -jar picard.jar FastqToSam \
FASTQ= $i \
FASTQ2= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_${BASH_REMATCH[3]}_R2_001.fastq.gz \
OUTPUT= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_unmapped.bam \
READ_GROUP_NAME= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]}.HGJFJBBXX.${BASH_REMATCH[3]} \
SAMPLE_NAME= ${BASH_REMATCH[2]} \
LIBRARY_NAME= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]} \
PLATFORM_UNIT= HGJFJBBXX.${BASH_REMATCH[3]} \
PLATFORM= illumina
done

# 1.2 MarkIlluminaAdapters
for i in *_unmapped.bam; do
  #F12149-L1_S1_unmapped.bam
    [[ $i =~ (.*)_unmapped.bam ]]
    java -Xmx8G -jar -jar ~/software/picard-tools-2.7.1/picard.jar MarkIlluminaAdapters \
  I=$i \
  O=${BASH_REMATCH[1]}_unmapped_xt.bam \
  M=${BASH_REMATCH[1]}_unmapped_xt_metrics.txt

##############################################
#2 Mapping the Fastq files to the reference #
##############################################

# 2.0 index Reference genome
bwa index GasAcuV1.fa

# 2.1  SamToFastq & bwa mem & MergeBamAlignment
java -Xmx16g -jar picard.jar SamToFastq \
INPUT=!file!_unmapped_xt.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 \
INTERLEAVE=true \
| \
bwa mem -M -t 16 -p GasAcuV1.fa /dev/stdin \
| \
java -Xmx16G -jar picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=!file!_unmapped.bam \
OUTPUT=!file!_clean.bam \
CREATE_INDEX=true \
R=GasAcuV1.fa \
MAX_INSERTIONS_OR_DELETIONS=-1 \
ATTRIBUTES_TO_RETAIN=XS \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
CLIP_ADAPTERS=false \
ADD_MATE_CIGAR=true
done

############################################
#3 BAM files - manipulating and statistics #
############################################

# 3.1  MarkDuplicates

java -Xmx100G -jar picard.jar MarkDuplicates \
      I=!file!_clean.bam \
      O=!file!_clean_dedup.bam \
      CREATE_INDEX=true \
      METRICS_FILE=!file!_clean_dedup_metrics.txt \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

# 3.2 statistics
# 3.2.1
echo 'picard'
java -Xmx16G -jar picard.jar CollectMultipleMetrics \
      I=!file!_clean_dedup.bam \
      O=!file!_picard_multiple_metrics \
      R=GasAcuV1.fa \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=QualityScoreDistribution \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectGcBiasMetrics

java -Xmx16G -jar picard.jar CollectRawWgsMetrics \
     I=!file!_clean_dedup.bam \
     O=!file!_picard_raw_wgs_metrics.txt \
     R=GasAcuV1.fa \
     INCLUDE_BQ_HISTOGRAM=true

java -Xmx16G -jar picard.jar CollectWgsMetrics \
     I=!file!_clean_dedup.bam \
     O=!file!_picard_wgs_metrics.txt \
     R=GasAcuV1.fa \
     INCLUDE_BQ_HISTOGRAM=true

java -Xmx16G -jar picard.jar CollectWgsMetricsWithNonZeroCoverage \
       I=!file!_clean_dedup.bam  \
       O=!file!_picard_collect_wgs_metrics.txt \
       CHART=!file!_picard_collect_wgs_metrics.pdf  \
      R=GasAcuV1.fa \
# 3.2.2
echo 'GATK'
 java -jar GenomeAnalysisTK-3.7.jar \
   -T DepthOfCoverage \
   -R GasAcuV1.fa \
   -o !file!_out \
   -I !file!_clean_dedup.bam

####################################################################
#4 Getting the coverage from  picard_collect_wgs_metrics.txt files #
####################################################################   
# import os
import pandas as pd
import numpy as np

# Define a function to extract the values from a single file
def extract_values(file_path):
    values = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("CATEGORY"):
                categories = line.strip().split('\t')
            elif line.startswith("WHOLE_GENOME"):
                data = line.strip().split('\t')
                values = {categories[i]: data[i] for i in range(len(categories))}
    return values

# Specify the directory where your files are located
directory = 'PATH/GEOMAR_WGS/'
data = []

# Loop through the files and extract values
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        file_path = os.path.join(directory, filename)
        values = extract_values(file_path)
        sample_id = filename.replace("picard_collect_wgs_metrics.txt", "").split("_S")[1].rstrip('_')  # Extract e.g. "S48" from the filename
        data.append({'ID': sample_id, 'File_Name': filename, 'MEAN_COVERAGE': float(values['MEAN_COVERAGE']), 'SD_COVERAGE': float(values['SD_COVERAGE']), 'MEDIAN_COVERAGE': float(values['MEDIAN_COVERAGE'])})

# Sort the DataFrame by the numeric ID in ascending order
df = pd.DataFrame(data)
df['ID'] = df['ID'].astype(int)
df = df.sort_values(by='ID', ascending=True)

# Calculate the average and standard deviation of MEAN_COVERAGE
average_mean_coverage = df['MEAN_COVERAGE'].mean()
sd_mean_coverage = df['MEAN_COVERAGE'].std()

# Save the DataFrame to a CSV file
df.to_csv('MEAN_coverage_from_picard_collect_wgs_metrics.csv', index=False)

print(f'Average MEAN_COVERAGE: {average_mean_coverage}')
print(f'Standard Deviation of MEAN_COVERAGE: {sd_mean_coverage}')
