###QC fastq files (paired-end seq reads) using bbtools

FEXT=fastq #Set .fastq extension for any output file

# IN1=${SRRID}_Red_Sea_metagenome_sample_${SAMPLEID}_1.fastq.gz
# IN2=${SRRID}_Red_Sea_metagenome_sample_${SAMPLEID}_2.fastq.gz

IN1=$1 #First fastq file
IN2=$2 #Second fastq file
SRRID=`echo $IN1 | awk -F "_" '{print $1}'` #Fastq file ID common to both, e.g. for a file named ERR3980584_1.fastq.gz this yields: ERR3980584

OUT1=${SRRID}_merged_reads.${FEXT}
OUTU1=${SRRID}_unmerged1_reads.${FEXT}
OUTU2=${SRRID}_unmerged2_reads.${FEXT}

#Merge paired-end reads
bbmerge-auto.sh in1=${IN1} in2=${IN2} out=${OUT1} outu1=${OUTU1} outu2=${OUTU2}
echo "Merging paired-end seq reads is done!"
echo ""

#Adapter-trimming
TRIM_IN=${OUT1}
TRIM_OUT=${SRRID}_Adap_Trim.${FEXT}
bbduk.sh in=${TRIM_IN} \
            out=${TRIM_OUT} \
            ref=$HOME/bbmap/resources/adapters.fa \
            ktrim=r k=27 hdist=1 tpe tbo -eoom
echo "Adapter trimming is done!"
echo ""


#Quality trimming: quality-trim to Q20 using the Phred algorithm
QTRIM_IN=${TRIM_OUT}
QTRIM_OUT=${SRRID}_Qual_Trim.${FEXT}
bbduk.sh in1=${QTRIM_IN} \
            out1=${QTRIM_OUT} \
            qtrim=r trimq=20 -eoom
echo "Quality trimming is done!"
echo ""


#Quality filtering: This will discard reads with average quality below 20.
QFILT_IN=${QTRIM_OUT}
QFILT_OUT=${SRRID}_Qual_Filt.${FEXT}
bbduk.sh in1=${QFILT_IN} \
            out1=${QFILT_OUT} \
            maq=20 -eoom
echo "Quality filtering is done!"
echo ""

#Delete raw reads and any intermmediate files: keep only the QC samples
rm ${IN1} ${IN2} ${OUTU1} ${OUTU2} ${OUT1} ${TRIM_OUT} ${QTRIM_OUT}
