#!/bin/bash
#PBS -o sailfish_log.out
#PBS -j oe
#PBS -q shared
#PBS -l advres=genomics
#PBS -N sailfF_run
#PBS -l nodes=1:ppn=12
#PBS -l walltime=5:00:00

cd $SCRATCH

echo 'Job_Started: '`date`

# module set up
module purge
module load sailfish/0.10.0

# ---------------- version ---------------------#
#echo "	picard version:" $PICARD_TOOLS_HOME
echo "	current node:" `hostname`
echo "	current location:" $SCRATCH
#echo "	BWA version:" `which bwa`
#echo "	samtools version:" `which samtools`
#echo "	java version:" `which java`
echo "	SailFish_Version:" `which sailfish` 
# --------------------------------------------- #

REF=/ref_transcript/Mus_musculus.GRCm38.cdna.all.ensembl.fa
OUTDIR=/sailfish

# -------------------------- Indexing ------------------------------ #
# Command Line Options:
#  -v [ --version ]           print version string
#  -h [ --help ]              produce help message
#  -t [ --transcripts ] arg   Transcript fasta file(s).
#  -m [ --tgmap ] arg         file that maps transcripts to genes
#  -k [ --kmerSize ] arg      Kmer size.
#  -o [ --out ] arg           Output stem [all files needed by Sailfish will be
#                             of the form stem.*].
#  -p [ --threads ] arg (=48) The number of threads to use concurrently.
#  -f [ --force ]
#  	sailfish index -t <ref_transcripts> -o <out_dir> -k <kmer_len> -p numThreads
# ------------------------------------------------------------------- #
mkdir ${OUTDIR}/index_dir
sailfish index -t $REF \
               -o ${OUTDIR}/index_dir \
			   -k 19 \
               -p 12

echo 'Indexing_Done: '`date` 

# ----------------------- Quantification ------------------------------ #
#   -i [ --index ] arg                    Sailfish index [output of the "Sailfish index" command
#   -l [ --libtype ] arg                  Format string describing the library type
#   -r [ --unmated_reads ] arg            List of files containing unmated reads of (e.g. single-end reads)
#   -1 [ --mates1 ] arg   -2 [ --mates2 ] arg  File containing the # mates #2 mates
#   --no_bias_correct                     turn off bias correction
#   -m [ --min_abundance ] arg (=0)       transcripts with an abundance (KPKM) lower than this value will be reported at zero.
#   -o [ --out ] arg                      Basename of file where estimates are written
#   -n [ --iterations ] arg (=1000)       number of iterations to run the optimzation
#   -d [ --delta ] arg (=0.0050000000000000001) consider the optimization to have converged if the relative change in the estimated abundance of all
#                                         transcripts is below this threshold
#   -p [ --threads ] arg (=48)            The number of threads to use when counting kmers
#   -f [ --force ]                        Force the counting phase to rerun, even if a count databse exists.
#   -a [ --polya ]                        polyA/polyT k-mers should be discarded

#  <libtype> "T=SE:S=A"   Library format { type:single end, relative orientation:none, strandedness:antisense }
# ---------------------------------------------------------------------------- #

for FASTQ_INPUT_FILE in $(ls /data/Sample*);do
	sampleN=$(basename $FASTQ_INPUT_FILE | cut -d. -f 1)
 
	mkdir ${OUTDIR}/${sampleN}

	sailfish quant -i ${OUTDIR}/index_dir \
				   -l SR \
				   -r $FASTQ_INPUT_FILE \
				   -o ${OUTDIR}/${sampleN} \
				   -p 12

done

echo 'Job_Done: '`date` 

exit