# Input parameter: a fasta (not fastq!) file.


export BLASTDB=/net/shendure/vol10/jnburton/extern/blast/ # directory containing "nt" database

PARAMS="-perc_identity 95 -evalue 1e-30 -word_size 50"

echo `date`: blastn
blastn -query $1 -db nt $PARAMS -out $1.blastn -outfmt "7 bitscore sskingdoms sblastnames sscinames salltitles"
echo `date`: done! output is in $1.blastn
