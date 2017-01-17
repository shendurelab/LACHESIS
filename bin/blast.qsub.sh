for CHUNK in {1..190}; do
    qsub blast.sh $CHUNK
done

# run ChunkFasta before running this
# then run this
# then, after blast is done, run the following command:
#for i in {1..190}; do cat assembly.$i.blast.out >> assembly.blast.out; done