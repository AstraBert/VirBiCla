for i in 16S_ribosomal_RNA 18S_fungal_sequences 28S_fungal_sequences SSU_eukaryote_rRNA LSU_prokaryote_rRNA ref_viruses_rep_genomes
do
    echo "Started with $i"
    bash ../shell/retrieve_blastdb.sh \
        -db $i \
        -o ${i}.fsa
    size=$(du -m ${i}.fsa)
    pigz -p 5 -d ${i}.fsa.gz
    seqs=$(grep -o ">" ${i}.fsa | wc -l)
    echo $seqs
    pigz -p 5 ${i}.fsa
    echo "Finished with $i: it is $size MB in size and has $seqs sequences"
done