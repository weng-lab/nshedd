for i in $(seq 5 2 25); do
    kmc -ci0 -cs1000000000 -cx1000000000000000000000 -k${i} ../fastqs/D3-B_cells-HC_R1_001.fastq.gz ${i}mers .
    kmc_tools transform ${i}mers dump ${i}mers.txt
done
