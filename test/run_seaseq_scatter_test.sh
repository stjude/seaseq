bsub -P watcher -q compbio -J seaseqscat -o seaseqscat-out -e seaseqscat-err -N ../seaseq_wrap_scatter.sh inputymlscatter.yml SEASEQSCATTER
