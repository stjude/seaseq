bsub -P watcher -q compbio -J seaseq -o seaseq-out -e seaseq-err -N ../seaseq_wrap.sh inputyml.yml SEASEQ

