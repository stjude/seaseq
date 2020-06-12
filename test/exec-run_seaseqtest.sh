bsub -P watcher -q compbio -J seaseqexec -o seaseqexec-out -e seaseqexec-err -N ../exec-seaseq_wrap.sh inputyml.yml EXECSEASEQ
