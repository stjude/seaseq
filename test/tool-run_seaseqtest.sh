bsub -P watcher -q compbio -J seaseqtool -o seaseqtool-out -e seaseqtool-err -N ../tool-seaseq_wrap.sh inputyml.yml TOOLSEASEQ
