#CHECK DOCUMENTATION FOR FULL USE: https://velocyto.org/velocyto.py/tutorial/index.html

#RUN SMARTSEQ
velocyto run_smartseq2 -o OUTPUT -m repeat_msk.gtf -e MyTissue plateX/*.bam mm10_annotation.gtf

#RUN CHROMIUM
velocyto run10x -m repeat_msk.gtf mypath/sample01 somepath/refdata-cellranger-mm10-1.2.0/genes/genes.gtf

#RUN DROPSEQ
./droptag -c ./configs/indrop_v1_2.xml ~/mydata/SRR5945694_2.fastq.gz ~/mydata/SRR5945694_1.fastq.gz
dropest -m -V -b \
          -o ~/mydata/SRR5945695/SRR5945695_dropEst \
          -g ~/cellranger/refdata-cellranger-mm10-1.2.0/genes/genes.gtf \
          -L eiEIBA \
          -c ~/mysource/dropEst/configs/config_desc.xml \
          ~/mydata/SRR5945695_1/SRR5945695_1Aligned.sortedByCoord.out.bam

velocyto run_dropest -o ~/mydata/SRR5945695_results -m rep_mask.gtf ~/mydata/SRR5945695_1/correct_SRR5945695_1Aligned.sortedByCoord.out.tagged.bam mm10_annotation.gtf
