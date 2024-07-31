#for file in KCIS02_Pre/*varbin.50k.txt
#do
#outdir=KCIS02_Pre_cnv/
#name=${file/"KCIS02_Pre/"/}
#name=${name/varbin.50k.txt/}
#Rscript cbs.R -i KCIS02_Pre -o $outdir -d $file -s $name > ${name}.cbs.out 2> ${name}.error.log
#done

for file in KCIS05_Pre/*varbin.50k.txt
do
outdir=KCIS05_Pre_cnv/
name=${file##*/}
name=${name/varbin.50k.txt/}
Rscript cbs.R -i KCIS05_Pre -o $outdir -d $file -s $name > log/${name}.cbs.out 2> log/${name}.error.log
done

for file in KCIS06_Pre/*varbin.50k.txt
do
outdir=KCIS06_Pre_cnv/
name=${file##*/}
name=${name/varbin.50k.txt/}
Rscript cbs.R -i KCIS06_Pre -o $outdir -d $file -s $name > log/${name}.cbs.out 2> log/${name}.error.log
done

for file in KCIS07_Pre/*varbin.50k.txt
do
outdir=KCIS07_Pre_cnv/
name=${file##*/}
name=${name/varbin.50k.txt/}
Rscript cbs.R -i KCIS07_Pre -o $outdir -d $file -s $name > log/${name}.cbs.out 2> log/${name}.error.log
done

for file in KCIS08_Pre/*varbin.50k.txt
do
outdir=KCIS08_Pre_cnv/
name=${file##*/}
name=${name/.varbin.50k.txt/}
Rscript cbs.R -i KCIS08_Pre -o $outdir -d $file -s $name > log/${name}.cbs.out 2> log/${name}.error.log
done

for file in KCIS09_Pre/*varbin.50k.txt
do
outdir=KCIS09_Pre_cnv/
name=${file##*/}
name=${name/.varbin.50k.txt/}
Rscript cbs.R -i KCIS09_Pre -o $outdir -d $file -s $name > log/${name}.cbs.out 2> log/${name}.error.log
done


