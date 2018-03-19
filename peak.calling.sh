# input: sinlge base level coverage 
# output: every window with its read depth, number of C or G and p values
module load samtools
module load bedtools
#!/bin/bash 
#20160811_5mC_Jump_Seq_2.4ng.mm9.umi_encoded_adaptor_removed_no_mismatch.sorted.dedup.bam.plus.cove.base.level.shiftonebase.gz
path1=/project2/xinhe/Shengtong/Shengtong_JumpSeq/5mC_2.4ng
path2=/project2/xinhe/Shengtong/Shengtong_JumpSeq/mouse.genome
cd $path1/
prefix1=20160811_5mC_Jump_Seq_2.4ng
prefix2=mm9.umi_encoded_adaptor_removed_no_mismatch.sorted.dedup.bam
################################ 
winlength=21
strand=plus
################# get windows of equal size  and its reads depth  chromosome by chromosome ##########
gzip -d $prefix1.$prefix2.$strand.cove.base.level.shiftonebase.gz
rm $prefix1.$strand.$winlength.bp.win.5prime.cove
for i in `cut -f1 $path2/mm9.genome`
do   rm $prefix1.$strand.chr # clear unnecessary files if any, that's a very important step!
     rm $prefix1.$strand.reads.in.window
     rm $prefix1.$strand.endpt
     rm $prefix1.$strand.com
     grep  "\<$i\>" $prefix1.$prefix2.$strand.cove.base.level.shiftonebase  > $prefix1.$strand.chr
     numbase=$(wc -l < $prefix1.chr)
     totalcove=$( awk '{sum+=$3} END {print sum}' $prefix1.$strand.chr)
     awk '{sum+=$3} NR% 21 ==0  {print sum; sum=0}' $prefix1.$strand.chr  > $prefix1.$strand.reads.in.window
     subcove=$( awk '{sum+=$1} END {print sum}'  $prefix1.$strand.reads.in.window)
     awk ' NR% 21 ==0  {print $1, "\t", $2}' $prefix1.$strand.chr  > $prefix1.$strand.endpt
     lastbase=$(tail -1 $prefix1.$strand.endpt | awk ' {print $2}')
     paste -d" "  $prefix1.$strand.endpt    $prefix1.$strand.reads.in.window > $prefix1.$strand.com
     echo numbase is $numbase 
     echo lastbase is $lastbase
     if (($numbase<=$lastbase)); then  # the length of chromosome is 21*N
     awk '{print $1,"\t", $2-20, "\t", $2, "\t", $3}' $prefix1.$strand.com >> $prefix1.$strand.$winlength.bp.win.5prime.cove
     else                                 # the length of chromosome is NOT 21*N, there is residual
      awk '{print $1,"\t", $2-20, "\t", $2, "\t", $3}' $prefix1.$strand.com >> $prefix1.$strand.$winlength.bp.win.5prime.cove
      covlastwin=$(($totalcove-$subcove))
      echo $i $(($lastbase+1)) $numbase $covlastwin >> $prefix1.$strand.$winlength.bp.win.5prime.cove # append the last window < 21bp
     fi
done 
rm $prefix1.$strand.chr # clear unnecessary files if any, that's a very important step!
rm $prefix1.$strand.reads.in.window
rm $prefix1.$strand.endpt
rm $prefix1.$strand.com
gzip $prefix1.$prefix2.$strand.cove.base.level.shiftonebase
######################### only use effective windows with at least one reads #########################
awk ' $4>0 {print $1, "\t", $2, "\t", $3, "\t", $4}' $prefix1.$strand.$winlength.bp.win.5prime.cove |awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' > $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed
rm $prefix1.$strand.$winlength.bp.win.5prime.cove
echo get windows with positive number of reads 
################## get sequence of effective windows ############################
bedtools getfasta -fi $path2/mm9.fasta -bed $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed -fo  $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.seq
echo get sequence of windows with positive reads
################### count C's (G's) in each effective window for minus (plus) strand ##################### 
awk 'NR%2==0'  $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.seq  | sed 's/[^(c,C)]//g' | awk '{ print length }' >  $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC 
echo count the number of Cs in every window
paste -d" " $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed   $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC |awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' > $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC.combine 
echo for each window, combine its coverage and Cs 
############ output reads depth and number of C (G) in every window to compute p value ##############
awk '{print $4, "\t", $5 }' $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC.combine > $prefix1.$strand.numberC.reads.txt
echo output number of reads and Cs to compute p values 
rm $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed
rm $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.seq
rm $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC
######################### calculate p value ##############
module load R 
Rscript $path1/$prefix1.$strand.pvalue.R
echo calculate p values 
######################## get effective windows and it's p value ###########
dos2unix $prefix1.$strand.pval.txt
paste -d" " $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC.combine  $prefix1.$strand.pval.txt | awk '$5>0 {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5, "\t", $6}'  >  $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC.combine.pvalue 
echo filter out windows without C and combine p value 
rm $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC.combine
rm $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.pvalue.txt
gzip $prefix1.$strand.$winlength.bp.win.5prime.pstv.cove.bed.numberC.combine.pvalue
echo done


