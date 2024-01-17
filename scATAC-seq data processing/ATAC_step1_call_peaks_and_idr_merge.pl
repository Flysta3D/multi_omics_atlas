#!/usr/bin/perl
###This file is used for apply peak calling and merging steps of ATAC batches.
###The data dir structure is like: /workdir/E0-2h(timeslot)/ATAC***(batch_id)/ATAC***(batch_id).fragments.tsv.gz
###The result will be stored in /workdir/E0-2h(timeslot)/merge-peak-idr/total.merged.narrowpeak
$workdir="/workdir";
opendir(WORKDIR,"$workdir");
while($timewindow = readdir WORKDIR){
    if($timewindow !=~ /E*h/){next;}
    opendir(DIR1,"$workdir/$timewindow");
    $i=0;
    $exist=0;
    ###call peak from fragment files by MACS2 and sort each batch MACS2 results
    while($file= readdir DIR1){
        if($file =~ /^ATAC/){
            $batch[$i]=$file;
            $i++;
            system("macs2 callpeak -t $workdir/$timewindow/$file/*fragments.tsv.gz -n $file --shift -75 --extsize 150 --nomodel -B --SPMR -g dm --outdir $workdir/$timewindow/$file");
            system("sort -k8,8nr $workdir/$timewindow/$file/$file\_peaks.narrowPeak > $workdir/$timewindow/$file/$file\_peaks.sorted.narrowPeak");
            ###You can remove files produced by MACS2 that is not necessary for further steps at here. Only *_peaks.sorted.narrowPeak is needed. 
            #system("rm $workdir/$timewindow/$file\_treat_pileup.bdg");
            #system("rm $workdir/$timewindow/$file\_summits.bed");
            #system("rm $workdir/$timewindow/$file\_peaks.xls");
            #system("rm $workdir/$timewindow/$file\_peaks.narrowPeak");
            #system("rm $workdir/$timewindow/$file\_control_lambda.bdg");
        }
        if($file eq "merge-peaks-idr"){
           $exist = 1;
        }
    }
    ###make dir if merge-peak-idr dir is not exist
    if($exist == 0){
        system("mkdir $workdir/$timewindow/merge-peaks-idr");
    }
    ###apply idr for each batch pairs in one timewindow
    for($i=0;$i<scalar(@batch);$i++){
        for($j=$i+1;$j<scalar(@batch);$j++){
            system("idr --samples $workdir/$timewindow/$batch[$i]/$batch[$i]\_peaks.sorted.narrowPeak $workdir/$timewindow/$batch[$j]/$batch[$j]\_peaks.sorted.narrowPeak --input-file-type narrowPeak --rank p.value --output-file $workdir/$timewindow/merge-peaks-idr/ $batch[$i]\_$batch[$j].merged.narrowPeak");
        }
    }
    ###sort idr files and merge together by bedtools
    opendir(DIR2,"$workdir/$timewindow/merge-peaks-idr");
    $i=0;
    while($file= readdir DIR2){
        if($file =~ /^ATAC/ && $file =~ /merged\.narrowPeak$/){
            $narrowPeaks[$i]=$file;
        }
    }
    foreach $narrowfile (@narrowPeaks){
        system("cat $workdir/$timewindow/merge-peaks-idr/$narrowfile >> $workdir/$timewindow/merge-peaks-idr/total.unmerged.narrowPeak");
        system("rm $workdir/$timewindow/merge-peaks-idr/$narrowfile");
    }
    system("sort -k1,1 -k2,2n $workdir/$timewindow/merge-peaks-idr/total.unmerged.narrowPeak > $workdir/$timewindow/merge-peaks-idr/total.unmerged.sorted.narrowPeak");
    system("rm $workdir/$timewindow/merge-peaks-idr/total.unmerged.narrowPeak");
    system("bedtools merge -i $workdir/$timewindow/merge-peaks-idr/total.unmerged.sorted.narrowPeak > $workdir/$timewindow/merge-peaks-idr/total.merged.narrowPeak");
    system("rm $workdir/$timewindow/merge-peaks-idr/total.unmerged.narrowPeak");

}
