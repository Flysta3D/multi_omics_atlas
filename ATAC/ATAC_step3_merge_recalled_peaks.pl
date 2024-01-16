#!/usr/bin/perl
###This file is used for peak merging of different ATAC time slots for further whole embedding.
###The data dir structure is like: /workdir/E0-2h(timeslot)/kept_peaks.bed
###The result will be stored in /workdir/total.new.merged.bed
$workdir="/workdir";
opendir(WORKDIR,"$workdir");
while($timewindow = readdir WORKDIR){
    if($timewindow !=~ /E*h/){next;}
    opendir(DIR1,"$workdir/$timewindow/analysis/03.peak_by_clusters_and_reclustering/");
    while($file= readdir DIR1){
        if($file == "kept_peaks.bed"){
            system("cat $workdir/$timewindow/analysis/03.peak_by_clusters_and_reclustering/kept_peaks.bed >> $workdir/total.new.unmerged.bed");
        }
    }
    system("sort -k1,1 -k2,2n $workdir/total.new.unmerged.bed > $workdir/total.new.unmerged.sorted.bed");
    system("rm  $workdir/total.new.unmerged.bed");
    system("bedtools merge -i $workdir/total.new.unmerged.sorted.bed > $workdir/total.new.merged.bed");
    system("rm $workdir/total.new.unmerged.sorted.bed");

}
