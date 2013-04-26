#!/usr/bin/perl
# programmer: zocean
# usage:
# input:
# output:

use LWP::Simple;

my $esearch = "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" . "db=nucleotide&usehistory=y&rettype=gb&retmode=text&seq_start=1&seq_stop=100&id=";
open(IN,$ARGV[0]) or die();
#my $id= "ref|NM_001146472.1|";
while(<IN>){
chomp();
my $id=$_;
my @p = split(/\|/,$id);
my $q=$p[1];
my $sign=0;
my $esearch_result = get($esearch . $q);
my @esearch_table=split("\n",$esearch_result);
foreach $item(@esearch_table){
    $sign=0;
    if($item=~/ORGANISM/){
        $item=~s/\s+ORGANISM\s+//;
        print "Id:$id\tOrganism:$item\n";
        $sign=1;}
}
if ($sign=0){print "Id:$id\tOrganism:NOT FOUND\n";}
}
