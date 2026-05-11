#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;

my $usage = "pal2nal.pl protein.aln cds.aln [options] > output\n";
my $outformat = "fasta";
my $nogap = 0;
my $nomismatch = 0;
my $nocdsstop = 0;
my $verbose = 0;
GetOptions(
    "output=s" => \$outformat,
    "nogap" => \$nogap,
    "nomismatch" => \$nomismatch,
    "nocdsstop" => \$nocdsstop,
    "verbose" => \$verbose,
) or die $usage;
die $usage unless @ARGV >= 2;
my ($prot_file, $cds_file) = @ARGV;

# 读取蛋白质比对
my %prot_seqs;
my $id = "";
open(PROT, $prot_file) or die "Cannot open $prot_file";
while(<PROT>) {
    chomp;
    if(/^>(\S+)/) { $id = $1; $prot_seqs{$id} = ""; }
    elsif($id) { $prot_seqs{$id} .= $_; }
}
close(PROT);

# 读取原始 CDS 序列
my %cds_seqs;
$id = "";
open(CDS, $cds_file) or die "Cannot open $cds_file";
while(<CDS>) {
    chomp;
    if(/^>(\S+)/) { $id = $1; $cds_seqs{$id} = ""; }
    elsif($id) { $cds_seqs{$id} .= $_; }
}
close(CDS);

# 检查 ID 匹配
my @ids = keys %prot_seqs;
for my $id (@ids) {
    die "No CDS for $id" unless exists $cds_seqs{$id};
}

# 生成密码子比对
for my $id (@ids) {
    my $prot = $prot_seqs{$id};
    my $cds = $cds_seqs{$id};
    my $codon_seq = "";
    my $prot_len = length($prot);
    my $cds_len = length($cds);
    die "$id: CDS length not multiple of 3" if $cds_len % 3 != 0;
    my $codon_num = $cds_len / 3;
    for(my $i = 0; $i < $prot_len; $i++) {
        my $aa = substr($prot, $i, 1);
        if($aa eq '-') {
            $codon_seq .= "---";
        } else {
            if($i >= $codon_num) { die "$id: protein longer than CDS"; }
            my $codon = substr($cds, $i*3, 3);
            $codon_seq .= $codon;
        }
    }
    if($outformat eq "fasta") {
        print ">$id\n";
        for(my $j=0; $j<length($codon_seq); $j+=60) {
            print substr($codon_seq, $j, 60), "\n";
        }
    } else {
        # 其他格式略
    }
}
