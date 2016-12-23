#!/usr/bin/perl -w
use strict;
use IO::File;
use miRnaQTL;

my $qtl_files="QTL_UHT.txt";

my $foo=SNPtools::miRnaQTL->new();

$foo->add_SNP(-SNP_file=>"miRNA-UCSC.txt");

$foo->add_annotation(-Annotation_file=>"pre-miRNA-gff.txt");

$foo->add_QTL(-QTL_files=>"$qtl_files");

$foo->add_outfile(-out_file=>"susu.txt");

$foo->find_miRnaQTL();
