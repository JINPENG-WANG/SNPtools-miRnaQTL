package SNPtools::miRnaQTL;

use 5.006;
use strict;
use warnings;

=head1 NAME

SNPtools::miRnaQTL - The great new SNPtools::miRnaQTL!

=head1 VERSION

Version 1.0

=cut

our $VERSION = '1.0';


=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use SNPtools::miRnaQTL;

    my $foo = SNPtools::miRnaQTL->new();
	
	$foo->add_SNP(-SNP_file=>"");

	$foo->add_annotation(-Annotation_file=>"");

	$foo->add_QTL(-QTL_files=>"");

	$foo->add_outfile(-out_file=>"");

	$foo->find_miRnaQTL();


=cut

use IO::File;
use vars qw (%fields);

%fields = (
	snp => undef,
	annotation => undef,
	QTL => undef,
	out_file => undef,
);


=head1 SUBROUTINES/METHODS

=head2 new

	$foo = SNPtools::miRnaQTL->new();

=cut

sub new {
	my $class = shift;
	my $self = {
		%fields,
	};
	
	bless ($self, $class);
	return $self;
}

=head2 add_SNP
	
	$foo = SNPtools::miRnaQTL->add_SNP(-SNP_file=>"");

=cut

sub add_SNP {
	my $self = shift;
	my %parameters = @_;
	for (keys %parameters) {
		if($_=~/^-SNP_file$/){
			$self->{snp}=$parameters{$_};
		}
		else{
			die "Unacceptable parameters, please read README!\n";
		}
	}
}




=head2 add_annotation

		$foo->add_annotation(-Annotation_file=>"");
		
=cut

sub add_annotation {
	my $self=shift;
	my %parameters= @_;
	for (keys %parameters) {
		if($_=~/^-Annotation_file$/){
			$self->{annotation}=$parameters{$_};
		}
		else{
			die "Unacceptable parameters, please read README!\n";
		}
	}
}


=head2 add_QTL

		$foo->add_QTL(-QTL_files=>"");

=cut


sub add_QTL {
	my $self=shift;
	my %parameters= @_;
	for (keys %parameters) {
		if($_=~/^-QTL_files$/){
			$self->{QTL}=$parameters{$_};
		}
		else{
			die "Unacceptable parameters, please read README!\n";
		}
	}

}

=head2 add_outfile

	$find -> add_outfile(-out_file=>"");

			Output the result into the out file.
=cut

sub add_outfile {
	my $self = shift;
	my %parameters= @_;
	for (keys %parameters) {
		if($_=~/^-out_file$/){
			$self->{out_file}=$parameters{$_};
		}
		else{
			die "Unacceptable parameters, please read README!\n";
		}
	}
}


=head2 find_miRnaQTL

		$foo->find_miRnaQTL();

=cut


sub find_miRnaQTL {
	my $self = shift;
	my $snp_file=$self->{snp};
	my $annotation_file=$self->{annotation};
	my $QTL_files=$self->{QTL};
	my $out_file=$self->{out_file};


	if($snp_file){
		print "The SNP file user provides is:\t", $snp_file,"\n";
	}else{
		die "No SNP file provided!\n";
	}
	if($annotation_file){
		print "The annotation file user provides is:\t", $annotation_file,"\n";
	}else{
		die "No annotation file provided!\n";
	}
	if($QTL_files){
		print "The QTL files user provides are:\t", $QTL_files,"\n";
	}else{
		die "No QTL files provided!\n";
	}
	if($out_file){
		print "The output file user provides is:\t", $out_file,"\n";
	}else{
		die "No output file provided!\n";
	}


	my $snp_fh=IO::File->new("$snp_file",'r');
	my $annotation_fh=IO::File->new("$annotation_file",'r');
	my $out_fh=IO::File->new(">$out_file");

	my @QTL_files=split /\s+/, $QTL_files;

	my %qtl;
	my @lables;
	my $lable_head="";
	for my $qtl_file(@QTL_files){
		my $qtl_lable;
		if($qtl_file=~/(.+\/)?QTL_(\w+)\.txt$/){
			$qtl_lable=$2;
			push @lables, $qtl_lable;
			$lable_head.="\tSymbol_$qtl_lable\tQTL_ID_$qtl_lable";
		}
		my $qtl_fh=IO::File->new("$qtl_file",'r');

		while(<$qtl_fh>){
			chomp;
			my $line=$_;
			my @eles=split /\t/,$line;
			my ($QTL_ID,$chr,$qtl_start,$qtl_end)=@eles[0,2,5,6];
			$qtl{$qtl_lable}{$chr}{$QTL_ID}{start}=$qtl_start;
			$qtl{$qtl_lable}{$chr}{$QTL_ID}{end}=$qtl_end;
		}
	}

	$out_fh->print("SNP_ID\trefNCBI\trefUCSC\tAllele\tChromosome\tchromStart\tchromEnd\tPre_miRNA_ID\tPre-miRNA\tStrand\tStart\tEnd\tIn pre-miRNA\tMiRNA\tIn mature sequence\tIn seed region\tMiRNA\tIn mature sequence\tIn seed region");
	$out_fh->print("$lable_head\n");
	my %gff1; my %gff2; my (%miR,%mir);
	my $mir_ID;


	while(<$annotation_fh>){
		chomp;
		my $line=$_;
		my @eles=split /\t/, $line;
		my ($chr,$start,$end,$strand,$ID,$name)=@eles[0,3,4,5,6,8];
		if($chr=~/chr(\w+)/){
			$chr=$1;
		}
		if($ID=~/ID=(.+)/){
			$ID=$1;
		}
		if($line=~/Derives_from=(.+)/){
			my $miR_ID=$ID;
			$gff2{$chr}{$mir_ID}{$miR_ID}{start}=$start;
			$gff2{$chr}{$mir_ID}{$miR_ID}{end}=$end;
			if($name=~/Name=(bta.+)/){
				my $miR_name=$1;
				$miR{$miR_ID}=$miR_name;
			}
		}else{
			if($name=~/Name=(bta.+)/){
				my $mir_name=$1;
				$mir_ID=$ID;
				$gff1{$chr}{$mir_ID}{strand}=$strand;
				$gff1{$chr}{$mir_ID}{start}=$start;
				$gff1{$chr}{$mir_ID}{end}=$end;
				$mir{$mir_ID}=$mir_name;
			}
		}
	}

	while(<$snp_fh>){
		chomp;
		my $line=$_;
		if($line=~/^rs/){
			my @eles=split /\t/, $line;
			my ($id,$chr,$snp_start,$snp_end,$strand,$ref_N,$ref_U,$allele)=@eles[0,1,2,3,4,5,6,7];
			if($chr=~/chr(\w+)/){
				$chr=$1;
			}
			my @QTL_info;
			

			for my $qtl_lable (@lables){
				my $qtl_id="-";
				if($qtl{$qtl_lable}{$chr}){
					for my $QTL_ID (keys $qtl{$qtl_lable}{$chr}){
						my $qtl_start=$qtl{$qtl_lable}{$chr}{$QTL_ID}{start};
						my $qtl_end=$qtl{$qtl_lable}{$chr}{$QTL_ID}{end};
						 if( $snp_start>=$qtl_start && $snp_start <=$qtl_end || $snp_end>=$qtl_start && $snp_end<=$qtl_end){
							$qtl_id=$QTL_ID;
							last;
						}
					}
				}
				push @QTL_info, $qtl_lable;
				push @QTL_info, $qtl_id;
			}
			my $qtls_info=join("\t", @QTL_info);

			#print "@QTL_info\n";
			my $flag=0;
			if($gff1{$chr}){
				for my $mir_ID (keys %{$gff1{$chr}}){
					my $pre_start=$gff1{$chr}{$mir_ID}{start};
					my $pre_end=$gff1{$chr}{$mir_ID}{end};
					my $pre_strand=$gff1{$chr}{$mir_ID}{strand};
					my $mir_name=$mir{$mir_ID};
					if($snp_start>=$pre_start && $snp_start<=$pre_end|| $snp_end>=$pre_start && $snp_end<=$pre_end){
						$flag++;
						my $in_pre_miRNA=0;
						if($pre_strand eq "+"){
							$in_pre_miRNA=$snp_end-$pre_start+1;
							my $info="$id\t$ref_N\t$ref_U\t$allele\t$chr\t$snp_start\t$snp_end\t$mir_ID\t$mir_name\t$pre_strand\t$pre_start\t$pre_end\t$in_pre_miRNA";
							for my $miR_ID(keys %{$gff2{$chr}{$mir_ID}}){
								my $mature_start=$gff2{$chr}{$mir_ID}{$miR_ID}{start};
								my $mature_end=$gff2{$chr}{$mir_ID}{$miR_ID}{end};
								my $miR_name=$miR{$miR_ID};
								my $in_mature="-";
								my $inseed="-";
								if($snp_start>=$mature_start && $snp_start<=$mature_end){
									$in_mature=$snp_start-$mature_start;
									#$info=$info."$miR_name\t$in_mature\t";
									if($in_mature<=8 && $in_mature>=1){
										$inseed="inseed";
										#$info=$info."$inseed\t";
									}
								}
								$info=$info."\t$miR_name\t$in_mature\t$inseed";
							}
							if(keys %{$gff2{$chr}{$mir_ID}}<2){
								$info=$info."\t-\t-\t-";
							}
							$out_fh->print("$info\t$qtls_info\n");
						}
						if($pre_strand eq "-"){
							$in_pre_miRNA=$pre_end-$snp_start;
							my $info="$id\t$ref_N\t$ref_U\t$allele\t$chr\t$snp_start\t$snp_end\t$mir_ID\t$mir_name\t$pre_strand\t$pre_start\t$pre_end\t$in_pre_miRNA";
							for my $miR_ID(keys %{$gff2{$chr}{$mir_ID}}) {
								my $mature_start=$gff2{$chr}{$mir_ID}{$miR_ID}{start};
								my $mature_end=$gff2{$chr}{$mir_ID}{$miR_ID}{end};
								my $miR_name=$miR{$miR_ID};
								my $in_mature="-";
								my $inseed="-";
								if($snp_start>=$mature_start && $snp_start<=$mature_end){
									$in_mature=$mature_end-$snp_start;
									#$info=$info."$miR_name\t$in_mature\t";
									if($in_mature<=8 && $in_mature>=1){
										$inseed="inseed";
										#$info=$info."$inseed\t";
									}
								}
								$info=$info."\t$miR_name\t$in_mature\t$inseed";
							}
							if(keys %{$gff2{$chr}{$mir_ID}}<2){
								$info=$info."\t-\t-\t-";
							}
							$out_fh->print("$info\t$qtls_info\n");
						}					
					}
				}
			}
			if($flag == 0){
				print("$line\n");
			}
		}
	}

}




=head1 AUTHOR

Jinming Huang, Jinpeng Wang, C<< <wangjinpeng0225 at 163.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-snptools-mirnaqtl at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=SNPtools-miRnaQTL>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc SNPtools::miRnaQTL


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=SNPtools-miRnaQTL>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/SNPtools-miRnaQTL>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/SNPtools-miRnaQTL>

=item * Search CPAN

L<http://search.cpan.org/dist/SNPtools-miRnaQTL/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2016 Jinming Huang, Jinpeng Wang.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

1; # End of SNPtools::miRnaQTL
