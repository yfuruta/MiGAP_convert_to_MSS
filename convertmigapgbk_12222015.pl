#!/usr/bin/perl
#convertmigapgbk.pl
#Dec 18 2015
#Convert MiGAP genbank file to NCBI genbank file by removing and/or adding information.
######input format######
#file_directory<tab>strain_name
########################
#

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;

if (@ARGV != 1){
	die "usage: $0 filelist";
}

my $filelist = $ARGV[0];

open (FILE, $filelist);


while (<FILE>){
	/(.*)\t(.*)/;
	my $file = $1;
	my $strain = $2;
	my $seqio_obj = Bio::SeqIO->new(-file => $file,
									-format => "genbank");
	my $trnacount = 1;
	my $rrnacount = 1;

	open (OUT, ">./converted_gbk/${strain}_12242015.gbk");
	open (OUT2, ">./converted_gbk/${strain}_DDBJformat_12242015.txt");

	my @gapstart;
	my @gapend;
	my $lastprintedgapnum = -1;

	while (my $seq_obj = $seqio_obj->next_seq){

		for my $feat_obj ($seq_obj->get_SeqFeatures){

			if ($feat_obj->primary_tag eq "source"){
				##Extract position information
				my $start = $feat_obj->location->start;
				my $end = $feat_obj->location->end;
				my $strand = $feat_obj->location->strand;
				my $position;
				if ($strand == 1){
					$position = $start."\.\.".$end;
				}elsif ($strand == -1){
					$position = "complement\(".$start."\.\.".$end."\)";
				}

				##Printing description for genbank
				print OUT "LOCUS       ",$strain," ",$end," bp    DNA circular\n";
				print OUT "DEFINITION  ",$strain,"\n";
				print OUT "FEATURES             Location/Qualifiers\n";
				print OUT "     source          ",$position,"\n";
				print OUT "                     \/organism\=\"Helicobacter pylori ",$strain,"\"\n";
				print OUT "                     \/mol_type\=\"genomic DNA\"\n";
				print OUT "                     \/strain\=\"",$strain,"\"\n";

				##Printing description for DDBJ format
				print OUT2 $strain,"\tsource\t",$position,"\torganism\tHelicobacter pylori ",$strain,"\n";
				print OUT2 "\t\t\tmol_type\tgenomic DNA\n";
				print OUT2 "\t\t\tstrain\t",$strain,"\n";
				print OUT2 "\t\t\tcountry\tJapan\n";
				print OUT2 "\t\t\thost\tHomo sapiens\n";
				print OUT2 "\t\t\tisolation_source\tHuman stomach\n";

				#Listing gaps
				my $genomeseq = $seq_obj->seq;
#				print $genomeseq,"\n";
				while ($genomeseq =~ /(N+)/g){
					if (length($1)>=10){ #Change here for gap length threshold.
						push (@gapstart,length($`) + 1);
						push (@gapend, length($`) + length($1));
					}
				}
				print $#gapstart,"\n";

			}elsif ($feat_obj->primary_tag eq "CDS"){
				my $gene = "default";
				my $locus_tag;
				my $transl_table;
				my $codon_start;
				my $product = "default";
				my $translation;

				##Extract position information
				my $start = $feat_obj->location->start;
				my $end = $feat_obj->location->end;
				my $strand = $feat_obj->location->strand;
				my $position;
				if ($strand == 1){
					$position = $start."\.\.".$end;
				}elsif ($strand == -1){
					$position = "complement\(".$start."\.\.".$end."\)";
				}

				##Extract locus tag value
				if ($feat_obj->has_tag("note")){
					for my $value2 ($feat_obj->get_tag_values("note")){
						if ($value2 =~ /gene_(\d+)/){
							my $num = sprintf "%04d",$1;
							if ($strain =~ /(.*)_plasmid/){
								$locus_tag = "HP".$1."p_".$num;
							}else{
								$locus_tag = "HP".$strain."_".$num;
							}
							last;
						}
					}
				}
				##Check overlap with gaps and the cds
				my $gapnum = 0;
				my @overlapstart;
				my @overlapend;
				my @overlappedgapnum;
				my $overlaplength = 0;
				while ($gapnum <= $#gapstart){
					if ($start > $gapend[$gapnum]){
						if ($gapnum > $lastprintedgapnum){
							&printgap($gapstart[$gapnum]."\.\.".$gapend[$gapnum]);
							$lastprintedgapnum = $gapnum;
						}
					}elsif ($start <= $gapend[$gapnum]){
						if ($end >= $gapstart[$gapnum]){
							my $x;
							my $y;
							push (@overlappedgapnum, $gapnum);
							if ($start + 30 <= $gapstart[$gapnum]){
								push (@overlapstart, $gapstart[$gapnum]);
								$y = $gapstart[$gapnum];
							}else {
								push (@overlapstart, $start);
								$y = $start;
							}
							if ($end - 30 >= $gapend[$gapnum]){
								push (@overlapend, $gapend[$gapnum]);
								$x = $gapend[$gapnum];
							}else {
								push (@overlapend, $end);
								$x = $end;
							}
							$overlaplength = $overlaplength + ($x-$y);
						}else{

							last;
						}
					}
					$gapnum++;
				}

				my @CDSgapposition;				
				#Remove CDS if more than 75% are gap.
				if ($overlaplength >= ($end-$start)*0.75){
					print $locus_tag,"\tAlmost gap! No print!\n";
					if ($overlappedgapnum[0]>$lastprintedgapnum){
						&printgap($gapstart[$overlappedgapnum[0]]."\.\.".$gapend[$overlappedgapnum[0]]);
						$lastprintedgapnum = $overlappedgapnum[0];
					}
					next;
				}
				#Extract nonoverlapped region in cds!!

				if ($#overlappedgapnum == -1){
					push (@CDSgapposition, "cds_".$start."\.\.".$end);
				}else{
					my $nonoverlapstart = $start;
					my $overlapcount = 0;
					while($overlapcount<=$#overlappedgapnum){
						if ($overlapstart[$overlapcount] > $start && $overlapcount==0){
							push (@CDSgapposition, "cds_".$start."\.\.>".($overlapstart[$overlapcount]-1));
						}
						push (@CDSgapposition, "gap_".$overlappedgapnum[$overlapcount]."_".$gapstart[$overlappedgapnum[$overlapcount]]."\.\.".$gapend[$overlappedgapnum[$overlapcount]]);
						if (($overlapcount+1) <= $#overlappedgapnum){
							push (@CDSgapposition, "cds_<".($overlapend[$overlapcount]+1)."\.\.>".($overlapstart[$overlapcount+1]-1));
						}else{
							if ($overlapend[$overlapcount]< $end){
								push (@CDSgapposition, "cds_<".($overlapend[$overlapcount]+1)."\.\.".$end);
							}
						}
						$overlapcount++;
					}
				}

				##Extract translation table value
				if ($feat_obj->has_tag("transl_table")){
					for my $value3 ($feat_obj->get_tag_values("transl_table")){
						$transl_table = $value3;
						last;
					}				
				}

				##Extract codon start value
				if ($feat_obj->has_tag("codon_start")){
					for my $value4 ($feat_obj->get_tag_values("codon_start")){
						$codon_start = $value4;
						last;
					}				
				}

				my $tempgene = "default";
				##Extract product value
				if ($feat_obj->has_tag("product")){
					for my $value5 ($feat_obj->get_tag_values("product")){
						if ($value5 =~ /^[a-zA-z0-9][a-zA-z0-9]$/){
							last;
						}
						$product = $value5;
						unless ($product =~ /[A-Z][A-Z][A-Z]/){
							$product = lcfirst($product);
						}
						if ($product =~ /(.*)\((.*)\)$/){
							$product = $1;
							$tempgene = $2;
							if ($tempgene !~ /^[a-z][a-z][a-z][A-Z]$/){
								$tempgene = "default";
							}
							print $product,"\t",$tempgene,"\n";
						}
						last;
					}				
				}

				##Extract gene tag value
				if ($feat_obj->has_tag("gene")){
					for my $value1 ($feat_obj->get_tag_values("gene")){
						if ($value1 =~ /_/ || $value1 =~ /^[a-zA-z0-9][a-zA-z0-9]$/){
							last;
						}
						$gene = $value1;
						last;
					}
				}elsif ($tempgene !~ /default/){
					$gene = $tempgene;
				}

				##Extract translation value
				if ($feat_obj->has_tag("translation")){
					for my $value6 ($feat_obj->get_tag_values("translation")){
						$translation = $value6;
						last;
					}				
				}

				##Printing description
				if ($#CDSgapposition == 0){
					&printcds($position, $gene, $locus_tag, $transl_table, $codon_start, $product, $translation);
				}else{
					foreach (@CDSgapposition){
						if (/cds_(.*)/){
							my $positionnp = $1;

							if ($strand == -1){
								$positionnp =~ /\d+\.\.\<?\>?(\d+)/;
								my $startnp = $1;
								if (($end-$startnp)%3 == 1){
									$codon_start = 3;
									print $locus_tag,"\t",$codon_start,"\n";
								}elsif (($end-$startnp)%3 == 2){
									$codon_start = 2;
									print $locus_tag,"\t",$codon_start,"\n";
								}
								$positionnp = "complement\(".$positionnp."\)";
							}else {
								$positionnp =~ /(\d+)\.\.\<?\>?\d+/;
								my $startnp = $1;
								if (($startnp-$start)%3 == 1){
									$codon_start = 3;
									print $locus_tag,"\t",$codon_start,"\n";
								}elsif (($startnp-$start)%3 == 2){
									$codon_start = 2;
									print $locus_tag,"\t",$codon_start,"\n";
								}

							}
							&printcdswotrans($positionnp, $gene, $locus_tag, $transl_table, $codon_start, $product);
						}elsif (/gap_(\d+)_(.*)/){
							if ($1 > $lastprintedgapnum){
								&printgap($2);
								$lastprintedgapnum = $1;
							}
						}
						print $_,"\n";
					}
				}

			}elsif ($feat_obj->primary_tag eq "tRNA"){
				my $product;
				my $locus_tag;

				##Extract position information
				my $start = $feat_obj->location->start;
				my $end = $feat_obj->location->end;
				my $strand = $feat_obj->location->strand;
				my $position;
				if ($strand == 1){
					$position = $start."\.\.".$end;
				}elsif ($strand == -1){
					$position = "complement\(".$start."\.\.".$end."\)";
				}

				##Extract product value
				if ($feat_obj->has_tag("product")){
					for my $value1 ($feat_obj->get_tag_values("product")){
						$product = $value1;
						last;
					}				
				}

				##Building locus tag
				my $num = sprintf "%02d", $trnacount;
				if ($strain =~ /(.*)_plasmid/){
					$locus_tag = "HP".$1."p_t".$num;
				}else{
					$locus_tag = "HP".$strain."_t".$num;
				}
				$trnacount++;

				##Printing description for genbank format
				print OUT "     tRNA            ",$position,"\n";
				print OUT "                     \/locus_tag\=\"",$locus_tag,"\"\n";
				if ($product !~ /default/){
					print OUT "                     \/product\=\"",$product,"\"\n";
				}else {
					print OUT "                     \/product\=\"hypothetical tRNA\"\n";
				}

				##Printing description for DDBJ format
				print OUT2 "\ttRNA\t",$position,"\tlocus_tag\t",$locus_tag,"\n";
				if ($product !~ /default/){
					print OUT2 "\t\t\tproduct\t",$product,"\n";
				}else {
					print OUT2 "\t\t\tproduct\thypothetical tRNA\n";
				}

			}elsif ($feat_obj->primary_tag eq "rRNA"){
				my $product;
				my $locus_tag;

				##Extract position information
				my $start = $feat_obj->location->start;
				my $end = $feat_obj->location->end;
				my $strand = $feat_obj->location->strand;
				my $position;
				if ($strand == 1){
					$position = $start."\.\.".$end;
				}elsif ($strand == -1){
					$position = "complement\(".$start."\.\.".$end."\)";
				}

				##Extract product value
				if ($feat_obj->has_tag("product")){
					for my $value1 ($feat_obj->get_tag_values("product")){
						$product = $value1;
						if ($product =~ /\-16S/){
							$product = "16S ribosomal RNA";
						}
						last;
					}				
				}

				##Building locus tag
				my $num = sprintf "%02d", $rrnacount;
				if ($strain =~ /(.*)_plasmid/){
					$locus_tag = "HP".$1."p_r".$num;
				}else{
					$locus_tag = "HP".$strain."_r".$num;
				}
				$rrnacount++;

				##Printing description
				print OUT "     rRNA            ",$position,"\n";
				print OUT "                     \/locus_tag\=\"",$locus_tag,"\"\n";
				if ($product !~ /default/){
					print OUT "                     \/product\=\"",$product,"\"\n";
				}else {
					print OUT "                     \/product\=\"hypothetical rRNA\"\n";
				}

				##Printing description for DDBJ format
				print OUT2 "\trRNA\t",$position,"\tlocus_tag\t",$locus_tag,"\n";
				if ($product !~ /default/){
					print OUT2 "\t\t\tproduct\t",$product,"\n";
				}else {
					print OUT2 "\t\t\tproduct\thypothetical rRNA\n";
				}

			}
		}
	}
	open (IN2, $file);
	my $seqflag = 0;
	while (<IN2>){
		if ($seqflag == 1){
			print OUT $_;
		}elsif ($_=~/ORIGIN/){
			$seqflag = 1;
			print OUT $_;
		}
	}
}


sub printcds{
	my $position = $_[0];
	my $gene = $_[1];
	my $locus_tag = $_[2];
	my $transl_table = $_[3];
	my $codon_start = $_[4];
	my $product = $_[5];
	my $translation = $_[6];

	print OUT "     CDS             ",$position,"\n";
	if ($gene !~ /default/){
		print OUT "                     \/gene\=\"",$gene,"\"\n";
	}
	print OUT "                     \/locus_tag\=\"",$locus_tag,"\"\n";
	print OUT "                     \/transl_table\=\"",$transl_table,"\"\n";
	print OUT "                     \/codon_start\=\"",$codon_start,"\"\n";
	if ($product !~ /default/){
		print OUT "                     \/product\=\"",$product,"\"\n";
	}else {
		print OUT "                     \/product\=\"hypothetical protein\"\n";
	}
	print OUT "                     \/translation\=\"",$translation,"\"\n";

	##Printing description for DDBJ format
	print OUT2 "\tCDS\t",$position,"\tlocus_tag\t",$locus_tag,"\n";
	if ($gene !~ /default/){
		print OUT2 "\t\t\tgene\t",$gene,"\n";
	}
	print OUT2 "\t\t\ttransl_table\t",$transl_table,"\n";
	print OUT2 "\t\t\tcodon_start\t",$codon_start,"\n";
	if ($product !~ /default/){
		print OUT2 "\t\t\tproduct\t",$product,"\n";
	}else{
		print OUT2 "\t\t\tproduct\thypothetical protein\n";
	}
	print OUT2 "\t\t\ttranslation\t",$translation,"\n";
}

sub printcdswotrans{
	my $position = $_[0];
	my $gene = $_[1];
	my $locus_tag = $_[2];
	my $transl_table = $_[3];
	my $codon_start = $_[4];
	my $product = $_[5];
#	my $translation = $_[6];
	
	print OUT "     CDS             ",$position,"\n";
	if ($gene !~ /default/){
		print OUT "                     \/gene\=\"",$gene,"\"\n";
	}
	print OUT "                     \/locus_tag\=\"",$locus_tag,"\"\n";
	print OUT "                     \/transl_table\=\"",$transl_table,"\"\n";
	print OUT "                     \/codon_start\=\"",$codon_start,"\"\n";
	if ($product !~ /default/){
		print OUT "                     \/product\=\"",$product,"\"\n";
	}else {
		print OUT "                     \/product\=\"hypothetical protein\"\n";
	}
#	print OUT "                     \/translation\=\"",$translation,"\"\n";

	##Printing description for DDBJ format
	print OUT2 "\tCDS\t",$position,"\tlocus_tag\t",$locus_tag,"\n";
	if ($gene !~ /default/){
		print OUT2 "\t\t\tgene\t",$gene,"\n";
	}
	print OUT2 "\t\t\ttransl_table\t",$transl_table,"\n";
	print OUT2 "\t\t\tcodon_start\t",$codon_start,"\n";
	if ($product !~ /default/){
		print OUT2 "\t\t\tproduct\t",$product,"\n";
	}else{
		print OUT2 "\t\t\tproduct\thypothetical protein\n";
	}
#	print OUT2 "\t\t\ttranslation\t",$translation,"\n";
}

sub printgap{
	my $gappos = $_[0];

	print OUT "     assembly_gap    ",$gappos,"\n";
	print OUT "                     \/gap_type\=\"within scaffold\"\n";
	print OUT "                     \/linkage_evidence\=\"paired-ends\"\n";

	print OUT2 "\tassembly_gap\t",$gappos,"\testimated_length\tknown\n";
	print OUT2 "\t\t\tgap_type\twithin scaffold\n";
	print OUT2 "\t\t\tlinkage_evidence\tpaired-ends\n";
}

