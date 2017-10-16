#! /usr/bin/perl -w
use strict;

my %info;
my %geneinfo;
open(IN,$ARGV[0]) or die;
while(<IN>){
	chomp;
	next if (/^#/);
	my ($gene,$transcript,$chr,$strand,$txstart,$txend,$cdsstart,$cdsend,$exonnumber,$exonstarts,$exonends) = (split /\t/);
	my @exonstarts = split /,/, $exonstarts;
	my @exonends = split /,/, $exonends;
	my $length = 0;
	for (my $i=0; $i<@exonstarts; $i++){
		my $exonlen = $exonends[$i] - $exonstarts[$i] ;
		$length += $exonlen;	
	}

	
	$info{$gene}{$transcript}{"length"} = $length;
	$info{$gene}{$transcript}{"exonnumber"} = $exonnumber; 
	

}
close IN;

open(OUT1, '>>', 'gencode.gtf.refFlat.geneinfo');
open(OUT2, '>>', 'gencode.gtf.refFlat.transcriptinfo');

print OUT1 "gene_id\tNtranscripts\ttranscriptlength\tNexon\n";
print OUT2 "gene_id\ttranscript\ttranscriptlength\tNexon\n";

foreach my $gene (keys %info)
{
	my $ntranscripts = 0;
	my $length = 0 ;
	my $exonnumber = 0;
	foreach my  $transcript (keys %{$info{$gene}})
	{
		$ntranscripts++;
		$length += $info{$gene}{$transcript}{"length"};	
		$exonnumber += $info{$gene}{$transcript}{"exonnumber"};
	print OUT2 $gene,"\t",$transcript,"\t",$length,"\t",$exonnumber,"\n";
	}
	print OUT1 $gene,"\t",$ntranscripts,"\t",$length/$ntranscripts,"\t",$exonnumber/$ntranscripts, "\n";
}

close OUT1;
close OUT2;

