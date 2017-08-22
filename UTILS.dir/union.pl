#!/usr/bin/env perl

use strict;
use Pod::Usage;
use Getopt::Long;

=head1 NAME

union.pl

=head1 SYNOPSIS

USAGE: union.pl --in <file> --seg_name <segment col num> --seg_start <start col num> --seg_end <end col num>

=head1 DESCRIPTION

 Segment operations - UNION
 Generate the union of all segments described in a file in the following format

 Input:
 seg1	4	5
 seg1	6	10
 seg1	7	12
 seg2	1	2	
 seg1	1	2
 seg2	3	2

 Output:
 seg1	4	12
 seg2	1	3
 seg1	1	2
 
 For testing purposes use as input file:
 /seq/aspergillus1/gustavo/devel/tab_file/file_for_tests/union.pl.txt

=head1 AUTHOR
Gustavo C. Cerqueira 2013

=cut


my $in;
my $segment_col;
my $start_col;
my $end_col;

GetOptions(	'in=s'			=> \$in,
			'seg_name=i'	=> \$segment_col,
			'seg_start=i' 	=> \$start_col,
			'seg_end=i'		=> \$end_col );
		
		

if( not defined($in) ){
   print STDERR "!!!! Parameter --in is required !!!!\n\n";  
   pod2usage(-verbose => 1 ,-exitval => 2);
} 

if( not defined($segment_col) ){
   print STDERR "!!!! Parameter --seg_name is required !!!!\n\n";  
   pod2usage(-verbose => 1 ,-exitval => 2);
} 

if( not defined($start_col) ){
   print STDERR "!!!! Parameter --seg_start is required !!!!\n\n";  
   pod2usage(-verbose => 1 ,-exitval => 2);
} 

if( not defined($end_col) ){
   print STDERR "!!!! Parameter --seg_end is required !!!!\n\n";  
   pod2usage(-verbose => 1 ,-exitval => 2);
} 


open IN, $in or die "Unable to open file $in\n";

my %segVec;
my %sortedSegVec;

while(<IN>){
	my $line = $_;
	chomp( $line );
	my @cols = split "\t", $line;
	
	my $seg_name = $cols[ $segment_col - 1 ];
	my $start    = $cols[ $start_col - 1 ];
	my $end      = $cols[ $end_col - 1 ];
	
	#print " $seg_name $start $end \n";
	
	if( $end > $start ){
		push @{$segVec{$seg_name}}, { type => "start", coord => $start };
		push @{$segVec{$seg_name}}, { type => "end", coord => $end };
	}else{
		push @{$segVec{$seg_name}}, { type => "start", coord => $end };
		push @{$segVec{$seg_name}}, { type => "end", coord => $start };
	}
	
}

close(IN);

#print $segVec{seg1}[0]{type} . "\n";
#print $segVec{seg1}[0]{coord} . "\n";

foreach my $seg_name (keys %segVec) {
	@{$sortedSegVec{$seg_name}} = sort {$a->{coord} <=> $b->{coord}} @{$segVec{$seg_name}};
	
	my $union_start = -1;
	my $collapsed_ranges = 0;
	 
	foreach my $item ( @{$sortedSegVec{$seg_name}} ){
		my $type  = $item->{type};
		my $coord = $item->{coord};
		
		if ( $type eq "start" ){
			#print "start: $coord\n";
			$union_start = $coord if ( $collapsed_ranges == 0 );				
			$collapsed_ranges++;
		}else{
			#print "end: $coord\n";
			$collapsed_ranges--;
			my $union_end = $coord;
			print "$seg_name\t$union_start\t$union_end\n"  if ( $collapsed_ranges == 0 );				
		}
	}
		
}  

#getc();
#print $sortedSegVec{seg1}[0]{type} . "\n";
#print $sortedSegVec{seg1}[0]{coord} . "\n";


# Evaluating union




	
