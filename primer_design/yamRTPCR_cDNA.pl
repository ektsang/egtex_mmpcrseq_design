#!/usr/bin/perl -w

=head1 NAME

yamPCR.pl: Yet another multiplex PCR primer design program

=head1 SYNOPSIS

yamPCR.pl --targetFile='./targetInfo.txt' [--annotation='UCSC_TableBrowser_gencode19'] [--maxAmpliconSize=400] [--minAmpliconSize=150] [--numberCandidatePrimers=20] [--maxPrimerTm=64] [--minPrimerTm=57] [--excludedPrimers='excluded-primers.txt'] [--includedPrimers='included-primers.txt']

Note: all arguments within the square brackets are optional.

=head1 AUTHORS - Kun Zhang

Email: kzhang@genetics.med.harvard.edu

=head1 DEPENDENCY

The following CPAN modules are required:
(1) DBI
(2) Bio::Seq
(3) Bio::SeqIO
(4) Bio::SearchIO
(5) Bio::Tools::BPlite
(6) Algorithm::Numerical::Shuffle qw /shuffle/
(7) FileHandle
(8) IPC::Open3
(9) Getopt::Long


Additional programs required:
(1) MIT Primer3 (the subroutine for Tm calculation was written in order to use the SantaLucia's NN parameters)
(2) A repeat library of the target organism
(3) NCBI Blast
(4) A formatted blast database 

=head1 LAST MODIFIED

08/06/2005

=cut

=head1 FUNCTIONS

=cut
use strict;

use lib qw(/usr/lib64/perl5);


use DBI;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Algorithm::Numerical::Shuffle qw /shuffle/;
use Sort::Array qw/Sort_Table/;
use File::Temp qw/ tempfile tempdir /;
use IPC::Open3;
use Getopt::Long;

my @bestPossiblePrimer;
my @candidateDividers;
my %primerSetCompatibleTable;
my %ID2Seq;

#====================================================
#----------Parameters need to be customized----------
#====================================================
my $BLAST_DIR = "/users/etsang/software/ncbi-blast-2.5.0+/bin/"; #TO BE CHANGED
## NOTE: $GENOMEDIR is expected to exist as an environment variable, set in bashrc, for instance
my $GENOME_DATA_DIR = $ENV{'GENOMEDIR'} . "/egtex_masked/"; #contains (masked) genome, each chr in a fasta file. (e.g., chr1.fa)
my $PRIMER3_EXE = "./primer3_KZ";
my $CHECKPD_EXE = "./checkPD";
my $GENOME_BLAST_DB = "hg19.gencode19.blast.db";
my $REP_LIB = "humrep_and_simple.fa";
my $NNParamFile = "NN_param.txt";
my $GENE_ANNOTATION_FULL = "UCSC_TableBrowser_gencode19"; # only used to get transcript to gene mapping; should contain the transcripts in the blast database
my $GENE_ANNOTATION = "UCSC_TableBrowser_gencode19"; # should be equal to or a subset of the full annotation
#====================================================

my $readLen = 75;
my $primerMaxLen = 25;
my $primerMinLen = 17;
my $primerOptLen = 21;
my $primerMaxTm = 64;
my $primerMinTm = 57;
my $primerMaxGC = 70;
my $primerMinGC = 30;
my $ampliconMaxLen = 400;
my $ampliconMinLen = 150;
my $PCRSizeLimit = 10000; # max amplicon length of interactions
my $primerMaxPolyX = 4;
my $numberCandidatePrimers = 20;
my $primerSelfAny = 8.0;
my $primerSelfEnd = 4.0;
my $specificity_threshold = 1.0001;
my %excludedPrimers;
my @includedPrimers;




#Some constants
my $R = 1.987;
my %deltaH;
my %deltaS;

=head2 main

 Title   : main
 Function: The main entry of this program
 Returns : None
 Args    : None

=cut

sub main(){

    if(scalar(@ARGV) < 1){
	print STDERR "yamPCR.pl: Yet another multiplex PCR primer design program.\n";
	print STDERR "USAGE:  yamPCR.pl --targetFile='./targetInfo.txt' [--annotation='UCSC_TableBrowser_gencode19'] [--maxAmpliconSize=400] [--minAmpliconSize=150] [--numberCandidatePrimers=20] [--maxPrimerTm=64] [--minPrimerTm=57] [--excludedPrimers='excluded-primers.txt'] [--includedPrimers='included-primers.txt']\n";
	print STDERR "Format of the target information file:\n";
	print STDERR "locusID  chromosome  chrom_start  chrom_end strand\n";
	print STDERR "Example:\n";
	print STDERR "geneA  chr3  1238943 1239305 +\n";
	print STDERR "geneB  chr5  572398 573100 -\n";
	return;
    }
    
    my $start_time = time();
    
    my $targetInfoFile;
    my $excludedPrimerFile;
    my $includedPrimerFile;
    my $result = GetOptions("targetFile=s" => \$targetInfoFile,
			    "annotation=s" => \$GENE_ANNOTATION,
			    "maxAmpliconSize=i" => \$ampliconMaxLen,
			    "minAmpliconSize=i" => \$ampliconMinLen,
			    "numberCandidatePrimers=i" => \$numberCandidatePrimers,
			    "maxPrimerTm=i" => \$primerMaxTm,
			    "minPrimerTm=i" => \$primerMinTm,
			    "excludedPrimers=s" => \$excludedPrimerFile,
			    "includedPrimers=s" => \$includedPrimerFile);
    # for testing purposes - Emily
    print STDERR "Gene annotation for target search: ", $GENE_ANNOTATION, "\n";
    
    readNNParam($NNParamFile);
    
    #Read a list of primers that should be excluded.
    readExcludedPrimers($excludedPrimerFile) if($excludedPrimerFile);

    #Read a list of primers that should be included.
    readIncludedPrimers($includedPrimerFile) if($includedPrimerFile);

    #Read the target informations
    my %ampliconInfo = parseTargetInfo($targetInfoFile);
    
    die("No valid target read. Aborted!\n") if(scalar(keys(%ampliconInfo)) < 1);
    
    #Get correspondence between transcripts and genes
    my %transcriptGeneMap = parseAnnotation($GENE_ANNOTATION_FULL);
    
    #Design candidate primers using the Primer3 program
    getCandidatePrimer(\%ampliconInfo, $numberCandidatePrimers, $ampliconMinLen, $ampliconMaxLen);
    
    #Map all primer sequences to the template sequence (human genome) by Blast search
    my %primerMappingInfo = mapPrimers(\%ampliconInfo, $GENOME_BLAST_DB);

    #Filter out primers with multiple amplicons, primer dimers, etc.
    primerFilter(\%primerMappingInfo, \%ampliconInfo, $numberCandidatePrimers, $ampliconMinLen, $ampliconMaxLen, \%transcriptGeneMap);

    #Search a set of multiplex primers
    getMultiplexPrimerSet(\%ampliconInfo, \%primerMappingInfo, $PCRSizeLimit, \%transcriptGeneMap);

    #Generate some simple reports
    makeReports(\%ampliconInfo);

    print STDERR "Total time used: ", int((time()-$start_time)/60), " minutes\n";
}



=head2 parseTargetInfo

 Title   : parseTargetInfo
 Function: To parse the target information file
 Returns : A hash table containing target informations.
 Args    : The input file name.

=cut

sub parseTargetInfo(){
    my $fileName = shift;
    open(INFILE, "$fileName") || die("Error in opening target file ($fileName)!\n");
    my %ampliconInfo;
    while(my $line = <INFILE>){
	chop($line);
	next if (!$line);
	my ($locus,$chr,$start,$end,$strand_information) = split(/[\t ]+/, $line);
	
	next if(!$locus || !$chr || !$start || !$end ||! $strand_information);
	next if($chr !~ /chr[0-9XYxy]+/); #rui 10/20 #chr[0-9XY]
	next if($end < $start);
	if($end-$start+1 > $ampliconMaxLen){
	    print STDERR "$locus: requested amplicon size exceeds maximal amplicon length. Skipped\n";
	    next;
	}
	$chr =lc($chr);		
	$ampliconInfo{$locus}->{'chr'} = $chr;
	$ampliconInfo{$locus}->{'request'}->{'start'} = $start;
	$ampliconInfo{$locus}->{'request'}->{'end'} = $end;
	$ampliconInfo{$locus}->{'request'}->{'strand'} = $strand_information;
	if ($strand_information eq "+" ) {
	    $ampliconInfo{$locus}->{'limit'}->{'start'} = $start-($readLen-5);
	    $ampliconInfo{$locus}->{'limit'}->{'end'} = $end+($ampliconMaxLen-$primerMinLen);
	}
	elsif ($strand_information eq "-" )
	{
	    $ampliconInfo{$locus}->{'limit'}->{'start'} = $start-($ampliconMaxLen-$primerMinLen);
	    $ampliconInfo{$locus}->{'limit'}->{'end'} = $end+($readLen-5);
	}
    }
    return %ampliconInfo;
}

=head2 parseAnnotation

 Title   : parseAnnotation
 Function: To parse the annotation file to extract transcript to gene correspondence
 Returns : A hash table mapping transcript names to gene names
 Args    : The input file name.

=cut

sub parseAnnotation(){
    my $fileName = shift;
    open(INFILE, "$fileName") || die("Error in opening target file ($fileName)!\n");
    my %transcriptGeneMap;
    while(my $line = <INFILE>){
	chop($line);
	next if (!$line);
	my @fields = split(/[\t ]+/, $line);
	$transcriptGeneMap{$fields[1]} = $fields[12]; # transcript: gene
    }
    return %transcriptGeneMap;
}


=head2 makeReports

 Title   : makeReports
 Function: To generate a report.
 Returns : None
 Args    : $h_ampliconInfo	The handle of %ampliconInfo hash table

=cut

sub makeReports(){
    my $h_ampliconInfo = shift;
    my %ampliconInfo = %{$h_ampliconInfo};
    my $totalAmpliconLength = 0;
    foreach my $amplicon (sort keys(%ampliconInfo)) {
	if ($ampliconInfo{$amplicon}->{'LP'}){
	    print $amplicon, "\t";
	    print $ampliconInfo{$amplicon}->{'LP'}, "\t",
	    shortOligoTm($ampliconInfo{$amplicon}->{'LP'}, 200, 1.5, 50, 0.4, 0.0, 50), "\t",
	    $ampliconInfo{$amplicon}->{'RP'}, "\t",
	    shortOligoTm($ampliconInfo{$amplicon}->{'RP'}, 200, 1.5, 50, 0.4, 0.0, 50), "\t";
	    my @words = split(/[:]/, $ampliconInfo{$amplicon}->{'amplicon'});
	    my $RPlength = length($ampliconInfo{$amplicon}->{'RP'});
	    $words[3] = int($words[3]) + $RPlength - 1;
	    print $words[1], ":", $words[2], "-", $words[3], "\t";
	    $totalAmpliconLength += $words[3]-$words[2];
	    print $words[3]-$words[2], "\t";
	    print "\n";
	}
    }
    print STDERR "Total length of region covered by these amplicons is $totalAmpliconLength (bp).\n";
}


=head2 readExcludedPrimers

 Title   : readExcludedPrimers
 Function: To read a list of primers that should be excluded from a file.
 Returns : None
 Args    : $fileName	The name of file containing excluded primers.

=cut

sub readExcludedPrimers(){
    my $fileName = shift;
    open(INFILE, "$fileName") || die("Failed to read the file ($fileName) containing excluded primers!\n");
    while(my $line = <INFILE>){
	chop($line);
	$line =~ s/[ \t]//g;
	$excludedPrimers{uc($line)} = 1;
    }
    close(INFILE);
}

=head2 readIncludedPrimers

 Title   : readIncludedPrimers
 Function: To read a list of primers that should be included from a file.
 Returns : None
 Args    : $fileName	The name of file containing included primers.

=cut

sub readIncludedPrimers(){
    my $fileName = shift;
    open(INFILE, "$fileName") || die("Failed to read the file ($fileName) containing included primers!\n");
    while(my $line = <INFILE>){
	chop($line);
	$line =~ s/[ \t]//g;
	$line =~ s/\r//g;
	print STDERR "$line\n";
	push(@includedPrimers, uc($line));
    }
    close(INFILE);
}


=head2 getCandidatePrimer

 Title   : getCandidatePrimer
 Function: To find a set of candidate primer using the MIT Primer3 program
 Returns : None
 Args    :  $h_ampliconInfo      a reference to the %ampliconInfo hash table
            $numCandidatePrimer  number of candidate primers wanted for each amplicon
            $minLen 	         maximal amplicon size
            $maxLen 	         minimal amplicon size

=cut

sub getCandidatePrimer(){
    my $h_ampliconInfo = shift;
    my $numCandidatePrimer = shift;
    my $minLen = shift;
    my $maxLen = shift;
    my %ampliconInfo = %{$h_ampliconInfo};

    # emily added variables to count drop outs
    my $numberNoTranscript = 0;
    my $numberNoPrimer = 0;
    my $numberTooShort = 0;

    my %EndAllowed;
    $EndAllowed{'A'} = 1;
    $EndAllowed{'T'} = 1;
    $EndAllowed{'C'} = 1;

    my $numPrimer3CandidatePrimer = $numCandidatePrimer*4;
    print STDERR "Searching candidate primers using Primer3 ... \n";

    #group amplicons based on chromosome
    my %chrSet;
    foreach my $ampliconID (keys(%ampliconInfo)){
	next if (!$ampliconInfo{$ampliconID}->{'chr'});
	push(@{$chrSet{$ampliconInfo{$ampliconID}->{'chr'}}}, $ampliconID);
    }

    foreach my $chr (keys(%chrSet)){
	my $templateFile = $GENOME_DATA_DIR. $chr . ".fa";
	my $SeqIN = Bio::SeqIO->new(-file => $templateFile , '-format' => 'Fasta');
	my $templateSeq = $SeqIN->next_seq();
	print STDERR "amplicon_ID\tstart_limit\tend_limit\trange\ttarget_start\ttarget_length\t[minLen maxLen]\n";
	foreach my $ampliconID (sort @{$chrSet{$chr}}){
	    my $start = $ampliconInfo{$ampliconID}->{'limit'}->{'start'};
	    $start = 0 if $start < 0;
	    my $end = $ampliconInfo{$ampliconID}->{'limit'}->{'end'};
	    $end =$templateSeq->length() if ($end > $templateSeq->length);
	    #TRY and make primers to seq within read length of LEFT of site
	    my $chromfound = 0;
	    my $position = int($ampliconInfo{$ampliconID}->{'request'}->{'start'});
	    my $startbuffer = int($ampliconInfo{$ampliconID}->{'request'}->{'start'}) - int($ampliconInfo{$ampliconID}->{'limit'}->{'start'});
	    my $endbuffer = int($ampliconInfo{$ampliconID}->{'limit'}->{'end'}) - int($ampliconInfo{$ampliconID}->{'request'}->{'start'});
	    my ($editpos, $startpos, $endpos);
	    my $geneseq = '';
	    open (GENEFILE, $GENE_ANNOTATION) or die "error opening gene annotation file: $!\n";
	    GENELOOP: while (<GENEFILE>) { #for each variant, loop through gene annotation file and check which exon each site is
		chomp;
		my @fieldsref = split;
		my ($chromref, $strand, $txstart, $txend) = ($fieldsref[2], $fieldsref[3], $fieldsref[4], $fieldsref[5]);
		$chromref = lc($chromref);
		if ($chr eq $chromref) {
		    $chromfound = 1;
		    my $exoncount = int($fieldsref[8]);
		    my @exonstarts = split(/\,/, $fieldsref[9]);
		    my @exonends = split(/\,/, $fieldsref[10]);
		    for (my $i = 0; $i < @exonstarts; $i++) { #for each exon, check if variant lies in it
			$geneseq = join '', $geneseq, $templateSeq->subseq(int($exonstarts[$i])+1, int($exonends[$i])); # accumulate gene seq
			$editpos = length($geneseq) - (int($exonends[$i])-$position) if (int($exonstarts[$i]) < $position && int($exonends[$i]) >= $position);  #exons are 0-indexed!!
		    }
		    if ($editpos) {
			$startpos = $editpos - $startbuffer;
			$endpos = $editpos + $endbuffer;
			$startpos = 0 if ($startpos < 0);
			$endpos = length($geneseq) + 1 if ($endpos > length($geneseq) + 1);
		    }
		    last GENELOOP if ($editpos);
		} 
		elsif ($chromfound == 1) { #if all genes on the chromosome of variant are checked, end search
		    last GENELOOP;
		}
	    }
	    if (!$editpos) {
		print STDERR "No transcript containing $ampliconID, will be ignored.\n";
		$numberNoTranscript ++;
		delete($$h_ampliconInfo{$ampliconID});
		next;
	    }

	    my ($h_fileIN, $name_fileIN) = tempfile();
	    my ($h_fileOUT, $name_fileOUT) = tempfile();
	    my @input;
	    push @input, "PRIMER_SEQUENCE_ID=", $templateSeq->id, ":$start - $end\n";
	    my $inputseq="";
	    $inputseq = substr $geneseq, $startpos, $endpos-$startpos;
	    # reverse complement if on reverse strand
	    if ($ampliconInfo{$ampliconID}->{'request'}->{'strand'} eq "-") {
		$inputseq = uc ($inputseq) ;
		$inputseq=reverse ($inputseq);
		$inputseq=~tr/ATGC/TACG/;
		
	    }
	    print STDERR "$inputseq\t$startpos\t$endpos\n";
	    push @input, "SEQUENCE=", $inputseq , "\n";
	    my $targetStart;
	    if ($ampliconInfo{$ampliconID}->{'request'}->{'strand'} eq "+") {
		$targetStart = $editpos-$startpos-3;
	    }
	    elsif ($ampliconInfo{$ampliconID}->{'request'}->{'strand'} eq "-") {
		$targetStart = $endpos-$editpos-3;
	    }			

	    # emily changed target length (was set to 5, which works well for snps but not indels)
	    my $targetLen = $ampliconInfo{$ampliconID}->{'request'}->{'end'} - $ampliconInfo{$ampliconID}->{'request'}->{'start'} + 5;
	    
	    # emily added for debugging
	    my $inputseqLength = length($inputseq);
	    if ($inputseqLength < $minLen) {
		print STDERR "$ampliconID input sequence for primer3 is shorter ($inputseqLength) than minimum amplicon length ($minLen).\n";
		$numberTooShort ++;
		delete($$h_ampliconInfo{$ampliconID});
		next;
	    }
	    print STDERR "$ampliconID\t$start\t$end\t", $end-$start, "\t$targetStart\t$targetLen\t[$minLen,$maxLen]\n";
	    push @input, "TARGET=$targetStart,$targetLen\n",
	    "PRIMER_ATC_CLAMP=2\n",
	    "PRIMER_SELF_END=$primerSelfEnd\n",
	    "PRIMER_SELF_ANY=$primerSelfAny\n",
	    "PRIMER_OPT_SIZE=$primerOptLen\n",
	    "PRIMER_MIN_SIZE=$primerMinLen\n",
	    "PRIMER_MAX_SIZE=$primerMaxLen\n",
	    "PRIMER_MIN_TM=$primerMinTm\n",
	    "PRIMER_OPT_TM=", ($primerMinTm+$primerMaxTm)/2, "\n",
	    "PRIMER_MAX_TM=$primerMaxTm\n",
	    "PRIMER_MIN_GC=$primerMinGC\n",
	    "PRIMER_MAX_GC=$primerMaxGC\n",
	    "PRIMER_MAX_POLY_X=$primerMaxPolyX\n",
	    "PRIMER_MAX_MISPRIMING=12.0\n",
	    "PRIMER_MAX_END_STABILITY=8.0\n",
	    "PRIMER_MISPRIMING_LIBRARY=$REP_LIB\n",
	    "PRIMER_NUM_RETURN=$numPrimer3CandidatePrimer\n",
	    "PRIMER_PRODUCT_SIZE_RANGE=", $minLen, "-" , $maxLen, "\n=\n";

	    my $cmd = $PRIMER3_EXE;
	    print $h_fileIN @input; 
	    $h_fileIN->close();
	    system("$cmd < $name_fileIN > $name_fileOUT");
	    my $LP = '';
	    my $RP = '';
	    while(my $line =$h_fileOUT->getline){	
		chop($line);
		if($line =~ /PRIMER_LEFT[_0-9]+SEQUENCE/){
		    my @words = split(/=/, $line);
		    $LP = uc($words[1]);
		}elsif($line =~ /PRIMER_RIGHT[_0-9]+SEQUENCE/){
		    my @words = split(/=/, $line);
		    $RP = uc($words[1]);
		    if($LP && !$excludedPrimers{uc($LP)}
		       && !$excludedPrimers{uc($RP)}
		       && $EndAllowed{substr($LP, length($LP)-1, 1)}      # 3' end of the primer must be A or T
		       && $EndAllowed{substr($RP, length($RP)-1, 1)}
		       && substr($LP, length($LP)-5, 5) =~ /[CG]/gi       # at least one C/G on the last five base pairs of 3' end
		       && substr($RP, length($RP)-5, 5) =~ /[CG]/gi){
			#&& substr($LP, length($LP)-2, 2) !~ /G/gi       # no G on the last two base pairs of 3' end
			#&& substr($RP, length($RP)-2, 2) !~ /G/gi){
			push(@{$ampliconInfo{$ampliconID}->{'cLP'}}, $LP);
			push(@{$ampliconInfo{$ampliconID}->{'cRP'}}, $RP);
		    }
		}elsif($line =~ /ERROR/){
		    print STDERR "Error message from Primer3:\n$line\n";
		}
	    }
	    unlink($name_fileIN);
	    unlink($name_fileOUT);

	    if(!$ampliconInfo{$ampliconID}->{'cLP'}){
		print STDERR "Primer3 couldn't find any primer for $ampliconID, will be ignored.\n";
		$numberNoPrimer ++;
		delete($$h_ampliconInfo{$ampliconID});
	    }
	}
    }
    
    
    foreach my $amp (keys(%{$h_ampliconInfo})){
	next if (!$$h_ampliconInfo{$amp}->{'cLP'});
	my $size = scalar(@{$$h_ampliconInfo{$amp}->{'cLP'}});
	my @LPs = @{$$h_ampliconInfo{$amp}->{'cLP'}};
	my @RPs = @{$$h_ampliconInfo{$amp}->{'cRP'}};
	my @primer_list;
	for(my $i=0; $i<$size; $i++){
	    my @primers = ($LPs[$i], $RPs[$i]);
	    push(@primers, max_dG_profile_mismatch(\@primers));
	    push(@primer_list, join(",", @primers));
	}
	
	my @sorted_primer_list = Sort_Table(
	    cols		=> '3',
	    field		=> '3',
	    sorting  	=> 'ascending',
	    structure	=> 'csv',
	    separator	=> ',',

	    data		=> \@primer_list
	    );

	delete($$h_ampliconInfo{$amp}->{'cLP'});
	delete($$h_ampliconInfo{$amp}->{'cRP'});
	$size = $numCandidatePrimer*2 if($size > $numCandidatePrimer*2);
	my @scores;
	for(my $i=0; $i<$size; $i++){
	    my ($LP, $RP, $score) = split(",", $sorted_primer_list[$i]);
	    push(@{$$h_ampliconInfo{$amp}->{'cLP'}}, $LP);
	    push(@{$$h_ampliconInfo{$amp}->{'cRP'}}, $RP);
	    push(@scores, $score);
	}	
	#print join(";", @scores), "\n";
	undef(@primer_list);			
    }
    
    # Emily added print statement to give number of amplicons for which primer design failed
    print STDERR "Number of amplicons for which primers were not designed because no transcript contained the site of interest: $numberNoTranscript.\n";
    print STDERR "Number of amplicons for which primers were not designed because primer3 could not find primers: $numberNoPrimer.\n";
    print STDERR "Number of amplicons for which primers were not designed because template sequence was shorter than amplicon min length: $numberTooShort.\n";
}


=head2 mapPrimers

 Title   : mapPrimers
 Function: To map all candidate primers to the template sequence (such as the human genome) by Blast
 Returns : None
 Args    : $h_ampliconInfo    a reference to the %ampliconInfo hash table
           $templateFile      the file name of the template sequence

=cut

sub mapPrimers(){
    my $h_ampliconInfo = shift;
    my $templateFile = shift;
    my %ampliconInfo = %{$h_ampliconInfo};
    my @primerList;
    my %primerMappingInfo;

    my $primerMappingInfoFile = "primer_mapping_info.txt";
    
    foreach my $id (keys(%ampliconInfo)){
	@primerList = @{$ampliconInfo{$id}->{'cLP'}};
	my $i = 0;
	foreach my $primer (@primerList){
	    next if (!$primer);
	    $primerMappingInfo{$primer}->{'tag'} = $id . "_cLP_" . $i ;
	    $ID2Seq{$primerMappingInfo{$primer}->{'tag'}} = $primer;
	    $i++;
	    $primerMappingInfo{$primer}->{'T95'} = shortOligoTm($primer, 200, 1.5, 50, 0.4, 0.0,  95);
	}
	@primerList =  @{$ampliconInfo{$id}->{'cRP'}};
	$i = 0;
	foreach my $primer (@primerList){
	    next if (!$primer);
	    $primerMappingInfo{$primer}->{'tag'} = $id . "_cRP_" . $i ;
	    $ID2Seq{$primerMappingInfo{$primer}->{'tag'}} = $primer;
	    $i++;
	    $primerMappingInfo{$primer}->{'T95'} = shortOligoTm($primer, 200, 1.5, 50, 0.4, 0.0, 95);
	}
	undef(@primerList);
    }

    my $i = 0;
    foreach my $includedPrimer (@includedPrimers){
	next if($includedPrimer !~ /[acgtACGT]/);
	$primerMappingInfo{$includedPrimer}->{'tag'} = "IncludedPrimers:" . $i ;
	$ID2Seq{$primerMappingInfo{$includedPrimer}->{'tag'}} = $includedPrimer;
	$i++;
	$primerMappingInfo{$includedPrimer}->{'T95'} = shortOligoTm($includedPrimer, 200, 1.5, 50, 0.4, 0.0, 95);
    }

    readPrimerMappingInfo($primerMappingInfoFile, \%primerMappingInfo);		
    
    my $unmapped_primer = 0;

    my $primerFile = sprintf("primers_%x.fa", rand(65536));
    open(PRIMER, ">$primerFile") || die("Error in opening temporary primer file!\n");
    foreach my $primer (keys(%primerMappingInfo)){
	next if($primer !~ /[a-zA-Z]/);	
	next if($primerMappingInfo{$primer}->{'targetPos'});	
	next if(!$primerMappingInfo{$primer}->{'tag'});
	print PRIMER ">", $primerMappingInfo{$primer}->{'tag'}, "\n";
	print PRIMER $primer, "\n\n";
	$unmapped_primer++;
    }
    close(PRIMER);

    my $blastOutFile = sprintf("blastout_%x.txt", rand(65536));
    print STDERR "Mapping primers ($unmapped_primer oligos) to the template using blast ...\n";
    # next line modified by ET 2016/11/23
    my $cmd = $BLAST_DIR . "blastn -db $templateFile -num_threads 4 -word_size 7 -evalue 10 -query $primerFile -out $blastOutFile";
    system($cmd);

    print STDERR "Parsing blast results ...\n";	
    
    my $blast_out = new Bio::SearchIO(-format => 'blast', -file => $blastOutFile);
    while(my $result = $blast_out->next_result()){
        my $primerID = $result->query_name();
        while(my $hit = $result->next_hit()){
	    my $match_seqid = $hit->name();
	    while(my $hsp = $hit->next_hsp()){
		my $strand = $hsp->strand('hit') > 0 ?  '+' : '-';
		my $match_percentage = $hsp->percent_identity();
		my ($match_query_start, $match_query_end) = $hsp->range('query');
		my ($match_target_start, $match_target_end) = $hsp->range('hit');
		my $match_len = $match_target_end - $match_target_start + 1;
		my $primerSeq  =  $ID2Seq{$primerID};
		my $homologyScore = $match_len * (1.0-(1.0-$match_percentage/100.0)*2.0);
		my $hit_string = $hsp->hit_string();
		my $query_string = $hsp->query_string();
		my $homology_string = $hsp->homology_string();
		my $hit_percentage_annealed;
		if($homology_string =~ / /){
		    $hit_percentage_annealed = $homology_string =~ /  / ||
			$hit_string =~ /-/ || 
			$query_string =~ /-/ ? #ignore sequences with more than one consecutive mismatch or an indel
			0 : oligoPercentageAnnealedWithMismatch($hit_string, $query_string, $homology_string, 200, 1.5, 50, 0.4, 0.0, $primerMinTm-5) - 0.000001;	
		}else{
		    $hit_percentage_annealed = oligoPercentageAnnealed($hit_string, 200, 1.5, 50, 0.4, 0.0, $primerMinTm-5) - 0.000001;							 
		}
		
		if($hit_percentage_annealed > 0.0001){
		    push(@{$primerMappingInfo{$primerSeq}->{'matches'}},
			 ($match_seqid . ":" . $match_target_start . ":" . $strand . ":" . $hit_percentage_annealed . ":" . $primerSeq));
		    
		}
		if(length($primerSeq) == $match_len && $match_percentage == 100){
		    if(!$primerMappingInfo{$primerSeq}->{'targetPos'}){
			$primerMappingInfo{$primerSeq}->{'targetPos'}	=
			    $match_seqid . ":" . $match_target_start . ":" . $strand;
			$primerMappingInfo{$primerSeq}->{'excluded'} = 0;
		    }
		}
	    }

    	}
    }

    print STDERR "Done!\n";
    unlink($primerFile);
    unlink($blastOutFile);

    #Save the blast results for future use
#    dumpPrimerMappingInfo($primerMappingInfoFile, \%primerMappingInfo);			

    return %primerMappingInfo;
}

sub readPrimerMappingInfo(){
    my $fileName = shift;
    my $h_primerMappingInfo = shift;
    if(open(INFILE, "$fileName")){ 
	while(my $line = <INFILE>){
	    chop($line);
	    my ($primer, $targetPos, $match) = split(/\t/, $line);
	    next if($primer !~ /[a-zA-Z]/);	
	    next if(!$$h_primerMappingInfo{$primer}->{'tag'});
	    $$h_primerMappingInfo{$primer}->{'targetPos'} = $targetPos;	
	    my @matches = split(/,/, $match);
	    foreach $match (@matches){
		push(@{$$h_primerMappingInfo{$primer}->{'matches'}}, $match);				
	    }
	}
	close(INFILE);
    }else{
	print STDERR "$fileName not found. Skipped!\n";
    }	
}

sub dumpPrimerMappingInfo(){
    my $fileName = shift;
    my $h_primerMappingInfo = shift;
    open(OUTFILE, ">$fileName") || die("Error in opening file $fileName.\n");
    foreach my $primer (keys(%{$h_primerMappingInfo})){
	print OUTFILE $primer, "\t", $$h_primerMappingInfo{$primer}->{'targetPos'}, "\t", join(",", @{$$h_primerMappingInfo{$primer}->{'matches'}}), "\n" if($$h_primerMappingInfo{$primer}->{'matches'} && $$h_primerMappingInfo{$primer}->{'targetPos'});
    }
    close(OUTFILE);
}


=head2 primerFilter

 Title   : primerFilter
 Function: A filter to remove primer pairs that 			
 		(1) may generate more than one amplicon
 		(2) belong or may interact with the list of excluded primers
 Returns : None
 Args    : $h_primerMappingInfo 	a reference to the %primerMappingInfo hash table
           $h_ampliconInfo 		a reference to the %ampliconInfo hash table
           $numCandidatePrimer		the number of candidate primer wanted
           $sizeLimit 			the size limit of PCR
           $h_transcriptGeneMap         a reference to the %transcriptGeneMap hash table

=cut

sub primerFilter(){
    my $h_primerMappingInfo = shift;
    my $h_ampliconInfo = shift;
    my $numCandidatePrimer = shift;
    my $ampliconMinLen = shift;
    my $ampliconMaxLen = shift;
    my $h_transcriptGeneMap = shift;


    my %ampliconInfo = %{$h_ampliconInfo};
    my %primerMappingInfo = %{$h_primerMappingInfo};

    print STDERR "Filtering candidate primers ... \n";

    my $numberMoreThanOneAmplicon = 0;
    my $numberNoAmplicon = 0;
    my $numberInteractWithIncludedPrimers = 0;
    my $numberOversizeAmplicon = 0;
    my $candidatePrimerPairs = 0;
    
    # Emily added counter for number of sites with primers before and after filtering
    my $numberAmpliconsBeforeFiltering = 0;
    my $numberAmpliconsRemovedByFiltering = 0;

    print STDERR "Number of candidate primer pairs before filtering: \n";
    foreach my $amplicon (keys(%ampliconInfo)) {
	my @cLP = @{$ampliconInfo{$amplicon}->{'cLP'}};
	# Emily added counter update
	$numberAmpliconsBeforeFiltering ++;
	print STDERR "$amplicon:", scalar(@cLP), "\n";
    }

    foreach my $amplicon (keys(%ampliconInfo)) {
	my @cLP = @{$ampliconInfo{$amplicon}->{'cLP'}};
	my @cRP = @{$ampliconInfo{$amplicon}->{'cRP'}};
	my @filtered_cLP;
	my @filtered_cRP;
	my $size = scalar(@cLP);
	for (my $i=0; $i< $size; $i++) {
	    my $LP = pop(@cLP);
	    my $RP = pop(@cRP);

	    next if ($primerMappingInfo{$LP}->{'excluded'} || $primerMappingInfo{$RP}->{'excluded'} );


	    last if(scalar(@filtered_cLP) >= $numCandidatePrimer);

	    my @primerSet;
	    push(@primerSet, ($LP, $RP));

	    #remove primer pairs that do not generate one amplicon
	    my @amplicons = checkAmplification(\@primerSet, $h_primerMappingInfo, $ampliconMaxLen*2, $h_transcriptGeneMap);
	    my $relative_mass = 0;
	    my $firstAmpliconSize = 0;
	    
	    for(my $i=0; $i<scalar(@amplicons); $i++){
		my ($mass,$chr,$chr_start,$chr_end) = split(/[:]/, $amplicons[$i]);
		$firstAmpliconSize = $chr_end-$chr_start+1 if($i==0);
		$relative_mass += $mass;
	    }
	    
	    if($relative_mass > $specificity_threshold){
		$numberMoreThanOneAmplicon ++;
		next;
	    }

	    if ($firstAmpliconSize == 0) {
		$numberNoAmplicon ++;
		next;
	    } elsif($firstAmpliconSize < $ampliconMinLen || $firstAmpliconSize > $ampliconMaxLen ){
		$numberOversizeAmplicon ++;
		print STDERR "Amplicon size outside boundaries [$ampliconMinLen,$ampliconMaxLen]. The amplicon has length: $firstAmpliconSize.\n";
		next;
	    }
	    
	    #check any possible interaction between the pair of candidate primers and all included primers.
	    if(scalar(@includedPrimers)){
		my $compatible = 1;
		foreach my $includedPrimer (@includedPrimers){
		    my @primerSet;
		    push(@primerSet, ($LP, $RP, $includedPrimer));
		    my @amplicons = checkAmplification(\@primerSet, $h_primerMappingInfo, $ampliconMaxLen*2, $h_transcriptGeneMap);
		    my $relative_mass = 0;
		    for(my $i=0; $i<scalar(@amplicons); $i++){
			my @fields = split(/:/, $amplicons[$i]);
			$relative_mass += $fields[0];
		    }
		    if($relative_mass > $specificity_threshold){
			$compatible = 0;
			$numberInteractWithIncludedPrimers ++;
			last;
		    }
		    undef(@primerSet);
		}
		next if (!$compatible);

		my @primerSet = ($LP, $RP, @includedPrimers);
		for(my $i=2; $i < scalar(@primerSet); $i++){
		    for(my $j = 0; $j < 2; $j++){
			next if($primerSet[$i] !~ /[acgtACGT]/ || $primerSet[$j] !~ /[acgtACGT]/);
			my $scores = `$CHECKPD_EXE $primerSet[$i] $primerSet[$j]`;
			my @fields = split(/[ \t]/, $scores);
			if($fields[0] > $primerSelfAny*100 || $fields[1] > $primerSelfEnd*100){
			    $compatible = 0;
			    $numberInteractWithIncludedPrimers ++;
			    last;
			}
		    }
		}
		next if (!$compatible);
	    }

	    push(@filtered_cLP, $LP);
	    push(@filtered_cRP, $RP);

	}
	delete($$h_ampliconInfo{$amplicon}->{'cLP'});
	delete($$h_ampliconInfo{$amplicon}->{'cRP'});
	push(@{$$h_ampliconInfo{$amplicon}->{'cLP'}}, @filtered_cLP);
	push(@{$$h_ampliconInfo{$amplicon}->{'cRP'}}, @filtered_cRP);
	undef(@filtered_cLP);
	undef(@filtered_cRP)
    }
    print STDERR "Number of candidate primers removed due to having having more than one amplicon: $numberMoreThanOneAmplicon.\n";
    print STDERR "Number of candidate primers removed due to having having no amplicon: $numberNoAmplicon.\n";
    print STDERR "Number of candidate primers removed due to having interaction with included primers: $numberInteractWithIncludedPrimers.\n";
    print STDERR "Number of candidate primers removed due to amplicon size: $numberOversizeAmplicon.\n";
    print STDERR "Good candidate primer pairs left after filtering:\n";
    foreach my $amplicon (sort keys(%{$h_ampliconInfo})) {
	my @cLP = @{$$h_ampliconInfo{$amplicon}->{'cLP'}};
	my @cRP = @{$$h_ampliconInfo{$amplicon}->{'cRP'}};
	if(scalar(@cLP) == 0){
	    # Emily added counter update
	    $numberAmpliconsRemovedByFiltering ++;
	    print STDERR "No good candidate primers could be found for $amplicon!\n";
	    delete($$h_ampliconInfo{$amplicon});
	}
	print STDERR $amplicon, ": ", scalar(@cLP), "\n";
    }
    # Emily added print statements about primers lost during filtering
    print STDERR "Number of amplicons with candidate primers before filtering: $numberAmpliconsBeforeFiltering.\n";
    print STDERR "Number of amplicons for which all candidate primers were filtered out: $numberAmpliconsRemovedByFiltering.\n";
    my $remaining = $numberAmpliconsBeforeFiltering - $numberAmpliconsRemovedByFiltering;
    print STDERR "Number of amplicons with candidate primers after filtering: $remaining.\n";

    print STDERR "Done!\n";
}

=head2 getMultiplexPrimerSet

 Title   : getMultiplexPrimerSet
 Function: To find a set of multiplexing primer pairs that do not interact with each other by exhaustive searching
 Returns : None
 Args    : $h_primerMappingInfo 	a reference to the %primerMappingInfo hash table
           $h_ampliconInfo 		a reference to the %ampliconInfo hash table
           $sizeLimit 			the size limit of PCR
           $h_transcriptGeneMap         a reference to the %transcriptGeneMap hash table

=cut

sub getMultiplexPrimerSet(){
    my $h_ampliconInfo = shift;
    my $h_primerMappingInfo = shift;
    my $sizeLimit = shift;
    my $h_transcriptGeneMap = shift;

    my %primerMappingInfo = %{$h_primerMappingInfo};
    my %ampliconInfo = %{$h_ampliconInfo};
    my @selectedPrimers;
    my $ampPos = 0;
    my %targetAmpliconTable;

    foreach my $ampliconID (keys(%ampliconInfo)) {
	$targetAmpliconTable{$ampliconID} = 1;
    }
    my @targetAmplicons = sort(keys(%targetAmpliconTable));

    print STDERR "Searching for a set of multiplex primers ... \n";
    while(scalar(@targetAmplicons) > 0 &&
	  !getCompatiblePrimers(\@selectedPrimers, \@targetAmplicons, 0, $h_ampliconInfo, $h_primerMappingInfo, $sizeLimit, $h_transcriptGeneMap)){

	print STDERR "No compatible primer set found for all amplicons.\n";
	print STDERR "Attempting to find a primer set for maximal number of amplicons ... \n";

	#if no good primer set can be found for all amplicons, remove the amplicon having the most bad interactions with others, and try again.
	my %interactionTable;
	foreach my $setA (keys(%primerSetCompatibleTable)) {
	    foreach my $setB (keys(%{$primerSetCompatibleTable{$setA}})) {
		next if($primerSetCompatibleTable{$setA}->{$setB} > 0 && $primerSetCompatibleTable{$setA}->{$setB} < $specificity_threshold*2);
		my @setAInfo = split(/:/, $setA);
		my @setBInfo = split(/:/, $setB);
		next if(!$targetAmpliconTable{$setAInfo[0]} || !$targetAmpliconTable{$setBInfo[0]});
		$interactionTable{$setAInfo[0]}->{$setAInfo[1]} ++;
		$interactionTable{$setBInfo[0]}->{$setBInfo[1]} ++;
	    }
	}

	#The max-min principle is used here. Basically, the score for each amplicon is the least number of bad interactions among all candidate primers
	#and the amplicon with maximal score is removed.
	#
	my $maxMinInteraction = 0;
	my $maxMinAmplicon;
	foreach my $ampliconID (keys(%interactionTable)) {
	    my $minInteraction = 1000;
	    foreach my $primerPos (keys(%{$interactionTable{$ampliconID}})) {
		$minInteraction = $interactionTable{$ampliconID}->{$primerPos}
		if($interactionTable{$ampliconID}->{$primerPos} < $minInteraction);
	    }
	    $interactionTable{$ampPos}->{'min'} = $minInteraction;
	    print STDERR "(", $ampliconID, ",", $minInteraction, ")\n";
	    if($maxMinInteraction < $minInteraction){
		$maxMinInteraction = $minInteraction;
		$maxMinAmplicon = $ampliconID;
	    }
	}
	print STDERR "Amplicon $maxMinAmplicon ejected\n";
	delete($targetAmpliconTable{$maxMinAmplicon}) if($maxMinAmplicon);
	@targetAmplicons = keys(%targetAmpliconTable);
	print STDERR "candidate amplicons: ", join(";", @targetAmplicons), "\n";
	undef(@selectedPrimers);
    }

    foreach my $item (@selectedPrimers) {
	my @fields = split(/:/, $item);
	my @cLP = @{$ampliconInfo{$fields[0]}->{'cLP'}};
	my @cRP = @{$ampliconInfo{$fields[0]}->{'cRP'}};
	$$h_ampliconInfo{$fields[0]}->{'LP'} = $cLP[$fields[1]];
	$$h_ampliconInfo{$fields[0]}->{'RP'} = $cRP[$fields[1]];
	my @primerSet = ($cLP[$fields[1]], $cRP[$fields[1]]);
	my @amplicon =  checkAmplification(\@primerSet, $h_primerMappingInfo, $sizeLimit, $h_transcriptGeneMap);
	$$h_ampliconInfo{$fields[0]}->{'amplicon'} = join(";", @amplicon);
    }
}

=head2 getCompatiblePrimers

 Title   : getCompatiblePrimers
 Function: The helper function of getMultiplexPrimerSet. This function conduct exhaustive searching by recursion.
 Returns : 1 if good primers found
 		0 if no good primers found
 Args    : $h_selectedPrimers 	        a reference to the @selectedPrimers list, containing all good primers found.
           $h_targetAmplicons 	        a reference to the @targetAmplicon list, containing all target amplicons.
           $sizeLimit 	                the size limit of PCR
           $ampPos 	                the position of the target amplicon within the @targetAmplicon list. Each time
                                        this function is called, it finds primers for this target amplicon that do not
                                        interact with good primers selected previously (for the amplicons in front of this one).
           $h_ampliconInfo 		a reference to the %ampliconInfo hash table
           $h_primerMappingInfo 	a reference to the %primerMappingInfo hash table
           $sizeLimit 			the size limit of PCR
           $h_transcriptGeneMap         a reference to the %transcriptGeneMap hash table

=cut

sub getCompatiblePrimers(){
    my $h_selectedPrimers = shift;
    my $h_targetAmplicons = shift;
    my $ampPos = shift;
    my $h_ampliconInfo = shift;
    my $h_primerMappingInfo = shift;
    my $sizeLimit = shift;
    my $h_transcriptGeneMap = shift;

    my %primerMappingInfo = %{$h_primerMappingInfo};
    my @targetAmplicons = @{$h_targetAmplicons};
    if($ampPos == scalar(@targetAmplicons)){
	return 1;
    }

    my %ampliconInfo = %{$h_ampliconInfo};

    my @cLP = @{$ampliconInfo{$targetAmplicons[$ampPos]}->{'cLP'}};
    my @cRP = @{$ampliconInfo{$targetAmplicons[$ampPos]}->{'cRP'}};

    for (my $i=0; $i < scalar(@cLP) ; $i++) {

	push(@{$h_selectedPrimers}, $targetAmplicons[$ampPos]. ":$i");
	my @selectedPrimers = @{$h_selectedPrimers};

	my $compatible = 1;

	for(my $j=0; $j < scalar(@selectedPrimers)-1; $j++){
	    if(!&checkTwoPrimerSets($selectedPrimers[$j],
				    $selectedPrimers[scalar(@selectedPrimers)-1],
				    $h_targetAmplicons,
				    $h_ampliconInfo,
				    $h_primerMappingInfo,
				    $sizeLimit, $h_transcriptGeneMap)){
		$compatible = 0;
		last;
	    }
	}

	if($compatible && &getCompatiblePrimers($h_selectedPrimers, $h_targetAmplicons, $ampPos+1, $h_ampliconInfo, $h_primerMappingInfo, $sizeLimit, $h_transcriptGeneMap) == 1){
	    return 1;
	}
	pop(@{$h_selectedPrimers});
    }
    return 0;
}

=head2 checkTwoPrimerSets

 Title   : checkTwoPrimerSets
 Function: To check whether there is bad interaction (that leads to more than two amplicons) between two primer pairs.
 Returns : 1 if the two primer pairs are compatible
 		0 if the two primer pairs are not compatible
 Args    : $setA	                the first primer set
           $setA	                the second primer set
           $h_targetAmplicons 	        a reference to the @targetAmplicon list, containing all target amplicons.
           $h_ampliconInfo 		a reference to the %ampliconInfo hash table
           $h_primerMappingInfo 	a reference to the %primerMappingInfo hash table
           $sizeLimit 			the size limit of PCR
           $h_transcriptGeneMap         a reference to the %transcriptGeneMap hash table

=cut

sub checkTwoPrimerSets(){
    my $setA = shift;
    my $setB = shift;
    my $h_targetAmplicons = shift;
    my $h_ampliconInfo = shift;
    my $h_primerMappingInfo = shift;
    my $sizeLimit = shift;
    my $h_transcriptGeneMap = shift;
    

    if (!$primerSetCompatibleTable{$setA}->{$setB}){
	my %ampliconInfo = %{$h_ampliconInfo};
	my @targetAmplicons = @{$h_targetAmplicons};
	my @setAInfo = split(/:/, $setA);

	my $ampliconA = $setAInfo[0];
	my @cLP = @{$ampliconInfo{$ampliconA}->{'cLP'}};
	my @cRP = @{$ampliconInfo{$ampliconA}->{'cRP'}};
	my @primerSet = ($cLP[$setAInfo[1]], $cRP[$setAInfo[1]]);
	if($setA ne $setB){
	    my @setBInfo = split(/:/, $setB);
	    my $ampliconB = $setBInfo[0];
	    @cLP = @{$ampliconInfo{$ampliconB}->{'cLP'}};
	    @cRP = @{$ampliconInfo{$ampliconB}->{'cRP'}};
	    push(@primerSet, ($cLP[$setBInfo[1]], $cRP[$setBInfo[1]]));
	}

	#check potential primer dimers
	for(my $i=1; $i < scalar(@primerSet); $i++){
	    for(my $j = 0; $j < $i; $j++){
		my $scores = `$CHECKPD_EXE $primerSet[$i] $primerSet[$j]`;
		my @fields = split(/[ \t]/, $scores);
		if($fields[0] > $primerSelfAny*100 || $fields[1] > $primerSelfEnd*100){
		    $primerSetCompatibleTable{$setA}->{$setB} = -1;
		    return 0;
		}
	    }
	}

	#search the number of amplicons
	my @amplicons = checkAmplification(\@primerSet, $h_primerMappingInfo, $sizeLimit, $h_transcriptGeneMap);
	my $relative_mass = 0;
	for(my $i=0; $i<scalar(@amplicons); $i++){
	    my @fields = split(/:/, $amplicons[$i]);
	    $relative_mass += $fields[0];
	}
	$primerSetCompatibleTable{$setA}->{$setB} = $relative_mass;
    }
    return ($primerSetCompatibleTable{$setA}->{$setB} < $specificity_threshold*2 
	    && $primerSetCompatibleTable{$setA}->{$setB} > 0 ) ? 1 : 0;
}

=head2 checkAmplification

 Title   : checkAmplification
 Function: To find amplicons for a set of primers
 Returns : None
 Args    : $h_Primers	            a reference to the @primers list, containing the sequence of primers
           $h_primerMappingInfo     a reference to the %primerMappingInfo hash table
           $sizeLimit 		    the size limit of PCR
           $h_transcriptGeneMap     a reference to the %transcriptGeneMap hash table

=cut

sub checkAmplification(){
    my $h_Primers = shift;
    my $h_primerMappingInfo = shift;
    my $sizeLimit = shift;
    my $h_transcriptGeneMap = shift;
    
    my @Primers = @{$h_Primers};
    my %primerMappingInfo = %{$h_primerMappingInfo};

    my @amplicons;
    my %MappingInfo;
    my %pos2Seq;
    foreach my $primer (@Primers){
	foreach my $match (@{$primerMappingInfo{$primer}->{'matches'}}) {
	    my @fields = split(/:/, $match);
	    push(@{$MappingInfo{$fields[0]}->{$fields[2]}}, $fields[1] + $fields[3]);
	    $pos2Seq{$fields[0]}{$fields[1]} = $fields[4];
	}
    }

    my %processedAmplicons;
    
    foreach my $chr (keys(%MappingInfo)){
	foreach my $forward_position (@{$MappingInfo{$chr}->{'+'}}){
	    foreach my $reverse_position (@{$MappingInfo{$chr}->{'-'}}){
		my $size = int($reverse_position) - int($forward_position)+1;
		if($size > 0 && $size < $sizeLimit){
		    my $yield = ($forward_position - int($forward_position)) * ($reverse_position - int($reverse_position));
		    if ($yield > 0.000001){
			my $FPseq = $pos2Seq{$chr}{int($forward_position)};
			my $RPseq = $pos2Seq{$chr}{int($reverse_position)};
			my $gene = $$h_transcriptGeneMap{$chr};
			my $ampliconString = sprintf("%s-%s-%d-%.4f", $FPseq, $RPseq, $size, $yield);
			my $found = 0;
			# only add amplicon if it's new (i.e. different gene, yield, length, FPseq, or RPseq)
			if (!defined $processedAmplicons{$gene}) {
			    $processedAmplicons{$gene} = [$ampliconString];
			    push(@amplicons,
				 sprintf("%6.5f:%s_%s:%d:%d:%s:%s", $yield, $gene, $chr, int($forward_position), int($reverse_position), $FPseq , $RPseq));
			} elsif (!grep {$_ eq $ampliconString} @{$processedAmplicons{$gene}}) {
			    push(@{$processedAmplicons{$gene}}, $ampliconString);
			    push(@amplicons,
				 sprintf("%6.5f:%s_%s:%d:%d:%s:%s", $yield, $gene, $chr, int($forward_position), int($reverse_position), $FPseq , $RPseq));
			}
		    }
		}
	    }
	}
    }
    # DEBUGGING PRINTS
#    foreach my $g (keys(%processedAmplicons)) {
#	print $g, " : ";
#	foreach my $amp (@{$processedAmplicons{$g}}){
#	    print $amp, ", ";
#	}
#	print "\n\n";
#    }

    return if(scalar(@amplicons) < 1);
    
    #sort all amplicons based on yields
    my @sorted_amplicons = Sort_Table(
	cols		=> '5',
	field		=> '1',
	sorting  	=> 'descending',
	structure	=> 'csv',
	separator	=> '\:',
	data		=> \@amplicons
	);

    return @sorted_amplicons;
}


=head2 shortOligoTm

 Title   : shortOligoTm
 Function: To calculate Tm for short oligonucleotides based on SantaLucia J. PNAS 95:1460-1465 (1998).
 		Tm is adjusted for MgCl2 and DMSO concentration based on von Ahsen et al. Clinical Chemistry 47: 1956-61 (2001)
 Returns : $Tm, the annealing temperature
 Args    : $seq                        the oligonucleotide sequence (should be in upper case)
           $C_primer                   the concentration of oligonucleotide in nM
           $C_Mg                       the concentration of magnesium in mM
           $C_MonovalentIon            the concentration of mono-vaent ions (Na+, Ka+) in mM
           $C_dNTP                     the concentration of dNTP in mM
           $percentage_DMSO            the concentration of DMSO in (v/v)%
           $percentage_annealed        the percentage of templates that anneal to primers

=cut

sub shortOligoTm(){
    my $seq = shift;
    my $C_primer = shift;   # nM
    my $C_Mg = shift;       # mM
    my $C_MonovalentIon = shift;    #mM
    my $C_dNTP = shift;	#mM
    my $percentage_DMSO = shift;
    my $percentage_annealed = shift; #percentage of templates that anneal to primers

    $seq =~ s/[ \t\n]+//g;
    $percentage_annealed = 50.0 if (!$percentage_annealed);
    $percentage_annealed /= 100.0;

    my $C_SodiumEquivalent = $C_MonovalentIon + 120 * sqrt($C_Mg-$C_dNTP);
    my $seqLength = length($seq);
    my $dH = $deltaH{'pm'}->{substr($seq, 0, 1)} + $deltaH{'pm'}->{substr($seq, $seqLength-1, 1)};
    my $dS = $deltaS{'pm'}->{substr($seq, 0, 1)} + $deltaS{'pm'}->{substr($seq, $seqLength-1, 1)};
    $seq = uc($seq);
    for(my $i = 0; $i < $seqLength - 1; $i ++){
	$dH += $deltaH{'pm'}->{substr($seq, $i, 2)};
	$dS += $deltaS{'pm'}->{substr($seq, $i, 2)};
	print STDERR "$seq $i ", substr($seq, $i, 2), "\n" if(!$deltaS{'pm'}->{substr($seq, $i, 2)});
    }
    $dS += 0.368 * $seqLength * log($C_SodiumEquivalent/1000.0);
    my $Tm = sprintf("%5.2f", ($dH * 1000) / ($dS + $R * (log($C_primer*(1-$percentage_annealed)/$percentage_annealed)-21.4164)) - 273.15 - 0.75*$percentage_DMSO);
    return $Tm;
}


=head2 oligoPercentageAnnealed

 Title   : oligoPercentageAnnealed
 Function: To calculate the percentage of annealed oligo at a given temperature.
 Returns : $Tm, the annealing temperature
 Args    : $seq                        the oligonucleotide sequence (should be in upper case)
           $C_primer                   the concentration of oligonucleotide in nM
           $C_Mg                       the concentration of magnesium in mM
           $C_MonovalentIon            the concentration of mono-vaent ions (Na+, Ka+) in mM
           $C_dNTP                     the concentration of dNTP in mM
           $percentage_DMSO            the concentration of DMSO in (v/v)%
           $temperature                the temperature of experiment

=cut

sub oligoPercentageAnnealed(){
    my $seq = shift;
    my $C_primer = shift;   # nM
    my $C_Mg = shift;       # mM
    my $C_MonovalentIon = shift;    #mM
    my $C_dNTP = shift;	#mM
    my $percentage_DMSO = shift;
    my $temperature = shift; #temperature
    $seq = uc($seq);

    my $C_SodiumEquivalent = $C_MonovalentIon + 120 * sqrt($C_Mg-$C_dNTP);
    my $seqLength = length($seq);
    my $dH = $deltaH{'pm'}->{substr($seq, 0, 1)} + $deltaH{'pm'}->{substr($seq, $seqLength-1, 1)};
    my $dS = $deltaS{'pm'}->{substr($seq, 0, 1)} + $deltaS{'pm'}->{substr($seq, $seqLength-1, 1)};
    for(my $i = 0; $i < $seqLength - 1; $i ++){
	$dH += $deltaH{'pm'}->{substr($seq, $i, 2)};
	$dS += $deltaS{'pm'}->{substr($seq, $i, 2)};
    }
    $dS += 0.368 * $seqLength * log($C_SodiumEquivalent/1000.0);

    my $temp = (($dH * 1000) /($temperature + 273.15 + 0.75*$percentage_DMSO)  - $dS)/$R;
    my $percentage_annealed = sprintf("%5.4f", 1.0/(exp($temp+21.4164)/$C_primer+1.0));
    return $percentage_annealed;
}


=head2 oligoPercentageAnnealedWithMismatch

 Title   : oligoPercentageAnnealedWithMismatch
 Function: To calculate the percentage of annealed oligo at a given temperature.
 Returns : $Tm, the annealing temperature
 Args    : $seq                        the oligonucleotide sequence (should be in upper case)
           $C_primer                   the concentration of oligonucleotide in nM
           $C_Mg                       the concentration of magnesium in mM
           $C_MonovalentIon            the concentration of mono-vaent ions (Na+, Ka+) in mM
           $C_dNTP                     the concentration of dNTP in mM
           $percentage_DMSO            the concentration of DMSO in (v/v)%
           $temperature                the temperature of experiment

=cut

sub oligoPercentageAnnealedWithMismatch(){
    my $probeSeq = shift;
    my $targetSeq = shift;
    my $homologySeq = shift;        
    my $C_primer = shift;   # nM
    my $C_Mg = shift;       # mM
    my $C_MonovalentIon = shift;    #mM
    my $C_dNTP = shift;	#mM
    my $percentage_DMSO = shift;
    my $temperature = shift; #temperature
    my $C_SodiumEquivalent = $C_MonovalentIon + 120 * sqrt($C_Mg-$C_dNTP);
    my $seqLength = length($probeSeq);
    $probeSeq = uc($probeSeq);
    $targetSeq = uc($targetSeq);
    # Emily added print statements for debugging
#    print "IN OLIGO PERCENTAGE ANNEALED WITH MISMATCH\n";
#    print "probe seq: $probeSeq.\n";
#    print "target seq: $targetSeq.\n";
#    print "delta s mm: $deltaS{'mm'}.\n";
#    print "delta h mm: $deltaH{'mm'}.\n";
    
    my %compBase=('A', 'T', 'T','A', 'C', 'G', 'G', 'C');
    my $compTargetSeq='';
    for(my $i=0; $i< length($targetSeq); $i++){
	$compTargetSeq = $compTargetSeq . $compBase{substr($targetSeq,$i,1)};
    }
    
    my $dH = $deltaH{'pm'}->{substr($probeSeq, 0, 1)} + $deltaH{'pm'}->{substr($probeSeq, $seqLength-1, 1)};
    my $dS = $deltaS{'pm'}->{substr($probeSeq, 0, 1)} + $deltaS{'pm'}->{substr($probeSeq, $seqLength-1, 1)};
    for(my $i = 0; $i < $seqLength - 1; $i ++){
	if(substr($homologySeq,$i+1,1) eq '|'){
	    $dH += $deltaH{'pm'}->{substr($probeSeq, $i, 2)};
	    $dS += $deltaS{'pm'}->{substr($probeSeq, $i, 2)};
	}else{
	    $dH += $deltaH{'mm'}->{substr($probeSeq, $i, 2)}->{substr($compTargetSeq, $i, 2)};
	    $dS += $deltaS{'mm'}->{substr($probeSeq, $i, 2)}->{substr($compTargetSeq, $i, 2)};
	    $i++;
	    my $seqA = substr($compTargetSeq, $i+1, 1).substr($compTargetSeq, $i, 1);
	    my $seqB = substr($probeSeq, $i+1, 1).substr($probeSeq, $i, 1);
	    $dH += $deltaH{'mm'}->{$seqA}->{$seqB};
	    if(!$deltaH{'mm'}->{$seqA}->{$seqB}){
		#print "$probeSeq\n$homologySeq\n$compTargetSeq\n$i,$seqA,$seqB\n";
	    }
	    $dS += $deltaS{'mm'}->{$seqA}->{$seqB};
	}
    }
    $dS += 0.368 * $seqLength * log($C_SodiumEquivalent/1000.0);
    
    my $temp = (($dH * 1000) /($temperature + 273.15 + 0.75*$percentage_DMSO)  - $dS)/$R;
    my $percentage_annealed = sprintf("%5.4f", 1.0/(exp($temp+21.4164)/$C_primer+1.0));
    return $percentage_annealed;
}

=head2 readNNParam

 Title   : readNNParam
 Function: To read Nearest-neighbor parameters from a file.
 Returns : none
 Args    : $NNParamFile  the name of the NN parameter file. 

=cut

sub readNNParam(){
    my $NNParamFile = shift;
    open(NNFILE, $NNParamFile) || die("Error in opening NN parameter file!");
    my $line = <NNFILE>;
    while($line = <NNFILE>){
	chop($line);
	my ($seqF, $seqR, $dH, $dS, $mismatch) = split(/[:\t]/, $line);
	if(!$mismatch){
	    $deltaH{'pm'}->{$seqF} = $dH;
	    $deltaS{'pm'}->{$seqF} = $dS;
	}else{
	    $deltaH{'mm'}->{$seqF}->{$seqR} = $dH;
	    $deltaS{'mm'}->{$seqF}->{$seqR} = $dS;			
	}		
    }
    close(NNFILE);
}

=head2 delta_G_profile

 Title   : delta_G_profile
 Function: To calculate delta_G at 37C for all 6-bp sub-sequences of an oligonucleotide sequence
 Returns : @dGs an array containing delta_G values of all sub-sequences (from 3'-end to 5'-end)
 Args    : $seq the oligonucleotide sequence (should be in upper case)
=cut

sub delta_G_profile(){
    my $seq = shift;
    my @dGs;
    $seq = uc($seq);
    my $seqLength = length($seq);
    for(my $i = $seqLength-6; $i>0; $i--){
	my $dH = $deltaH{'pm'}->{substr($seq, $i, 1)} + $deltaH{'pm'}->{substr($seq, $i+5, 1)};
	my $dS = $deltaS{'pm'}->{substr($seq, $i, 1)} + $deltaS{'pm'}->{substr($seq, $i+5, 1)};         
	for(my $j = 0; $j < 5; $j ++){
	    $dH += $deltaH{'pm'}->{substr($seq, $i+$j, 2)};
	    $dS += $deltaS{'pm'}->{substr($seq, $i+$j, 2)};
	}
	my $dG = sprintf("%12.2f", $dH*1000-$dS*310.15);
	push(@dGs, $dG);
    }
    return @dGs;
}

=head2 delta_G_profile_comp

 Title   : delta_G_profile_comp
 Function: To compare the delta_G profiles of two oligonucleotide sequences
 Returns : $score
 Args    : $h_dGs_A         the handle of the first delta_G profile
           $h_dGs_B         the handle of the second delta_G profile
=cut

sub delta_G_profile_comp(){
    my $dGs_A_coded = shift;
    my $dGs_B_coded = shift;
    my @dGs_A = split(/,/, $dGs_A_coded);
    my @dGs_B = split(/,/, $dGs_B_coded);
    my $size_A = scalar(@dGs_A);
    my $size_B = scalar(@dGs_B);
    my $min_size = $size_A < $size_B ? $size_A : $size_B;
    my $score = 0;
    my $weight = 1.0;
    for(my $i = 0; $i< $min_size; $i++){
	my $ddG = ($dGs_A[$i] - $dGs_B[$i]);
	$score += $ddG*$ddG*$weight;
	$weight*=0.95;
    }
    $score = sprintf("%12.2f", $score);
    return $score;
}


=head2 max_dG_profile_mismatch

 Title   : max_dG_profile_mismatch
 Function: To identify the maximal differential delta_G profile score among a set of oligonucleotides
 Returns : $score
 Args    : @h_oligos   A handle of the list containing sequences of a set of oligonucleotides
 
=cut

sub max_dG_profile_mismatch(){
    my $h_oligos = shift;
    my %oligo_dG_profiles;
    foreach my $oligo (@{$h_oligos}){
	my @dGs = &delta_G_profile($oligo);
	$oligo_dG_profiles{$oligo} = join(",", @dGs);
    }
    my $max_score = 0;
    my $size = scalar(@{$h_oligos});
    for(my $i=0; $i<$size-1; $i++){
	for(my $j=$i+1; $j<$size; $j++){
	    my $score = &delta_G_profile_comp($oligo_dG_profiles{$$h_oligos[$i]}, $oligo_dG_profiles{$$h_oligos[$j]});
	    $max_score = $score if($score > $max_score);
	}
    }
    return $max_score;
}



main();
