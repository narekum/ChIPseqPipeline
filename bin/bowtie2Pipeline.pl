#! /usr/bin/perl

####################################################
##  Narendra Kumar, PhD                            #
##  Epigenetics Unit                               #
##  Institute of Cancer Sciences                   #
##  University of Glasgow, UK                      #
##  narekum@gmail.com,narendra.kumar@glasgow.ac.uk #
####################################################

BEGIN {
        use Cwd qw(realpath cwd);
        use File::Basename;
        our ($fn, $dir) = fileparse(realpath($0));
}
$dir=~s/\/$//;

#use strict ;
use Getopt::Long;
use File::Copy 'move';

print '
  _____ _     _____ _____                         _____ _____ _____  ______ _      _____ _   _ ______
 / ____| |   |_   _|  __ \                       |  __ \_   _|  __ \|  ____| |    |_   _| \ | |  ____|
| |    | |__   | | | |__) |_____ ___  ___  __ _  | |__) || | | |__) | |__  | |      | | |  \| | |__
| |    | ._ \  | | |  ___/______/ __|/ _ \/ _` | |  ___/ | | |  ___/|  __| | |      | | | . ` |  __|
| |____| | | |_| |_| |          \__ \  __/ (_| | | |    _| |_| |    | |____| |____ _| |_| |\  | |____
 \_____|_| |_|_____|_|          |___/\___|\__, | |_|   |_____|_|    |______|______|_____|_| \_|______|
                                             | | 
                                             |_| Bowtie2, Fastqc, samtools

';

my $cmd=$0." ".join(" ",@ARGV); ### command line copy

my $mapdir;
my $species;
my $quals;
my $threads=2;
my $dontzip="y";
my $extra;
my $insertsize=150;

my $start_time ;
my $help;
my $time_tag=$start_time=time;

GetOptions ('h|help'=>\$help,                     # --help         : print this help
            "d|dir=s" => \$mapdir,                # --dir          : directory where sequences are kept  (required)
            "s|species=s" => \$species,           # --species      : species e.g. hg19  (required)
            "q|quals=s" => \$quals,               # --quals        : quals (required)
            "t|threads=s" => \$threads,           # --threads      : number of processors to be used (default 2)
            "e|extra=s" => \$extra,               # --extra        : extra bowtie options enclosed in quote ""
            "z|zip=s" => \$dontzip,               # --zip          : zip the fastq files after the run y=yes; n=no (default y)
            "i|insertsize=i" => \$insertsize,     # --insertsize   : insertsize, length of the reads (default 150)
           )
or &PRINTHELP($fn);


if(defined $help || !$mapdir || !$species || !$quals ) {
        &PRINTHELP($fn);
        exit;
}

my @options=( "$cmd",
           "--dir            $mapdir",
           "--species        $species",
           "--quals          $quals",
           "--threads        $threads",
           "--extra          \"$extra\"",
           "--zip            $dontzip",
           "--insertsize     $insertsize",
);

print "** Running $fn with the following options\n\n" ;
print "      $_\n" foreach @options ;
print "\n";

###############################################################################
#Set the values of following variables to the their respective locations of 
#your system
my $home="/home/naren";
my $bowtie="$home/local/builds/bowtie2-2.1.0/bowtie2";
my $index="/mnt/WDMyBook2TB/library/indexes/genomes/bowtie2-indexes/$species";
my $fastqc="$home/local/builds/FastQC/fastqc";
my $compressor="pigz";
my $samtools="samtools";
my $bamToBed="bamToBed";
my $fetchChromSizes="fetchChromSizes";
my $bedGraphToBigWig="bedGraphToBigWig";
my $contaminant_list="$home/local/builds/FastQC/Contaminants/contaminant_list.txt";
###############################################################################

my $bedToBedGraph="python $dir/bedToBedGraph.py";
my $extendbed="python $dir/extendbed.py";
my $countBowtieReads="$dir/countBowtieReads";
my $mapsuffix=".$species.sam";

################################################################################
## COPY SEQUENCES INTO SEQUENCE FOLDER
## CREATE ONE IF NOT EXIST
if ( -d "$mapdir/sequences" ) {
        print  "$mapdir/sequences data folder already exists - Not moving to sequences folder\n";
} else {
        mkdir "$mapdir/sequences" ;
        opendir(my $dh, $mapdir) || die "Can't opendir $mapdir: $!";
        my @fastqfiles = grep { /fastq.gz$|fastq$/ && -f "$mapdir/$_" } readdir($dh);
    	closedir $dh;
	foreach (@fastqfiles){
		print "moving $mapdir/$_ to $mapdir/sequences\n";
        	move ("$mapdir/$_", "$mapdir/sequences/") or die "move failed: $!";
	}
}

################################################################################
## QUALITY CHECK BY FASTQC 
## if the quality folder already exists we can skip, otherwise call FastQC

if ( -d "$mapdir/quality" ) {
        print  "$mapdir/quality data folder already exists - Not checking quality\n";
} else {
	mkdir "$mapdir/quality" ;
	system ("$fastqc --nogroup --threads $threads --contaminants $contaminant_list -f fastq $mapdir/sequences/*.fastq $mapdir/sequences/*.fastq.gz");
	system ("mv $mapdir/sequences/*fastqc* $mapdir/quality/");
}

################################################################################
## Do the alignment
## if sam directory exists - do nothing 

if ( -d "$mapdir/sam" ) {
	print  "$mapdir/sam data folder already exists - Not aligning\n";
} else {
	# if bam data folder exists - do nothing 
	if ( -d "$mapdir/bam" ) {
		print  "$mapdir/bam data folder already exists - Not aligning\n";
	} else {

		# unzip the files 
		print "unzipping the fastq files\n";
		system ("gunzip $mapdir/sequences/*.gz");
		print "unzipping - DONE\n";
		system ("$bowtie --version > $mapdir/bowtie.log");

		mkdir "$mapdir/sam";
		#open LOG, ">>$mapdir/bowtie.log";

		my @fastqfiles = READDIR("$mapdir/sequences", '.fastq$');
		foreach (@fastqfiles){
			my $outfile= $_ . $mapsuffix ;
			my $bowtiecmd="$bowtie $extra --$quals -p $threads -t -x $index -U $mapdir/sequences/$_ -S $mapdir/sequences/$outfile --un $mapdir/sequences/$_.unaligned";
			#print LOG "$bowtiecmd\n";
			system ("echo $bowtiecmd >> $mapdir/bowtie.log");
			print "***\n$bowtiecmd\n";
			system ("$bowtiecmd >> $mapdir/bowtie.log");
			move ("$mapdir/sequences/$outfile", "$mapdir/sam/");
		}

		mkdir "$mapdir/unaligned";
		system ("mv $mapdir/sequences/*.fastq.unaligned $mapdir/unaligned/");
		#system ("mv $mapdir/sequences/*.fastq.max $mapdir/unaligned/");

                if  ($dontzip eq "y" ) {
                        # compress to save space
                        system ("$compressor -r $mapdir/sequences/*");
                        system ("$compressor -r $mapdir/unaligned/*");
		}
	}
}

################################################################################
# make bams and index them

if ( -d "$mapdir/bam" ) {
	print "BAM data folder already exists - Not creating or indexing BAMs\n";
} else {
	mkdir "$mapdir/bam";
	my @samfiles=READDIR("$mapdir/sam", '.sam$');
	foreach (@samfiles){
		$basename=$_;
		$basename=~s/.sam//;
		print "$_ $basename\n";
		system ("$samtools view -bS $mapdir/sam/$_ > $mapdir/sam/$basename.bam");
		system ("$samtools sort $mapdir/sam/$basename.bam $mapdir/sam/$basename.sorted");
		system ("rm -f $mapdir/sam/$basename.bam");
		move ("$mapdir/sam/$basename.sorted.bam","$mapdir/sam/$basename.bam" );
		system ("$samtools index $mapdir/sam/$basename.bam");
		move ("$mapdir/sam/$basename.bam", "$mapdir/bam/");
		move ("$mapdir/sam/$basename.bam.bai", "$mapdir/bam/");
	}
	system ("rm -rf $mapdir/sam/");
}

################################################################################
# convert to bed files

if ( -d "$mapdir/bed" ){
	print "BED data folder already exists - Not converting to BED\n";
} else {
        mkdir "$mapdir/bed";
	my @bamfiles=READDIR("$mapdir/bam", '.bam$');
	foreach (@bamfiles){
		print "$_\n";
		my $awkcmd='BEGIN{OFS="\t"}{$4=".";$5="0";print}';
		#print "$awkcmd \n";
		system ("$bamToBed -i $mapdir/bam/$_ | awk '$awkcmd' | sort -k1,1 -k2,2n -k3,3n -k6,6 > $mapdir/bam/$_.bed");
		system ("uniq $mapdir/bam/$_.bed $mapdir/bam/$_.uniq.bed");
		move ("$mapdir/bam/$_.bed", "$mapdir/bed/");
		move ("$mapdir/bam/$_.uniq.bed", "$mapdir/bed/");
	}

}

################################################################################
# make bigWig files

if ( -d "$mapdir/bigWig" ) {
        print  "bigWig data folder already exists - Not making bigWig files\n";
} else {
	if ( ! -f "$dir/../chrmSizes/chrmSizes.$species" ){
		print "$dir/../chrmSizes/chrmSizes.$species not found\n--fetching now\n";
		system ("$fetchChromSizes $species > $dir/../chrmSizes/chrmSizes.$species");
	}	

	# extend the bed file to the insertsize
	if (-d "$mapdir/extendedbed"){
		print "Extended BED data folder already exists - not making extended beds\n";
	} else {
		mkdir "$mapdir/extendedbed";
		system ("$extendbed $mapdir/bed/ $species $insertsize");
		system ("mv $mapdir/bed/*.extended $mapdir/extendedbed/");
	}

	mkdir "$mapdir/bigWig";
	my @extended=READDIR("$mapdir/extendedbed",'.extended$');
	foreach (@extended){
		print "$_\n";
		system ("$bedToBedGraph --bed $mapdir/extendedbed/$_ --bedGraph $mapdir/extendedbed/$_.bedGraph --cores $threads --normalise --keepDupes --chrmSorted");
		system ("$bedGraphToBigWig $mapdir/extendedbed/$_.bedGraph $dir/../chrmSizes/chrmSizes.$species $mapdir/extendedbed/$_.bigWig");
		move ("$mapdir/extendedbed/$_.bigWig","$mapdir/bigWig");
	}

	system ("rm -rf $mapdir/extendedbed");
	
}

system ("$countBowtieReads $mapdir");

sub READDIR {
	my ($dir,$pattern) = @_;
	opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
        my @files = grep { /$pattern/ && -f "$dir/$_" } readdir($dh);
        closedir $dh;
	return (@files);

}

sub PRINTHELP {
        my ($fn)=@_;

        print <<DES;
Usage: $fn < --dir MAPDIRECTORY --species SPECIES --quals QUALS > [Options]

Options:
  -h, --help       Show this help message.
  -d, --dir        directory where sequences are kept  (required)
  -s, --species    species e.g. hg19  (required)
  -q, --quals  	   quals (required)
                   available quals: phred33-quals
                                    phred64-quals
                                    solexa-quals
                                    solexa1.3-quals
                                    integer-quals
  -t, --threads    number of processors to be used (default 2)
  -e, --extra      extra bowtie options enclosed in quote ""
  -z, --zip        zip the fastq files after the run y=yes; n=no (default y)
  -i, --insertsize extend the reads to this length (default 150)

Examples:
\$ perl $fn --dir mapdirectory --species hg19 --quals phred33-quals
\\or\\
\$ perl $fn  -d mapdirectory -s hg19 -q phred33-quals -t 8 -z n -i 150

DES
exit;
}
