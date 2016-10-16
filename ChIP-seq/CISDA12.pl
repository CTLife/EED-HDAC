#!/usr/bin/env perl5
use  strict;
use  warnings;
use  v5.18;
## perl5 version >= 5.18,   you can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## Help Infromation
my $HELP_g = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use CISDA (ChIP-Seq Data Analyzer), version 0.7.2, 2016-04-17.      
        CISDA is a Pipeline for Single-end and Paired-end ChIP-Seq Data Analysis by Integrating Lots of Softwares.

        Step 6: Convert BAM to BIGWIG by using deepTools.
        Usage:  
               perl  CISDA8.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]   [-fragLen n]
        For instance: 
               perl  CISDA8.pl    -in 6-FinalBAM/4_Subread     -out 9-deepTools/4_Subread    -fragLen 200  >> CISDA8.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -in  inputDir       inputDir is the name of your input folder that contains your BAM files, the suffix of the BAM files must be ".bam".  (no default)
        -out outDir         outDir is the name of your output folder that contains running results (Bigwig format) of this step.  (no default)
        -fragLen n          n is fragment length for single-end reads.  (no default)        
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';

## Version Information
my $version_g = "  The Sixth Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.2, 2016-04-17.";

## Keys and Values
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "help");         }       ## when the number of command argumants is odd. 
my %args = @ARGV;

## Initialize  Variables
my $input_g   = '6-FinalBAM/4_Subread';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '9-deepTools/4_Subread';     ## This is only an initialization  value or suggesting value, not default value.
my $fragLen_g = 200;                         ## This is only an initialization  value or suggesting value, not default value.

## Available Arguments
my $available = "  -version    -help    -in    -out     -fragLen     ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA8.pl  -help". ';
    print "\n\n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version_g\n";   exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP_g\n";      exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in' };       }else{say   "\n -in  is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'};       }else{say   "\n -out is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-fragLen' }   )     { $fragLen_g= $args{'-fragLen'};   }else{say   "\n -fragLen is required.\n\n";    say  "\n$HELP_g\n\n";  exit 0; }      

## Conditions
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$fragLen_g =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";

## Print Command Arguments to Standard Output
print  "\n\n
        ################ Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
                Fragment length for single-end reads:  $fragLen_g
        ###############################################################  
\n\n";
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say    "\n\n\n\n\n\n##################################################################################################";
say    "Running......";
sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}
my $output2_g = "$output_g/QC_Results";
&myMakeDir($output_g);
&myMakeDir($output2_g);
opendir(my $DH_input_g, $input_g)  ||  die;     
my @inputFiles_g = readdir($DH_input_g);
my $pattern_g  = "[-.0-9A-Za-z]+";
my $numCores_g = 4;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step......" ;
sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output2_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output2_g/VersionsOfSoftwares.txt   2>&1");
}
&printVersion("bamCoverage            --version");
&printVersion("bamCompare             --version");
&printVersion("bamPEFragmentSize      --version");
&printVersion("bigwigCompare          --version");
&printVersion("computeGCBias          --version");
&printVersion("computeMatrix          --version");
&printVersion("correctGCBias          --version");
&printVersion("estimateScaleFactor    --version");
&printVersion("multiBamSummary        --version");
&printVersion("multiBigwigSummary     --version");
&printVersion("plotCorrelation        --version");
&printVersion("plotCoverage           --version");
&printVersion("plotFingerprint        --version");
&printVersion("plotHeatmap            --version");
&printVersion("plotPCA                --version");
&printVersion("plotProfile            --version");
###################################################################################################################################################################################################





###################################################################################################################################################################################################                 
{
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the input file names ......";
my @groupFiles = ();     
my $fileNameBool = 1;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {   
        next unless $inputFiles_g[$i] =~ m/\.bam$/;  
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] !~ m/^QC_Results$/;
        next unless $inputFiles_g[$i] !~ m/^unpaired/;
        next unless $inputFiles_g[$i] !~ m/^removed_/;
        say   "\t......$inputFiles_g[$i]" ; 
        my $temp = $inputFiles_g[$i]; 
        $groupFiles[++$#groupFiles] = $inputFiles_g[$i];  
        $temp =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/   or  die   "wrong-1: ## $temp ##";
        $temp =~ m/_(Rep[1-9])\.bam$/  or    die   "wrong-2: ## $temp ##";
        if($temp !~ m/^((\d+)_($pattern_g)_(Rep[1-9]))(_[1-2])?\.bam$/) {
             $fileNameBool = 0;
        }
}
if($fileNameBool == 1)  { say    "\n\t\tAll the file names are passed.\n";  }
@groupFiles   = sort(@groupFiles);
my $numGroup  = 0;
my $noteGroup = 0;
for ( my $i=0; $i<=$#groupFiles; $i++ ) { 
    $groupFiles[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])/  or  die;
    my $n1 = $1;
    $n1>=1  or  die;
    if($noteGroup != $n1) {say "\n\t\tGroup $n1:";  $numGroup++; }
    say  "\t\t\t$groupFiles[$i]";
    $noteGroup = $n1; 
}
say  "\n\t\tThere are $numGroup groups.";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting BAM files in input folder ......";
my @BAMfiles_g = ();
open(seqFiles_FH, ">", "$output2_g/BAM-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bam$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    next unless $inputFiles_g[$i] !~ m/^removed_/;
    say    "\t......$inputFiles_g[$i]"; 
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.bam$/  or  die;  
    $BAMfiles_g[$#BAMfiles_g+1] =  $inputFiles_g[$i];
    say        "\t\t\t\tBAM file:  $inputFiles_g[$i]\n";
    say   seqFiles_FH  "BAM file:  $inputFiles_g[$i]\n";

}
say   seqFiles_FH  "\n\n\n\n\n";  
say   seqFiles_FH  "All BAM files:@BAMfiles_g\n\n\n";
say        "\t\t\t\tAll BAM files:@BAMfiles_g\n\n";
my $num1 = $#BAMfiles_g + 1;
say seqFiles_FH   "\nThere are $num1 BAM files.\n";
say         "\t\t\t\tThere are $num1 BAM files.\n";
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $BigWig1  = "$output_g/1-Extend-E7";
my $BigWig2  = "$output_g/2-Extend-1X";
my $BigWig3  = "$output_g/3-Extend-RPKM";
&myMakeDir($BigWig1);
&myMakeDir($BigWig2);
&myMakeDir($BigWig3);

my $EffectiveGenomeSize = 2150570000;   ## for mouse
my $norFact = (10**7)*$fragLen_g; 
## There are three normalization methods: 10^7, 1X, RPKM.  The second method (1X) is the most reasonable method.

say   "\n\n\n\n\n\n##################################################################################################";
say   "Convert BAM to BigWig  ......";
for (my $i=0; $i<=$#BAMfiles_g; $i++) {
    my $temp = $BAMfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$BAMfiles_g[$i]";
    system(" bamCoverage   --numberOfProcessors $numCores_g    --bam $input_g/$temp.bam     --outFileName $BigWig1/$temp.bw     --outFileFormat bigwig       --normalizeTo1x  $norFact                --extendReads $fragLen_g      --binSize 10   --ignoreDuplicates  ");                
    system(" bamCoverage   --numberOfProcessors $numCores_g    --bam $input_g/$temp.bam     --outFileName $BigWig2/$temp.bw     --outFileFormat bigwig       --normalizeTo1x  $EffectiveGenomeSize    --extendReads $fragLen_g      --binSize 10   --ignoreDuplicates  ");                    
    system(" bamCoverage   --numberOfProcessors $numCores_g    --bam $input_g/$temp.bam     --outFileName $BigWig3/$temp.bw     --outFileFormat bigwig       --normalizeUsingRPKM                     --extendReads $fragLen_g      --binSize 10   --ignoreDuplicates  ");                                             
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $PE_FragmentSize  = "$output_g/4-PE-FragmentSize";
&myMakeDir($PE_FragmentSize);
say   "\n\n\n\n\n\n##################################################################################################";
say   "Estimate Fragment Size  ......";
for (my $i=0; $i<=$#BAMfiles_g; $i++) {
    my $temp = $BAMfiles_g[$i]; 
    $temp =~ s/\.bam$//  ||  die; 
    say   "\t......$BAMfiles_g[$i]";
    system("  bamPEFragmentSize   --histogram $PE_FragmentSize/$temp.png  --numberOfProcessors $numCores_g    --plotTitle $temp    --binSize 1000  $input_g/$temp.bam ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/5-bigwigCompare";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bamfile  $input_g/$temp.bam   --ratio  log2   --binSize 100   --numberOfProcessors $numCores_g     ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/6-computeGCBias";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/7-computeMatrix";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/8-correctGCBias";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/9-estimateScaleFactor";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/10-multiBigwigSummary";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/11-plotCorrelation";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/12-plotCoverage";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/13-plotFingerprint";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################






###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/14-plotHeatmap";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/15-plotPCA";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################






###################################################################################################################################################################################################
{
my $computeMatrix  = "$output_g/16-plotProfile";   
&myMakeDir($computeMatrix);
say   "\n\n\n\n\n\n##################################################################################################";
say   "bigwigCompare  ......";
for (my $i=0; $i<$#BAMfiles_g; $i++) {
    my $temp1 = $BAMfiles_g[$i]; 
    my $temp2 = $BAMfiles_g[$i+1]; 
    $temp1 =~ s/\.bam$//  ||  die; 
    $temp2 =~ s/\.bam$//  ||  die; 
    system("bigwigCompare   --bigwig1  $output_g/1-Extend-E7/$temp1.bw   --effectiveGenomeSize $EffectiveGenomeSize   --genome mm10.2bit  --fragmentLength $fragLen_g   --numberOfProcessors $numCores_g  --GCbiasFrequenciesFile 1.txt   --plotFileFormat  svg  --biasPlot 1.svg   --regionSize 300 ");                
}
}
###################################################################################################################################################################################################



###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
