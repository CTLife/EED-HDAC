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

        Step 11: Call peak and Compute reads density by using DANPOS2.
        Usage:  
               perl  CISDA11.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]   
        For instance: 
               perl  CISDA11.pl    -in 6-FinalBAM/4_Subread     -out 11-DANPOS2/4_Subread   >> CISDA11.runLog  2>&1
     
        -------------------------------------------------------------------------------------------------------------------
        Optional arguments:
        -v, --version        Show version number of this program and exit.

        -h, --help           Show this help message and exit.

        Required arguments:
        -in  inputDir       inputDir is the name of your input folder that contains your BAM files, the suffix of the BAM files must be ".bam".  (no default)

        -out outDir         outDir is the name of your output folder that contains running results of this step.  (no default)
        ------------------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines such as RASDA, MESDA and HISDA,
        please see the web site: https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ He lab, yongp@outlook.com, Academy for Advanced Interdisciplinary Studies 
        and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.    
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------  
';

## Version Information
my $version_g = "  The Eleventh Step of CISDA (ChIP-Seq Data Analyzer), version 0.7.2, 2016-04-17.";

## Keys and Values
if ($#ARGV   == -1) { print  "\n$HELP_g\n\n";  exit 0; }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0) { @ARGV = (@ARGV, "help");         }       ## when the number of command argumants is odd. 
my %args = @ARGV;

## Initialize  Variables
my $input_g   = '6-FinalBAM/4_Subread';      ## This is only an initialization  value or suggesting value, not default value.
my $output_g  = '11-DANPOS2/4_Subread';        ## This is only an initialization  value or suggesting value, not default value.

## Available Arguments
my $available = "  -version    -help    -in    -out    ";
my $boole_g = 0;
while( my ($key, $value) = each %args ) {
    if($available !~ m/\s$key\s/) {print  "    Cann't recognize $key !!\n";  $boole_g = 1; }
}
if($boole_g == 1) {
    print "\n    The Command Line Arguments are wrong!\n";
    print   '    Please see help message by using "perl  CISDA11.pl  -help". ';
    print "\n\n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version_g\n";   exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP_g\n";      exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in' };       }else{say   "\n -in  is required.\n";          say  "\n$HELP_g\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'};       }else{say   "\n -out is required.\n";          say  "\n$HELP_g\n";    exit 0; }

## Conditions
$input_g   =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";
$output_g  =~ m/^\S+$/   ||  die   "\n$HELP_g\n\n";

## Print Command Arguments to Standard Output
print  "\n\n
        ################ Arguments ###############################
                Input  folder:  $input_g
                Output folder:  $output_g
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
my $output2_g  = "$output_g/QC_Results";
my $Dtriple_g  = "$output_g/1_Dtriple";
my $Stat_g     = "$output_g/2_Stat";
my $WIQ_g      = "$output_g/3_WIQ";
my $Profile_g  = "$output_g/4_Profile";
&myMakeDir($output2_g);
&myMakeDir($Dtriple_g);
&myMakeDir($Stat_g);
&myMakeDir($WIQ_g);    
&myMakeDir($Profile_g);   
opendir(my $DH_input_g, $input_g)  ||  die;     
my @inputFiles_g = readdir($DH_input_g);
my $pattern_g  = "[-.0-9A-Za-z]+";
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
&printVersion("danpos.py           -h");
&printVersion("danpos.py  dtriple  -h");
&printVersion("danpos.py  stat     -h");
&printVersion("danpos.py  wiq      -h");
&printVersion("danpos.py  profile  -h");
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
my $num1_g = $#BAMfiles_g + 1;
say seqFiles_FH   "\nThere are $num1_g BAM files.\n";
say         "\t\t\t\tThere are $num1_g BAM files.\n";
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $NGSinput_g = "";
my $bool_g = 0;
for(my $i=0; $i<=$#inputFiles_g; $i++) {
    next unless $inputFiles_g[$i] =~ m/\.bam$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    next unless $inputFiles_g[$i] !~ m/^removed_/;  
    if($inputFiles_g[$i] =~ m/^\d+_Input/i) {
        $NGSinput_g = "$input_g/$inputFiles_g[$i]";  
        $bool_g++;
    }  
}
if ($bool_g != 1) { print  "#$bool_g#\n";  die;}
open(runFH, ">", "$output2_g/runLog.txt")  or die;
print   runFH   "\nNGS input file: $NGSinput_g\n\n";
print   "        \nNGS input file: $NGSinput_g\n\n";
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $filesDANPOS = "";
for(my $i=0; $i<=$#BAMfiles_g; $i++) {
    $filesDANPOS = $filesDANPOS."$input_g/$BAMfiles_g[$i],";
}
$filesDANPOS =~ s/\,$//  or  die;
system("danpos.py  dtriple   $filesDANPOS   --paired 0   --pheight 1e-3  --height 0   --testcut 0   --out $Dtriple_g      --save 1   --bg $NGSinput_g      --edge 1    --count 10000000    --span 10  --clonalcut 1e-10    --frsz 200   --extend 100     >>$Dtriple_g/runLog.txt   2>&1");                                  
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n";





## END
