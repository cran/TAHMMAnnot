#License:  This file is part of ChIPmix for TAG Project.

#    ChIPmix for TAG Project is free software: you can redistribute it and/ or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    ChIPmix for TAG Project is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with ChIPmix for TAG Project.  If not, see <http://www.gnu.org/licenses/>.


###############################################################################################
#Pre-process_normalisation.pl         08-10-08
#
#Aim :   To read the arguments and run the program normalisation_data.pl
#-----------------------------------------------------------------------------------------------
#Main function : - To read the arguments and run the program normalisation_data.pl.               
#-----------------------------------------------------------------------------------------------
#Arguments : The arguments are written in the command line : perl normalisation_data_Rep1.pl 2 3 fileOUT fileIN
#		 - nbpuce = number of chip in the experiment for one biological replicate (maximum 2 chips).
#		 - nblame = number of arrays in one chip (maximum 3 arrays).
#		 - fileOUT = name of the output data file which is the file with 6 columns
#                  named SONDE, PUCE, LAME, DYE, TRAITEMENT, VALEUR (SONDE=identifier of the probe).
#		 - fileIN = name of the input data files, they correspond to the pair files provided 
#		            by NimbleGen: one file by chip, by array and by treatment for one biological replicate.
#-----------------------------------------------------------------------------------------------
#Output :        
#------------------------------------------------------------------------------------------------
#Dependencies : 
#------------------------------------------------------------------------------------------------
#Authors : Caroline Bérard: caroline.berard@agroparistech.fr
#          Marie-Laure Martin-Magniette: marie_laure.martin@agroparistech.fr
#----------------------------------------------------------------------------------------------
#
#Date : 08/10/2008
#-----------------------------------------------------------------------------------------------
#License : GPL
#          (c) 2008 Institut National de la Recherche Agronomique
#-----------------------------------------------------------------------------------------------
#Warning : 
################################################################################################

#!/usr/bin/perl -w


$nbpuce = $ARGV[0];

$nblame = $ARGV[1];

$ficout = $ARGV[2];
$path = $ARGV[3];

my $fic1 = $ARGV[4];
my $fic2 = $ARGV[5];
my $fic3 = $ARGV[6];
my $fic4 = $ARGV[7];
my $fic5 = $ARGV[8];
my $fic6 = $ARGV[9];

my $fic7 = $ARGV[10];
my $fic8 = $ARGV[11];
my $fic9 = $ARGV[12];
my $fic10 = $ARGV[13];
my $fic11 = $ARGV[14];
my $fic12 = $ARGV[15];



my $arguments = "$fic1 $fic2 $fic3 $fic4 $fic5 $fic6 $fic7 $fic8 $fic9 $fic10 $fic11 $fic12";

system("perl $path/exec/normalisation_data.pl $nbpuce $nblame $ficout $arguments");
