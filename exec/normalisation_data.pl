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
#normalisation_data.pl         08-10-08
#
#Aim :   To transform the data pair.txt from NimbleGen into the input data file for Nomalisation.R.               
#-----------------------------------------------------------------------------------------------
#Main function : - Create a file with a special form used to do the normalisation by ANOVA.               
#-----------------------------------------------------------------------------------------------
#Arguments : The arguments are written in the command line : perl normalisation_data_Rep1.pl 2 3 fileOUT fileIN
# 		(In the ChIPmix for TAG Project, this program is called from the Pre-process_normalisation.pl 
#		 program so the arguments are also coming directly from this program).
#		 - nbpuce = number of chip in the experiment.
#		 - nblame = number of arrays in one chip.
#		 - fileOUT = name of the output data file which is the file with 6 columns
#                  named SONDE, PUCE, LAME, DYE, TRAITEMENT, VALEUR (SONDE=identifier of the probe).
#		 - fileIN = name of the input data file, it corresponds to the pair files provided 
#		            by NimbleGen by treatment.
#-----------------------------------------------------------------------------------------------
#Output :        - 1 file: The file with 6 columns which will use by the program Normalisation.R
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

my %sonde;

$nbpuce = $ARGV[0];

$nblame = $ARGV[1];


$arg=3;

for($p=1; $p<=$nbpuce; $p++)
{
	for($l=1; $l<=$nblame; $l++)
	{
		$fic1 = $ARGV[$arg];
						
		open (FIC1, $fic1) || die("Couldn't open $fic1 : $!");
		while(<FIC1>)
		{
			@tab = split (/\s+/,$_); #Separe tous les mots de la ligne et les met dans un tab
			if ($tab[1] =~ /BLOCK1/) # On prend que les sondes 'normales'
			{
				#l'ID est le num de l'exp + le num du chr + la position
				$ID = $tab[0].','.$tab[3];
				#la 1ère case du tab de hash est le numéro de la puce
				$sonde{"$ID"}[0] = $p;
				#la 2ème case du tab de hash est le numéro de la lame
				$sonde{"$ID"}[1] = $l;
				#la 3ème case du tab de hash est le dye
				@tab2 = split (/_/,$tab[0]); 
				$sonde{"$ID"}[2] = $tab2[1];	
				#la 4ème case du tab de hash est le traitement IP ou INPUT
				$sonde{"$ID"}[3] = "Int1";
				#la 5ème case du tab de hash est la valeur du signal
				$sonde{"$ID"}[4] = $tab[8];
			}

		}
		close(FIC1);

		$fic2 = $ARGV[($arg+1)];

		open (FIC2, $fic2) || die("Couldn't open $fic2 : $!");
		while(<FIC2>)
		{
			@tab = split (/\s+/,$_); #Separe tous les mots de la ligne et les met dans un tab
			if ($tab[1] =~ /BLOCK1/) # On prend que les sondes 'normales'
			{
				#l'ID est le num de l'exp + le num du chr + la position
				$ID = $tab[0].','.$tab[3];
				#la 1ère case du tab de hash est le numéro de la puce
				$sonde{"$ID"}[0] = $p;
				#la 2ème case du tab de hash est le numéro de la lame
				$sonde{"$ID"}[1] = $l;
				#la 3ème case du tab de hash est le dye
				@tab2 = split (/_/,$tab[0]); 
				$sonde{"$ID"}[2] = $tab2[1];	
				#la 4ème case du tab de hash est le traitement IP ou INPUT
				$sonde{"$ID"}[3] = "Int2";
				#la 5ème case du tab de hash est la valeur du signal
				$sonde{"$ID"}[4] = $tab[8];
			}

		}
		close(FIC2);
		$arg = $arg + 2;

	}
}


$nbchr = 0;
my %test;
open (FIC2, $fic2) || die("Couldn't open $fic2 : $!");
while(<FIC2>)
{
	@tab = split (/\s+/,$_); #Separe tous les mots de la ligne et les met dans un tab
	if ($tab[1] =~ /BLOCK1/) # On prend que les sondes 'normales'
	{
		@tab2 = split(/\:/, $tab[2]);
		@tab3 = split('r', $tab2[0]);
		push (@chr, $tab3[1]) unless(exists $test{$tab2[0]});
		$nbchr = $nbchr + 1 unless(exists $test{$tab2[0]});
		$test{$tab2[0]}=1;
	}
}
close(FIC2);

#Le fichier Informations.txt est codé en dur, on ne peut pas lui changer son nom
open (OUT2, "> Informations.txt") || die("Couldn't create the output file : $!");

print OUT2 ($nbpuce."\t");
print OUT2 ($nbchr."\t");
for($i=0; $i<=$#chr; $i++)
{
	print OUT2 ($chr[$i]."\t");
}
print OUT2 ("\n");

#fichier que l'on va creer
$fic = $ARGV[2];
open (OUT, "> $fic") || die("Couldn't create $fic : $!");

#Affichage dans le fichier
print OUT ("ID\t"."CHIP\t"."ARRAY\t"."DYE\t"."TREATMENT\t"."VALUE\n");
foreach $val (keys (%sonde)){
	print OUT ($val."\t".$sonde{"$val"}[0]."\t".$sonde{"$val"}[1]."\t".$sonde{"$val"}[2]."\t".$sonde{"$val"}[3]."\t".$sonde{"$val"}[4]."\n");
	
}





