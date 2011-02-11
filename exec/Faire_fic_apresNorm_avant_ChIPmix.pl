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
#Faire_fic_apresNorm_avant_ChIPmix.pl         08-10-08
#
#Aim :   To transform the data after normalisation into the input data file for ChIPmix.R.
#-----------------------------------------------------------------------------------------------
#Main function : - Create files with a special form (3 columns) which will use to do the ChIPmix.R program.                  
#-----------------------------------------------------------------------------------------------
#Arguments : The arguments are written in the command line : perl Faire_fic_apresNorm_avant_ChIPmix.pl fileIN MoyenneDyeswap
#	       - fileIN = name of the input data file, it corresponds to the normalized data.
#	       - MoyenneDyeswap = average over the 2 chips of a dye-swap. If TRUE, intensities 
#                                  are averaged on the dye-swap. Otherwise, we are working by chip separately.  
#                                    (default value = TRUE)
#-----------------------------------------------------------------------------------------------
#Output :      - One file by chromosome with 3 columns named ID, INPUT et IP (ID=identifier of the probe,
#		 INPUT = log2(INPUT), IP = log2(IP)). 
#		 If MoyenneDyeswap = FALSE, we have one file by chip and by chromosome.
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

$fic1 = $ARGV[0];

$dye = $ARGV[1];


open (FIC2, "Informations.txt"); 
while(<FIC2>)
{
	@tab = split (/\s+/,$_); #Separe tous les mots de la ligne et les met dans un tab
	$nbpuce = $tab[0];
	$nbchr = $tab[1];
	for($i=0;$i<$nbchr;$i++)
	{
		$nomchr[$i] = $tab[$i+2];
		if($tab[$i+2] =~ /[0-9]/)
		{
			$chr[$i] = '0'.$tab[$i+2];
		}
		else
		{
			$chr[$i] = $tab[$i+2];
		}
	}
}


open (FIC1, $fic1); 
@name = split (/.txt/,$fic1); #Separe tous les mots de la ligne et les met dans un tab

#je fais comme si y avait au max 2 puces ??
my %P1IP;
my %P1INPUT;
my %P2IP;
my %P2INPUT;

while(<FIC1>)
{
	@tab = split (/\s+/,$_); #Separe tous les mots de la ligne et les met dans un tab
	if ($tab[0] =~ /CHR/) #eviter la ligne d'entête
	{
		@tab2 = split (/\,/,$tab[0]);
		@tab3 = split(/\"/,$tab2[1]);
		$ID = $tab3[0];
		if($tab[4] =~ /Int2/)
		{
			if ($tab[1] =~ /1/)
			{
				$P1IP{"$ID"} = $tab[5];
			}
			else
			{
				$P2IP{"$ID"} = $tab[5];
			}
		}
		else
		{
			if ($tab[1] =~ /1/)
			{
				$P1INPUT{"$ID"} = $tab[5];
			}
			else
			{
				$P2INPUT{"$ID"} = $tab[5];
			}
		}
	}
}


if ($dye =~ /T/)
{
	for($c=0;$c<$nbchr;$c++)
	{
		$fic = "FIC_chr.$c";
		$name = "MoyDye_CHR".$chr[$c].".txt";
		open ($fic, "> $name");
		print $fic ("ID\t"."IS1\t"."IS2\n");  #Int1=INPUT
	}

	foreach $val (keys (%P1IP)){

		$moyIP = ($P1IP{"$val"} + $P2IP{"$val"}) /2;
		$moyINPUT = ($P1INPUT{"$val"} + $P2INPUT{"$val"}) /2;


		for($c=0;$c<$nbchr;$c++)
		{
			$num = 'CHR'.$chr[$c];
			if($val =~ /$num/)
			{
				$fic = "FIC_chr.$c";
				print $fic ($val."\t".$moyINPUT."\t".$moyIP."\n");
			}
		}

	}

}
else
{
	for($p=1;$p<=$nbpuce;$p++)
	{
		for($c=0;$c<$nbchr;$c++)
		{
			$fic = "FIC_puce.$p._chr.$c";
			$name = "Puce".$p."_CHR".$chr[$c].".txt";
			open ($fic, "> $name");
			print $fic ("ID\t"."IS1\t"."IS2\n");
		}
	}

	
	foreach $val (keys (%P1IP)){
	
		for($p=1;$p<=$nbpuce;$p++)
		{
			for($c=0;$c<=$nbchr;$c++)
			{
				$num = 'CHR'.$chr[$c];
				if($val =~ /$num/)
				{
					$fic = "FIC_puce.$p._chr.$c";
					if($p == 1)
					{
						print $fic ($val."\t".$P1INPUT{"$val"}."\t".$P1IP{"$val"}."\n");
					}
					elsif($p == 2)
					{
						print $fic ($val."\t".$P2INPUT{"$val"}."\t".$P2IP{"$val"}."\n");
					}
				}
			}
		}

	}
}


