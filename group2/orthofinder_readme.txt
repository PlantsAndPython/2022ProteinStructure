#this tutorial follows the steps laid out by https://davidemms.github.io/menu/tutorials.html

#make directory to do this in
mkdir ~/orthofinder_tutorial
cd ~/orthofinder_tutorial

#download orthofinder
wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder_glibc-2.15.tar.gz
  #glibc is tied to the version of OS that the computer is running, the HPCC uses glibc2.17 so the latest orthofinder which requires glibc2.25 won't work

#extract the package and move to the operating directory
tar xzvf OrthoFinder_glibc-2.15.tar.gz
cd OrthoFinder/

#to test if orthofinder is working properly, have orthofinder produce it's help command
 ./orthofinder -h
#if it worked, time to get the data for your 

#get peptide sequences, preferably in fasta or fa form
wget [ftp link]
#for this sample ortholog finder I selected a combination of the top 25 genomes from pdb (https://www.rcsb.org/stats/distribution-source-organism-natural) as well as a number of established model organisms as availability allowed from ensembl.com
#the taxa used for this are:
	#Anopheles gambiae
	#Apis mellifera
	#Arabidopsis thaliana
	#Azotobacter vinelandii
	#Bos taurus
	#Brassica oleracea
	#Caenorhabditis elegans
	#Danio rerio
	#Drosophila melanogaster
	#Equus caballus
	#Escherichia coli
	#Gallus gallus
	#Haloarcula marismortui
	#Homo sapiens
	#Hordeum vulgare
	#Mus musculus
	#Oryctolagus cuniculus
	#Oryza sativa
	#Ovis aries
	#Rattus norvegicus
	#Saccharomyces cerevisiae
	#Staphylococcus aureus
	#Thermosynechococcus elongatus
	#Thermus thermophilus
	#Triticum aestivum
	#Zea mays
#a note on the natural organism distribution for PDB, the organisms most represented by structures for PDB have a number of organisms that lack genomes as the structures that fill PDB are not only technically limited but also research limited, meaning for example, a lot of plant molecular biologists study arabidopsis, the first plant that is most commonly represented by PDB structures is Thaumatococcus daniellii for repeated structure conformations of it's namesake protein thaumatin.

#currently, the method for selecting taxa is by manually querying the subsetted PDB files for the source organisms and then searching the genomes. Due to the decentralized nature of genome databases with more niche genomes, it is doubtful that this orthofinder process will be easily scalable 


#now that you downloaded your genome peptide sequences, unzip them
gunzip *.gz

#extract longest script for each gene to save time and accuracy for orthofinder
for f in *fa ; do python ~/Nick/orthofinder/OrthoFinder/tools/primary_transcript.py $f ; done
for f in *fasta ; do python ~/Nick/orthofinder/OrthoFinder/tools/primary_transcript.py $f ; done

#run orthofinder
~/orthofinder/OrthoFinder/orthofinder -f primary_transcripts/
#this takes around 2 hours to run

#so the output files of a orthofinder can be viewed by
cd primary_transcripts/OrthoFinder
ls

#as a quality control, check that both the run report says that at least 80% of the genes were assigned to orthogroup in Comparative_Genomics_Statistics/Statistics_Overall.tsv
#as a second qc, check the species tree that all species are properly sorted in ./Species_Tree/SpeciesTree_rooted.txt
#you can fit it in http://etetoolkit.org/treeview/

#orthologues can be viewed in the orthologues subdirectory in each species which compares by species pairwise
Orthologues/Orthologues_Drosophila_melanogaster/Drosophila_melanogaster__v__Homo_sapiens.tsv

#for what is more relevant to the project, the orthogroups directory will have all the total orthogroups listed in a tab deliminated file 
Orthogroups/Orthogroups.tsv


