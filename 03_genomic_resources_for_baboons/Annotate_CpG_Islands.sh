# CpG Island Annotation: Emboss: http://emboss.open-bio.org/html/use/apbs06.html#GroupsAppsTableNucleiccpgislandsR6
	module load emboss; module load postgresql; cpgplot /data/tunglab/shared/genomes/panubis1/Panubis1.0.fa -outfeat panubis1.emboss.newcpgreport.gff -outfile panubis1.emboss.out -noplot #default settings 
	sed -i '1,5d' panubis1.emboss.newcpgreport.gff
	 
	awk '{print $1,$4,$5}' panubis1.emboss.newcpgreport.gff > tmp 
  sed -e 's/ /\t/g' tmp > panubis1.emboss.cpg.bed ; rm tmp