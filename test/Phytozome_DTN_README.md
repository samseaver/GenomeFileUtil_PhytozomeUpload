# Loading Phytozome Reference Genomes via GFU

## File Types (using Arabidopsis thaliana TAIR10 as an example)

* Required Files:
    * Assembly in FASTA format i.e. Athaliana_167_TAIR9.fa.gz
    * Gene models in GFF3 format i.e. Athaliana_167_TAIR10.gene.gff3.gz
* Optional Files:
    * Functions in defline file (optional) i.e. Athaliana_167_TAIR10.defline.txt
* Gene names: synonym file i.e. Athaliana_167_TAIR10.synonym.txt
* Predicted annotations: annotation_info i.e. Athaliana_167_TAIR10.annotation_info.txt

/homes/chicago/seaver/genomes/Phytozome/PhytozomeV13/Athaliana/TAIR10/annotation
/homes/chicago/seaver/genomes/Phytozome/PhytozomeV13/Athaliana/TAIR10/annotation/Athaliana_167_TAIR10.annotation_info.txt
/homes/chicago/seaver/genomes/Phytozome/PhytozomeV13/Athaliana/TAIR10/annotation/Athaliana_167_TAIR10.defline.txt
/homes/chicago/seaver/genomes/Phytozome/PhytozomeV13/Athaliana/TAIR10/annotation/Athaliana_167_TAIR10.synonym.txt
/homes/chicago/seaver/genomes/Phytozome/PhytozomeV13/Athaliana/TAIR10/annotation/Athaliana_167_TAIR10.gene.gff3.gz




You need to mount the phytozome files, as transferred to the dtn, as a separate volume in run_test.sh:
-v /homes/chicago/seaver/genomes:/kb/module/genomes
