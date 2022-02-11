#!/usr/bin/env python

import shutil
import glob
import gzip
import os

genomes_dir = "/homes/chicago/seaver/genomes/Phytozome/PhytozomeV13"

search_path = os.path.join(genomes_dir,"*","*","annotation","*.gene.gff3.gz")
mislinked_genomes_list = list()
for gff_file in glob.glob(search_path):
    if("Csativus" in gff_file):
        continue

    gff_line_list=list()
    gene_ids=list()
    mislinked_genome = False
    with gzip.open(gff_file,'rt') as gzfh:
        for current_line in gzfh.readlines():
            current_line=current_line.strip('\r\n')

            if current_line.isspace() or current_line == "" or current_line.startswith("#"):
                gff_line_list.append(current_line)
                continue

            # Split line
            try:
                (contig_id, source_id, feature_type, start, end,
                 score, strand, phase, attributes) = current_line.split('\t')
            except ValueError:
                raise ValueError(f"unable to parse {current_line}")

            # Populating with attribute key-value pair
            # This is where the feature id is from
            attribute_dict = dict()
            for attribute in attributes.split(";"):
                attribute = attribute.strip()

                # Sometimes empty string
                if not attribute:
                    continue

                # Use of 1 to limit split as '=' character can also be made available later
                # Sometimes lack of "=", assume spaces instead
                if "=" in attribute:
                    key, value = attribute.split("=", 1)

                elif " " in attribute:
                    key, value = attribute.split(" ", 1)

                else:
                    logging.debug(f'Unable to parse {attribute}')
                    continue

                attribute_dict[key]=value

            if(feature_type == "gene"):
                gene_ids.append(attribute_dict['Name'])

            if(feature_type == 'mRNA'):
                #
                # Because of how Phytozome organizes IDs and Names in the GFF file
                # We need to preserve the IDs for maintaining links between the
                # various feature/mrna/cds dictionaries when running GFU *but*
                # The actual final identifiers used in the genome object are the
                # 'Names'.
                #
                # However, if the 'Name' matches between genes and mrnas, as it can
                # happen, then the wrong 'parent' is identified in the GFU code here:
    # https://github.com/kbaseapps/GenomeFileUtil/blob/master/lib/GenomeFileUtil/core/FastaGFFToGenome.py#L763-L767
                # So, here we try to find if the Name matches, and we append '.mRNA'
                # so to allow the GFU code to find the right parent
                #
                if(attribute_dict['Name'] in gene_ids):
                    mislinked_genome=True
                    mislinked_genomes_list.append(gff_file)
                    
                    # Edit mRNA Name
                    attribute_dict['Name']=attribute_dict['Name']+".mRNA"

                    # Rebuild attributes
                    attributes = ";".join(f"{key}={val}" for key, val in attribute_dict.items())

            new_line = '\t'.join([contig_id, source_id, feature_type, start, end,
                                  score, strand, phase, attributes])
            gff_line_list.append(new_line)
    gzfh.close()

    if(mislinked_genome is True):
        print("Moving and rewriting gff_file :"+gff_file)
 
        shutil.move(gff_file,gff_file+".Original")

        with gzip.open(gff_file,'wt') as gzfh:        
            gzfh.write("\n".join(gff_line_list))
        gzfh.close()
