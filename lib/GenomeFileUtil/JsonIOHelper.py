import os
import json
from DataFileUtil.DataFileUtilClient import DataFileUtil


def download_genome_to_json_files(token, genome_ref, target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    file_name_to_data_map = {}
    dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=token, service_ver='dev')
    genome_data = dfu.get_objects({'object_refs': [genome_ref]})['data'][0]
    genome_obj = genome_data['data']
    genome_meta = genome_data['info'][10]
    file_name_to_data_map["genome.json"] = genome_obj
    file_name_to_data_map["genome.meta.json"] = genome_meta
    if 'genbank_handle_ref' in genome_obj:
        gbk_file_name = "genome.gbk"
        dfu.shock_to_file({'handle_id': genome_obj['genbank_handle_ref'], 
                           'file_path': os.path.join(target_dir, gbk_file_name)})
        genome_obj['genbank_handle_ref'] = gbk_file_name
    if 'contigset_ref' in genome_obj:
        contigset_data = dfu.get_objects({'object_refs': [genome_obj['contigset_ref']
                                                          ]})['data'][0]
        contigset_obj = contigset_data['data']
        contigset_meta = contigset_data['info'][10]
        file_name_to_data_map["contigset.json"] = contigset_obj
        file_name_to_data_map["contigset.meta.json"] = contigset_meta
        genome_obj['contigset_ref'] = "contigset.json"
    elif 'assembly_ref' in genome_obj:
        assembly_data = dfu.get_objects({'object_refs': [genome_obj['assembly_ref']
                                                         ]})['data'][0]
        assembly_obj = assembly_data['data']
        assembly_meta = assembly_data['info'][10]
        file_name_to_data_map["assembly.json"] = assembly_obj
        file_name_to_data_map["assembly.meta.json"] = assembly_meta
        genome_obj['assembly_ref'] = "assembly.json"
        fasta_handle_ref = assembly_obj['fasta_handle_ref']
        fasta_file_name = "assembly.fa"
        dfu.shock_to_file({'handle_id': fasta_handle_ref, 
                           'file_path': os.path.join(target_dir, fasta_file_name)})
        assembly_obj['fasta_handle_ref'] = fasta_file_name
        assembly_obj['external_source_id'] = fasta_file_name
        if 'taxon_ref' in assembly_obj:
            taxon_obj = dfu.get_objects({'object_refs': [assembly_obj['taxon_ref']
                                                             ]})['data'][0]['data']
            file_name_to_data_map["taxon.json"] = taxon_obj
            assembly_obj['taxon_ref'] = "taxon.json"
            if 'taxon_ref' in genome_obj:
                genome_obj['taxon_ref'] = "taxon.json"
            taxon_obj['parent_taxon_ref'] = ""
    for target_file_name in file_name_to_data_map:
        with open(os.path.join(target_dir, target_file_name), 'w') as f:
            json.dump(file_name_to_data_map[target_file_name], f)


def compare_genome_json_files(dir1, dir2):
    files_to_compare = {'assembly.fa': False, 'genome.gbk': False,
                        'assembly.json': True, 'assembly.meta.json': True, 
                        'genome.json': True, 'genome.meta.json': True, 
                        'taxon.json': True}
    dif_files = {}
    for file_name in files_to_compare:
        data1 = ""
        data2 = ""
        if os.path.exists(os.path.join(dir1, file_name)):
            with open(os.path.join(dir1, file_name), 'r') as f:
                data1 = f.read()
        if os.path.exists(os.path.join(dir2, file_name)):
            with open(os.path.join(dir2, file_name), 'r') as f:
                data2 = f.read()
        if files_to_compare[file_name]:
            data1 = json.dumps(json.loads(data1), indent=4, sort_keys=True)
            data2 = json.dumps(json.loads(data2), indent=4, sort_keys=True)
        if data1 != data2:
            dif_files[file_name] = True
            lines1 = data1.split('\n')
            lines2 = data2.split('\n')
            if len(lines1) != len(lines2):
                print("Difference in [" + file_name + "]: number of lines " + 
                      str(len(lines1)) + " != " + str(len(lines2)))
            line_count = max(len(lines1), len(lines2))
            for pos in range(0, line_count):
                line1 = "" if pos >= len(lines1) else lines1[pos]
                line2 = "" if pos >= len(lines2) else lines2[pos]
                if line1 != line2:
                    print("Difference in [" + file_name + "], line " + str(pos + 1) + ":")
                    print("    1 - [" + line1 + "]")
                    print("    2 - [" + line2 + "]")
                    break
    return dif_files
