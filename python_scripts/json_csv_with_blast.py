import json
import csv
import os
import argparse

def above_two(data_records):
    index=[]
    for i in range(len(data_records)):
        modules = data_records[i]["modules"].keys()
        if len(modules) > 2:
            # print(i, len(modules))
            index.append(i)
        else:
            continue
    return index

# Reads all of the json files in the dir
def list_json_files(directory):
    json_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.json'):
                json_files.append(os.path.join(root, file))
    return json_files

def process_files(file_paths, csv_filename):
    
    # Check if the output file already exists
    file_exists = os.path.isfile(csv_filename)
    
    with open(csv_filename, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        if not file_exists:
            writer.writerow(["Bin_ID", "Contig_ID", "Gene_Coordinates", "Ctg_Number", "Ctg_Tag_Id", "Mibig_biosynthetic_gene_cluster", "Mibig_Product", "Blast_name", "Blast_genecluster","Blast_annotation","Blast_perc_coverage", "Blast_perc_ident", "Blast_blastscore", "Blast_evalue", "Blast_locus_tag" ])
        
        for file_path in file_paths:

            with open(file_path, 'r') as files:
                data = json.load(files)
                index_ = above_two(data["records"])
                
                for i in index_:
                    bin_id = data["input_file"]
                    contig_id = data["records"][i]["modules"]["antismash.modules.clusterblast"]["record_id"]
                    blast_info = data["records"][i]["modules"]["antismash.modules.clusterblast"]["general"]["results"][0]["ranking"][1][1]["pairings"]
                    mibig_entries = data["records"][i]["modules"]["antismash.modules.clusterblast"]["knowncluster"]["mibig_entries"]["1"]

                    # Dict for each ctg_tagid
                    blast_ = {}
                    for blast in blast_info:
                        ctg_name = blast[0].split("|")[4]
                        blast_["ctg_cordinates"] = blast[0].split("|")[2]
                        blast_name, blast_genecluster, blast_start, blast_end, blast_strand, blast_annotation, blast_perc_ident, blast_blastscore, blast_perc_coverage, blast_evalue, blast_locus_tag  = blast[2].values()
                        blast_[ctg_name] = {"blast_name" : blast_name, 
                                            "blast_genecluster": blast_genecluster,
                                            "blast_annotation": blast_annotation,
                                            "blast_perc_coverage":blast_perc_coverage,
                                            "blast_perc_ident": blast_perc_ident,
                                            "blast_blastscore": blast_blastscore,
                                            "blast_evalue": blast_evalue,
                                            "blast_locus_tag": blast_locus_tag,
                                            "cordinates": blast[0].split("|")[2] }

                    blast_keys = list(blast_.keys())[1:]
                    mibig_entries_keys = [ ctg_tagid for ctg_tagid in mibig_entries.keys()]
        

                    for ix, vl in  enumerate(blast_keys):
                        ctg_id = vl.split("_")[0]
                        ctg_tagid = vl
                        ctg_cordinates = blast_[vl]["cordinates"]
                        blast_name = blast_[vl]["blast_name"]
                        blast_genecluster= blast_[vl]["blast_genecluster"]
                        blast_annotation = blast_[vl]["blast_annotation"]
                        blast_perc_coverage = blast_[vl]["blast_perc_coverage"]
                        blast_perc_ident = blast_[vl]["blast_perc_ident"]
                        blast_blastscore = blast_[vl]["blast_blastscore"]
                        blast_evalue = blast_[vl]["blast_evalue"]
                        blast_locus_tag = blast_[vl]["blast_locus_tag"]
                        biosynthetic_gene_cluster_id=None 
                        mibig_entries_product=None

                        if vl not in mibig_entries_keys:
                            writer.writerow([
                                    bin_id, 
                                    contig_id,
                                    ctg_cordinates, 
                                    ctg_id, 
                                    ctg_tagid, 
                                    biosynthetic_gene_cluster_id, 
                                    mibig_entries_product, 
                                    blast_name,
                                    blast_genecluster,
                                    blast_annotation,
                                    blast_perc_coverage,
                                    blast_perc_ident,
                                    blast_blastscore,
                                    blast_evalue,
                                    blast_locus_tag
                                ])



                    for ctg_tagid, vl in mibig_entries.items():
                        ctg_id = ctg_tagid.split("_")[0]
                        mibig_entry = vl[0]
                        biosynthetic_gene_cluster_id = mibig_entry[2]
                        mibig_entries_product = mibig_entry[4]
                        blast_info_ = blast_.get(ctg_tagid)

                        if blast_info_:
                            ctg_cordinates = blast_info_["cordinates"]
                            blast_name = blast_info_["blast_name"]
                            blast_genecluster= blast_info_["blast_genecluster"]
                            blast_annotation = blast_info_["blast_annotation"]
                            blast_perc_coverage = blast_info_["blast_perc_coverage"]
                            blast_perc_ident = blast_info_["blast_perc_ident"]
                            blast_blastscore = blast_info_["blast_blastscore"]
                            blast_evalue = blast_info_["blast_evalue"]
                            blast_locus_tag = blast_info_["blast_locus_tag"]

                        else:
                            ctg_cordinates = None
                            blast_name = None
                            blast_genecluster= None
                            blast_annotation = None
                            blast_perc_coverage = None
                            blast_perc_ident = None
                            blast_blastscore = None
                            blast_evalue = None
                            blast_locus_tag = None
                            
                        writer.writerow([
                            bin_id, 
                            contig_id,
                            ctg_cordinates, 
                            ctg_id, 
                            ctg_tagid, 
                            biosynthetic_gene_cluster_id, 
                            mibig_entries_product, 
                            blast_name,
                            blast_genecluster,
                            blast_annotation,
                            blast_perc_coverage,
                            blast_perc_ident,
                            blast_blastscore,
                            blast_evalue,
                            blast_locus_tag
                        ])

def main():
    parser = argparse.ArgumentParser(description="Process JSON files and output combined results to a CSV file. Please note dont run with the same output filename it will add the same results in the file.")
    parser.add_argument('--result_dir',"-r", help='JSON files to process')
    parser.add_argument('-o', '--output', default='combined_output.csv', help='Output CSV file name (default: combined_output.csv)')
    
    args = parser.parse_args()
    
    list_jsons = list_json_files(args.result_dir)

    process_files(list_jsons, args.output)
    print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main()
