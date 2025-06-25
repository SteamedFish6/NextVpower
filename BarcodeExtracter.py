#Barcode Extracter for NextVpower, Extracting Barcodes from Nextstrain Phylogenetic Trees
#Author: Zhenyu Guo

## Usage:
## 1. Clone nextstrain data repository:
## git clone https://github.com/nextstrain/nextclade_data.git /path/of/NextClade_resource #example path
## 2. Run BarcodeExtracter.py to dump barcodes files:
## python BarcodeExtracter.py -i /path/of/NextClade_resource/nextclade_data/data/nextstrain -o /path/of/NextVpower/resource

## Note: measles tree's file name is different from other virus, the file-collecting function will detect it separately.

import os
import shutil
import json
import pandas as pd
import argparse

Default_Resource_Path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resource")
def _getExtracterArgs():
    parser = argparse.ArgumentParser(description='''.''')
    parser.add_argument('-i', "--input", type=str, help="[Dir] path of folder that stores folders of virus", required=True)
    parser.add_argument('-o', "--output", type=str, help="[Dir] path of output barcodes and reference sequences (default: ./resource)", default="resource")
    parser.add_argument("--sepchar", type=str, help="[Char] a character that separates lineages meta-information (default: |)", default="|")
    return parser.parse_args()

Use_Info_Dict = {'sars-cov-2': ['pango', 'alias', 'clade', 'clade_who'],
                 'flu': ['name', 'pango', 'clade'],
                 }
Use_Info_All = ['name', 'clade']

class TreeNode:
    # "children":[], `children` is a list contains node(s), in `tree.json` the key is 'branches'
    def __init__(self, tree_dict: dict, gene_list: list, vir_name='sars-cov-2', phylo_level=1):
        try: #self.attr_dict
            name = tree_dict["name"]
            attribute = tree_dict["node_attrs"]
        except:
            self.attr_dict = {'name':'0', 'clade':'None'}
        else:
            if vir_name == 'sars-cov-2':
                clade = attribute["clade_nextstrain"]["value"]
                alias = attribute["partiallyAliased"]["value"]
                pango = attribute["Nextclade_pango"]["value"]
                try:
                    clade_who = attribute["clade_who"]["value"]
                except:
                    clade_who = "None"
                self.attr_dict = {'name':name, 'alias':alias, 'pango':pango, 'clade':clade, 'clade_who':clade_who}
            elif vir_name == 'flu':
                clade = attribute["clade_membership"]["value"]
                try:
                    pango = attribute["proposedSubclade"]["value"]
                except:
                    pango = "unassigned"
                self.attr_dict = {'name':name, 'pango':pango, 'clade':clade}
            else:
                # alias = "None"
                # pango = "None"
                # clade = "None"
                # clade_who = "None"
                clade = attribute["clade_membership"]["value"]
                self.attr_dict = {'name':name, 'clade':clade}
        
        self.mut_nuc = tree_dict["branch_attrs"]["mutations"].get("nuc", [])
        
        mut_aa_dict = {}
        for gene in gene_list:
            mut_aa_dict[gene] = tree_dict["branch_attrs"]["mutations"].get(gene, [])
        self.mut_aa_dict = mut_aa_dict
        
        # self.branches
        branches = tree_dict.get("children", [])
        self.branches = branches
        
        self.phylo_level = phylo_level
        

class Tree():

    def __init__(self, tree_dict: dict, gff_dict: dict, vir_name='sars-cov-2'):
        self.tree = tree_dict
        self.vir_name = vir_name
        self.gff_dict = gff_dict.copy()
        self.gene_list = list(gff_dict.keys())
        self.gene_list.remove("nuc")
        self.root = TreeNode(tree_dict, self.gene_list, vir_name=self.vir_name, phylo_level=1)
    
    def getVarAnno_levelTraversal(self, info_sepr='|', use_info=None):
        gene_dict = self.gff_dict
        # nuc_dict = gene_dict.pop("nuc")
        del gene_dict["nuc"]
        # nuc_mut_list = []
        for gene in gene_dict:
            gene_dict[gene]["aa_mut"] = []
        
        # count nuc_mut, aa_mut, record them to df
        record_df = pd.DataFrame(columns=['Base Changes', 'Pos', 'Region', 'Gene ID', 'AAC', 'Lineage', 'Phylo Level'])
        record_df.index.name = 'pos'
        queue = [self.root] #queue for level traversal
        
        res = {}
        abs_mut_queue = [self.root.mut_nuc] #record "absolute path", sync with `queue`
        
        i = 0
        used_nuc_aa_mut_pos_dict = {} #key:(nuc_mut, aa_mut), value:nuc_pos
        
        while queue:
            current = queue.pop(0) #queue out
            lineage = current.attr_dict["name"]
            phylo_level = current.phylo_level
            
            # nuc_pos_mut_dict = {int(site.strip("ACGTN-")): site for site in current.mut_nuc}
            nuc_pos_mut_dict = {int(site[1:-1]): site for site in current.mut_nuc}
            
            for gene in current.mut_aa_dict:
                if "segments" in gene_dict[gene]:
                    start = gene_dict[gene]["segments"][0]["start"]
                else:
                    start = gene_dict[gene]["start"]
                region = gene_dict[gene]["type"]
                
                for aa_site in current.mut_aa_dict[gene]:
                    aa_pos = int(aa_site[1:-1])
                    nuc_codon_pos_list = [3*aa_pos+start-3, 3*aa_pos+start-2, 3*aa_pos+start-1]
                    
                    for nuc_codon_pos in nuc_codon_pos_list:
                        if nuc_codon_pos in nuc_pos_mut_dict:
                            nuc_site = nuc_pos_mut_dict.pop(nuc_codon_pos)
                            nuc_aa_tuple = (nuc_site, aa_site)
                            if nuc_aa_tuple not in used_nuc_aa_mut_pos_dict:
                                used_nuc_aa_mut_pos_dict[nuc_aa_tuple] = phylo_level
                                record_df.loc[i] = [nuc_site, nuc_codon_pos, region, gene, aa_site, lineage, phylo_level]
                                i += 1
            
            for nuc_pos in nuc_pos_mut_dict:
                nuc_site = nuc_pos_mut_dict[nuc_pos]
                nuc_aa_tuple = (nuc_site, '_')
                if nuc_aa_tuple not in used_nuc_aa_mut_pos_dict:
                    used_nuc_aa_mut_pos_dict[nuc_aa_tuple] = phylo_level
                    record_df.loc[i] = [nuc_site, nuc_pos, '-', '-', '-', lineage, phylo_level]
                    i += 1
            
            ## collect result
            abs_mut = abs_mut_queue.pop(0)
            current_abs_mut = list(set(abs_mut + current.mut_nuc))
            ### remove `use_name[:4] == "NODE"`, only include nodes with mutation site(s)
            if current.mut_nuc:
                if current.attr_dict["name"] != '0':
                    use_name = info_sepr.join([current.attr_dict[info] for info in use_info]) if use_info else current.attr_dict['name']
                    # use_name = "{}_{}".format(use_name, phylo_level) ## todo: add `phylo_level` to lineage name
                    if use_name[:4] != "NODE":
                        res[use_name] = current_abs_mut
            
            if current.branches:
                next_phylo_level = phylo_level + 1
                for branch in current.branches:
                    queue.append(TreeNode(branch, self.gene_list, vir_name=self.vir_name, phylo_level=next_phylo_level)) #queue in
                    abs_mut_queue.append(current_abs_mut)
        
        record_df.sort_values(by=["Pos", "Phylo Level"], inplace=True)
        
        return res, record_df


def collectTreeAndFastaAndGff(path: str, target_fname='tree.json', fasta_fname='reference.fasta', gff_fname='genome_annotation.gff3') -> dict:
    '''Collect tree.json and reference.fasta
    '''
    collect_dict_fname = {}
    for tup in os.walk(path):
        sub_path, foldername, fname_list = tup
        # for fname in fname_list:
        #     if fname == target_fname:
        #         collect_dict_fname[sub_path] = [os.path.join(sub_path, fname), os.path.join(sub_path, fasta_fname), os.path.join(sub_path, gff_fname)]
        if target_fname in fname_list and fasta_fname in fname_list and gff_fname in fname_list:
            collect_dict_fname[sub_path] = [os.path.join(sub_path, target_fname), os.path.join(sub_path, fasta_fname), os.path.join(sub_path, gff_fname)]
        
    return collect_dict_fname

def convertDict2DF(res_dict: dict) -> pd.DataFrame:
    all_sites = []
    for lineage in res_dict:
        all_sites = list(set(all_sites + res_dict[lineage]))
    
    # df = pd.DataFrame(columns=sorted(all_sites), index=sorted(list(res_dict.keys())), dtype=int) #cols: mut_sites; rows: lineages
    df = pd.DataFrame(columns=all_sites, index=list(res_dict.keys()), dtype=int) #cols: mut_sites; rows: lineages
    for lineage in res_dict:
        mut_site_series = pd.Series(1, index=res_dict[lineage], dtype=int) ##todo: 1 -> level
        df.loc[lineage] = mut_site_series
        
    df.fillna(0, inplace=True)
    df.index.name = "Lineages"
    # df = df.astype("int8")
    return df

# def readGff(fname: str) -> pd.DataFrame:
#     new_df = pd.DataFrame(columns=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
#     i = 0
#     with open(fname, 'r') as gff:
#         for line in gff:
#             if line[0] not in ('#', '\n', ''):
#                 line_list = line.split()
#                 new_df.loc[i] = line_list
#                 i += 1
#     new_df['start'] = new_df['start'].astype(int)
#     new_df['end'] = new_df['end'].astype(int)
#     return new_df


if __name__ == "__main__":
    params = _getExtracterArgs()
    Fpath = params.input
    Outpath = params.output
    Var_Anno_Outpath = params.output
    Info_Sepr = params.sepchar
    
    File_Collection_Dict = collectTreeAndFastaAndGff(Fpath)
    File_Collection_Dict_measles = collectTreeAndFastaAndGff(Fpath, "measles_nextclade.json", "measles_reference_N450.fasta", "measles_reference_N450.gff3")
    File_Collection_Dict_prrsv1 = collectTreeAndFastaAndGff(Fpath, "tree.json", "lelystad_orf5.fasta", "lelystad_orf5.gff3")
    File_Collection_Dict.update(File_Collection_Dict_measles)
    File_Collection_Dict.update(File_Collection_Dict_prrsv1)
    
    for Sub_Path in File_Collection_Dict:
        Tree_Fname, Fasta_Fname, Gff_Fname = File_Collection_Dict[Sub_Path]
        Rel_Path = os.path.relpath(Sub_Path, Fpath)
        DirSplit_List = Rel_Path.split(os.sep)
        Vir_Name = DirSplit_List[0]
        Save_Name = '_'.join(DirSplit_List)
        Use_Info = Use_Info_Dict.get(Vir_Name, Use_Info_All)
        Barcode_OutName = os.path.join(Outpath, "clade_barcodes_{}.csv".format(Save_Name)) #replace subfolder name, / _
        Fasta_OutName = os.path.join(Outpath, "clade_refseq_{}.fasta".format(Save_Name))
        Gff_Outname = os.path.join(Outpath, "clade_gff_{}.gff".format(Save_Name))
        Var_Anno_OutName = os.path.join(Var_Anno_Outpath, "clade_varAnno_{}.tsv".format(Save_Name))
        print(Rel_Path, end='\t')
        
        with open(Tree_Fname, 'r') as f:
            Data = json.load(f)
        Clade_Tree = Tree(Data["tree"], Data["meta"]["genome_annotations"], vir_name=Vir_Name)
        Res, Var_Anno_df = Clade_Tree.getVarAnno_levelTraversal(use_info=Use_Info, info_sepr=Info_Sepr)
        print(len(list(Res.keys())))
        
        Df = convertDict2DF(Res)
        Df.to_csv(Barcode_OutName)
        shutil.copyfile(Fasta_Fname, Fasta_OutName)
        shutil.copyfile(Gff_Fname, Gff_Outname)
        Var_Anno_df.to_csv(Var_Anno_OutName, sep='\t', index=False)
