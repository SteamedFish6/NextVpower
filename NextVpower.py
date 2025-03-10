#Next V-Power, adapted from "Virus Phylogenetic Resolver of Wastwater-based Epidemiology (V-Power)"
#Author: Zhenyu Guo

import os
import numpy as np
import pandas as pd
import cvxpy as cp
# from scipy.optimize import minimize, LinearConstraint
import argparse

def _getVpowerArgs():
    parser = argparse.ArgumentParser(description='''Next V-Power, adapted from Virus Phylogenetic Resolver of Wastwater-based Epidemiology (V-Power).''')
    group1 = parser.add_argument_group("Demixing arguments for the demix solver")
    group1.add_argument('-i', "--input", type=str, help="[File/Dir] path of input sample table file, .vcf file or .vcf files under a folder", required=True)
    group1.add_argument('-o', "--output", type=str, help="[File] path of output table file (default: ./demix_result.tsv)", default="demix_result.tsv")
    group1.add_argument('-b', "--barcode", type=str, help="[File] specify a usher_barcodes.csv as input barcode matrix (default: ./resource/usher_barcodes.csv)", default="resource/usher_barcodes.csv")
    group1.add_argument('-l', "--maxlineages", type=int, help="[Int] maximum number of demixing lineages (default: 100)", default=100)
    
    group2 = parser.add_argument_group("Sample Processing arguments for the [--input] sample handler")
    # group2.add_argument('-v', "--vcfs", action='store_true', help="[Flag] parse *.vcf files under input folder")
    group2.add_argument('-r', "--minrate", type=float, help="[Float] filter mutation sites with mutation rate lower than setting threshold in sample vectors (default: 0.0)", default=0.0)
    group2.add_argument('-d', "--mindepth", type=float, help="[Float] filter mutation sites with depth lower than setting threshold in *.vcf files (default: 0.0)", default=0.0)
    group2.add_argument('-a', "--ann_file", type=str, help="[File] specify a var_anno.tsv as input variation annotation table (default: ./resource/var_anno.tsv)", default="resource/var_anno.tsv")
    
    group3 = parser.add_argument_group("Barcode Processing arguments for the barcode filter")
    group3.add_argument('-n', "--nsites", type=int, help="[Int] filter lineages with fewer than [Int] mutation sites (default: 20)", default=20)
    group3.add_argument('-k', "--klineages", type=int, help="[Int] retain \"key\" mutation sites present in more than [Int] lineages (default: 200)", default=200)
    group3.add_argument("--barfilter2", action='store_true', help="[Flag] use old barcode filter to handle barcode matrix")
    group3.add_argument("--merge", action='store_true', help="[--barfilter2 only][Flag] merge lineages with completely identical mutation sites in the barcode matrix")
    
    group4 = parser.add_argument_group("Middle file output arguments for middle processes")
    group4.add_argument("--ann_outpath", type=str, help="[Dir] if not None, add save the annotated *.vcf table files under a folder (optional)", default=None)
    group4.add_argument("--vcsample", type=str, help="[File] save the sample table file converted from *.vcf files (optional)", default=None)
    group4.add_argument("--fsample", type=str, help="[File] save the filtered sample table file (optional)", default=None)
    group4.add_argument("--fbarcode", type=str, help="[File] save the filtered barcode matrix file (optional)", default=None)
    group4.add_argument("--potentials", type=str, help="[File] save potential sites not recorded in barcode but present in samples (optional)", default=None)
    
    parser.add_argument('--version', action='version', version="NextVpower_v0.13")
    return parser.parse_args()


def collectFile(path: str, isfolder=False, ftype='vcf', fullpath=False) -> list:
    '''Collect files/folders under a path according to extends.
    '''
    file_list = os.listdir(path)
    collect_list = []
    for file_name in file_list:
        full_pathname = os.path.join(path, file_name)
        if isfolder == True:
            if os.path.isdir(full_pathname):
                if fullpath == False:
                    collect_list.append(file_name)
                else:
                    collect_list.append(full_pathname)
        else:
            if file_name[-(len(ftype)+1):] == '.{}'.format(ftype):
                if fullpath == False:
                    collect_list.append(file_name)
                else:
                    collect_list.append(full_pathname)
    collect_list = sorted(collect_list)
    return collect_list

def readBarcode(fname: str) -> pd.DataFrame:
    df = pd.read_csv(fname, sep=',', index_col=0)
    df = df.T
    df.index.name = "Base Changes"
    return df

def readAnnoDF(fname: str) -> pd.DataFrame:
    df = pd.read_csv(fname, sep='\t', index_col=0)
    df.index.name = "Base Changes"
    return df

def readSampleTable(fname: str) -> pd.DataFrame:
    df = pd.read_csv(fname, sep='\t', index_col=0)
    df.index.name = "Base Changes"
    return df


def convertVcf2DF(vcfname: str) -> pd.DataFrame:
    '''Read VCF file created by freebayes and convert it to VCF DataFrame.
    '''
    col_index = ["Base Changes", "Position", "Type", "REF Base", "ALT Base", "Depth", "MAF"]
    df = pd.DataFrame(columns=col_index)
    df.set_index("Base Changes", inplace=True)
    vcffile = open(vcfname, 'r')
    for line in vcffile:
        if line[0] != '#':
            line_compo = line[:-1].split('\t')
            if len(line_compo) != 0:
                pos_str = line_compo[1]
                ref = line_compo[3]
                alt_ls = line_compo[4].split(',')
                info_list = line_compo[7].split(';')
                ao_ls = info_list[1].split('=')[1].split(',')
                ro = eval(info_list[3].split('=')[1])
                #dp = eval(info_list[2].split('=')[1])
                # dp = ao + ro #dp = sum(ao_ls) + ro
                typ_ls = info_list[4].split('=')[1].split(',')
                
                for alt, ao_str, typ in zip(alt_ls, ao_ls, typ_ls):
                    ao = eval(ao_str)
                    maf = ao / (ao + ro) * 100
                    
                    if typ == 'mnp' or typ == 'complex':
                        ref_alt_len = min(len(ref), len(alt))
                        for i in range(ref_alt_len):
                            s_ref = ref[i]
                            s_alt = alt[i]
                            s_pos = str(eval(pos_str) + i)
                            s_base_change = s_ref + s_pos + s_alt
                            depth = "{}:{} {}:{}".format(s_alt, ao, s_ref, ro)
                            df.loc[s_base_change] = pd.Series({"Position": s_pos, "Type": typ, "REF Base": s_ref, "ALT Base": s_alt, "Depth": depth, "MAF": maf})
                    
                    ##fix an issue: convert mutiple INDEL to single INDEL
                    else:
                        if typ == 'del':
                            alt_len = len(alt)
                            if alt_len > 1:
                                s_alt = alt[-1]
                                s_ref = ref[alt_len-1:]
                                s_pos = str(eval(pos_str) + alt_len - 1)
                            else:
                                s_ref = ref
                                s_alt = alt
                                s_pos = pos_str
                        elif typ == 'ins':
                            ref_len = len(ref)
                            if ref_len > 1:
                                s_ref = ref[-1]
                                s_alt = alt[ref_len-1:]
                                s_pos = str(eval(pos_str) + ref_len - 1)
                            else:
                                s_ref = ref
                                s_alt = alt
                                s_pos = pos_str
                        else:
                            s_ref = ref
                            s_alt = alt
                            s_pos = pos_str
                            
                        s_base_change = s_ref + s_pos + s_alt
                        depth = "{}:{} {}:{}".format(alt, ao, ref, ro)
                        df.loc[s_base_change] = pd.Series({"Position": s_pos, "Type": typ, "REF Base": s_ref, "ALT Base": s_alt, "Depth": depth, "MAF": maf})
    df["MAF"] = df["MAF"].astype("float64") #change dtype of column in pd.DataFrame
    return df

def BackPasteAnno(vcf_df: pd.DataFrame, anno_df: pd.DataFrame, outname: str=None) -> pd.DataFrame:
    '''Back paste annotion to VCF DataFrame.
    '''
    using_index = list(set(vcf_df.index) & set(anno_df.index))
    vcf_df["Region"] = anno_df.loc[using_index, "Region"]
    vcf_df["Gene ID"] = anno_df.loc[using_index, "Gene ID"]
    vcf_df["AAC"] = anno_df.loc[using_index, "AAC"]
    vcf_df["Effect"] = anno_df.loc[using_index, "Effect"]
    vcf_df.fillna('unknown', inplace=True)
    if outname != None:
        vcf_df.to_csv(outname, sep='\t')
    return vcf_df

def FilterVcfDF(var_df: pd.DataFrame, min_depth=0.0, outname=None) -> pd.DataFrame:
    '''Filter mutation sites with low Depth in VCF DataFrame.
    '''
    if min_depth > 0:
        out_var_df = pd.DataFrame(columns=var_df.columns)
        for idx in var_df.index:
            depth_contents = var_df.loc[idx, 'Depth'].split(' ')
            depth = eval(depth_contents[0].split(':')[1]) + eval(depth_contents[1].split(':')[1])
            if depth >= min_depth:
                out_var_df.loc[idx] =var_df.loc[idx]
    else:
        out_var_df = var_df
    if outname != None:
        out_var_df.to_csv(outname, sep='\t')
    return out_var_df

def CollectSampleVar(var_df_dict: dict, outname=None) -> pd.DataFrame:
    '''Collect 'MAF' column in VCF DataFrame of all samples, merge them, and fillna.
    Output is Sample DataFrame.
    '''
    all_index = []
    for spname in var_df_dict:
        all_index = all_index + list(var_df_dict[spname].index)
    using_index = sorted(list(set(all_index)))
    sample_list = list(var_df_dict.keys())
    sp_df = pd.DataFrame(columns=sample_list, index=using_index)
    sp_df.index.name = "Base Changes"
    for spname in var_df_dict:
        var_df = var_df_dict[spname]
        sp_df[spname] = var_df['MAF']
    sp_df.fillna(0.0, inplace=True)
    sp_df = sp_df / 100
    if outname != None:
        sp_df.to_csv(outname, sep='\t')
    return sp_df


def FilterSPDF_by_rate(sp_df: pd.DataFrame, min_rate=0.0, outname=None) -> pd.DataFrame:
    '''Filter mutation sites with low mutation rate in Sample DataFrame.
    '''
    if min_rate > 0:
        out_sp_df = pd.DataFrame(columns=sp_df.columns)
        for site in sp_df.index:
            rate_contents = sp_df.loc[site]
            if rate_contents.max() >= min_rate: #确保该位点至少有一个样品rate超过阈值
                out_sp_df.loc[site] = rate_contents
    else:
        out_sp_df = sp_df
    if outname != None:
        out_sp_df.to_csv(outname, sep='\t')
    return out_sp_df

def FilterBarcode_old(barcode_df: pd.DataFrame, sp_df: pd.DataFrame, lineage_num: int, remove_duplicates=True, outname=None) -> pd.DataFrame:
    '''Filter Barcode DataFrame.
    '''
    #column: lineages, row: mutation sites
    #Remove mutation sites not recorded in sp_df.
    using_index = list(set(barcode_df.index) & set(sp_df.index))
    out_barcode_df = barcode_df.loc[using_index]
    
    #Remove all-zero columns.
    using_columns = []
    for col in out_barcode_df.columns:
        if sum(out_barcode_df[col]) != 0:
            using_columns.append(col)
    out_barcode_df = out_barcode_df[using_columns]
    
    if remove_duplicates == True: #Remove identical columns.
        new_df = out_barcode_df[using_columns[0]].to_frame()
        new_df = new_df.T #column: mutation sites, row: lineages
        new_df.set_index(pd.Index(['Group0']), inplace=True)
        group_dict = {'Group0': using_columns[0]}
        i = 1
        for lineage in using_columns[1:]:
            lineage_col = out_barcode_df[lineage]
            need_added = True
            for group_name in new_df.index: #Search identical lineages.
                align_col = new_df.loc[group_name] - lineage_col
                if align_col.max() == 0 and align_col.min() == 0:
                    need_added = False #If identical lineage exists, no need to add it to barcode_df.
                    group_dict[group_name] = "{};{}".format(group_dict[group_name], lineage) #Merge identical lineage name.
                    break
            if need_added == True:
                new_group_name = "Group{}".format(i)
                new_df.loc[new_group_name] = lineage_col
                group_dict[new_group_name] = lineage
                i += 1
        new_df.set_index(pd.Index(list(group_dict.values())), inplace=True)
        out_barcode_df = new_df.T
        out_barcode_df.index.name = "Lineages"
    
    #Sort by num of mutation sites.
    using_columns_dict = {}
    for col in out_barcode_df.columns:
        using_columns_dict[col] = sum(out_barcode_df[col])
    using_columns = sorted(using_columns_dict, key=using_columns_dict.get, reverse=True)
    
    #Reserve lineage with more mutation sites.
    if len(using_columns) > lineage_num:
        using_columns = using_columns[:lineage_num]
    out_barcode_df = out_barcode_df[using_columns]
    
    #Output.
    if outname != None:
        out_barcode_df.to_csv(outname, sep='\t')
    return out_barcode_df

def FilterBarcode(barcode_df: pd.DataFrame, sp_df: pd.DataFrame, key_lineage_num=200, min_sites=20, outname=None) -> pd.DataFrame:
    '''Filter Barcode DataFrame. This barcode filter was authored by Kun Yang and adapted by Zhenyu Guo.
    '''
    #column: lineages, row: mutation sites
    site_num = barcode_df.sum(axis=0)
    barcode_df = barcode_df.loc[:, site_num >= min_sites] #Reserve columns with sum>20
    
    lineage_num = barcode_df.sum(axis=1)
    key_sites = list(lineage_num.loc[lineage_num >= key_lineage_num].index)
    sample_sites = list(sp_df.index)
    reserved_sites = set(key_sites + sample_sites) #Take the union.
    junk_sites = set(barcode_df.index) - reserved_sites #Sites to be discarded
    junk_lineages = set() #Lineages to be discarded
    for j_site in junk_sites:
        row = barcode_df.loc[j_site]
        j_lineages = set(row.loc[row == 1].index) #Find lineages defined by junk_sites.
        junk_lineages.update(j_lineages) #Record lineages to be discarded.
    reserved_lineages = set(barcode_df.columns) - junk_lineages #Lineages to be reserved
    barcode_df = barcode_df[list(reserved_lineages)]
    
    #Reserve mutation sites recorded in sp_df.
    using_index = list(set(barcode_df.index) & set(sp_df.index))
    barcode_df = barcode_df.loc[using_index]
    
    #Output.
    if outname != None:
        barcode_df.to_csv(outname, sep='\t')
    return barcode_df

def FilterSPDF(barcode_df: pd.DataFrame, sp_df: pd.DataFrame, outname=None, outname_p=None) -> pd.DataFrame:
    '''Filter mutation sites not recorded in Barcode DataFrame, in Sample DataFrame.
    '''
    using_index = barcode_df.index
    out_sp_df = sp_df.loc[using_index]
    out_sp_df.reindex(using_index)
    
    if outname != None:
        out_sp_df.to_csv(outname, sep='\t')
    
    if outname_p != None: #Get unrecorded sites (with rate of each sample) and save them.
        potential_sites = list(set(sp_df.index) - set(using_index))
        pot_sp_df = sp_df.loc[potential_sites]
        # pot_sp_df.reindex(potential_sites)
        pot_sp_df.to_csv(outname_p, sep='\t')
    return out_sp_df


def SolveDemix(sample_df: pd.DataFrame, barcode_df: pd.DataFrame, max_lineage_num=100, min_result=0.0):
    '''Demixing via minimize solver.
    '''
    mutation_matrix = barcode_df.to_numpy()
    sp_var_matrix = sample_df.to_numpy()
    mutation_matrix = mutation_matrix[:, 0: max_lineage_num]
    lineage_num = np.size(mutation_matrix, 1)
    sample_num = np.size(sp_var_matrix, 1)
    result_all = np.zeros((lineage_num, sample_num))
    
    for sample_idx in range(sample_num):
        sample_name = sample_df.columns[sample_idx]
        x = cp.Variable(mutation_matrix.shape[1])
        cost = cp.norm(mutation_matrix @ x - sp_var_matrix[:, sample_idx], 1)
        constraints = [sum(x) == 1, x >= 0]
        prob = cp.Problem(cp.Minimize(cost), constraints)
        solver = cp.CLARABEL
        
        try:
            prob.solve(verbose=False, solver=solver)
        except cp.error.SolverError:
            print("Failed to demix sammle '{}', most likely due to low sequencing coverage.".format(sample_name))
            # raise ValueError("Failed to demix sammle '{}', most likely due to low sequencing coverage.".format(sample_name))
        else:
            result_all[:, sample_idx] = x.value
    
    result_all[result_all < min_result] = 0
    lineage_name_index = barcode_df.columns[:lineage_num]
    out_df = pd.DataFrame(result_all, index=lineage_name_index, columns=sample_df.columns)
    out_df['sum'] = out_df.sum(axis=1)
    out_df.sort_values(by='sum', ascending=False, inplace=True)
    return out_df

# def SolveDemix(sample_df: pd.DataFrame, barcode_df: pd.DataFrame, max_lineage_num=100, method='SLSQP', tol=1e-8, maxiter=50000) -> pd.DataFrame:
#     '''Demixing via minimize solver.
#     '''
#     mutation_matrix = barcode_df.to_numpy()
#     sp_var_matrix = sample_df.to_numpy()
#     mutation_matrix = mutation_matrix[:, 0: max_lineage_num]
#     lineage_num = np.size(mutation_matrix, 1)
#     sample_num = np.size(sp_var_matrix, 1)
#     result_all = np.zeros((lineage_num, sample_num))
#     for sample_idx in range(sample_num):
#         start_point = (np.array([1 for i in range(lineage_num)]).T / lineage_num)
#         arguments = (mutation_matrix, sp_var_matrix[:, sample_idx])
#         # jac: None; hess, hessp: None
#         bounds = [(0, 1) for i in range(lineage_num)]
#         Aeq = np.array([1 for i in range(lineage_num)])
#         beq = 1
#         constraints = LinearConstraint(Aeq, beq, beq)
#         options = {'maxiter': maxiter, 'disp': False}
#         res = minimize(fun=Preva, x0=start_point, args=arguments, method=method, bounds=bounds, constraints=constraints, tol=tol, options=options)
#         if res.success:
#             result_all[:, sample_idx] = res.x
#     lineage_name_index = barcode_df.columns[:lineage_num]
#     out_df = pd.DataFrame(result_all, index=lineage_name_index, columns=sample_df.columns)
#     out_df['sum'] = out_df.sum(axis=1)
#     out_df.sort_values(by='sum', ascending=False, inplace=True)
#     return out_df

# def Preva(X, *args):
#     mm, pp = args
#     return sum(abs(mm.dot(X) - pp))


###tools: Handle raw AnnoDF
# def readAnnoDF(fname: str, outname=None) -> pd.DataFrame:
#     raw_df = pd.read_csv(fname, sep='\t', usecols=['Genome Position', 'Gene Region', 'Annotation Type', 'Base changes:Virus Number', 'Protein.Position.Amino Acid Change'])
#     out_df = pd.DataFrame(columns=["Base Changes", "Region", "Gene ID", "AAC", "Effect"])
#     out_df.set_index("Base Changes", inplace=True)
#     # region_name_dict = {'5\'UTR': 'five_prime_UTR'}
#     for idx in raw_df.index:
#         base_change_rawls = raw_df.loc[idx, 'Base changes:Virus Number'].split(':')[0].split('>')
#         base_change = "{}{}{}".format(base_change_rawls[0], raw_df.loc[idx, 'Genome Position'], base_change_rawls[1])
#         gene_id = raw_df.loc[idx, 'Gene Region']
#         #region = gene_id if gene_id in ('5\'UTR', '3\'UTR') else 'CDS'
#         aac_raw = raw_df.loc[idx, 'Protein.Position.Amino Acid Change'].split('.')[-1]
#         aac_match = re.match(r'(\d+)(.*)>(.*)', aac_raw)
#         if not aac_match:
#             aac = '-'
#         else:
#             aac = "p.{}{}{}".format(aac_match.group(2), aac_match.group(1), aac_match.group(3))
#             if aac[2] == '-':
#                 aac_match = re.match(r'(\d+)-(\d+)(.*)>(.*)', aac_raw)
#                 aac = "p.{}{}_{}{}del".format(aac_match.group(3)[0], aac_match.group(1), aac_match.group(3)[-1], aac_match.group(2))
#         effect = raw_df.loc[idx, 'Annotation Type']
#         out_df.loc[base_change] = pd.Series({"Region": region, "Gene ID": gene_id, "AAC": aac, "Effect": effect})
#     out_df.to_csv(outname, sep='\t')
#     return out_df


if __name__ == "__main__":
    params = _getVpowerArgs()
    fpath = params.input
    
    if params.ann_outpath != None:
        AnnoDF = readAnnoDF(params.ann_file)
    
    # if params.vcfs:
    if os.path.isdir(fpath):
        print("Parsing vcf files...", end='')
        Var_DF_Dict = {}
        for fname in collectFile(fpath, fullpath=True, ftype='vcf'):
            VcfDF = convertVcf2DF(fname)
            spname = os.path.basename(fname).split('.', 1)[0]
            Var_DF_Dict[spname] = FilterVcfDF(VcfDF, min_depth=params.mindepth)
            if params.ann_outpath != None:
                ann_outname = os.path.join(params.ann_outpath, os.path.basename(fname).split('.')[0]+".ann.tsv")
                AnnTab = BackPasteAnno(VcfDF, AnnoDF, ann_outname)
        SP_df_raw = CollectSampleVar(Var_DF_Dict, outname=params.vcsample)
        print("")
    elif fpath[-4:] == '.vcf':
        print("Parsing vcf files...", end='')
        Var_DF_Dict = {}
        VcfDF = convertVcf2DF(fpath)
        spname = os.path.basename(fpath).split('.', 1)[0]
        Var_DF_Dict[spname] = FilterVcfDF(VcfDF, min_depth=params.mindepth)
        if params.ann_outpath != None:
            ann_outname = os.path.join(params.ann_outpath, os.path.basename(fpath).split('.')[0]+".ann.tsv")
            AnnTab = BackPasteAnno(VcfDF, AnnoDF, ann_outname)
        SP_df_raw = CollectSampleVar(Var_DF_Dict, outname=params.vcsample)
        print("")
    else:
        SP_df_raw = FilterSPDF_by_rate(readSampleTable(fpath), min_rate=params.minrate)
    
    print("Handling barcode matrix and sample matrix...", end='')
    if params.barfilter2:
        Barcode_df = FilterBarcode_old(readBarcode(params.barcode), SP_df_raw, lineage_num=params.maxlineages, remove_duplicates=params.merge, outname=params.fbarcode)
    else:
        Barcode_df = FilterBarcode(readBarcode(params.barcode), SP_df_raw, key_lineage_num=params.klineages, min_sites=params.nsites, outname=params.fbarcode)
    SP_df = FilterSPDF(Barcode_df, SP_df_raw, outname=params.fsample, outname_p=params.potentials)
    
    print("\nDemixing...", end=' ')
    Out_DF = SolveDemix(SP_df, Barcode_df, max_lineage_num=params.maxlineages)
    
    Out_DF.to_csv(params.output, sep='\t')
    print("Done.")
    
    