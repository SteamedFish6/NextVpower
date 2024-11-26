# NextVpower
Next V-Power, adapted from "Virus Phylogenetic Resolver of Wastwater-based Epidemiology (V-Power).
=======
Next V-Power is a tool for SARS-CoV-2 lineage demixing via amplicon sequencing data.

Requirements:
---------------
 - python3
 - numpy
 - pandas
 - scipy

This tool was developed on Windows, and tested on Linux.

Installation:
---------------
1.Clone this repository to your local directory. 
2.Install requiments in a python environment. Skip if you have got these requirments.
```sh
pip install numpy pandas scipy
```
3.Change directory to where you cloned this repository to.
```sh
cd /path/of/NextVPower
```
4.Extract `usher_barcodes.zip` and `var_anno.zip` under current folder.
Output files should be `usher_barcodes.csv` and `var_anno.csv`.
```sh
unzip usher_barcodes.zip
```
```sh
unzip var_anno.zip
```
5.Run python command to check installation and see help.
```sh
python NextVpower.py -h
```

> [!NOTE]
> `usher_barcodes.csv` was generated via [Freyja](https://github.com/andersen-lab/Freyja) and this file was copied from [Freyja repository](https://github.com/andersen-lab/Freyja/blob/main/freyja/data/usher_barcodes.csv).
>
> `var_anno.csv` was downloaded from [NGDC: RCoV19 - Variation Annotation](https://ngdc.cncb.ac.cn/ncov/variation/annotation) and I made some format convertion on it.

Usage: 
---------------

1. Demix from input sample table file, and save result in result.tsv:
```sh
python NextVpower.py -i PP_raw_example.tsv -o demix_result_example.tsv
```
2. Demix from input *.vcf files under a folder, retain top 200 lineages, and save result in result.tsv:
```sh
python NextVpower.py -i vcf_example -v -l example -o demix_result_vcf_example_top200.tsv
```
3. Merge lineages with completely identical mutation sites when Demixing:
```sh
python NextVpower.py -i PP_raw_example.tsv -m -o demix_result_example_merged.tsv
```
4. Add annotation to *.vcf files according to variation annotation table and demix:
```sh
python NextVpower.py -i vcf_example -v -o demix_result_vcf_example.tsv --ann_outpath ann_tab_example
```
5. Demix from input *.vcf files under a folder, save result in result.tsv, and save middle data to files:
```sh
python NextVpower.py -i vcf_example -v -o demix_result_vcf_example.tsv --vcsample PP_raw_example.tsv --fbarcode MMFF_example.tsv --fsample PPFF_example.tsv
```

Please see detailed usage by typing `python NextVpower.py -h` or in [source code](NextVpower.py)

>Example data was published in [Langjun Tang, Zhenyu Guo, Xiaoyi Lu, Junqiao Zhao, Yonghong Li, Kun Yang,
Wastewater multiplex PCR amplicon sequencing revealed community transmission of SARS-CoV-2 lineages during the outbreak of infection in Chinese Mainland,
*Heliyon*, Volume 10, Issue 15, 2024, e35332.](https://doi.org/10.1016/j.heliyon.2024.e35332)


Publications
------------
This project was not published yet —— but you can still have a try on your SARS-CoV-2 amplicon data.

License
-------
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

