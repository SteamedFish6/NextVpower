# NextVpower
Next V-Power, adapted from "Virus Phylogenetic Resolver of Wastewater-based Epidemiology (V-Power).
=======
Next V-Power is a tool for multi-virus lineage demixing from amplicon sequencing data, for wastewater and other mixed samples.

> [!IMPORTANT]
> (Updated on 2025/6/25) Now we have generated barcode matrixes of multi-virus extracted from [NextClade phylogenetic tree](https://github.com/nextstrain/nextclade_data/tree/release). To specify a barcode matrix, please use parameter `-b` to select a file.
>
> Apart from demixing sars-cov-2 with usher_barcodes, you can try demixing other kind of viruses with our new barcode matrix, including IAV, IBV, hmpv, measles, mpox, rsv, dengue, prrsv1 and marburgvirus.
>
> Parameters `-n` and `-k` may vary for different kind of viruses.
>
> Our barcode extracter tool available (see `BarcodeExtracter.py`), then you can update barcode matrixes from NextClades source on time!

Basic flow sheet:
---------------
![image](https://github.com/user-attachments/assets/8471c6cc-c55c-4156-81a8-8c379aaa2b2e)

Note: Red numbers indicate adjustable parameters. Green numbers may change according to your sequencing data and the updating of the barcode matrix A.


Requirements:
---------------
 - python3
 - numpy
 - pandas
 - cvxpy

This tool was developed on Windows, and tested on Linux.

Installation:
---------------
1.Clone this repository to your local directory. 

2.Install requiments in a python environment. Skip if you have reached these requirments.
```sh
pip install numpy pandas cvxpy
```
3.Change directory to where you cloned this repository to.
```sh
cd /path/of/NextVpower #example path
```
4.Uncompress `resource.tar.gz` under current directory.
Output folder should contains `usher_barcodes.csv`, `var_anno.csv`, and other barcodes files like `clade_barcodes_flu_h1n1pdm_ha_CY121680.csv`. 
```sh
tar -zxvf resource.tar.gz
```
 - For Windows, please uncompress the tar.gz file to the folder where they are current located.

5.Run python command to check installation and see help.
```sh
python NextVpower.py -h
```


> [!NOTE]
> `usher_barcodes.csv` was generated via [Freyja](https://github.com/andersen-lab/Freyja) and this file was copied from [Freyja repository](https://github.com/andersen-lab/Freyja/blob/main/freyja/data/usher_barcodes.csv).
>
> `var_anno.csv` was downloaded from [NGDC: RCoV19 - Variation Annotation](https://ngdc.cncb.ac.cn/ncov/variation/annotation) and we made some format convertion on it (Only for sars-cov-2 annotation).
>
> Unlike Freyja, NextVpower takes into account INDEL-type mutations when processing input data, which may provide more accurate lineage demixing results. Also, NextVpower can detect the potential mutation sites of each sample.

> [!NOTE]
> `Vpower.m` and `Vpower2.m` are MATLAB scripts of V-Power.
>
> `NextVpower.py` is a standalone program, which has integrated the function of `Vpower.m` and `Vpower2.m` without using MATLAB engine.

Usage: 
---------------

1. Demix from input sample mutation table file, and save result to result.tsv:
```sh
python NextVpower.py -i PP_raw_example.tsv -o demix_result_example.tsv
```
2. Demix from input *.vcf files under a folder, filter mutation sites with mutation rate lower than 0.1, filter mutation sites with depth lower than 10, and save result to result.tsv:
```sh
python NextVpower.py -i vcf_example -r 0.1 -d 10 -o demix_result_vcf_example.tsv
```
3. Set barcode filter criteria, filter lineages with fewer than 0.02 of total "key" mutation sites (default 0.01), retain "key" mutation sites present in more than 0.1 of total lineages (default 0.05):
```sh
python NextVpower.py -i PP_raw_example.tsv -n 0.02 -k 0.1 -o demix_result_example_300_30.tsv
```
4. Add annotation to sample mutation table according to variation annotation table and demix:
```sh
python NextVpower.py -i vcf_example -o demix_result_vcf_example.tsv --ann_vcsample Annotated_PP_raw_example.tsv
```
5. Demix from input *.vcf files under a folder, save result in result.tsv, and save middle data to files:
```sh
python NextVpower.py -i vcf_example -o demix_result_vcf_example.tsv --vcsample PP_raw_example.tsv --fbarcode MMFF_example.tsv --fsample PPFF_example.tsv --potentials potential_sites_example.tsv
```

Please see detailed usage by running `python NextVpower.py -h` or in [source code](NextVpower.py)

>Example data was published in [Langjun Tang, Zhenyu Guo, Xiaoyi Lu, Junqiao Zhao, Yonghong Li, Kun Yang,
Wastewater multiplex PCR amplicon sequencing revealed community transmission of SARS-CoV-2 lineages during the outbreak of infection in Chinese Mainland,
*Heliyon*, Volume 10, Issue 15, 2024, e35332.](https://doi.org/10.1016/j.heliyon.2024.e35332)

> [!TIP]
> How to create VCF files?
> 
> - 1.Align sample.reads.fq to reference.fa by running `bwa mem`, `samtools sort` and `samtools view`, to generate sample.bam.
> 
> - 2.Use variant calling tools ([freebayes](https://github.com/freebayes/freebayes) is recommended) to generate sample.vcf from sample.bam.

Here's an example of demixing using reads file:
```sh
bwa index reference.fa
samtools faidx reference.fa
bwa mem -t 4 reference.fa sample.reads.fq | samtools view -bS --threads 4 | samtools sort --threads 4 -o sample.bam
samtools index sample.bam
freebayes --haplotype-length 0 --pooled-continuous -m 20 -p 1 --strict-vcf -q 13 --min-coverage 10 -F 0.4 -f reference.fa sample.bam | \
    bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT | \
    bcftools norm -f reference.fa -s -m +any -Ov > sample.vcf
python NextVpower.py -i sample.vcf -b virus_barcode.csv -o sample.demix_result.tsv
```


Publications
------------
This project was not published yet —— but you can still have a try on your virus amplicon sequencing data.

License
-------
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

