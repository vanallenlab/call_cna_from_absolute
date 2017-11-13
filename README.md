# CallCNAFromAbsolute
Given gene annotated segtab files from ABSOLUTE generate amplifications, LOH, and deletions based on process from Brastianos, Carter et al Cancer Discovery 2015. Implementation by David Liu.

Installation
------------
Click the green Clone or Download button above, copy the presented web URL, and then paste it into your terminal, preceded by the command 'git clone'.
  `git clone https://github.com/vanallenlab/call_cna_from_absolute.git`

Parameters
----------
The following parameters can be provided to an AbsMafAnalyzer object:
* `input_dir`: The path to directory containing the input (gene-level annotated segment) *.annotated files.

* `output_dir`: Optional. The path to the desired output directory. If not provided, output files will by default be placed in the input_dir provided.

* `build_for_bands`: Optional. Set this parameter to "hg19" in order to add cytoband information. In addition to creating a file summarizing the copy
number alterations on a cytoband level using hg19 coordinates, this will add two extra columns to the output files, "band" and "arm" (q.31.1 and q might, for example).

Usage
-----
* Without providing output directory: `python CallCNAFromAbsolute.py geneLevelAnnotated/ --build_for_bands hg19`
* Providing output directory: `python CallCNAFromAbsolute.py geneLevelAnnotated/ --output_dir resultsFolder/ --build_for_bands hg19`

Output
------
Example output row:

|genes	|chr|	start|	start_gene|	start_exon|	end	|end_gene	|segment_end_exon	|Num_Probes	|sample	|modal_total_cn	|expected_total_cn	|rescaled.cn.a1	|rescaled.cn.a2	|focality_1	|focality_2	|called_CNA1	|called_CNA2	|fr_below_1	|fr_above_1	|fr_below_2	|fr_above_2	|band	|arm|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
 |AC007435.1|	2|	172967066|	DLX2|	0-|	179666966|	TTN|	|	618|	OCSCC-OC033-TP-NB-SM-F3R7J-SM-F3R8D|	5|	4.99996|	0|	5|	1|	0.987910135|	del|	amp|	0|	0.980218414|	0.958944192|	0.012089865|	q31.1|	q|

Example cytoband summary file .tsv:

```
amp  11q  11q13.2
amp	11q	11q13.2
amp	13q	13q12.2
amp	15q	15q15.1
amp	4q	4q35.1
amp	7q	7q11.23
amp	9q	9q34.3
del	13q	13q21.1
del	17p	17p12
del	17p	17p13.1
del	17p	17p13.3
del	19q	19q13.42
del	3p	3p21.31
del	3p	3p26.3
del	4p	4p14
del	4p	4p16.3
del	6q	6q12
del	6q	6q13
del	7q	7q22.1
del	9p	9p24.3
del	9q	9q34.3
high amp	11q	11q13.3
```
