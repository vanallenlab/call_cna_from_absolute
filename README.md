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

* `build_for_bands`: Optional. To create summary .tsv files tabulating the cytobands for all amplifications and deletions using hg19 coordinates, set this parameter to "hg19"

Usage
-----
* Without providing output directory: `python CallCNAFromAbsolute.py geneLevelAnnotated/ --build_for_bands hg19`
* Providing output directory: `python CallCNAFromAbsolute.py geneLevelAnnotated/ --output_dir resultsFolder/ --build_for_bandss hg19`
