# input file format

All input file were tab-separate.

## transcripts

|ENSEMBL | barcode1 | fluorescence_probe_binding_region | barcode2 | SYMBOL | fluorescence_label|
| :------: | :------: | :------: | :------: | :------: | :------: |
|ENSMUSG00000004891|A|TGTCTATGTTTACAGCGGGC|C|Nes|DO1-G|


## junction
|AlternativeSplice | barcode1;fluorescence_probe_binding_region;barcode2 | barcode1;fluorescence_probe_binding_region;barcode2 |
| :------: | :------: | :------: |
|1:4886832-4889456+_1:4886832-4889459+|T;CAGTGAATGCGAGTCCGTCT;G|T;CAGTGAATGCGAGTCCGTCT;G|

## circRNA
|BackSplice | barcode1 | fluorescence_probe_binding_region | barcode2|
| :------: | :------: | :------: |:------: |
|1:10058824-10059823+|T|CAGTGAATGCGAGTCCGTCT|G|