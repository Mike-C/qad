# qad
Experimental BWT(qlfc) based compressor

**Version v1.13**
* Removed libdivsufsort, instead parallel (OpenMP) bucket-sort with reduced memory for sorting stage, ~1.35N is used, while a bit slower, it does not handle highly redundant inputs fast as of yet (no special cases), so sorting can really slow down on files containing repetitive runs of characters. Currently there is up-to 8-character deep radix-sort, creating buckets for quick sort, really small buckets < 8 characters, are sorted by insertion-sort.
* Fast mode 'c2' and binary modes 'cb'/'cfb'/'c2b' not maintained in this release, only high compression 'cf' and default compression 'c' are. If you use any of unsupported switches, unsupported methods will be ignored.
* Improved predictions for context selection for compressing tree branches.
* Small files and files with small alphabets are not handled very well, should be fixed in some of later releases. Also archives now contain small amount of redundant data in their footer, so compression is not so suitable for smaller files.
* Console output now contains time of comp./deco.
* Backward compatibility is not maintained.
