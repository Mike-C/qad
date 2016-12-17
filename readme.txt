QAD, Qlfc based BWT compressor v1.13, (c) Petr Kucera 2014 - 2016

You are using this software at your own risk. Author does not bear responsibility for any potential damage concerning the use. Note, due to experimental nature and to lack of thorough tessting, decompression may fail on some files.

Usage: qad 'c'/'cf'/'c2'/'d' input-file output-file ('c' for default compress, 'cf' default mode with experimental filter possibly suiting text data, 'c2' for fast compress, 'd' to unpack). If you want to output raw binary data, with no entropy coding, use either a 'cb'/'cfb'/'c2b' switch. 'v' print curent program version.


Version notes:

QAD v1.13
* Removed libdivsufsort, parallel bucket-sort with reduced memory for sorting part, ~1.35N used instead, while a bit slower, it does not handle
  highly redundant inputs fast yet, so it can really slow down on files containing repetitive runs of same characters.
  Currently there is 3-character deep radix-sort, creating buckets for quick sort, really small buckets < 8 characters, are sorted by insertion-sort.
* Fast mode 'c2' and binary modes 'cb'/'cfb'/'c2b' not maintained in this release, only high compression 'cf' and medium compression 'c' are.
  If you use any of unsupported switches, unsupported methods will be ignored.
* Improved predictions for context selection for compressing tree branches.
* Small files and files with small alphabets not handled very well, should be fixed in some of later releases.
  Also archives now contain small amount of redundant data in their footer, so are not so suitable for small files.
* Console output now contains time for comp/dec.
* Backward compatibility not maintained.

QAD v1.12
* New swith 'b' added. Combined with other compression modes (i.e. 'cb'/'cfb'/'c2b'), program outputs raw binary data without entropy coding.
* Bug fix for program parameters parsing, if the name of compressed file started with '2' or 'f' characters before, wrong compression mode might have been selected.
* Ranks 2+ have been split into two channels, according to first rank, similar to as what was done before with RLE (or the zero rank), this follows that first rank split. Their subcontexts (in these new extra channels) are choosed acording to up to four previous ranks.

QAD v1.11
* Context selection was slightly improved, now also previous character is taken into account.
* New parameter v shows current program version.
* Small reduction of memory requirements on decoder part.
* various bug-fixes, i.e. header/footer parsing, fixed state for escape characters, compressing files with small alphabets
* Added shell script for testing. Script 'q' creates log file q.log, where you can watch online compression results i.e. with tail -f q.log. Script takes same parameters as program for compression, i.e: ./qad c bib bib.qad

QAD v1.10
* Secondary context selection added on tree-based RLE, improved context selection. 

QAD v1.09
* Added secondary context selection for ranks, according ranks( their metrics ) just to follow. 

QAD v1.08
* Slight tune of RC to match better the latest changes.
* 64bit binaries should be from now on included as well, not just a 32bit emulated code.

QAD v1.07
* Constant value removed from rank context selection formula. Its been replaced by length of last continuous RLE segment, while RLE context selection formula now use previous rank, instead of actual one.

QAD v1.06
* Ranks were split into two channels, like has been done with RLE previously.
* Decoder now uses less memory

QAD v1.05
* 8 passes in default mode reduced to 5, speeding up compression and decompression.
* Second rle channel( in default mode as well ) now use mostly the same context as first one, to benefit from latest changes.

QAD v1.04
* 2 & 3 context state mix. Bit 1 in context 2 is more favorized if there is 1 in context 3, rather than 0.
* Improved context selection for the RLE

QAD v1.03
* Fix for buffer overflow(as to number of subcontexts) on higly repetitive data
* Buffer rotation for default mode decode has been reduced for only one rle channel, it should be removed in future
* There is still a risk of buffer underflow for RC(mostly on 64bit) which needs to be fixed.

QAD v1.02
* Memory req. changes on encoder side, BWT computation buffer is now reused 
* Number of subcontexts reduced ten times, previous number took quite a lot of memory, there is a danger of overflow on high repetitive data, fix needed.

QAD v1.01
* A quick-fix for fast mode decode
* Added bit history into RC for the last 7 bits, weights are yet not set to suit contexts.

QAD v1.0
* Improved and simplified default compression (c) mode and filtered mode (cf).
* Removed rank context selection according to RLE (small compression gains, new context selection of RLE according ranks uses those states instead).
* To solve memory management issues from previous version, decoder algorithm has been simplifies. Decoder now temporary uses more memory and is slower on default compression or cf mode.
* Both default and filtered mode are not backward compatible with previous version.

QAD v0.9a
* The previous change turned to be highly input data sensitive. Due to very rough freq. distribution approximation, It has been stripped of the default compression mode and added as an experimental filter under cf archiver switch option, otherwise archive compressed on default mode compressed with version 0.8 should be compatible with this version too. Fast mode archives from 0.7, onwards.

QAD v0.8a
* Removed constant (on default mode), for RLE context selection according to rank( >= TREESIZE(9) before),  this should now approximate the actual rle frequency distribution and should improve compression ratio on default mode.
* This change also amends archive footer (adding an extra integer), and breaking backward compatibility with previous versions.

QAD v0.7a
* Added libdivsufsort 2.0.1 (MIT Licence), memory management increases
* 32bit and emulated 64bit binaries now included both for Linux and Windows

QAD v0.6a
* Slightely improves RLE, you can turn on/off DelayRleBitMantisaCtx and RLEtoRankPrediction at compile time
* At compile time, you can also define, if you want to build 32bit/64bit by Build32bit / Build64bit. Or try other output aligment than bit or byte, yet this is only experimental
* From now on, there is a source code included, licence changes to GPL 2.0s

QAD v0.5a
* Utilizes a modification of Subbotin Range Coder instead of Arithmetic Coder.
* Slightely worse compression ratio with faster execution, even though the Win 32bit executable was compiled with 64bit emulation only, which slows execution down. Linux binary has not been included this time, because 64bit code emulation crashes gcc. For 32bit version a 32bit RC with a bit-aligned encode output will be used.

QAD v0.4a
* Improved RLE, added new compression mode, currently unoptimized. To switch between modes use c/c2. (c - new mode, c2 faster compression mode, the previous one).
* Several contexts packed together, to improve memory management further
* A new package header 

QAD v0.3a
* Fix for footer parsing issue, when alphabet size is < 10. For these, the structure has been changed. Breaking the backward compatibility with previous versions, as of yet footer is not compressed.

QAD v0.2a
* Added RLE, which provide a slight compression boost and speed-up

QAD v0.1a
* There is no context mixing or RLE involved, affecting compression ratio.
* There is no BWT transform in the demo, future versions should include either a divsufsort or a simple quicksort and change licence to FLOSS.
* Header/Footer is not packed, the reverse compatibility will possibly not be kept in future versions
* Slowest part is possibly entropy coder (AC in this case), RC should be used instead 
* Binary is distributed both as crosscompiled Windows version and Linux version
* For some reason Windows version tends to be almost as twice as fast on Linux under Wine than Linux binary.