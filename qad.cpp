/*
 * qad.cpp - base class providing RPC like mechanisms
 *
 * (c) 2014 - 2016 Petr Kucera, kuc406.moxo.cz <quadro2/at/email.cz>
 *
 * The BWT-QLFC based compressor, utilizing a cyclic tree and RLE
 * compression of either branches or states. Final entropy coder,
 * either 32-bit or 64-bit is based on a Subbotin's Range Coder
 * implementation.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program (see COPYING); if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301 USA.
 *
 */

#include <time.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <unistd.h>
//#include <math.h>
#include <stdlib.h>
//#include <iomanip>
#include <cstring>
//#include <fcntl.h>

//#define Build32bit // for 32-bit RC with bit-aligned output

#define Build64bit // for 64-bit RC with byte-aligned output

#define VERSION 1.13

#define DelayRleBitMantisaCtx // affects both compression methods 'c' and 'c2'
//#define RLEtprankPrediction	// affects compression method 'c' only, this is implemented also vice-versa,
								// as one of changes between compression methof c and c2
								// note since version 1.0, RLEtprankPrediction is not used
//#define DEBUG


#ifdef Build64bit
 // Byte-aligned 64bit RC output limists
 #define TOP 0x0100000000000000
 #define ALI 8 // output aligment in bits
#endif

#ifdef Build32bit
 // Bit-aligned 32bit RC output limists
 #define TOP 0x80000000
 #define ALI 1
#endif

// nibble bassed, 32 bit
//#define TOP 0x10000000

// 2D array -> 1D recomputation
#define p( ctx, mode, c ) ctx + ( mode * CTXs ) + ( ( c ) * SubCTXs * CTXs )

// decoder read single bit for bit-aligned input
#ifdef Build32bit
 #define GetRCb( c, p ) ( c << 1 ) + \
	( ( ( unsigned char )inBuf[ ctx ][ 3 + ( p >> 3 ) ] >> ( 7 - ( p & 7 ) ) ) & 1 )
#else
 #define GetRCb( c, p ) ( c << 1 ) + \
	( ( ( unsigned char )inBuf[ ctx ][ 7 + ( p >> 3 ) ] >> ( 7 - ( p & 7 ) ) ) & 1 )
#endif

#define MAXCHAR 0xFF

using namespace std;

unsigned long size;

unsigned int bwtIndex;	// BWT Index
bool bwtI = false;

const int READ_SEG = 16384;
const int BUFPOS = 8;
const int TREESIZE = 9;	// cyclic tree to encode mtf/qlfc ranks

// maximum treshold for (rdiv = char / TREESIZE), treeRLE method
const int TREETRESH = 15; // set 0 to ignore, default 15

//RC Variables
const int alphabet = 1;	// size of alphabet
const int CTXs = 12;	// number of contexts
const int SubCTXs = 430000; // /subcontexts

struct read64 {			// struct to help parse large numbers
#ifdef Build64bit		// -not so effective
	unsigned int  b1;
#endif
	unsigned short b2;
	unsigned char b3;
	unsigned char b4;
} ;

#ifdef Build64bit
unsigned long long LR[ CTXs ]; // left + range
unsigned long long RCLow[ CTXs ] = { 0 }, RCRange[ CTXs ];
#endif
#ifdef Build32bit
unsigned int LR[ CTXs ];
unsigned int RCLow[ CTXs ] = { 0 }, RCRange[ CTXs ];
#endif
read64 * readLow[ CTXs ], * readLr[ CTXs ]; // struct pointer for parsing large numbers

unsigned int * fq;	// symbol cumulative prob buffer (1D, will be mapped into 2D with p( macro)
unsigned int * sRt;	// symbol rank table (for cum. prob. table re-scalle RC hashing)
unsigned int pos[ CTXs ] = { 0 }, isCTXUsed[ CTXs ] = { 0 }; // possition in buffer, context usage identifier
unsigned int RCbufP[ CTXs ] = { 0 }; // buffer position

// DeRC
#ifdef Build64bit
unsigned long long code[ CTXs ]; // input code for RC decoder
#endif
#ifdef Build32bit
unsigned int code[ CTXs ];
#endif
unsigned long pozInBufBin[ CTXs ] = { BUFPOS };
char * inBuf[ CTXs ];				// input buffer

//RC Variables ends

char * outBuf[ CTXs ], storedChars[ MAXCHAR + 1 ]; //variety buffers
long b[ MAXCHAR + 2 ], i, p0, p1, treshold, pivot = -1, aC = 0, trueMax = 0, sumLen;

int inFileOffst, pOffset, tt[ MAXCHAR + 1 ], sf[ MAXCHAR + 3 ], sf2[ MAXCHAR + 3 ], 
	sc[ MAXCHAR + 3 ], pendPos[ MAXCHAR + 3 ], order[ MAXCHAR + 3 ];	// some sf- symbol freq. are not needed

unsigned long rleSeq = 0;	// rleSeq - length of RLE sequence
unsigned long r2[ 8 ];		// approximate RLE distribution
long rank = 0, rdiv = 0;	// rdiv - rank modulo tree size
unsigned long RLEfreq[ MAXCHAR + 1 ], RLEtop[ 4 ], RLEdistEst = 0, RNKcnt = 0, RNKcnt2 = 0;

unsigned int * bq;

unsigned char * data, * bQlfc, * bRLE, stripFilters = 0, globalMode = 0; // compression mode, 0 = default, 1 = fast
unsigned char binOut = 0, bits[ MAXCHAR + 1 ], bitH[ MAXCHAR + 1 ] = { 0 };  // bit history (7 bits);
unsigned char lbit[ CTXs ] = { 0 };

ofstream * outF;

clock_t start = clock();

/*
 Precompute symbol freq. table with assumption that frequency of
 higher ranks will gradually decrease after qlfc/mtf transform
*/
void FillInAritfq( void ) {
	int t = TREESIZE, f[ 25258 /*7355*/ ];
	//#pragma omp parallel for
	for( int i = 0; i <= 14670; i++ ) f[ i ] = ( 14670 + 1 - i ) << 1;
	//#pragma omp parallel for
	for( int i = 0; i <= 14670; i++ ) {
		int b = 0;
		int m = ( i ) % TREESIZE;
		if( m >= 8 ) b = 1, fq[ p( 0, i / t, b ) ] += f[ i ];
		else  b = 0, fq[ p( 0, i / t, b ) ] += f[ i ], fq[ p( 0, i / t, 1 ) ] += f[ i ];
		if( m <= 2 ) {
			if( CTXs > 1 ) b = 1, fq[ p( 1, i / t, b ) ] += f[ i ];
			if( m == 0 ) {
				if( CTXs > 2 ) b = 0, fq[ p( 2, i / t, b ) ] += f[ i ], fq[ p( 2, i / t, 1 ) ] += f[ i ];
			} else {
				if( CTXs > 2 ) b = 1, fq[ p( 2, i / t, b ) ] += f[ i ];
				if( m == 1 ) {
					if( CTXs > 3 ) b = 0, fq[ p( 3, i / t, b ) ] += f[ i ], fq[ p( 3, i / t, 1 ) ] += f[ i ];
				} else {
					if( CTXs > 3 ) b = 1, fq[ p( 3, i / t, b ) ] += f[ i ];
				}
			}
		} else {
			if( CTXs > 1 ) b = 0, fq[ p( 1, i / t, b ) ] += f[ i ],  fq[ p( 1, i / t, 1 ) ] += f[ i ];
			if( m <= 4 ) {
				if( CTXs > 4 ) b = 0, fq[ p( 4, i / t, b ) ] += f[ i ], fq[ p( 4, i / t, 1 ) ] += f[ i ];
				if( m == 3 ) {
					if( CTXs > 6 ) b = 1, fq[ p( 6, i / t, b ) ] += f[ i ];
				} else {
					if( CTXs > 6 ) b = 0, fq[ p( 6, i / t, b ) ] += f[ i ], fq[ p( 6, i / t, 1 ) ] += f[ i ];
				}
			} else {
				if( CTXs > 4 ) b = 1, fq[ p( 4, i / t, b ) ] += f[ i ];
				if( m <= 6 ) {
					if( CTXs > 5 ) b = 1, fq[ p( 5, i / t, b ) ] += f[ i ];
					if( CTXs > 8 ) b = 1, fq[ p( 8, i / t, b ) ] += f[ i ];
					if( CTXs > 10) b = 1, fq[ p( 10, i / t, b ) ] += f[ i ];
					if( m == 5 ) {
						if( CTXs > 7 ) b = 1, fq[ p( 7, i / t, b ) ] += f[ i ];
						if( CTXs > 9 ) b = 1, fq[ p( 9, i / t, b ) ] += f[ i ];
					} else {
						if( CTXs > 7 ) b = 0, fq[ p( 7, i / t, b ) ] += f[ i ], fq[ p( 7, i / t, 1 ) ] += f[ i ];
						if( CTXs > 9 ) b = 0, fq[ p( 9, i / t, b ) ] += f[ i ], fq[ p( 9, i / t, 1 ) ] += f[ i ];
					}
				} else {
					if( CTXs > 5 ) b = 0, fq[ p( 5, i / t, b ) ] += f[ i ], fq[ p( 5, i / t, 1 ) ] += f[ i ];
					if( CTXs > 8 ) b = 0, fq[ p( 8, i / t, b ) ] += f[ i ], fq[ p( 8, i / t, 1 ) ] += f[ i ];
					if( CTXs > 10) b = 0, fq[ p( 10, i / t, b ) ] += f[ i ], fq[ p( 10, i / t, 1 ) ] += f[ i ];
				}
			}
		}
	}
}

/*
	Initialize range decoder, parse compressed file
*/
int DeBinRCInit( char * deFile ) {

  ifstream file( deFile, ios::in | ios::binary | ios::ate );
  unsigned long citO[ CTXs ], pom = 0;
  int ctxMark = 16383, ctxCount = 0, ctxO = 0;  // read used contexts

  long sh;
  if( file.is_open() )
  {
	size = ( unsigned long ) file.tellg();
	sh = size;
	file.seekg( sh - ( sizeof( char ) + sizeof( char ) ), ios::beg );
	file.read( ( char *) &( trueMax ), sizeof( char ) );
	file.read( ( char *) &( globalMode ), sizeof( char ) );
	RLEdistEst = globalMode >> 4;
	binOut = ( globalMode >> 1 ) & 1;
	globalMode &= 0x01;
	int sTail = ( 2 * sizeof( char ) ) + sizeof( char );
	if( globalMode != 1 ) {
		sTail += sizeof( unsigned int ) + sizeof( unsigned int );
		file.seekg( sh - ( sTail - sizeof( char ) ), ios::beg );
		file.read( ( char *) &( RNKcnt ), sizeof( unsigned int ) );
		file.read( ( char *) &( RNKcnt2 ), sizeof( unsigned int ) );
	}
	if( trueMax < 10 || binOut ) {
		ctxO = 2;
		file.seekg( sh - sizeof( char ) - sTail, ios::beg );
		file.read( ( char *) &( ctxMark ), 2 );
	}
	if( globalMode == 1 ) ctxMark &= 2047;

	for( i = 0; i < CTXs; i++ ) if( ( ( ctxMark >> i ) & 1 ) == 1 ) ctxCount++;

	int sub = 0, pairs = 0;
	if( trueMax > 0 && trueMax < MAXCHAR - 1 ) {
		file.seekg( sh - sTail - ctxO, ios::beg );
		file.read( ( char *) &( pairs ), sizeof( char ) );

		if( pairs > 0 ) {
			file.seekg( sh - ( ( ( pairs * 2 ) ) * sizeof( char ) ) - sTail - ctxO, ios::beg );
			for( i = 0; i < pairs; i++ ) {
				file.read( ( char *) &( p0 ), sizeof( char ) );
				file.read( ( char *) &( p1 ), sizeof( char ) );
				if( p1 > p0 )
					for( int r = p0; r <= p1; r++ ) storedChars[ sub++ ] = r;
			}
		}
	}

	sTail = ( globalMode != 1 ? sizeof( int )<<1 : 0 );
	int remCnt = 0, dex = sTail + ( ( trueMax > 0 && trueMax < MAXCHAR - 1) ? 2 : 1 );
	if( trueMax <= 127 ) remCnt = trueMax + 1; else remCnt = MAXCHAR + 1 - ( trueMax + 1 );

	if( remCnt - sub > 0 ) {
		remCnt -= sub;
		file.seekg( sh - (( dex + ( pairs << 1 ) + remCnt ) * sizeof( char ) ) - ctxO - 1, ios::beg );
		for( i = 0; i < remCnt; i++ ) {
			file.read( ( char *) &storedChars[ sub++ ], sizeof( char ) );
		}
		remCnt += pairs << 1;
	} else remCnt = pairs << 1;

	file.seekg( sh - ( ( dex + remCnt ) * sizeof( char ) ) - ( ctxCount * sizeof( int ) ) -
												 3 * sizeof( int ) - ctxO - 1, ios::beg );
	file.read( ( char *) &( bwtIndex ),		sizeof( int ) );
	file.read( ( char *) &( sf2[ 0 ] ), 2 * sizeof( int ) );
	file.read( ( char *) &( pos[ 0 ] ), /*CTXs*/ctxCount * sizeof( int ) );

	for( int ctx = 0; ctx < CTXs; ctx++ ) {
		if( ( ( ctxMark >> ctx) & 1 ) == 1 ) citO[ ctx ] = pos[ pom++ ]; else citO[ ctx ] = 0;
	}

	sf2[ MAXCHAR + 1 ] = sf2[ 1 ];
	sumLen = sf2[ 1 ];

	for( i = 1; i < MAXCHAR + 1; i++ )  sf2[ i ] = -1;

	remCnt = sub;
	if( remCnt == 0 )
		for( i = 1; i < MAXCHAR + 1; i++ ) sf2[ i ] = 1; // all symbols will be used
	if( trueMax <= 127 )
		for( i = 0; i < remCnt; i++ ) sf2[ ( unsigned char ) storedChars[ i ] + 1 ] = 1;
	if( trueMax > 127 && trueMax < MAXCHAR ) {
		for( i = 1; i < MAXCHAR + 1; i++ ) sf2[ i ] = 1; // all symbols will be used
		for( i = 0; i < remCnt; i++ ) sf2[ ( unsigned char )storedChars[ i ] + 1 ] = -1; // mask unused symbols
	}

	file.seekg( 0, ios::beg );
  }

  for( int ctx = 0; ctx < CTXs; ctx++ )
  {
	if( citO[ ctx ] == 0 ) continue;

	inBuf[ ctx ] = (char*) malloc ( citO[ ctx ] + 8 );
	if( inBuf[ ctx ] == NULL ) cout << "\n Not enough free memory", exit( 1 );

	file.read( ( char *) &inBuf[ ctx ][ 0 ], citO[ ctx ] );

#ifdef Build64bit
	for( int x = 0; x < 8; x++ ) {
#endif
#ifdef Build32bit
	for( int x = 0; x < 4; x++ ) {
#endif
	 code[ ctx ] <<= 8;
	 code[ ctx ] += ( unsigned char )inBuf[ ctx ][ x ];
	}
  }

  return 0;
}


void Bin( int ctx, unsigned char cC, bool finish = false ) {
	isCTXUsed[ ctx ] = 1;

	// if finish == TRUE, write out compressed buffer
	if( finish ) outF->write( ( char *) &outBuf[ ctx ][ 0 ], pos[ ctx ] );

	outBuf[ ctx ][ pos[ ctx ] ] = ( outBuf[ ctx ][ pos[ ctx ] ] << 1 ) + cC;

	if( ( ++RCbufP[ ctx ] & 7 ) == 0 ) pos[ ctx ]++;
}

//inline
bool DeBin( int ctx ) {

#ifdef Build64bit
	bool bin = ( ( code[ ctx ] >> 63 ) );
#else
	bool bin = ( code[ ctx ] >> 31 );
#endif
	code[ ctx ] = GetRCb( code[ ctx ], pozInBufBin[ ctx ] ), pozInBufBin[ ctx ]++;

	return bin;
}


#define STATEMIX_23 \
	if( ctx == 2 && cC ) fq[ p( 3, sub, 1 ) ] += 1, fq[ p( 3, sub, 0 ) ] += 1; \
	if( ctx == 10 && sub == 0 && cC ) fq[ p( 10, 1, 1 ) ] += 1, fq[ p( 10, 1, 0 ) ] += 1; \
	if( ctx == 3 ) { if( cC == 0 ) { if( fq[ p( 2, sub, 1 ) ] > 10 ) fq[ p( 2, sub, 1 ) ] -= 10; \
		} else fq[ p( 2, sub, 1 ) ] += 24; \
	}


int eff, ef[ CTXs ] = { 100 }, ef2[ CTXs ] = { 100 };
	// update character cum. prob. according bit history (the last 7 bits), or without it
#define UPDATE_PROBABILITY( b ) \
	if( ef[ b ] > 7900 ) \
		ef[ b ] >>= 1, ef2[ b ] >>= 1; \
	eff = ( ef[ b ] + 1 ) / ( ef2[ b ] + 1 ); \
	eff = 100 + ( eff > 8 ?( eff > 10 ? 59 : 29 ): 8 ); \
	if( ( ctx == 11 || ( ctx >= 2 && ctx <= 8 ) ) && sub <= 2 ) { \
		unsigned char id = ( ( ctx - 2 ) << 2 ) + sub, bitC = bits[ bitH[ id ] = ( bitH[ id ] << 1 ) + cC ]; \
		if( cC == 0 ) { \
			if( bitC <= 1 ) id = eff - ( 64 - ( bitC ? 48 : 0 ) ), fq[ ft[ 0 ] ] += id, fq[ ft[ 1 ] ] += id; \
			else fq[ ft[ 0 ] ] += eff + 8, fq[ ft[ 1 ] ] += eff + 8; \
		} else { \
			if( bitC >= 5 ) fq[ ft[ cC ] ] += eff - ( 16 - ( ( 7 - bitC ) << 1 ) ); \
			else fq[ ft[ 1 ] ] += eff + 2; \
		} \
	} else { \
		fq[ ft[ cC ] ] += eff; \
		if( cC == 0 ) fq[ ft[ 1 ] ] += eff; \
	} \
	ft[ 3 ] = p( ctx, sub, 3), sRt[ ft[ 0 ] ]++;

	// context-adaptive hash for freq. table re-scale, handled via &3
	// possibly little help only for context 2,(or below) which is quite redundand
	// edit: added context 11 (it's derived from context 2)
	// -don't go above 2^14 in total cumulative freq, the bottom limit
#define FREQUENCY_TABLE_RESCALE \
	if( fq[ ft[ 1 ] ] >= 13900 || ((sRt[ ft[ 0 ]] >= ( sRt[ ft[ 3 ] ] >> 2 ) + 30 && \
										( ctx == 2 || ctx == 11 ) )  ) ) \
	{ \
		if( ctx == 2 || ctx == 11 ) { \
			ft[ 2 ] = p( ctx, sub, 2 ), sRt[ ft[ 1 ] ]++; \
			sRt[ ft[ sRt[ ft[ 1 ] ] & 3 ] ] = sRt[ ft[ 0 ] ], sRt[ ft[ 0 ] ] = 0; \
		} \
		fq[ ft[ 1 ] ] -= fq[ ft[ 0 ] ] >> 1, fq[ ft[ 0 ] ] -= fq[ ft[ 0 ]] >> 1; \
		fq[ ft[ 1 ] ] -= ( fq[ ft[ 1 ] ] - fq[ ft[ 0 ] ] ) >> 1; \
	};

/*
 Range decoder, based on Subbotin's RC, ft-pointer into symbol cum.prob.
 table fq, which is 2D array, stored 1D. ctx-context, sub-subcontext
 (sharing the same Low, Range & code-space), cC-coded character
*/
#ifdef DEBUG
long binRCcnt[ CTXs ] = { 0 }, bitRLECcnt[ CTXs ] = { 0 }, treeRLEcnt = 0;
#endif
unsigned int DeBinRC( int ctx, int sub, unsigned char cC = 0 ) {
#ifdef DEBUG
	binRCcnt[ ctx ]++;
#endif
	if( binOut ) return DeBin( ctx );
	unsigned long ft[ 4 ], total = fq[ ft[ 1 ] = p( ctx, sub, 1 ) ];
	ft[ 0 ] = p( ctx, sub, 0 );

	RCRange[ ctx ] /= total;
	cC = ( ( ( code[ ctx ] - RCLow[ ctx ] ) / ( RCRange[ ctx ] ) ) >= fq[ ft[ 0 ] ] );
	RCLow[ ctx ] += ( cC == 0 ? 0 : fq[ ft[ 0 ] ] ) * RCRange[ ctx ];
	RCRange[ ctx ] *= ( unsigned long ) ( cC == 0 ? fq[ ft[ 0 ] ] : total - fq[ ft[ 0 ] ] );

	STATEMIX_23

	UPDATE_PROBABILITY( ctx )

	FREQUENCY_TABLE_RESCALE
	ef[ ctx ]++;

	LR[ ctx ] = RCLow[ ctx ] + RCRange[ ctx ];
	while( (
	//( readLow[ ctx ]->b4 == readLr[ ctx ]->b4 ) ) || //comparision instead xor, seems not work yet
		( RCLow[ ctx ] ^ ( LR[ ctx ] ) ) < TOP ) ||
			( RCRange[ ctx ] <= 0xFFFF && ( (RCRange[ ctx ] = -RCLow[ ctx ] ), 1 ) ) ) {
			ef2[ ctx ]++;
#ifdef Build64bit
// byte-aligned input
		code[ ctx ] = ( code[ ctx ] << ALI ) + ( unsigned char ) inBuf[ ctx ][ pozInBufBin[ ctx ]++ ];
#else
// bit-aligned input
		code[ ctx ] = GetRCb(code[ ctx ], pozInBufBin[ ctx ] ), pozInBufBin[ ctx ]++;
#endif
		RCLow[ ctx ] <<= ALI, RCRange[ ctx ] <<= ALI, LR[ ctx ] <<= ALI; // RCLow[ctx] + RCRange[ctx];
	}
	return cC;
}


/*
 Range coder, based on Subbotin's RC, ft-pointer into symbol freq.
 table fq, which is 2D array, stored 1D. ctx-context, sub-subcontext
 (sharing the same Low, Range & code-space), cC-coded character
*/
void BinRC( int ctx, int sub, unsigned char cC, bool finish = false ) {
#ifdef DEBUG
	binRCcnt[ ctx ]++;
#endif
	if( binOut ) { Bin( ctx, cC, finish ); return; }
	unsigned long ft[ 4 ], total = fq[ ft[ 1 ] = p( ctx, sub, 1 ) ]; isCTXUsed[ ctx ] = 1;
	ft[ 0 ] = p( ctx, sub, 0 );

	RCRange[ ctx ] /= total;
	RCLow[ ctx ] += ( cC ? fq[ ft[ 0] ] : 0 ) * RCRange[ ctx ];
	RCRange[ ctx ] *= ( cC ? total - fq[ ft[ 0 ] ] : fq[ ft[ 0 ] ] );

	STATEMIX_23

	UPDATE_PROBABILITY( ctx )

	FREQUENCY_TABLE_RESCALE
	ef[ ctx ]++;

	// if finish == TRUE, write out compressed buffer
	if( finish ) outF->write( ( char *) &outBuf[ ctx ][ 0 ], pos[ ctx ] );

	LR[ ctx ] = RCLow[ ctx ] + RCRange[ ctx ];
	while( (
	//( readLow[ ctx ]->b4 == readLr[ ctx ]->b4 ) ) || // comparision instead ^, somehow buggy
		( RCLow[ ctx ] ^ ( LR[ ctx ] ) ) < TOP ) ||
			( RCRange[ ctx ] <= 0xFFFF && ( ( RCRange[ ctx ] = -RCLow[ ctx ] ), 1 ) ) ) {
				ef2[ ctx ]++;
#ifdef Build64bit
		outBuf[ ctx ][ pos[ ctx ]++ ] = ( RCLow[ ctx ] >> 32 ) >> 24; // >> 56, byte aligned output
		//outBuf[ ctx ][ pos[ ctx ]++ ] = readLow[ ctx ]->b4; // byte aligned output (64bit emulation)
#else
// (seem larger aligments, up to 1Byte, are less error prone than 1bit with a 32bit Low)
		if( ( ++RCbufP[ ctx ] & ( ( 8 / ALI ) - 1 ) ) == 1 )
			outBuf[ ctx ][ pos[ ctx ] ] = readLow[ ctx ]->b4 >> ( 8 - ALI );
		else outBuf[ ctx ][ pos[ ctx ] ] =
			( outBuf[ ctx ][ pos[ ctx ] ] << ALI ) + ( readLow[ ctx ]->b4 >> ( 8 - ALI) ),
			( ( RCbufP[ctx] & ( ( 8 / ALI ) - 1) ) != 0 ) ? 1 : pos[ ctx ]++;
#endif
		RCLow[ ctx ] <<= ALI, RCRange[ ctx ] <<= ALI, LR[ ctx ] <<= ALI; //RCLow[ctx] + RCRange[ctx];
	}
	lbit[ ctx ] = cC;
}


long rleDelayedMantisaCTX[ CTXs ] = { 0 }, rleAvgCTX[ CTXs ][2] = { {1},{1} }, s00[ CTXs ] = { 0 };

unsigned long deBitRLE( int b, int levl, int d ) {
#ifdef DEBUG
	bitRLECcnt[ b ]++;
#endif
	unsigned long l = 0, sl = 0, ret = 0;

	// decode mantisa
	while( DeBinRC( 0 + b, rleDelayedMantisaCTX[ b ] + levl + ( ( s00[ b ] & 0xE ) > 0 ) + l++ ) == 0 ) ;
	s00[ b ] = ( ( s00[ b ] << 1 ) + ( l <= 1 ) ) & 0xF;

#ifdef DelayRleBitMantisaCtx
	if( b != 11 ) // delayed context back for mantisa, may hurt compression on small files
		rleDelayedMantisaCTX[ b ] = l << 7;
#endif
	// decode rest of RLE
	for( unsigned int t = 0; t < l; t++ ) ret |= DeBinRC( 0 + b, 40 + levl + ( l << 7 ) + sl++ ) << t;

	return ( ret | ( 1 << l ) ) - 1;
}


// beware, mantisa and rest of RLE is packed into same context to save memory
// it may hurt compression if they mix by an accident
void bitRLE( int b, int p, int levl ) {
#ifdef DEBUG
	bitRLECcnt[ b ]++;
#endif
	/*unsigned*/ long seq = r2[ p ] + 1, l = 0, sl = 0, m;

	if( seq <= 1 ) return; // sequence was of zero length

	m = seq;	 // first encode mantisa
	while( m > 1 ) BinRC( 0 + b, rleDelayedMantisaCTX[ b ] + levl + ( ( s00[ b ] & 0xE ) > 0 ) + l++,
		( ( ( m >>= 1 ) > 1 ) ? 0 : 1 ), false );
	s00[ b ] = ( ( s00[ b ] << 1 ) + ( l <= 1 /*seq <= 4*/ ) ) & 0xF;
	l <<= 7;	 // mantisa will chose context of RLE-bit, according to its length l
	while( seq > 1 ) BinRC( 0 + b, 40 + levl + ( l ) + sl++, (seq & 1), false ), seq >>= 1;

#ifdef DelayRleBitMantisaCtx
	if( b != 11 ) // delayed context back for mantisa, may hurt compression on small files
		rleDelayedMantisaCTX[ b ] = l;
#endif
}

bool decodeStart = true;
long forwr = 0, d1 = 0, d2 = 0, r24 = 0, a3 = 0, a1b = -1;
long ncdiv = -1, scdiv = -1, firstRun = 1, nc2 = -1, c4 = 0;
long nc4 = 0, nc8 = 8, ls0 = -1, o2 = 0, o3[ 0xFF ] = { 0 };
#define INIT_DeTreeRLE firstRun = 1, nc2 = ncdiv = scdiv = -1, c4 = 0, nc4 = 0, \
		nc8 = 8, ls0 = -1, decodeStart = true, forwr = 0, r24 = 0, o2 = 0, a3 = 0; \
		for(int s = 0; s <= 100; s++) o3[ s ] = 0;

void deTreeRLE_R25_Init( ) {
	a1b = 2;
	if( a3 <= 2 ) {
		a1b = ( DeBinRC( 10, 0 ) == 0 ), a3 = -1;
	} else a3 = 0;
	if( a1b > 0 )
		r2[ 5 ] = deBitRLE( 10, 0, 8 ) + ( a3 < 0 );
	else r2[ 5 ] = 1;
	a3 += r2[ 5 ]; // - ( a1b != 2 );
}

void deTreeRLEInit( ) {
	d2 = deBitRLE( 0, 0, 8);
	r2[ 4 ] = 2;
	if( d2 != 1 )
		r2[ 4 ] = 1 + ( ( 1 - DeBinRC( 1, 0 ) ) << 1 );
	r24 = r2[ 4 ];
	if( r2[ 4 ] > 1 )
		r2[ 4 ] += deBitRLE( 1, 0, 8 ) - 2;
}

int deTreeRLE( int mode, int rdivInit = 0 )
{
#ifdef DEBUG
	treeRLEcnt++;
#endif
	long rt = 0, rnk = 0, c2 = 0, b4 = 0, b8, rdiv0, rdiv0F;
	rdiv = 0;

	if( firstRun ) {
		if( r2[ 4 ] == 1 )
		{
			rdiv = decodeStart ? rdivInit : forwr;
			forwr = deBitRLE( 0, 0, 8 );
			r2[ 4 ] = 2;
			if( ( rdiv != 1 && forwr != 1 ) || r24 == 1 )
				r2[ 4 ] = 1 + ( ( 1 - DeBinRC( 1, 0 ) ) << 1 );
			r24 = r2[ 4 ];
			if( r2[ 4 ] > 1 )
				r2[ 4 ] += deBitRLE( 1, 0, 8 ) - 2;
			decodeStart = false;
		} else r2[ 4 ]--;
	}
	firstRun = 0;

	if( ncdiv != -1 ) rdiv = ncdiv;
	if( scdiv == -1 ) {
		if( r2[ 4 ] == 1 ) {
			ncdiv = decodeStart ? rdivInit : forwr;
			forwr = deBitRLE( 0, 0, 8 );
			r2[ 4 ] = 2;
			if( ( ncdiv != 1 && forwr != 1 ) || r24 ==1 )
				r2[ 4 ] = 1 + ( ( 1 - DeBinRC( 1, 0 ) ) << 1 );
			r24 = r2[ 4 ];
			if( r2[ 4 ] > 1 )
				r2[ 4 ] += deBitRLE( 1, 0, 8 ) - 2;
			decodeStart = false;
		}
		else r2[ 4 ]--, ncdiv = 0;
	} else ncdiv = scdiv;

	if( r2[ 4 ] == 1 ) {
		scdiv = decodeStart ? rdivInit : forwr;
		forwr = deBitRLE( 0, 0, 8 );
		r2[ 4 ] = 2;
		if( ( scdiv != 1 && forwr != 1 ) || r24 == 1 )
			r2[ 4 ] = 1 + ( ( 1 - DeBinRC( 1, 0 ) ) << 1 );
		r24 = r2[ 4 ];
		if( r2[ 4 ] > 1 )
		r2[ 4 ] += deBitRLE( 1, 0, 8 ) - 2;
		decodeStart = false;
	} else r2[ 4 ]--, scdiv = 0;

	ls0 -= ( rnk += TREESIZE * rdiv );
	rdiv0 = rdiv; // base value
	rdiv0F = rdiv & 0xF;

	int rdivTresh = rdiv, ncdTresh = ncdiv;
	if( TREETRESH && rdivTresh > TREETRESH ) rdivTresh = TREETRESH;
	if( TREETRESH && ncdTresh > TREETRESH ) ncdTresh = TREETRESH;

	long r = ( rdivTresh << 1 ) + ( ( ncdiv > rdiv ) | ( ls0 >= TREESIZE ) ),
		 ncd = ( ncdTresh << 1 ) + ( ( scdiv > ncdiv ) | ( rdiv > ncdiv ) );
	if( ( ncd ) >= SubCTXs - 1 ) ncd = SubCTXs - 2; // buffer overflow prottection
	if( ( r ) >= SubCTXs - 1 ) r = SubCTXs - 2;

	if( ( r <<= 1 ) >= SubCTXs - 1 ) r = SubCTXs - 2;
	r += ( ( o2 & 0x1E ) != 0 );
	if( ( ncd <<= 1 ) >= SubCTXs - 1 ) ncd = SubCTXs - 2;
	ncd += ( ( ( ( ( o2 << 1 ) + ( rdiv > 0 ) ) & 0x1F ) & 0x1E ) != 0 );

	if( mode ) {
		if( nc2 == -1 ) {
			c2 = ( rdiv == 0 && r2[ 5 ] == 1 ) ||
				  ( rdiv > 0 && DeBinRC( 9, r ) );
			if( rdiv == 0 ) {
				if( c2 ) {
					deTreeRLE_R25_Init( );
				} else r2[ 5 ]--;
			}
		} else c2 = nc2;
		nc2 = ( ncdiv == 0 && r2[ 5 ] == 1 ) ||
				( ncdiv > 0 && DeBinRC( 9, ncd ) );
		if( ncdiv == 0 ) {
			if( nc2 ) {
				deTreeRLE_R25_Init( );
			} else r2[ 5 ]--;
		}
	} else { // mode == 0
		if( nc2 == -1 ) {
			c2 = DeBinRC( 9, r );
		} else c2 = nc2;
		nc2 = DeBinRC( 9, ncd );
	}

	if( c2 )
	{ // rank > 2
		if( c4 == 0 && ( b4 = DeBinRC( 4, r ) ) && ( b8 = DeBinRC( binOut ? 4 : 8, r ) ) ) { };
		if( c4 && ( b4 = nc4 ) ) b8 = nc8;
		if( nc2 && ncdiv == rdiv ) {
			c4 = 1, nc8 = 0, nc4 = DeBinRC( 4, ncd );
			if( nc4 ) nc8 = DeBinRC( binOut ? 4 : 8, ncd );
		} else c4 = 0, nc4 = nc8 = ( ncdiv > rdiv ); // no direct ctx read
		if( TREETRESH && rdiv > TREETRESH ) rdiv = TREETRESH;
		if( ( rdiv <<= 1 ) >= SubCTXs - 1 ) rdiv = SubCTXs - 2; // buffer overflow prottection
		if( b4 == 0 )
			rnk += 3 + DeBinRC( 6, rdiv + ( nc4 | ( ls0 >= 4 ) | ( ( o3[ 60 ] & 0xE ) > 0 ) ) );
		else {
			if( b8 ) rnk += 7 + DeBinRC( 5, r );
			else rnk += 5 + DeBinRC( binOut ? 6 : 7, rdiv + ( nc8 | ( ls0 >= 6 ) |
				( ( o3[ 80 ] & 0xE ) > 0 ) ) );
		}
	} // <= 2
	else {
		if( TREETRESH && rdiv > TREETRESH ) rdiv = TREETRESH;
		if( ( rdiv <<= 1 ) >= SubCTXs - 1 ) rdiv = SubCTXs - 2; // overflow protection
		rdiv += ( ( nc2 & ( ncdiv == rdiv0 ) ) | ( ncdiv > rdiv0 ) | ( ls0 >= 2 ) |
			( ( o3[ rdiv0F ] & 0xE ) > 0 ) );
		if( DeBinRC( 2, rdiv ) ) rnk += 1 + DeBinRC( 3, rdiv );
	}
	o2 = ( ( o2 << 1 ) + ( rdiv0 > 0 ) ) & 0x1F;
	o3[ rdiv0F ] = ( ( o3[ rdiv0F ] << 1 ) + ( ( rt = rnk % TREESIZE ) >= 2 ) ) & 0xF;
	o3[ 60 ] = ( ( o3[ 60 ] << 1 ) + ( rt >= 5 ) ) & 0xF;
	o3[ 80 ] = ( ( o3[ 80 ] << 1 ) + ( rt >= 7 ) ) & 0xF;
	return ls0 = rnk;
}

long pp1 = -1, pz1, pz2, pr1 = -1, ps1 = -1, st1, d0 = -1, ordiv, rdivR2;
void treeRLE( long rnk, int mode )
{
	pp1 = ps1, pz2 = pz1 = rnk, ps1 = rnk = pr1, pr1 = pz1;
	if( rnk < 0 ) return;
#ifdef DEBUG
	treeRLEcnt++;
#endif
	if( ( ordiv = rdiv = rnk / TREESIZE) > 0 )
	{
		r2[ 3 ] = rdiv;
		if( TREETRESH && rdiv > TREETRESH ) rdiv = TREETRESH;

		bitRLE( 0, 3, 0);
		if( ( d2 != 1 && r2[ 3 ] != 1 ) || d0 == 0 )
			BinRC( 1, 0, r2[ 4 ] == 1, false ), r2[ 4 ]--;
		d0 = r2[ 4 ];
		if( r2[ 4 ] > 0 )
			bitRLE( 1, 4, 0 );
		d2 = r2[ 3 ];
		r2[ 4 ] = 1;
		if( ( rdiv <<= 1 ) >= SubCTXs - 1 ) rdiv = SubCTXs - 2; // buffer overflow prottection
		rnk %= TREESIZE;
	} else r2[ 4 ]++;

	pp1 -= ps1 - rnk, pz2 -= ps1 - rnk;
	rdiv += ( st1 = ( ( pz2 >= TREESIZE ) | ( pp1 >= TREESIZE ) ) );

	rdivR2 = ( ( rdiv << 1 ) >= SubCTXs - 1 ) ? SubCTXs - 2 : ( rdiv << 1 ); // buffer overflow prottection
	rdivR2 += ( ( o2 & 0x1E ) != 0 ); // 1E - previous rdiv used by pp1 already

	if( rnk <= 2 ) {
		if( mode == 0 || rdiv > 1 )
			BinRC( 9, rdivR2, 0, false ); // rdiv! ending footer for rc9
		else
			r2[ 5 ]++;
		rdiv = rdiv - st1 + ( ( pz2 >= 3 ) | ( pp1 >= 2 ) | ( ( o3[ ordiv & 0xF ] & 0xE ) > 0 ) );
		if( rnk == 0 ) BinRC( 2, rdiv, 0, false );
		else
			BinRC( 2, rdiv, 1, false ), BinRC( 3, rdiv, rnk == 2, false );
	} else {
		if( mode == 0 || rdiv > 1 )
			BinRC( 9, rdivR2, 1, false ); // rdiv! ending footer for rc9
		else {
			if( a3 <= 2 )
				BinRC( 10, 0, r2[ 5 ] == 1, false ), r2[ 5 ]--;
			a3 = r2[ 5 ];
			if( r2[ 5 ] > 0 )
			bitRLE( 10, 5, 0 );
			r2[ 5 ] = 1;
		}
		if( rnk <= 4 )
			BinRC( 4, rdivR2, 0, false ),
			BinRC( 6, rdiv - st1 + ( ( pz2 >= 5 ) | ( pp1 >= 4 ) |
				( ( o3[ 60 ] & 0xE ) != 0 ) ), rnk == 4, false ); // 5
		else {
			BinRC( 4, rdivR2, 1, false );
			if( rnk <= 6 )
				BinRC( binOut ? 4 : 8, rdivR2, 0, false ),
				BinRC( binOut ? 6 : 7, rdiv - st1 + ( ( pz2 >= 7 ) |
					( pp1 >= 6 ) | ( ( o3[ 80 ] & 0xE ) != 0 ) ), rnk == 6, false ); // 7
			else
				BinRC( binOut ? 4 : 8, rdivR2, 1, false ), BinRC( 5, rdivR2, rnk == 8, false );
		}
	}
	o2 = ( ( o2 << 1 ) + ( rdiv > 1 ) ) & 0x1F;
	o3[ ordiv & 0xF ] = ( ( o3[ ordiv & 0xF ] << 1 ) + ( rnk >= 2 ) ) & 0xF;
	o3[ 60 ] = ( ( o3[ 60 ] << 1 ) + ( rnk >= 5 ) ) & 0xF;
	o3[ 80 ] = ( ( o3[ 80 ] << 1 ) + ( rnk >= 7 ) ) & 0xF;
}

#ifdef DEBUG
int pass = 0;
void stats_debug_pass( ) {
	cout<<"\n--------------------------------------";
	cout<<"\n-------------- pass " << ( pass++ ) << " ----------------";
	cout<<"\n--------------------------------------";
	for(int f=0; f<CTXs; f++)
		cout<<"\n binRCcnt["<<f<<"]: "<< binRCcnt[f];
	cout<<"\n--------------------------------------";
	for(int f=0; f<CTXs; f++)
		cout<<"\n bitRLECcnt["<<f<<"]: "<< bitRLECcnt[f];
	cout<<"\n--------------------------------------";
	cout<<"\n treeRLEcnt: "<< treeRLEcnt;
}
#endif



int decode( char * deFile ) {

	DeBinRCInit( deFile );

	int rankCount = ( sumLen == 0 ? 0 : sf2[ 0 ] ), ouBrlePos = 0; sf2[ 0 ] = 0;
	unsigned char *ouBrle, *ouBrank;
	ouBrle = ( unsigned char *) malloc( sumLen + ( trueMax << 4 ) + 2 );
	if( ouBrle == NULL ) cout << "\n Not enough free memory", exit( 1 );
	memset( ouBrle, 0, sumLen + ( trueMax << 4 ) + 2 );
	ouBrank = ( unsigned char *) malloc( rankCount + trueMax + 2 );
	if( ouBrank == NULL ) cout << "\n Not enough free memory", exit( 1 );
	memset( ouBrank, 0, rankCount + trueMax + 2 );

	if( rankCount != 0 ) {
		if( ( int )( unsigned char ) globalMode == 0 ) {
			unsigned long bulk, upRankCount = 0;

			if( RNKcnt ) {

				unsigned long nextp = rankCount + trueMax + 1 - RNKcnt, nextp2 = 0, lrank = 0, rd0cnt = 0;

				// read rdiv for characters and rle ( rank / TREESIZE )
				// init
				deTreeRLEInit();
				deTreeRLE_R25_Init( );
				for( unsigned long ouBrnkPos = nextp, it = 0; it < (RNKcnt+0); it++) {
					rdiv = binOut ? deTreeRLE( 1, d2 ) + 1 : deTreeRLE( 1, d2 );
					ouBrank[ ++ouBrnkPos ] = rdiv + 1;
				}
				INIT_DeTreeRLE
#ifdef DEBUG
				stats_debug_pass( );
#endif
				deTreeRLEInit(); if( binOut ) r2[ 5 ] = deBitRLE( 10, 0, 8 );
				unsigned long ouBrnkPs = 0;
				for( unsigned long it = 0, lxtr = 0; it < RNKcnt+1; it++) {
					bulk = 2;
					++nextp;
					if( ( ( lrank != 2 && ouBrank[ nextp + 0 ] != 2 ) || lxtr == 1 ) && it < RNKcnt ) {
						lxtr = DeBinRC( 10, 0 );
						bulk = 1 + ( ( lxtr == 0 ) << 1 );
					} else lxtr = 0;
					lrank = ouBrank[ nextp + 0 ];
					if( bulk > 1 ) {
						bulk += deTreeRLE( binOut, d2 ) - 1;
					}

					nextp2 += bulk;

					if( it == RNKcnt ) break;
					upRankCount++;
					while( --bulk > 0 ) ouBrank[ouBrnkPs++] = 0;
					ouBrank[ ouBrnkPs++ ] = ouBrank[ nextp ];
				}
				if( bulk > 9 ) deTreeRLE( binOut, d2 );

				// zero last characters
				memset( ouBrank + ouBrnkPs, 0, rankCount + trueMax + 2 - ouBrnkPs );
				INIT_DeTreeRLE
#ifdef DEBUG
				stats_debug_pass( );
#endif
				RNKcnt = nextp2 - 1;

				long a1d = 0, a1c = -1; // init
				// add sign if rank is larger than 2 or not, when rdiv == 0
				deTreeRLEInit(); if( binOut ) r2[ 5 ] = deBitRLE( 10, 0, 8 );
				unsigned long skip = 0;
				for( unsigned long it = 0, ouBrnkPos = 0; ( it < RNKcnt - upRankCount ); ) {
					a1c = 2;
					if( a1d <= 2 && rd0cnt < RNKcnt2 ) {
						a1c = DeBinRC( 10, 0 );
					}
					if( a1c > 0 )
						it += ( skip = binOut ? deTreeRLE( 1, d2 ) + 1 + ( a1c == 1 ) :
							deTreeRLE( 0, d2 ) + ( a1c == 1 ) );
					else it += ( skip = 1 );
					a1d = skip - ( a1c != 2 );

					if( /*it > RNKcnt - upRankCount ||*/ rd0cnt == RNKcnt2 ) break;
					rd0cnt++;
					while( skip > 0 ) if( ouBrank[ ouBrnkPos++ ] == 0 ) skip--;
					ouBrank[ ouBrnkPos - 1 ] = 1; // mark all symbols larger than 2 and smaller than TREESIZE
				}
				if( skip >= 9 ) deTreeRLE( 0, d2 );
				INIT_DeTreeRLE
#ifdef DEBUG
				stats_debug_pass( );
#endif
				//  compute % TREESIZE for ranks, skip all zeros following escape-symbols(2)
				int ppc = -1, nC2 = -1, nextChar, ncdiv, nextSC, nC4, nC8,
					c4 = 0, b4 = 0, b8 = 0, o0 = 0, rdiv2;
				for( /*unsigned*/ long ouBrnkPos = 0, c2 = 0, l2 = 0, rdiv0, rdiv0F, ncd2;
					( unsigned )ouBrnkPos < RNKcnt; )
				{
					long rank = 0;
					if( ouBrank[ ouBrnkPos ] >= 2 ) {
						rdiv0 = rdiv = ouBrank[ ouBrnkPos ] - 1;
						ppc -= ( rank = rdiv * TREESIZE );
						if( ( rdiv <<= 1 ) >= SubCTXs - 1 ) rdiv = SubCTXs - 2; // buffer overflow prottection
						if( nC2 == -1 ) { // read first character
							l2 = ouBrank[ ouBrnkPos + ( ( unsigned )ouBrnkPos + 1 < RNKcnt ) ],
								l2 = l2 ? l2 - 1 : 0;
							int rdivR2 = rdiv + ( ( l2 > rdiv0 ) | ( ppc >= TREESIZE ) );
							if( ( rdivR2 <<= 1 ) >= SubCTXs - 1 ) rdivR2 = SubCTXs - 2; // buffer overflow prottection
							rdivR2 += ( o0 != 0 );
							c2 = DeBinRC( 9, rdivR2 );
						} else c2 = nC2;
					} else rdiv0 = rdiv = 0, c2 = ouBrank[ ouBrnkPos ];
					rdiv0F = rdiv0 & 0xF;

					nextChar = ( ( unsigned )ouBrnkPos + 1 < RNKcnt );
					nextSC = ( ( unsigned )ouBrnkPos + 2 < RNKcnt ),
					l2 = ouBrank[ ouBrnkPos + nextChar + nextSC ], l2 = l2 ? l2 - 1 : 0;
					if( ( ncdiv = ouBrank[ ouBrnkPos + nextChar ] - 1 ) >= 1 ) {
						ncd2 = ( ( ncdiv << 1 ) >= SubCTXs - 1 ) ? SubCTXs - 2 : ncdiv << 1;
						int ncdivR2 = ncd2 + ( ( l2 > ncdiv ) | ( rdiv0 > ncdiv ) );
						if( ( ncdivR2 <<= 1 ) >= SubCTXs - 1 ) ncdivR2 = SubCTXs - 2; // buffer overflow prottection
						ncdivR2 += ( ( ( ( o0 << 1 ) + ( rdiv0 > 0 ) ) & 0xFF ) != 0 );
						nC2 = nextChar ?
							DeBinRC( 9, ncdivR2 ) : nC2,
						rdiv += ( c2 & ( ( ncdiv >= ouBrank[ ouBrnkPos ] ) | ( ppc >= TREESIZE ) ) );
					} else ncd2 = ncdiv = 0, nC2 = ouBrank[ ouBrnkPos + nextChar ],
						rdiv += ( c2 & ( ppc >= TREESIZE ) );
					rdiv2 = ( ( rdiv << 1 ) >= SubCTXs - 1 ) ? SubCTXs - 2 : rdiv << 1;
					rdiv2 += ( o0 != 0 );

					if( c2 ) { // >= 3
						if( c4 == 0 && ( b4 = DeBinRC( 4, rdiv2 ) ) && ( b8 = DeBinRC( 8, rdiv2 ) ) ) { }; //rdiv2
						if( c4 && ( b4 = nC4 ) ) b8 = nC8;
						if( nC2 && ncdiv == ouBrank[ ouBrnkPos ] - 1 ) {
							c4 = 1, nC8 = 0;
							ncdiv = ncd2 + ( ( l2 > ncdiv ) | ( rdiv0 > ncdiv ) );
							if( ( ncdiv <<= 1 ) >= SubCTXs - 1 ) ncdiv = SubCTXs - 2; // buffer overflow prottection
							ncdiv += ( ( ( ( o0 << 1 ) + ( rdiv0 > 0 ) ) & 0xFF ) != 0 );
							nC4 = nextChar ? DeBinRC( 4, ncdiv ) : 0;
							if( nC4 ) nC8 = nextChar ? DeBinRC( 8, ncdiv ) : nC8;
						  // no direct ctx read
						} else c4 = 0, nC4 = nC8 = ( ncdiv > ouBrank[ ouBrnkPos ] - 1 );
						if( b4 == 0 ) rank += DeBinRC( 6, ( rdiv ^ ( rdiv & 1 ) ) +
							( nC4 | ( ppc >= 4 ) | ( o3[ 20 + rdiv0F ] > 0 ) ) ) + 3;
						else rank += ( b8 ?
							DeBinRC( 5, rdiv2 ) + 7 : //rdiv2
							DeBinRC( 7, ( rdiv ^ ( rdiv & 1 ) ) + ( nC8 | ( ppc >= 6 ) |
								( o3[ 40 + rdiv0F ] > 0 ) ) ) + 5 );
					} else { // for ranks <= 2
						rdiv += ( ( ppc >= 2 ) |
							( ( nC2 & ( ncdiv == rdiv0 ) ) | ( ncdiv > rdiv0 ) ) );
						if( DeBinRC( 2, rdiv ) ) rank += DeBinRC( 3, rdiv ) + 1;
					}
					o0 = ( ( o0 << 1 ) + ( rdiv0 > 0 ) ) & 0xFF;
					o3[ 20 + rdiv0F ] = ( ( o3[ 20 + rdiv0F ] << 1 ) + ( ( rank % TREESIZE ) >= 5 ) ) & 0xFF;
					o3[ 40 + rdiv0F ] = ( ( o3[ 40 + rdiv0F ] << 1 ) + ( ( rank % TREESIZE ) >= 7 ) ) & 0xFF;
					ouBrank[ ouBrnkPos++ ] = ppc = rank;
				}
#ifdef DEBUG
				stats_debug_pass( );
#endif
			} // RNKcnt != 0

			//	join rle and character buffers ouBrle, ouBrank,
			//	if symbol >= RLEdistEst read from first rle buffer, otherwise reconstruct second

			unsigned long uFP = 1, ouBrnkPos = 0; while( sf2[ uFP ] == -1 ) uFP++;
			unsigned long rnkSeq = 0, rnkSeq2 = 0, rleSecCTX = 0, rnkSecCTX = 0, prank = 0;
			unsigned long rnk0Seq = 0, rnk1Seq = 0, rnkSecCTX2 = 0, prank1 = 0;
			ouBrlePos = 0;

			int lr0 = 0, lr1 = 0, lr2 = 0, lr3 = 0, lr4 = 0;
			deTreeRLEInit();
			if( binOut ) r2[ 0 ] = deBitRLE( 7, 0, 8 ), r2[ 5 ] = deBitRLE( 10, 0, 8 );
			for( unsigned int i = 0; i < ( unsigned )rankCount; i++ )
			{
				rdiv = 0, rnkSeq++;
				// read first rle channel
				if( prank + rnkSeq < rleSecCTX ) {
					if( ( rdiv = deTreeRLE( binOut, d2 ) ) ) rnkSeq = -1;
				} else {	// read second rle channel
					if( ( rdiv = binOut ? ( r2[ 0 ] <= 1 ) : DeBinRC( binOut ? 8 : 10, 0 + ( lr0 > 0 ) ) )
						&& DeBinRC( binOut ? 8 : 10, 2 + ( lr1 > 0 ) ) )
						rdiv = 2 + deTreeRLE( binOut, d2 );
					if( binOut && --r2[ 0 ] <= 0 ) r2[ 0 ] = deBitRLE( 7, 0, 8 );
					lr0 = ( ( lr0 << 1 ) + ( rdiv == 0 ) ) & 0xF;
					lr1 = ( ( lr1 << 1 ) + ( rdiv >= 2 ) ) & 0xF;
					rnkSeq = 0, rleSecCTX = ( ( unsigned )rdiv <= RLEdistEst ? rdiv : RLEdistEst );
				}

				//  recover rank buffer according rank channels
				rnkSeq2++;	// read first rank channel
				if( prank >= ( unsigned )TREESIZE || ( rdiv >> 1 ) + rnkSeq2 < rnkSecCTX ) {
					if( ( prank = ouBrank[ ouBrnkPos++ ] ) ) rnkSeq2 = -1;
				} else {	// read second rank channel
					if( ( prank = DeBinRC( 11, 0 + ( lr2 > 0 ) ) ) ) {

						rnk1Seq++;
						if( prank1 >= ( unsigned )TREESIZE || ( rnk0Seq >> 1 ) + rnk1Seq < rnkSecCTX2 ) {
							if( DeBinRC( 11, 2 + ( lr4 > 0 ) ) ) prank = 2 + ouBrank[ ouBrnkPos++ ];
							prank1 = prank - 1;
							if( ( prank1 > 0 ) ) rnk1Seq = -1;
						} else {	// read second rank channel
							if( ( prank1 = DeBinRC( 11, 4 + ( lr3 > 0 ) ) ) && DeBinRC( 11, 6 + ( lr3 > 0 ) ) )
								prank1 = 2 + ouBrank[ ouBrnkPos++ ];
							prank = prank1 + 1;
							lr3 = ( ( lr3 << 1 ) + ( prank >= 2 ) ) & 0xF;
							rnk1Seq = 0, rnkSecCTX2 = ( prank - 1 <= RLEdistEst ? prank - 1 : RLEdistEst );
						}
						rnk0Seq = 0;

					} else rnk0Seq++;
					lr4 = ( ( lr4 << 1 ) + ( prank >= 2 ) ) & 0xF;
					lr2 = ( ( lr2 << 1 ) + ( prank >= 1 ) ) & 0xF;
					rnkSeq2 = 0, rnkSecCTX = ( prank <= RLEdistEst ? prank : RLEdistEst );
				}
				if( rdiv != 0 )
					for( unsigned int n = 0; n < ( unsigned )rdiv; n++ ) ouBrle[ ouBrlePos++ ] = 0;

				if( prank == ( unsigned )trueMax ) {
					sf2[ uFP++ ] = ouBrlePos;
					while( sf2[ uFP ] == -1 ) uFP++;
				} else ouBrle[ ouBrlePos++ ] = prank + 1;
			}
			memset( ouBrle + ouBrlePos, 0, sumLen - ouBrlePos );
			sf2[ uFP ] = sumLen;
#ifdef DEBUG
				stats_debug_pass( );
#endif
		} else if( ( int )globalMode == 1 ) { // " fast " mode decompression
			deTreeRLEInit(); // initial states
			deTreeRLE_R25_Init( );
			long rank = 0; long globalPos, uFP, nextC, rleSum = 0; bool uniq = false;
			for( int m = 0; m < 2; m++ ) {
				ls0 = - 1; // Init DeTreeRLE
				globalPos = 0; uFP = 1; while( sf2[ uFP ] == -1 ) uFP++; nextC = sf2[ uFP ];
				for( int k = 0; k < rankCount; k++ ) {
					rank = deTreeRLE( 1, d2 );

					if( m == 0 ) {
						if( rank != trueMax ) ouBrank[ globalPos++ ] = ( char ) rank + 1;
						if( rank == trueMax || k + 1 == rankCount )
						{
							int started = uFP;
							while( sf2[ uFP ] == -1 || ( uFP == started && uniq ) ) uFP++;
							uniq = true, sf2[ uFP ] = k; // globalPos, dont forget subtract escapes
						}
					}
					if( m == 1 ) {
						rleSum += rank;
						if( ouBrlePos + rank < sumLen ) ouBrlePos += rank;
						if( k >= nextC || k + 1 == rankCount ) { // is there next rank? add RLE sequence length
							sf2[ uFP++ ] += rleSum;
							while( sf2[ uFP ] == -1 ) uFP++;
							nextC = sf2[ uFP ];
							if( k + 1 == rankCount ) ouBrle[ ouBrlePos++ ] = (char)ouBrank[ globalPos++ ];
						} else ouBrle[ ouBrlePos++ ] = ( char )ouBrank[ globalPos++ ];
					}
					rank = 0;
				} // for k
			} // for m
		} // else if globalMode
	} // for cond rankCount != 0

	int norm = 0;
	for( int t = 1; t <= MAXCHAR + 1; t++ ) {
		if( sf2[ t ]!= -1 ) { if( globalMode == 1 ) sf2[ t ] -= norm++; } else sf2[ t ] = sf2[ t - 1 ];
		if( norm == trueMax + 1 ) sf2[ t ] =  sumLen;
	}

	// free unused buffers
	for( int ctx = 0; ctx < CTXs; ctx++ ) {
		free( inBuf[ ctx ] );
	}

	free( ouBrank ), free( sRt ), free( fq );
	if( rankCount ) { // qlfc detransform

		bQlfc = (unsigned char *) malloc( sumLen * 4 );
		if( bQlfc == NULL ) cout << "\n Not enough free memory", exit( 1 );
		long sCount = -1; unsigned long outZ = 0;
		for( int it = 0; it < MAXCHAR + 1; it++ )
		{
			sf[ it ] = sf2[ it ];
			if( sf2[ it ] >= sf2[ it + 1 ] ) pendPos[ it ] = MAXCHAR + 2; // unused chars = MAXCHAR + 2
			else {
				pendPos[ it ] = (unsigned char)ouBrle[ sf2[ it ] ];
				order[ (unsigned char)ouBrle[ sf2[ it ] ] ] = it; // order of symbols..
				sCount++;
			}
		}
		unsigned int ov[ MAXCHAR + 1 ] = { 0 };
		unsigned int * op[ MAXCHAR + 1 ];
		for( int it = 0; it < MAXCHAR + 1; it++ ) op[ it ] = ( unsigned int *) malloc( 4 );

		int in = b[ MAXCHAR + 1 ] = order[ 0 ]; // prepare MTF table start. position
		for( int i = 0; i < sCount; i++, in = b[ in ] ) b[ in ] = order[ i + 1 ];
		for( int it = 0; it < sumLen; it++ )
		{
			if( 0 != ( pendPos[ p0 = p1 = b[ MAXCHAR + 1 ] ] ) ) {
				long rank;
				for( rank = 1, p0 = b[ p1 ]; ; rank++, p1 = p0, p0 = b[ p1 ] )
					if( rank == pendPos[ p0 ] ) break;
				b[ p1 ] = b[ p0 ];
				b[ p0 ] = b[ MAXCHAR + 1 ];
				b[ MAXCHAR + 1 ] = p0; // actual decoded character
			}
			unsigned int relP = sf[ p0 ] - sf2[ p0 ];
			if( relP - ( ov[ p0 ] << 24 ) >= ( 1 << 24 ) ) {
				op[ p0 ] = ( unsigned int *) realloc( op[ p0 ], ++ov[ p0 ] << 2 );
				op[ p0 ][ ov[ p0 ] - 1 ] = outZ >> 2;
			}
			( ( unsigned int *) &bQlfc[ outZ ] )[ 0 ] = ( relP ); // init inverse BWT table, backward &0x00FFFFFF
			outZ += 3; //3; 2;
			if( sf[ p0 ] + 1 < sf2[ p0 + 1 ] )  // move to next pending value..
				pendPos[ p0 ] = ( unsigned char ) ouBrle[ ++sf[ p0 ] ];
			else pendPos[ p0 ] = MAXCHAR + 2;
			bQlfc[ outZ++ ] = p0; // putc
			if( it + 1 == sumLen ) {
				bq = ( unsigned int *) bQlfc;
				outZ >>= 2;
				for( unsigned int td = 0, k = 0; td < outZ; td++, k = 0 ) {
					unsigned char putc = ouBrle[ outZ - td - 1 ] = bQlfc[ ( bwtIndex * 4 ) + 3 ];
					for( ; k < ov[ putc ] && op[ putc ][ k ] <= bwtIndex; k++ ) ;
					bwtIndex = sf2[ putc ] + ( bq[ bwtIndex ] & 0x00FFFFFF ) + ( k << 24 );
					// bwtIndex = sf2[ putc ] + ( (unsigned int *) &bQlfc[ bwtIndex * 5 ] )[ 0 ];
				}
				//inverse_bw_transform( bQlfc, bQlfc, NULL, sumLen, bwtIndex );
				outF->write( (char *) &ouBrle[ 0 ], outZ );
				free( bQlfc );
				for( int c = 0; c < MAXCHAR + 1; c++ ) {
					free( op[ c ] );
				}
				outZ = 0;
			}
		}
	} else if( sumLen ) { // seq. of the same characters only
		bQlfc = (unsigned char *)malloc( ( sumLen + 1 ) );
		if( bQlfc == NULL ) cout << "\n Not enough free memory", exit( 1 );
		for( int i = 0; i < sumLen; i++ ) bQlfc[ i ] = storedChars[ 0 ];
		outF->write( (char *) &bQlfc[ 0 ], sumLen );
		free( bQlfc );
	}

	cout << "Done.\n\n Decoded into: 		" << ( sumLen ) << " character(s).\n";
	cout << " Decoding took:		"<< ((float)(clock() - start)/CLOCKS_PER_SEC ) << " sec\n";

	free( ouBrle );
	outF->close();
	return 0;
}


int cmpBranch( unsigned int b1, unsigned int b2, unsigned int depth = 2 ) {
	b1 += depth, b2 += depth; // skip (3 or 4 byte) radix sort processed data
	while( ( ( unsigned int *) ( data + b1 ) )[ 0 ] ==
	( (unsigned int *) ( data + b2 ) )[ 0 ] ) {
		b1 += 4, b2 += 4;
		if( b1 >= size + 3 ) b1 -= size;
		if( b2 >= size + 3 ) b2 -= size;
	}
	unsigned int p1 = ( ( unsigned int *) ( data + b1 ) )[ 0 ],
				 p2 = ( ( unsigned int *) ( data + b2 ) )[ 0 ];
	return __builtin_bswap32( p1 ) < __builtin_bswap32( p2 ); // __bswap_32, __constant_cpu_to_be32 endian conversion
}

void insertionSort( unsigned int * idx, unsigned int l, unsigned int p,
														unsigned int depth = 2 ) {
	unsigned int j, swp;
	for( unsigned int i = l + 1; i <= p; i++ )
	{
		j = i;
		while( j > l && cmpBranch( idx[ j ], idx[ j - 1 ], depth ) ) {
			swp = idx[ j - 1 ]; // top-down addition probably faster
			idx[ j - 1 ] = idx[ j ];
			idx[ j ] = swp;
			j--;
		}
	}
}

void quickSort( unsigned int * idx, unsigned int l, unsigned int p, unsigned int de = 2 ) {
	if( p - l <= 7 ) { insertionSort( idx, l, p ); return; }
	unsigned int lf = l, ri = p, mid = idx[ ( l + p ) >> 1 ], swp;
	do {
		while( ( idx[ lf ] != mid ) &&
			( 1 == cmpBranch( idx[ lf ], mid, de /*depth of already sorted data*/ ) ) ) lf++;
		while( ( idx[ ri ] != mid ) &&
			( 0 == cmpBranch( idx[ ri ], mid, de ) ) ) ri--; // { if( p - l == 1 ) return; ri--; }
		if( (lf <= ri) ) { swp = idx[ lf ];
			idx[ lf ] = idx[ ri ];
			idx[ ri ] = swp; lf += 1, ri -= 1; }
	} while( lf < ri );
  	if( l < ri ) quickSort( idx, l, ri );
  	if( lf < p ) quickSort( idx, lf, p );
}

#define RADIX_TBL_SIZE 256 * 256 * sizeof( unsigned int )

void rBuckQSort( unsigned int freq, unsigned int * prevIndex, unsigned int depth = 4, unsigned int abs = 0 ) {

	unsigned int * newIndex = (unsigned int *) malloc( freq * sizeof( unsigned int ) ),
				 * cnt = (unsigned int *) malloc( RADIX_TBL_SIZE );

	memset( cnt, 0, RADIX_TBL_SIZE );
	for ( unsigned int r = 0; r < freq; r++ )
		cnt[ __builtin_bswap16( ( (unsigned short *) &data[ prevIndex[ r ] + depth ])[ 0 ] ) ]++;
	int ttlfreq = 0, lfreq;
	for ( unsigned int r = 0; r < RADIX_TBL_SIZE / sizeof( unsigned int ); r++ )
		lfreq = cnt[ r ], cnt[ r ] = ttlfreq, ttlfreq += lfreq;
	for ( unsigned int r = 0; r < freq; r++ )
		newIndex[ cnt[ __builtin_bswap16( ( ( unsigned short *)
		&data[ prevIndex[ r ] + depth ] )[ 0 ] ) ]++ ] = prevIndex[ r ];
	lfreq = 0;
	unsigned int nfreq;
	for ( unsigned int r = 0; r < RADIX_TBL_SIZE / sizeof( unsigned int ); r++ ) {
		if( ( nfreq = cnt[ r ] - lfreq ) ) {
			if( nfreq > 8192 && depth <= 6 ) // stop recursion at given depth 2, 4, etc.
				rBuckQSort( nfreq, newIndex + lfreq, depth + 2, abs + lfreq );
			else {
				if( nfreq > 1 )
					quickSort( newIndex, lfreq, cnt[ r ] - 1, depth + 2 );
				for( unsigned int j = lfreq; j < cnt[ r ]; j++ ) {
					if( newIndex[ j ] ) bq[ abs + j ] = data[ newIndex[ j ] - 1 ];
					else bwtIndex += abs + j, bwtI = 1, bq[ abs + j ] = data[ size - 1 ];
				}
			}
		}
		lfreq = cnt[ r ];
	}
	free( newIndex ), free( cnt );

}

#define BUCKET_QSORT( coreN, from, to ) \
	maxSeg[ coreN ] >>= 8; \
	for( unsigned int b = from; b < to; b += sizeof( unsigned int ) ) { \
		if( ( freq = ( ( ( unsigned int *) &rdxSrt02[ b ] )[ 0 ] ) ) ) { \
			freq -= lstFreq; \
			if( freq > ( maxSeg[ coreN ] ) ) { \
				rBuckQSort( freq, bq + lstFreq, 2, lstFreq ); \
			} else { \
				if( freq > 1 ) \
					quickSort( bq, lstFreq, lstFreq + freq - 1 ); \
				for(unsigned int r = 0, lstFreqr; r < freq; r++) { \
					if( bq[ lstFreqr =  lstFreq + r ] ) bq[ lstFreqr ] = data[ bq[ lstFreqr ] - 1 ]; \
					else bwtIndex += lstFreqr, bwtI = 1, bq[ lstFreqr ] = data[ size - 1 ]; \
				} \
			} \
			lstFreq = freq + lstFreq; \
		} \
	} \

int main( int argc, char** argv ) {
	std::ios_base::sync_with_stdio(false);
	if( argc >= 2 && argv[ 1 ][ 0 ] == 'v' ) {
		cout << VERSION <<
#ifdef Build32bit
		", 32-bit";
#else
		", 64-bit";
#endif
		return 0;
	}
	if( argc != 4 || ( ( argv[ 1 ][ 0 ] != 'c' ) && ( argv[ 1 ][ 0 ] != 'd') ) ) {
		cout << "\n QAD, Qlfc based BWT compressor v" << VERSION << ", (c) Petr Kucera 2014 - 2016\n";
		cout << "\n You are using this software at your own risk! Author does not bear";
		cout << "\n any responsibility for any potential damage concerning its use.";
		cout << "\n\n Usage: * qad 'c'/'cf'/'c2'/'d' input-file output-file ( 'c' - default compress,";
		cout << "\n	'cf' - default mode with experimental filter possibly suiting text data,";
		cout << "\n	'c2' - fast compress, 'd' - to unpack )\n";
		cout << "\n	* If you want to output raw binary data with no entropy coding,";
		cout << "\n	use either a 'cb'/'cfb'/'c2b' switch.\n";
		cout << "\n	* 'v' - print curent program version.\n\n";
		return -1;
	}

//	if( ( argv[ 1 ][ 0 ] == 'c' ) &&
//		( strlen( argv[ 1 ] ) >= 2 && argv[ 1 ][ 1 ] == '2' ) ) globalMode = 1;

	if( argv[ 1 ][ 1 ] != 'f' ) stripFilters = 1;

//	if( argv[ 1 ][ 1 ] == 'b' ||
//		( strlen( argv[ 1 ] ) == 3 && argv[ 1 ][ 2 ] == 'b' ) ) binOut = 1;

	for( int ctx = 0; ctx <= CTXs; ctx++ ) {
		pozInBufBin[ ctx ] =  BUFPOS;
		RCRange[ ctx ] = TOP - 1;
		RCLow[ ctx ] = 0;
		readLow[ ctx ] = ( read64 *) &RCLow[ ctx ];
		readLr[ ctx ]  = ( read64 *) &LR[ ctx ];
	}

	for( i = 0; i <= MAXCHAR; i++ ) // lookup table for computation of the number of 1s in bit history
		bits[ i ] = ( i & 1 ) + ( ( i >> 1 ) & 1 ) + ( ( i >> 2 ) & 1 ) +
		( ( i >> 3 ) & 1 ) + ( ( i >> 4 ) & 1 ) + ( ( i >> 5 ) & 1 ) + ( ( i >> 6 ) & 1 );

	sRt = (unsigned int *) malloc( CTXs * SubCTXs * ( 4 ) * sizeof( unsigned int ) );
	fq = (unsigned int *) malloc( CTXs * SubCTXs * ( alphabet + 1 ) * sizeof( unsigned int ) );

	#pragma omp parallel for
	for( int ctx = 0; ctx < CTXs; ctx++ )
		for( int mode = 0; mode < SubCTXs; mode++ )
			for( int it = 0; it <= alphabet; it++ )
	{ fq [ p( ctx, mode, it ) ] = it + 1; /*sRt[ctx][mode][it&3] = 0;*/ } // set freq. counters

	r2[ 4 ] = r2[ 5 ] = r2[ 0 ] = r2[ 1 ] = r2[ 2 ] = 1; // initialize rle counters

	outF = new ofstream( argv[ 3 ], ios::binary );
	FillInAritfq();

	if( argv[ 1 ][ 0 ] == 'd' ) {

		cout << "\n Decompressing (" << argv[ 2 ] << " -> " << argv[ 3 ] << ").. ";

		decode( argv[ 2 ] );
		//DeBinRCInit( argv[ 2 ] );

		return 0;
	}

	FILE * file = fopen( argv[ 2 ], "rb" );

  if( /*file.is_open()*/ file != NULL )
  {
	fseek( file, 0, SEEK_END );
	size = ( unsigned long ) ftell( file );
	rewind( file ); // file.seekg( 0, ios::beg );

	cout << "\n Compressing (" << argv[2] << " -> " << argv[3] << ").. ";

	data = (unsigned char*) malloc( size + 17 );
	if( data == NULL ) cout << "\n Not enough free memory", exit( 1 );

	size = ( unsigned long ) fread( ( char *) data, 1, size, file );

	for (long l = 0; l <= 16; l++) { // trailingcharacters for beter alignment
		data[ size + l ] = data[ l ];
	}
	unsigned char * rdxSrt02 = (unsigned char*) malloc( RADIX_TBL_SIZE );
	memset( rdxSrt02, 0, RADIX_TBL_SIZE );
	unsigned int * rdxSrt = ( unsigned int *) rdxSrt02;

	#pragma omp parallel for
	for( unsigned int t = 0; t < size; t++ ) {
		#pragma omp atomic
		rdxSrt[ __builtin_bswap16( ( (unsigned short *) &data[ t ] )[ 0 ] ) ]++;
	}

	// prepare paralel sections borders and statistics for the 3 byte radix sort
	unsigned int freq, total = 0, procs = omp_get_num_procs(),
		parSec[ ( procs + 1 ) * 10 ][ 2 ], maxSeg[ procs * 10 ];
	memset( parSec, -1, ( procs + 1 ) * 2 * sizeof( unsigned int ) * 10 );
	memset( maxSeg, 0, procs * sizeof( unsigned int ) * 10 );
	for( unsigned int j = 0; j < RADIX_TBL_SIZE; j += 4 ) {
		if( ( freq = ( (unsigned int *) &rdxSrt02[ j ] )[ 0 ] ) ) {
			( (unsigned int *) &rdxSrt02[ j ] )[ 0 ] = total;
			for( int k = ( procs * 10 ) - 1; k > 0; k-- )
				if( total > ( ( size / ( procs * 10 ) ) * k ) ) {
					if( total < parSec[ k ][ 1 ] ) parSec[ k ][ 1 ] = total, parSec[ k ][ 0 ] = j;
					if( freq > maxSeg[ k ] ) maxSeg[ k ] = freq;
				}
			if( total <= ( size / ( procs * 10 ) ) && freq > maxSeg[ 0 ] ) maxSeg[ 0 ] = freq;
			total += freq;
		}
	}

	// failed to split symbol distribution on sections for paralel sorting, 1 thread not optimal
	parSec[ procs * 10 ][ 1 ] = size;
	parSec[ procs * 10 ][ 0 ] = RADIX_TBL_SIZE, parSec[ 0 ][ 0 ] = parSec[ 0 ][ 1 ] = 0;
	if( procs > 1 ) {
		if( ( signed )parSec[ procs - 1 ][ 1 ] == -1 )
			parSec[ procs - 1 ][ 1 ] = size, parSec[ procs - 1 ][ 0 ] = RADIX_TBL_SIZE;
		for( unsigned int s = 1; s <= procs - 1; s++ ) // for all sections
			if( ( signed )parSec[ s ][ 1 ] == -1 || parSec[ 1 ][ 0 ] == RADIX_TBL_SIZE )
				parSec[ s ][ 1 ] = parSec[ procs - 1 ][ 1 ], parSec[ s ][ 0 ] = RADIX_TBL_SIZE;
	}

	unsigned int rollStat;
	FILE * outTmpF = fopen( strcat( argv[ 3 ], ".tmp" ), "w+b" );
	bQlfc = ( unsigned char *) malloc ( 1 );
	for( int rr = 0; rr < 10; rr++ )
	{	//skip unused chars
		if( ( parSec[ ( ( rr + 1) * procs ) + 1 ][ 1 ] - parSec[ rr * procs ][ 1 ] ) == 0 ) continue;
		rollStat = parSec[ rr * procs ][ 1 ];
		unsigned int mx = parSec[ (rr + 1) * procs ][ 1 ] - parSec[ rr * procs ][ 1 ],
			fr = parSec[ rr * procs ][ 0 ] >> 2, to = parSec[ ( rr + 1 ) * procs ][ 0 ] >> 2;

		bq = ( unsigned int *) ( bQlfc = ( unsigned char *) realloc( bQlfc, mx << 2 ) );
		// fill rank sort indexes on positions in out. buffer
		for( unsigned int t = 0; t < size; t++ ) {
			unsigned short pos = __builtin_bswap16( ( ( unsigned short *) &data[ t ] )[ 0 ] );
			if( pos >= fr && pos < to )
				bq[ rdxSrt[ pos ]++ - rollStat ] = t;
		}

		for( unsigned int t = parSec[ rr * procs ][ 0 ], r = parSec[ ( rr + 1 ) * procs ][ 0 ];
						  t < r; t += 4 ) {
			if( ( ( unsigned int *) &rdxSrt02[ t ] )[ 0 ] ) {
				( ( unsigned int *) &rdxSrt02[ t ] )[ 0 ] -= rollStat;
			}
		}

		for( unsigned int t = 0; t < procs; t++ )
			parSec[ ( rr * procs ) + t ][ 1 ] -= rollStat;

		// quick sorts running in paralel on sections of unsorted data after 3-byte radix sort
		#pragma omp parallel
		{
			for( unsigned int r = rr * procs; r < ( rr + 1 ) * procs; r++ )
				#pragma omp sections nowait
				{
					#pragma omp section
					{
						unsigned int lstFreq = parSec[ r ][ 1 ], freq;
						//if( parSec[ r + 1 ][ 1 ] > parSec[ r ][ 1 ] )
						BUCKET_QSORT( r, parSec[ r ][ 0 ], parSec[ r + 1 ][ 0 ] )
					}
				}
		}

		for( unsigned int r = 0; r < mx; r++ )
			bQlfc[ r ] = ( unsigned char )bq[ r ];
		fwrite( (char*) &bQlfc[ 0 ], 1, mx, outTmpF );

		if( !bwtI ) bwtIndex += mx;
	}

	bQlfc = ( unsigned char *) realloc( bQlfc, /*( size < 100 ? 100 :*/ size * 4 /*)*/ );
	if( bQlfc == NULL ) cout << "\n Not enough free memory", exit( 1 );

	for( int oBufs = 0; oBufs < CTXs; oBufs++ )
		outBuf[ oBufs ] = (char *) &bQlfc[ size + ( ( ( size < 100 ? 100 : size ) * 3 * oBufs ) / CTXs) ];
#ifdef DEBUG
	cout << "\n  Sorting in:		"<< ((float)(clock() - start)/CLOCKS_PER_SEC ) << " sec. ";
#endif
	for( i = 0; i <= MAXCHAR + 1; sf2[ i ] = sc[ i ] = pendPos[ i ] = sf[ i ] = 0, i++ );
	for( i = 0; i <= MAXCHAR; tt[ i++ ] = 0, ( i > 0 ? b[ i - 1 ] = i : 1 ) );
	b[ MAXCHAR ] = b[ MAXCHAR + 1 ] = 0;

	unsigned int chunk = size / 100;
	rewind( outTmpF );

	for( unsigned int it = 0; it < size; ) {
		chunk = ( unsigned long ) fread( (char *) data, 1, chunk, outTmpF );
		for( unsigned int it2 = 0; it2 < chunk; it2++ ) {
			aC = data[ it2 ];
			sf[ aC + 1 ]++, sf2[ aC + 1 ]++, sc[ aC ]++;
		}
		it += chunk;
	}

	rewind( outTmpF );
	for( int it = 1; it <= MAXCHAR + 2; it++ ) {  // prepare sfuency table
		sf[ it ] += sf[ it - 1 ];
		sf2[ it ] = sf[ it ];
	}

	for( unsigned int it = 0; it < size; ) { // qlfc loop
		chunk = ( unsigned long ) fread( (char *) data, 1, chunk, outTmpF );
		for( unsigned int it2 = 0; it2 < chunk; it2++ ) {
			aC = data[ it2 ];
			inFileOffst++; long rank = 0;

			if( ( p1 = b[ MAXCHAR + 1 ] ) != aC ) {
				if( tt[ aC ] < pOffset ) {
					rank += treshold++;
					p1 = pivot;
				} else if( tt[ aC ] == pOffset ) {
					pOffset = inFileOffst;
					pivot = aC;
					treshold = 0;
				}

				while( true ) {  // passing symbol list
					if( ( ( p0 = b[ p1 ] ) ^ aC ) == 0 )
						{ rank++; b[ p1 ] = b[ aC ]; break; }
					if( ( ( p1 = b[ p0 ] ) ^ aC ) == 0 )
						{ rank += 2; b[ p0 ] = b[ aC ]; break; }
					if( ( ( p0 = b[ p1 ] ) ^ aC ) == 0 )
						{ rank += 3; b[ p1 ] = b[ aC ]; break; }
					if( ( ( p1 = b[ p0 ] ) ^ aC ) == 0 )
						{ rank += 4; b[ p0 ] = b[ aC ]; break; }
					rank += 4;
				}
				b[ aC ] = b[ MAXCHAR + 1 ];
				b[ MAXCHAR + 1 ] = aC;
				tt[ aC ] = inFileOffst;
			}
			if( pendPos[ aC ] == 0 ) // first characters
			{
				pendPos[ aC ]++;
				rank = pendPos[ MAXCHAR + 1 ]++;
				if( rank > trueMax ) trueMax = rank;
			}

			if( globalMode == 0 ) {
				if( rank == 0 ) RLEfreq[aC]++;
				else {
					if( RLEfreq[ aC ] < 4 && RLEfreq[ aC ] ) RLEtop[ RLEfreq[ aC ] ]++;
					if( RLEfreq[ aC ] == 0 ) RLEtop[ 0 ]++;
					RLEfreq[ aC ] = 0;
				}
			}

			bQlfc[ sf[ aC ]++ ] = rank;
		}
		it += chunk;
	} // qlfc loop
	rewind( outTmpF );
	fwrite( ( char *) bQlfc, 1, size, outTmpF );

	if( globalMode == 0 ) { // if we are in default mode, prepare RLE dist. estimates
		for( int i = 0; i < 4; i++ )
			if( RLEtop[ i ] == 0 ) {
				RLEdistEst = TREESIZE + 1;
				break;
			}

		if( RLEtop[ 1 ] && RLEtop[ 2 ] && RLEtop[ 3 ] ) { // beware division by zero
			RLEdistEst = ( RLEtop[ 0 ] / RLEtop[ 1 ] )  - ( ( RLEtop[ 1 ] / RLEtop[ 2 ] ) *
							 ( RLEtop[ 2 ] / RLEtop[ 3 ] ) );
			if( RLEdistEst > ( unsigned )TREESIZE  || RLEtop[ 0 ] / RLEtop[ 1 ] <= 2 ||
			  RLEtop[ 0 ] / RLEtop[ 1 ] >= 6 ) RLEdistEst = TREESIZE; // present top decoder lim., TREESIZE
		}
		if( RLEdistEst < 4 ) RLEdistEst = 4;  // this is not ideal for actual context selection top limit yet,
	}

	if( stripFilters ) RLEdistEst = TREESIZE;

	sf2[ 0 ] = 0; // clean rank counter
	long pp = -1, pz = -1, pr, ps = -1, st, pz2, rdivo;
	long poctyE[ 2 ] = { 0 }, a0 = -1, a1 = 0, o1 = 0, rdiv2;
	unsigned int ra0, rleSecCTX, rleSecCTX2, rnkSeq, prank;

	// currently processing previous rank, instead of actual, first line
#define COMPUTE_CUSTOM_TREE_RLE \
	pp = ps, pz2 = pz = rank, ps = rank = pr, pr = pz; \
	if( ( signed )rank == -1 ) { rleSeq = 0; continue; } \
	if( ( rdivo = rdiv = rank / TREESIZE ) > 0 ) \
	{ \
		if( l == 0 ) treeRLE( binOut ? rdiv - 1 : rdiv /*-1*/, 1 ), poctyE[ 0 ]++; \
		if( l == 1 ) { \
			if( ( d1 != 1 && rdiv != 1 ) || a0 == 0 ) \
				BinRC( 10, 0, r2[ 2 ] == 1, false ), r2[ 2 ]--; \
			a0 = r2[ 2 ], d1 = rdiv; \
			if( r2[ 2 ] > 0 ) \
				treeRLE( binOut ? r2[ 2 ] - 1 : r2[ 2 ] - 1, binOut ); \
			r2[ 2 ] = 1; \
		} \
		if( l > 2 ) rank %= TREESIZE; \
	} else if( l == 1 ) r2[ 2 ]++; \
	if( l == 2 && rdiv == 0 ) { \
		if( rank <= 2 ) r2[ 1 ]++; \
		else { \
			if( a1 <= 2 ) \
				BinRC( 10, 0, r2[ 1 ] > 1, false ), r2[ 1 ]--; \
			a1 = r2[ 1 ]; \
			if( r2[ 1 ] > 0 ) \
				treeRLE( binOut ? r2[ 1 ] - 1 : r2[ 1 ], binOut ); \
			r2[ 1 ] = 1; poctyE[ 1 ]++; \
		} \
	} \
	if( l == 3 ) { \
		rdivo &= 0xF; \
		pp -= ps - rank, pz2 -= ps - rank; \
		if( ( rdiv <<= 1 ) >= SubCTXs - 1 ) rdiv = SubCTXs - 2; \
		rdiv += ( st = ( ( pz2 >= TREESIZE ) | ( pp >= TREESIZE ) ) ); \
		rdiv2 = ( ( rdiv << 1 ) >= SubCTXs - 1 ) ? SubCTXs - 2 : ( rdiv << 1 ); \
		rdiv2 += ( o1 != 0 ); \
		if( rank <= 2 ) { \
			if( rdiv > 1 ) BinRC( 9, rdiv2, 0, false ); \
			rdiv += ( ( pz2 >= 3 ) | ( pp >= 2 ) ) - st; \
			if( rank == 0 ) BinRC( 2, rdiv, 0, false ); \
			else BinRC( 2, rdiv, 1, false ), BinRC( 3, rdiv, rank == 2, false ); \
		} else { \
			if( rdiv > 1 ) BinRC( 9, rdiv2, 1, false ); \
			if( rank <= 4) BinRC( 4, rdiv2, 0, false ), \
				BinRC( 6, rdiv - st + ( ( pz2 >= 5 ) | \
					( pp >= 4 ) | ( o3[ 20 + rdivo ] != 0 ) ), rank == 4, false ); \
			else { \
				BinRC( 4, rdiv2, 1, false ); \
				if( rank <= 6) BinRC( 8, rdiv2, 0, false ), \
					BinRC( 7, rdiv - st + ( ( pz2 >= 7 ) | \
						( pp >= 6 ) | ( o3[ 40 + rdivo ] != 0 ) ), rank == 6, false ); \
				else BinRC( 8, rdiv2, 1, false ), BinRC( 5, rdiv2, rank == 8, false ); \
			} \
		} \
		o1 = ( ( o1 << 1 ) + ( rdiv > 1 ) ) & 0xFF; \
		o3[ 20 + rdivo ] = ( ( o3[ 20 + rdivo ] << 1 ) + ( rank >= 5 ) ) & 0xFF; \
		o3[ 40 + rdivo ] = ( ( o3[ 40 + rdivo ] << 1 ) + ( rank >= 7 ) ) & 0xFF; \
	} // if l == 3

	for( int l = 0; l <= ( globalMode ? 1 : 4 ); l++ ) {
		 // initial values for each pass
		int lr0 = 0, lr1 = 0, lr = 0, lr4 = 0; d0 = ps1 = ps = -1;
		d1 = 0, d2 = 0, pr1 = pr = -1, ra0 = 1; while( !sf2[ ra0 ] ) ra0++;
		rleSecCTX2 = rleSecCTX = rleSeq = rnkSeq = prank = o2 = a3 = 0;

		for(int s = 0; s <= 100; s++) o3[ s ] = 0;
		unsigned int rnk0Seq = 0, rnk1Seq = 0, prank1 = 0;
		rewind( outTmpF );
		chunk = size / 100;
		for( unsigned int it = 0; it < size; /*it++*/ ) {
			chunk = ( unsigned long ) fread( (char *) data, 1, chunk, outTmpF );

			for( unsigned int it2 = 0; it2 < chunk; it2++ ) {
			long rank = data[ it2 ];

			if( it + it2 == ( unsigned )sf2[ ra0 ] )
				{ ra0++; while( sf2[ ra0 ] == sf2[ ra0 - 1 ] ) ra0++; it2--; rank = trueMax + 1; }

			if( rank > 0 ) {
				rank -= 1, rdiv = 0, rnkSeq++;
				if( l == 0 ) sf2[ 0 ]++;

				if( globalMode == 1 ) { // fast compression
					if( l == 0 ) treeRLE( rank, 1 );
					if( l == 1 ) treeRLE( rleSeq, 1 );
				} else {
					if( l != 4 ) {
						int skip = 0;
						if( prank < ( unsigned )TREESIZE && ( rleSeq >> 1 ) + rnkSeq >= rleSecCTX )
						{ // if x + rnkSeq >= rleSecCTX put rank into secondary channel
							if( rank >= 1 ) {
								if( l == 0 ) BinRC( 11, 0 + ( lr > 0 ), 1, false );
								rnk1Seq++;

								if( prank1 < ( unsigned )TREESIZE && ( rnk0Seq >> 1 ) + rnk1Seq >= rleSecCTX2 )
								{
									if( rank >= 2 ) {
										if( l == 0 ) BinRC( 11, 4 + ( lr4 > 0 ), 1, false );
										if( rank == 2 ) { if( l == 0 ) BinRC( 11, 6 + ( lr4 > 0 ),
											0, false ); skip = 1; }
										else{ if( l == 0 ) BinRC( 11, 6 + ( lr4 > 0 ), 1, false ),
											RNKcnt++; skip = 3; }
									} else { if( l == 0 ) BinRC( 11, 4 + ( lr4 > 0 ), 0, false ); skip = 1; }
									if( l == 0 ) lr4 = ( ( lr4 << 1 ) + ( rank >= 2 ) ) & 0xF;
									rnk1Seq = 0, rleSecCTX2 = ( ( unsigned )rank - 1 <= RLEdistEst ? rank - 1 : RLEdistEst );
								} else {
									rnk1Seq = ( rank - 1 > 0 ) ? -1 : rnk1Seq;
									if( rank == 1 ) { if( l == 0 ) BinRC( 11, 2 + ( lr1 > 0 ), 0, false ); skip = 1; }
									else{ if( l == 0 ) BinRC( 11, 2 + ( lr1 > 0 ), 1, false ), RNKcnt++; skip = 2; }
								}
								prank1 = rank - 1 /*+ ( rank >= 1 )*/;

								rnk0Seq=0;
							} else { if( l == 0 ) BinRC( 11, 0 + ( lr > 0 ), 0, false ); skip = 1; rnk0Seq++; }
							if( l == 0 ) lr1 = ( ( lr1 << 1 ) + ( rank >= 2 ) ) & 0xF;
							if( l == 0 ) lr = ( ( lr << 1 ) + ( rank >= 1 ) ) & 0xF;
							rnkSeq = 0, rleSecCTX = ( ( unsigned )rank <= RLEdistEst ? rank : RLEdistEst );
						} else rnkSeq = rank ? -1 : rnkSeq, RNKcnt = ( ( !l ) ? RNKcnt + 1: RNKcnt );
						prank = rank; //+ ( rank >= 1 );
						rleSeq = 0;
						if( skip == 1 ) continue;
						rank -= skip;

					} else {
						if( prank + rnkSeq >= rleSecCTX )
						{ // if previous rank + rnkSeq >= rleSecCTX put RLE to secondary channel
							if( rleSeq >= 1 ) {
								if( binOut ) bitRLE( 7, 0, 0 ), r2[ 0 ] = 1;
								else BinRC( binOut ? 8 : 10, 0 + ( lr0 > 0 ), 1, false );
								if( rleSeq == 1 ) BinRC( binOut ? 8 : 10, 2 + ( lr > 0 ), 0, false );
								else BinRC( binOut ? 8 : 10, 2 + ( lr > 0 ), 1, false ),
									treeRLE( rleSeq - 2, binOut );
							} else { if( binOut ) r2[ 0 ]++;
									else BinRC( binOut ? 8 : 10, 0 + ( lr0 > 0 ), 0, false ); }
							lr0 = ( ( lr0 << 1 ) + ( rleSeq == 0 ) ) & 0xF;
							lr = ( ( lr << 1 ) + ( rleSeq >= 2 ) ) & 0xF;
							rnkSeq = 0, rleSecCTX = ( rleSeq <= RLEdistEst ? rleSeq : RLEdistEst );
						} else {
							treeRLE( rleSeq, binOut ); rnkSeq = rleSeq ? -1 : rnkSeq;
						}
						rleSeq = 0, prank = rank; continue;
					}
					COMPUTE_CUSTOM_TREE_RLE
				} // if globalmode

				rleSeq = 0;
			} else rleSeq++;
			} // for it2
			it += chunk;
		} // for it
		// write out last characters
		if( !globalMode ) {
			if( l != 4 ) { long rank = pz; COMPUTE_CUSTOM_TREE_RLE }
			if( l == 0 || l == 4 ) treeRLE( TREESIZE, l == 0 || binOut );
			// flush RLE counters to be coded at the end of the block
			// this is buggy and redundand, lenghts can be recomputed
			if( l <= 4 && l != 3 ) {
				if( l == 2 && r2[ 1 ] >  1 ) {
					treeRLE( r2[ 1 ], binOut );
					if( r2[ 1 ] >= 9 ) treeRLE( 0, binOut );
				}
				if( l == 4 && r2[ 2 ] >= 1 ) {
					treeRLE( r2[ 2 ], binOut );
					//if( r2[ 2 ] >= 9 ) treeRLE( 0, binOut );
				}
				if( l == 1 && r2[ 2 ] >= 1 ) {
					treeRLE( r2[ 2 ] - 1, binOut );
					if( r2[ 2 ] > 9 ) treeRLE( 0, binOut );
				}

				if( l == 2 || l == 1 ) treeRLE( TREESIZE, binOut );
				if( binOut && ( l == 1 || l == 2 || l == 4 ) ) bitRLE( 10, 5, 0 );

				if( r2[ 5 ] >= 1 && ( l == 0 /*|| && l == 4*/ ) ) {
					if( a3 <= 2 ) BinRC( 10, 0, r2[ 5 ] == 1, false ), r2[ 5 ]--;
					if( r2[ 5 ] > 0 ) bitRLE( 10, 5, 0 );
				}
				if( binOut && l == 4 ) bitRLE( 7, 0, 0 );

				if( r2[ 4 ] >= 1 ) {
					r2[ 3 ] = 1, bitRLE( 0, 3, 0 );
					if( d0 == 0 ) BinRC( 1, 0, r2[ 4 ] == 1, false ), r2[ 4 ]--;
					d0 = r2[ 4 ];
					if( r2[ 4 ] > 0 ) bitRLE( 1, 4, 0 );
				}
				if( l != 0 || poctyE[ 0 ] ) { // upranks?
					// mark end of sequence symbols states
					r2[ 3 ] = 1;	 // mark rdiv
					bitRLE( 0, 3, 0 ), bitRLE( 0, 3, 0 );
					r2[ 4 ] = 1;	 // above TREESIZE
					if( d0 == 0 ) BinRC( 1, 0, r2[ 4 ] == 1, false ), r2[ 4 ]--;
					d0 = r2[ 4 ];
					if( r2[ 4 ] > 0 ) bitRLE( 1, 4, 0 );
					r2[ 4 ] = 1; // reset r2[ 4 ]
					if( d0 == 0 ) BinRC( 1, 0, r2[ 4 ] == 1, false ), r2[ 4 ]--;
					d0 = r2[ 4 ];
					if( r2[ 4 ] > 0 ) bitRLE( 1, 4, 0 ); // rdiv >> 1 = rdiv of previous character
					BinRC( 9, ( ( ( ( TREESIZE / TREESIZE ) << 1 ) +
						( ( rdiv >> 1 ) > 1 ) ) << 1 ) + ( ( o2 & 0x1E ) != 0 ), 0, false );
				}
				// finalize rle
				r2[ 1 ] = r2[ 2 ] = 1;
				r2[ 4 ] = r2[ 5 ] = 1;
				rleSeq = 0;
			}
		} else treeRLE( pz1, 1 );
#ifdef DEBUG
		stats_debug_pass( );
#endif
	} // for l block

	if( globalMode == 1 ) { // special endings..
		bitRLE( 1, 4, 0 );	// finalize RLE
		bitRLE( 10, 5, 0 );
	}

	// flush encoders
	for( int k = 0; k < CTXs; k++ ) {
		if( isCTXUsed[ k ] != 0 ) {
#ifdef DEBUG
			cout << "\n Finalizing CTX: " << k << " at " << pos[ k ] << " bytes";
#endif
			long len = pos[ k ];
			rdiv = rdiv > 3 ? 3 : 0;
			while( true ) { // 2 STOP bytes for each ctx probably not enough
				BinRC( k, rdiv, !lbit[ k ], false );
				if( pos[ k ] - len >= 3 ) { BinRC( k, rdiv, lbit[ k ], true ); break; }
				BinRC( k, rdiv, lbit[ k ], false );
				if( pos[ k ] - len >= 3 ) { BinRC( k, rdiv, !lbit[ k ], true ); break; }
			}
		}
	}
	outF->write( ( char *) &bwtIndex, sizeof( int ) );
	outF->write( ( char *) &sf2[ 0 ],sizeof( int ) );
	outF->write( ( char *) &sf2[ MAXCHAR + 1 ], sizeof( int ) ); // extra size

	int ctxMark = 0, pairs = 0;
	for( int ctx = 0; ctx < CTXs; ctx++ ) {	// rc. pointers (same buffer)
		if( isCTXUsed[ ctx ] != 0 ) {
			ctxMark |= 1 << ctx;
			outF->write( ( char *) &pos[ ctx ], sizeof( int ) );
		}
	}

	if( trueMax < MAXCHAR ) {	// append used / unused characters
		if( trueMax <= 127 ) {
			for( int ctx = 0; ctx <= MAXCHAR; ctx++ ) {
				if( ( sc[ ctx ] != 0 ) && ( ctx == 0 || ctx == MAXCHAR ||
						sc[ ctx - 1 ] == 0 || sc[ ctx + 1] == 0 ) ) {
					if( ( ctx == 0 && sc[ ctx + 1 ] != 0 ) || ( ctx > 0 && ctx < MAXCHAR &&
						sc[ ctx + 1 ] != 0 && sc[ ctx - 1 ] == 0 ) ) pairs++;
					if( ( ctx == 0 && sc[ ctx + 1 ] == 0) || ( ctx == MAXCHAR && sc[ ctx - 1 ] == 0 ) ||
						( ctx > 0 && ctx < MAXCHAR && sc[ ctx + 1 ] == 0 && sc[ ctx - 1 ] == 0 ) ) {
						outF->write( (char*) &ctx, sizeof( char ) );
					}
				}
			}
			for( int ctx = 0; ctx <= MAXCHAR; ctx++ ) {	// print intervals only..
				if( ( sc[ ctx ] != 0 ) && ( ctx == 0 || ctx == MAXCHAR ||
					 ( sc[ ctx - 1 ] == 0 && sc[ ctx + 1 ] != 0 ) ||
					 ( sc[ ctx + 1 ] == 0 && sc[ ctx - 1 ] != 0 ) ) ) {
					if( ( ctx == 0 && sc[ ctx + 1 ] == 0 ) || ( ctx == MAXCHAR && sc[ ctx - 1 ] == 0) ) ; else {
						outF->write( (char*) &ctx, sizeof( char ) );
					}
				}
			}
		} else { // trueMax > 127, reverse appending logic
			for( int ctx = 0; ctx <= MAXCHAR; ctx++ ) {
				if( ( sc[ ctx ] == 0 ) && ( ctx == 0 || ctx == MAXCHAR ||
					sc[ ctx - 1 ] != 0 || sc[ ctx + 1 ] != 0 ) ){
					if( ( ctx == 0 && sc[ ctx + 1 ] == 0) || ( ctx > 0 &&
						ctx < MAXCHAR && sc[ ctx + 1] == 0 && sc[ ctx - 1 ] != 0 ) ) pairs++;
					if( ( ctx == 0 && sc[ ctx + 1 ] != 0 ) || ( ctx == MAXCHAR && sc[ ctx - 1] != 0 ) ||
						 ( ctx > 0 && ctx < MAXCHAR && sc[ ctx + 1 ] != 0 && sc[ ctx - 1 ] != 0 ) ) {
						outF->write( (char*) &ctx, sizeof( char ) );
					}
				}
			}
			for( int ctx = 0; ctx <= MAXCHAR; ctx++ ) {	// print intervals only..
				if( ( sc[ ctx ] == 0 ) && ( ctx == 0 || ctx == MAXCHAR || ( sc[ ctx -1 ] != 0 &&
					sc[ ctx + 1 ] == 0 ) || ( sc[ ctx + 1 ] != 0 && sc[ ctx - 1 ] == 0 ) ) ) {
					if( ( ctx == 0 && sc[ ctx + 1 ] != 0 ) || ( ctx == MAXCHAR && sc[ ctx - 1 ] != 0 ) ) ; else {
						outF->write( (char*) &ctx, sizeof( char ) );
					}
				}
			}
		}
	}
	if( trueMax > 0 && trueMax < MAXCHAR - 1 ) outF->write( (char*) &pairs, sizeof( char ) );

	if( trueMax < 10 || binOut ) outF->write( (char*) &ctxMark, sizeof( char ) * 2 ); // mark used contexts

	// secondary RLE context if rank >= RLEdistEst, in default mode
	if( globalMode == 0 ) outF->write( (char*) &poctyE[ 0 ], sizeof( unsigned int ) );
	if( globalMode == 0 ) outF->write( (char*) &poctyE[ 1 ], sizeof( unsigned int ) );
	outF->write( (char*) &trueMax, sizeof( char ) ); // append number of used characters through input file

	globalMode += ( binOut ) << 1; // Save compressed or raw binary data
	globalMode += ( RLEdistEst ) << 4; // Save RLE distribution estimate
	outF->write( (char*) &globalMode, sizeof( char ) ); // append number of used characters through input file

	outF->close();

	argv[ 3 ][ strlen( argv[ 3 ] ) - 4 ] = 0; // remove .tmp
	ifstream dfiles( argv[ 3 ], ios::ate );
	dfiles.seekg( 0, std::ios::end );
	cout << "Done.\n\n  Input Size:		" << size << " character(s)";
	cout << "\n  Compressed into:	" << dfiles.tellg() << "; " <<
		( ( dfiles.tellg() / ( ( float ) ( size ? size : 1 ) ) ) * 8 ) << " bpc\n";
	cout << "  Compression took:	"<< ((float)(clock() - start)/CLOCKS_PER_SEC ) << " sec\n";
	remove( strcat( argv[ 3 ], ".tmp" ) );
	free( sRt );
	free( fq );
	free( bQlfc );
	free( data );
  }

  return 0;
}