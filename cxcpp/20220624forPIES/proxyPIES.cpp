#include <stdlib.h>
#include <stdio.h>
#include <string.h> // For memset?!?!

#define REAL_TYPE double
#define REAL_TYPE_FORMAT_STR "%lf"
#define msg( ... ) { char msgStr[32000]; snprintf( msgStr, 30000, __VA_ARGS__ ); printf( "[%s.%d] %s\n", __FILE__, __LINE__, msgStr ); }

int main( void ) {
	const char *const strInputFilename = "fort.300";
	const char *const strOutputFilename = "fort.100";
	const int sizeX = 800;
	bool verbo = false;
	//
	//
	//
	if (verbo) msg( "Attempting to allocate memory for vecX..." );
	REAL_TYPE *vecX = NULL;
	vecX = (REAL_TYPE *) malloc( sizeof(REAL_TYPE)*sizeX );
	if ( NULL == vecX ) {
		msg( "ERROR: Memory allocation for vecX failed." );
		return __LINE__;
	}
	memset( vecX, 0, sizeof(REAL_TYPE)*sizeX );
	if (verbo) msg( "Successfully allocated memory for vecX." );
	//
	const int sizeF = sizeX;
	if (verbo) msg( "Attempting to allocate memory for vecF..." );
	REAL_TYPE *vecF = NULL;
	vecF = (REAL_TYPE *) malloc( sizeof(REAL_TYPE)*sizeF );
	if ( NULL == vecF ) {
		msg( "ERROR: Memory allocation for vecF failed." );
		return __LINE__;
	}
	memset( vecF, 0, sizeof(REAL_TYPE)*sizeF );
	if (verbo) msg( "Successfully allocated memory for vecF." );
	//
	//
	//
	{
		if (verbo) msg( "Attempting to open input file..." );
		FILE *fPtr = fopen( strInputFilename, "r" );
		if ( NULL == fPtr ) {
			msg( "ERROR: Failed to open input file." );
			return __LINE__;
		}
		if (verbo) msg( "Successfully opened input file." );
		//
		if (verbo) msg( "Reading data from input file..." );
		for ( int n = 0; n < sizeX; n++ ) {
			int numItemsRead = fscanf( fPtr, REAL_TYPE_FORMAT_STR, &(vecX[n]) );
			if ( 1 != numItemsRead ) {
				if ( 0 == n ) {
					msg( "ERROR: Failed to read first expected value from input file." );
					return __LINE__;
				}
				msg( "ERROR: Failed to read next expected value from input file." );
				return __LINE__;
			}
			if (verbo) msg( "Read vecX(%d) = " REAL_TYPE_FORMAT_STR, n, vecX[n] );
		}
		{
			REAL_TYPE foo = 0;
			int retVal = fscanf( fPtr, REAL_TYPE_FORMAT_STR, &foo );
			if ( 0 < retVal ) {
				msg( "ERROR: Input file contains more data than expected." );
				return __LINE__;
			}
		}
		if (verbo) msg( "Successfully read %d values from input file.", sizeX );
		if (verbo) msg( "Closing input file." );
		fclose( fPtr );
	}
	//
	if (verbo) msg( "Calculating output..." );
	for ( int n = 0; n < sizeF; n++ ) {
		float x = ((double)n)*((double)sizeF)/(double(sizeX));
		int m = (int)( x + 0.5 );
		vecF[m] = (1.0+x)*vecX[n];
		msg( "Set vecF(%d) = " REAL_TYPE_FORMAT_STR, m, vecF[m] );
	}
	if (verbo) msg( "Finished calculating output." );
	//
	//
	{
		if (verbo) msg( "Attempting to open output file..." );
		FILE *fPtr = fopen( strOutputFilename, "w" );
		if ( NULL == fPtr ) {
			msg( "ERROR: Failed to open output file." );
			return __LINE__;
		}
		if (verbo) msg( "Successfully opened output file." );
		//
		if (verbo) msg( "Writing data to output file..." );
		for ( int n = 0; n < sizeF; n++ ) {
			int numCharsWritten = fprintf( fPtr, "%0.16lf\n", vecF[n] );
			if ( 0 >= numCharsWritten ) {
				if ( 0 == n ) {
					msg( "ERROR: Failed to write first value to output file." );
					return __LINE__;
				}
				msg( "ERROR: Failed to write next value to output file." );
				return __LINE__;
			}
		}
		//
		if (verbo) msg( "Successfully wrote %d values to output file.", sizeF );
		if (verbo) msg( "Closing output file." );
		fclose( fPtr );
	}
	//
	if (verbo) msg( "Freeing vecF." );
	free( vecF );
	if (verbo) msg( "Freeing vecX." );
	free( vecX );
return 0;
}
