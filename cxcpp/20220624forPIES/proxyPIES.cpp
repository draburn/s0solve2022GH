#include <stdlib.h>
#include <stdio.h>
#include <string.h> // For memset?!?!

#define REAL_TYPE double
#define REAL_TYPE_FORMAT_STR "%lf"
//#define msg( ... ) { char msgStr[32000]; snprintf( msgStr, 30000, __VA_ARGS__ ); printf( "[%s.%d] %s\n", __FILE__, __LINE__, msgStr ); }
#define msg( ... )

int main( void ) {
	msg( "Hello world!" );
	const char *const strInputFilename = "x.txt";
	const char *const strOutputFilename = "f.txt";
	const int sizeX = 800;
	//
	//
	//
	msg( "Attempting to allocate memory for vecX..." );
	REAL_TYPE *vecX = NULL;
	vecX = (REAL_TYPE *) malloc( sizeof(REAL_TYPE)*sizeX );
	if ( NULL == vecX ) {
		msg( "Memory allocation for vecX failed." );
		return __LINE__;
	}
	memset( vecX, 0, sizeof(REAL_TYPE)*sizeX );
	msg( "Successfully allocated memory for vecX." );
	//
	const int sizeF = sizeX;
	msg( "Attempting to allocate memory for vecF..." );
	REAL_TYPE *vecF = NULL;
	vecF = (REAL_TYPE *) malloc( sizeof(REAL_TYPE)*sizeF );
	if ( NULL == vecF ) {
		msg( "Memory allocation for vecF failed." );
		return __LINE__;
	}
	memset( vecF, 0, sizeof(REAL_TYPE)*sizeF );
	msg( "Successfully allocated memory for vecF." );
	//
	//
	//
	{
		msg( "Attempting to open input file..." );
		FILE *fPtr = fopen( strInputFilename, "r" );
		if ( NULL == fPtr ) {
			msg( "Failed to open input file." );
			return __LINE__;
		}
		msg( "Successfully open input file." );
		//
		msg( "Reading data from input file..." );
		for ( int n = 0; n < sizeX; n++ ) {
			int numItemsRead = fscanf( fPtr, REAL_TYPE_FORMAT_STR, &(vecX[n]) );
			if ( 1 != numItemsRead ) {
				if ( 0 == n ) {
					msg( "Failed to read first expected value from input file." );
					return __LINE__;
				}
				msg( "Failed to read next expected value from input file." );
				return __LINE__;
			}
			msg( "Read vecX(%d) = " REAL_TYPE_FORMAT_STR, n, vecX[n] );
		}
		{
			REAL_TYPE foo = 0;
			int retVal = fscanf( fPtr, REAL_TYPE_FORMAT_STR, &foo );
			if ( 0 < retVal ) {
				msg( "Input file contains more data than expected." );
				return __LINE__;
			}
		}
		msg( "Successfully read %d values from input file.", sizeX );
		msg( "Closing input file." );
		fclose( fPtr );
	}
	//
	msg( "Calculating output..." );
	for ( int n = 0; n < sizeF; n++ ) {
		float x = ((double)n)*((double)sizeF)/(double(sizeX));
		int m = (int)( x + 0.5 );
		vecF[m] = (1.0+x)*vecX[n];
		msg( "Set vecF(%d) = " REAL_TYPE_FORMAT_STR, m, vecF[m] );
	}
	msg( "Finished calculating output." );//
	//
	{
		msg( "Attempting to open output file..." );
		FILE *fPtr = fopen( strOutputFilename, "w" );
		if ( NULL == fPtr ) {
			msg( "Failed to open output file." );
			return __LINE__;
		}
		msg( "Successfully opened output file." );
		//
		msg( "Writing data to output file..." );
		for ( int n = 0; n < sizeF; n++ ) {
			int numCharsWritten = fprintf( fPtr, "%0.16lf\n", vecF[n] );
			if ( 0 >= numCharsWritten ) {
				if ( 0 == n ) {
					msg( "Failed to write first value to output file." );
					return __LINE__;
				}
				msg( "Failed to write next value to output file." );
				return __LINE__;
			}
		}
		//
		msg( "Successfully wrote %d values to output file.", sizeF );
		msg( "Closing output file." );
		fclose( fPtr );
	}
	//
	msg( "Freeing vecF." );
	free( vecF );
	msg( "Freeing vecX." );
	free( vecX );
return 0;
}
