#include <stdlib.h>
#include <stdio.h>

#define REAL_TYPE double
#define REAL_TYPE_FORMAT_STR "%lf"
#define msg( str ) printf( "[%s.%d] %s\n", __FILE__, __LINE__, str )

int main( void ) {
	int retCode = -1;
	int sz = 800;
	msg( "Hello world!" );
	FILE *fPtr = fopen( "x.txt", "r" );
	if ( NULL == fPtr ) {
		msg( "Failed to open input file." );
		return __LINE__;
	}
	//
	REAL_TYPE *vecX = NULL;
	vecX = (REAL_TYPE *) malloc( sizeof(REAL_TYPE)*sz );
	if ( NULL == vecX ) {
		msg( "Memory allocation for vecX failed." );
		return __LINE__;
	}
	//
	for ( int n = 0; n < sz; n++ ) {
		int retVal = fscanf( fPtr, REAL_TYPE_FORMAT_STR, &(vecX[n]) );
		//const char *const tempStr = sprintf( "Read vecX(%s) = %s\n", "%d", REAL_TYPE_FORMAT_STR );
		//printf( tempStr, n, vecX[n] );
		if ( 1 != retVal ) {
			msg( "Failed to read next expected value from input file." );
			return __LINE__;
		}
		//const char foo = sprintf( "Read vecX(%s) = %s.\n", "%d", REAL_TYPE_FORMAT_STR );
		//printf( foo, n, vecX[n] );
		//const char foo2 = sprintf( foo, n, vecX[n] );
		//msg( foo2 );
		printf( "Read vecX(%d) = %lf\n", n, vecX[n] );
	}
	{
		REAL_TYPE foo = 0;
		int retVal = fscanf( fPtr, REAL_TYPE_FORMAT_STR, &foo );
		if ( 0 < retVal ) {
			msg( "Input file contains more data than expected." );
			return __LINE__;
		}
	}
	msg( "Successfully read number of expected values from input file." );
	msg( "Closing file." );
	fclose( fPtr );
	//
	msg( "Freeing vecX." );
	free( vecX );
return 0;
}
