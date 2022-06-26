#include <stdlib.h>
#include <stdio.h>
#include <string.h> // For memset?!?!

#define REAL_TYPE double
#define REAL_TYPE_FORMAT_STR "%lf"
#define msg( ... ) { char msgStr[32000]; snprintf( msgStr, 30000, __VA_ARGS__ ); printf( "[%s.%d] %s\n", __FILE__, __LINE__, msgStr ); }
#define MAX_SIZE 32000

int main( void ) {
	int n0 = 0;
	int n1 = 0;
	int sizeX = 0;
	double epsX = 0.0;
	char fnameX0[MAX_SIZE];
	char fnameXN[MAX_SIZE];
	char fnameFN[MAX_SIZE];
	REAL_TYPE vecX0[MAX_SIZE];
	//
	{
		char fnameParam[MAX_SIZE];
		strncpy( fnameParam, "jevalDriver_param.txt", MAX_SIZE-1 );
		FILE *fPtr = fopen( fnameParam, "r" );
		if ( NULL == fPtr ) {
			msg( "ERROR: Failed to open parameter file \"%s\".", fnameParam );
			return __LINE__;
		}
		msg( "Opened paramter file \"%s\".", fnameParam );
		//
		{
			int numItemsRead = fscanf( fPtr, "%d", &n0 );
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read n0 from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
			if ( 0 >= n0 ) {
				msg( "ERROR: n0 is non-positive (%d).", n0 );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read n0 = %d.", n0 );
		//
		{
			int numItemsRead = fscanf( fPtr, "%d", &n1 );
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read n1 from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
			if ( 0 >= n0 ) {
				msg( "ERROR: n1 is non-positive (%d).", n1 );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read n1 = %d.", n1 );
		//
		{
			int numItemsRead = fscanf( fPtr, "%lg", &epsX );
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read epsX from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
			if ( 0.0 >= epsX ) {
				msg( "ERROR: epsX is non-positive (%le).", epsX );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read epsX = %le.", epsX );
		//
		{
			int numItemsRead = fscanf( fPtr, "%d", &sizeX );
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read sizeX from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
			if ( 0 >= n0 ) {
				msg( "ERROR: sizeX is non-positive (%d).", sizeX );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read sizeX = %d.", sizeX );
		//
		//
		{
			
			char formatStr [100];
			sprintf( formatStr, "%%%ds", MAX_SIZE-2 );
			int numItemsRead = fscanf( fPtr, formatStr, fnameX0 ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read fnameX0 from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read fnameX0 = \"%s\".", fnameX0 );
		//
		{
			
			char formatStr [100];
			sprintf( formatStr, "%%%ds", MAX_SIZE-2 );
			int numItemsRead = fscanf( fPtr, formatStr, fnameXN ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read fnameXN from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read fnameXN = \"%s\".", fnameX0 );
		//
		{
			
			char formatStr [100];
			sprintf( formatStr, "%%%ds", MAX_SIZE-2 );
			int numItemsRead = fscanf( fPtr, formatStr, fnameFN ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read fnameFN from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read fnameFN = \"%s\".", fnameFN );
		//
		fclose( fPtr );
	}
	msg( "Finished reading paramter file." );
	//
return 0;
}
