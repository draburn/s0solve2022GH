#include <stdlib.h>
#include <stdio.h>
#include <string.h> // For memset?!?!

#define msg( ... ) { char msgStr[32000]; snprintf( msgStr, 30000, __VA_ARGS__ ); printf( "[%s.%d] %s\n", __FILE__, __LINE__, msgStr ); }
#define MAX_FNAME_LENGTH 32000
#define MAX_SIZEX 32000

int main( void ) {
	int n0 = 0;
	int n1 = 0;
	int sizeX = 0;
	double epsX = 0.0;
	char fnameX0[MAX_FNAME_LENGTH];
	char fnameXN[MAX_FNAME_LENGTH];
	char fnameFN[MAX_FNAME_LENGTH];
	double vecX0[MAX_SIZEX];
	char cmdFeval[MAX_FNAME_LENGTH];
	char fnameOutputDir[MAX_FNAME_LENGTH];
	char cmdMkdir[MAX_FNAME_LENGTH];
	char cmdMv[MAX_FNAME_LENGTH];
	bool allowOverwrite = false;
	//
	{
		char fnameParam[MAX_FNAME_LENGTH];
		strncpy( fnameParam, "jevalDriver_param.txt", MAX_FNAME_LENGTH-1 );
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
			sprintf( formatStr, "%%%ds", MAX_FNAME_LENGTH-2 );
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
			sprintf( formatStr, "%%%ds", MAX_FNAME_LENGTH-2 );
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
			sprintf( formatStr, "%%%ds", MAX_FNAME_LENGTH-2 );
			int numItemsRead = fscanf( fPtr, formatStr, fnameFN ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read fnameFN from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read fnameFN = \"%s\".", fnameFN );
		//
		{
			char formatStr [100];
			sprintf( formatStr, "%%%ds", MAX_FNAME_LENGTH-2 );
			int numItemsRead = fscanf( fPtr, formatStr, cmdFeval ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read cmdFeval from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read cmdFeval = \"%s\".", cmdFeval );
		//
		{
			char formatStr [100];
			sprintf( formatStr, "%%%ds", MAX_FNAME_LENGTH-2 );
			int numItemsRead = fscanf( fPtr, formatStr, fnameOutputDir ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read fnameOutput from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read fnameOutputDir = \"%s\".", fnameOutputDir );
		//
		{
			char formatStr [100];
			sprintf( formatStr, "%%%ds", MAX_FNAME_LENGTH-2 );
			int numItemsRead = fscanf( fPtr, formatStr, cmdMkdir ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read fnameOutput from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read cmdMkdir = \"%s\".", cmdMkdir );
		//
		{
			char formatStr [100];
			sprintf( formatStr, "%%%ds", MAX_FNAME_LENGTH-2 );
			int numItemsRead = fscanf( fPtr, formatStr, cmdMv ); // Is this safe?
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read cmdMv from parameter file \"%s\".", numItemsRead, fnameParam );
				fclose( fPtr );
				return __LINE__;
			}
		}
		msg( "Read cmdMv = \"%s\".", cmdMv );
		//
		fclose( fPtr );
	}
	msg( "Finished reading paramter file." );
	//
	//
	//
	{
		FILE *fPtr = fopen( fnameX0, "r" );
		if ( NULL == fPtr ) {
			msg( "ERROR: Failed to open x0 file \"%s\".", fnameX0 );
			return __LINE__;
		}
		msg( "Opened x0 file \"%s\".", fnameX0 );
		//
		for ( int n = 0; n < sizeX; n++ ) {
			int numItemsRead = fscanf( fPtr, "%lg", &(vecX0[n]) );
			if ( 1 != numItemsRead ) {
				msg( "ERROR: fscanf() returned %d trying to read vecX0[%d] from file \"%s\".", numItemsRead, n, fnameX0 );
				fclose( fPtr );
				return __LINE__;
			}
		}
		{
			double foo = 0.0;
			int numItemsRead = fscanf( fPtr, "%lg", &foo );
			if ( 1 == numItemsRead ) {
				msg( "ERROR: vecX0 file \"%s\" contains more than %d items.", fnameX0, sizeX );
				fclose( fPtr );
				return __LINE__;
			}
		}
		//
		fclose( fPtr );
	}
	msg( "Finished reading x0 file." );
	//
	//
	{
		char cmdTemp[MAX_FNAME_LENGTH];
		snprintf( cmdTemp, MAX_FNAME_LENGTH, "%s %s", cmdMkdir, fnameOutputDir );
		msg( "Calling command \"%s\".", cmdTemp );
		system( cmdTemp );
	}
	msg( "Created output directory \"%s\" (or it already exists).", fnameOutputDir );
	//
	//
	for ( int n = n0; n <= n1; n++ ) {
		if (!allowOverwrite) {
			FILE *fPtr = fopen( fnameXN, "r" );
			if ( NULL != fPtr ) {
				msg( "ERROR: File \"%s\" already exists.", fnameXN );
				fclose( fPtr );
				return __LINE__;
			}
		}
		FILE *fPtr = fopen( fnameXN, "w" );
		if ( NULL == fPtr ) {
			msg( "ERROR: Failed to start writing to file \"%s\".", fnameXN );
			fclose( fPtr );
			return __LINE__;
		}
		msg( "Opened x file \"%s\".", fnameXN );
		//
		for ( int m = 0; m < sizeX; m++ ) {
			double x = vecX0[m];
			if ( n-1 == m ) {
				x += epsX;
			}
			int numCharsWritten = fprintf( fPtr, "%0.16le\n", x );
			if ( 0 >= numCharsWritten ) {
				msg( "Failed to write next value to x file %d.", n );
				return __LINE__;
			}
		}
		msg( "Successfully wrote %d values to x file %d \"%s\".", sizeX, n, fnameXN );
		fclose( fPtr );
		//
		msg( "Executing feval command \"%s\".", cmdFeval );
		{
			int retVal = system( cmdFeval );
			if ( 0 != retVal ) {
				msg( "WARNING: FEval command \"%s\" returned %d.", cmdFeval, retVal );
			}
		}
		//
		msg( "Moving files..." );
		{
			char fnameTemp[MAX_FNAME_LENGTH];
			snprintf( fnameTemp, MAX_FNAME_LENGTH, "%s/x%06d.txt", fnameOutputDir, n );
			if (!allowOverwrite) {
				FILE *fPtr = fopen( fnameTemp, "r" );
				if ( NULL != fPtr ) {
					msg( "ERROR: File \"%s\" already exists.", fnameTemp );
					fclose( fPtr );
					return __LINE__;
				}
			}
			char cmdTemp[MAX_FNAME_LENGTH];
			snprintf( cmdTemp, MAX_FNAME_LENGTH, "%s %s %s", cmdMv, fnameXN, fnameTemp );
			msg( "Calling command \"%s\".", cmdTemp );
			int retVal = system( cmdTemp );
			if ( 0 != retVal ) {
				msg( "ERROR: X move command \"%s\" returned %d.", cmdTemp, retVal );
				return __LINE__;
			}
		}
		{
			char fnameTemp[MAX_FNAME_LENGTH];
			snprintf( fnameTemp, MAX_FNAME_LENGTH, "%s/f%06d.txt", fnameOutputDir, n );
			if (!allowOverwrite) {
				FILE *fPtr = fopen( fnameTemp, "r" );
				if ( NULL != fPtr ) {
					msg( "ERROR: File \"%s\" already exists.", fnameTemp );
					fclose( fPtr );
					return __LINE__;
				}
			}
			char cmdTemp[MAX_FNAME_LENGTH];
			snprintf( cmdTemp, MAX_FNAME_LENGTH, "%s %s %s", cmdMv, fnameFN, fnameTemp );
			msg( "Calling command \"%s\".", cmdTemp );
			int retVal = system( cmdTemp );
			if ( 0 != retVal ) {
				msg( "ERROR: F move command \"%s\" returned %d.", cmdTemp, retVal );
				return __LINE__;
			}
		}
	}
	//
return 0;
}
