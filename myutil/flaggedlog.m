function retCode = flaggedlog( fileName="", lineNum=0, str="", fname="FLAGGEDLOG.txt"  );
	msg( fileName, lineNum, str );
	fid = fopen( fname, "a" );
	if ( -1 == fid )
		if ( 0 == nargout )
			error( str );
		else
			retCode = 1;
			return;
		endif
	endif
	warning( str );
	clockNow = clock();
	fwrite( fid, sprintf( "flog() called from file \"%s\" line %d.\n", fileName, lineNum ) );
	fwrite( fid, sprintf( "  Currently time is %02d:%02d:%02d.%03d.\n", ...
	  clockNow(4), clockNow(5), floor(clockNow(6)), floor(1000*(clockNow(6) - floor(clockNow(6)))) )  );
	fwrite( fid, sprintf( "  Message: \"%s\".\n", str ) );
	fwrite( fid, "\n" );
	fclose( fid );
retCode = 0;
return;
endfunction
