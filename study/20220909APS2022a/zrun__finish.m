msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, sprintf( "Finished suite \"%s\".", suiteName ) );
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
setprngstatedat(backup_prngStateDat); % May be redundant.
%
msg( __FILE__, __LINE__, sprintf( "Elapsed time is %gs.", time()-startTime ) );