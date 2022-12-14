%  Script...
%    commondefs
%  Overview...
%    Part of my common module.
%    This script does the following...
%      1. Defines verbosity level (VERBLEV__) constants;
%      2. Defines return code (RETCODE__) constants; and,
%      3. Defines msg_threshold-family function handles
%         such as msg_error() and msg_progress().
%    See msg_thresh and msg for more information.

% Verblev likely doesn't apply to built-in error handling;
%  "quiet" may be equivalent to "error".
% "Flagged" is for cases flagged during development: "Can this ever happen?".
% "Notify" is for undesirable situations where likely nothing is wrong.
% "Main" should be limited to one or two lines of output.
% These levels are intended to be broadly suitable;
%  more detailed control should probably be done in situ,
%  perhaps by locally using verblevs between these values.
VERBLEV__QUIET = 0;
VERBLEV__ERROR = 100;
VERBLEV__FLAGGED = 200;
VERBLEV__WARN = 300;
VERBLEV__NOTIFY = 400;
VERBLEV__MAIN = 500;
VERBLEV__PROGRESS = 600;
VERBLEV__COPIOUS = 700;
VERBLEV__UNLIMITED = 800;


% How broadly "success" is defined depends on the program.
% Using built-in error handling, the error return codes are probably moot.
RETCODE__SUCCESS = 0;
RETCODE__IMPOSED_STOP = 100;
RETCODE__ALGORITHM_BREAKDOWN = 200;
RETCODE__UNSPECIFIC_ERROR = 1000;
RETCODE__BAD_INPUT = 1100;
RETCODE__BAD_DEPENDENCY = 1200;
RETCODE__INTERNAL_INCONSISTENCY = 1300;
RETCODE__NOT_SET = 9999;

msg_error = @( verbLev, fileName, lineNum, msgStr )( ...
  msg_thresh( VERBLEV__ERROR, verbLev, fileName, lineNum, msgStr) );
msg_flagged = @( verbLev, fileName, lineNum, msgStr )( ...
  msg_thresh( VERBLEV__FLAGGED, verbLev, fileName, lineNum, msgStr) );
msg_warn = @( verbLev, fileName, lineNum, msgStr )( ...
  msg_thresh( VERBLEV__WARN, verbLev, fileName, lineNum, msgStr) );
msg_notify = @( verbLev, fileName, lineNum, msgStr )( ...
  msg_thresh( VERBLEV__NOTIFY, verbLev, fileName, lineNum, msgStr) );
msg_main = @( verbLev, fileName, lineNum, msgStr )( ...
  msg_thresh( VERBLEV__MAIN, verbLev, fileName, lineNum, msgStr) );
msg_progress = @( verbLev, fileName, lineNum, msgStr )( ...
  msg_thresh( VERBLEV__PROGRESS, verbLev, fileName, lineNum, msgStr) );
msg_copious = @( verbLev, fileName, lineNum, msgStr )( ...
  msg_thresh( VERBLEV__COPIOUS, verbLev, fileName, lineNum, msgStr) );

STEPTYPE__NEWTON = 0;
STEPTYPE__PICARD = 100;
STEPTYPE__PICARD_SCALED = 101;
STEPTYPE__GRADDIR = 200;
STEPTYPE__GRADDIR_SCALED = 201;
STEPTYPE__LEVCURVE = 300;
STEPTYPE__LEVCURVE_SCALED = 301;
STEPTYPE__GRADCURVE = 400;
STEPTYPE__GRADCURVE_SCALED = 401;
STEPTYPE__SPECIFIED_VECTOR = 800;

% ValLev: Validation level, qualitative amount of time (and other resources)
%  to spend on validation checks.
%  _ZERO: None. (May not be fully supported.)
%  _LOW: Allow up to estimated +3% run time at most.
%  _MEDIUM: Allow up to estimated +30% run time typically.
%  _HIGH: Allow up to estimated x3 run time for some cases.
%  _VERY_HIGH: Allow up to estimated x100 run time for some cases.
%  _UNLIMITED: Do any and all validations.
% AKA ValdLev.
VALLEV__ZERO = 0;
VALLEV__LOW = 100;
VALLEV__MEDIUM = 200;
VALLEV__HIGH = 300;
VALLEV__VERY_HIGH = 350;
VALLEV__UNLIMITED = 400;
VALDLEV__ZERO = 0;
VALDLEV__LOW = 100;
VALDLEV__MEDIUM = 200;
VALDLEV__HIGH = 300;
VALDLEV__VERY_HIGH = 350;
VALDLEV__UNLIMITED = 400;
%
eps025 = eps^0.25;
eps050 = eps^0.50;
eps075 = eps^0.75;
eps100 = eps;
eps125 = eps^1.25;
eps150 = eps^1.50;
eps175 = eps^1.75;
eps200 = eps^2;
eps250 = eps^2.50;
eps300 = eps^3;

% QualLev: Quality level, qualitative threshold for what constitutes "success".
%  _ZERO: No minimal quality level, likely useful only for testing.
%  _MEDIUM: Intended typical quality level.
%  _EXTREME: May break code, likely useful only for testing.
QUALLEV__ZERO = 0;
QUALLEV__VERY_LOW = 100;
QUALLEV__LOW = 200;
QUALLEV__MEDIUM = 300;
QUALLEV__HIGH = 400;
QUALLEV__VERY_HIGH = 500;
QUALLEV__EXTREME = 600;

% InvestLev: Investment level, qualitative amount of time (and other resources)
%  to spend tyring to achieve "success".
%  _ZERO: No minimal level of work, likely useful only for testing.
%  _MEDIUM: Intended typical effort level.
%  _EXTREME: May break code, likely useful only for testing.
INVESTLEV__ZERO = 0;
INVESTLEV__VERY_LOW = 100;
INVESTLEV__LOW = 200;
INVESTLEV__MEDIUM = 300;
INVESTLEV__HIGH = 400;
INVESTLEV__VERY_HIGH = 500;
INVESTLEV__EXTREME = 600;

%!test
%!	commondefs
%!	thisFile = "test commondefs";
%!	verbLev = VERBLEV__MAIN;
%!	disp( "" );
%!	disp( "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
%!	msg_error( verbLev, thisFile, __LINE__, "This should be displayed." );
%!	msg_main( verbLev, thisFile, __LINE__, "This should also be displayed." );
%!	msg_progress( verbLev, thisFile, __LINE__, "THIS SHOULD NOT BE DISPLAYED! ***" );
%!	disp("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
%!	disp( "*** Is the above text correct?  ***" );
%!	disp( "" );
