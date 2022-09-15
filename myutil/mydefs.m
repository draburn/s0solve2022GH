% Verblev likely doesn't apply to built-in error and warning messages.
%  _SILENT: No messages. (May not be fully supported.)
%  _ERROR: Only fatal error messages.
%  _FLAGGED: For flagging stuff in development.
%  _WARN: Include warnings.
%  _NOTIFY: Include information of significant interest, like _FLAGGED in production.
%  _MAIN: Should be limited to one or two lines of the main results.
%  _INFO: Include a bit of additional information of interest beyond _MAIN.
%  _PROGRESS: Include progress information to a varying extent.
%  _COPIOUS: Include all progress information and worthwhile internal stuff.
%  _UNLIMITED: Anything and everything.
VERBLEV__SILENT = 0;
VERBLEV__ERROR = 100;
VERBLEV__FLAGGED = 200;
VERBLEV__WARNING = 300;
VERBLEV__NOTE = 400;
VERBLEV__NOTICE = 400;
VERBLEV__MAIN = 500;
VERBLEV__INFO = 550;
VERBLEV__PROG = 600;
VERBLEV__PROGRESS = 600;
VERBLEV__DETAILS = 650;
VERBLEV__DETAILED = 650;
VERBLEV__COPIOUS = 700;
VERBLEV__UNLIMITED = 800;
%
% These are verb verblevs...
VERBLEV__WARN = 300;
VERBLEV__NOTIFY = 400;
VERBLEV__INFORM = 550;



% ValdLev: Validation level, qualitative amount of time (and other resources)
%  to spend on validation checks.
%  _ZERO: None. (May not be fully supported.)
%  _LOW: Allow up to estimated +3% run time at most.
%  _MEDIUM: Allow up to estimated +30% run time typically.
%  _HIGH: Allow up to estimated x3 run time for some cases.
%  _VERY_HIGH: Allow up to estimated x100 run time for some cases.
%  _UNLIMITED: Do any and all validations.
VALDLEV__ZERO = 0;
VALDLEV__LOW = 100;
VALDLEV__MEDIUM = 200;
VALDLEV__HIGH = 300;
VALDLEV__VERY_HIGH = 350;
VALDLEV__UNLIMITED = 400;



% How broadly "success" is defined depends on the program.
% Using built-in error handling, the error return codes are probably moot.
RETCODE__SUCCESS = 0;
RETCODE__IMPOSED_STOP = 100;
RETCODE__ALGORITHM_BREAKDOWN = 200;
RETCODE__UNSPECIFIC_ERROR = 1000;
RETCODE__BAD_INPUT = 1100;
RETCODE__BAD_DEPENDENCY = 1200;
RETCODE__INTERNAL_INCONSISTENCY = 1300;
RETCODE__NUMERICAL_ISSUE = 1400;
RETCODE__NOT_SET = 9999;



% 2022-09-14.
% Let's make this groot-specific, and use friendly strings.
STR_GROOT_FLAG__UNSET = "UNSET";
STR_GROOT_FLAG__CNVG = "cnvg";
STR_GROOT_FLAG__STOP = "STOP";
STR_GROOT_FLAG__STALL = "STALL";
