% Verblev likely doesn't apply to built-in error and warning messages.
%  _SILENT: No messages. (May not be fully supported.)
%  _ERROR: only fatal error messages.
%  _FLAGGED: for flagging stuff in development.
%  _WARN: include warnings.
%  _NOTIFY: include information of interest, like _FLAGGED in production.
%  _MAIN: should be limited to one or two lines of the main results.
%  _PROGRESS: include progress information to a varying extent.
%  _COPIOUS: include all progress information and worthwhile internal stuff.
%  _UNLIMITED: anything and everything.
VERBLEV__SILENT = 0;
VERBLEV__ERROR = 100;
VERBLEV__FLAGGED = 200;
VERBLEV__WARN = 300;
VERBLEV__NOTIFY = 400;
VERBLEV__MAIN = 500;
VERBLEV__PROGRESS = 600;
VERBLEV__COPIOUS = 700;
VERBLEV__UNLIMITED = 800;

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
VALDLEV__VERY_HIGH = 400;
VALDLEV__UNLIMITED = 900;
