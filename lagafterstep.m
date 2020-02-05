function [lagtime, x_atlagtime] = lagtime_afterstep(t, x, mark, tstep)
%
% For a quantity x (e.g., soil thickness at a point on a hillslope)
% responding to a step change in some other forcing, compute the time it
% takes for x to pass a certain threshold during its transition from its
% initial steady state to its final steady state.  The trickiest part of
% this is accounting for quantities that approach their final steady state
% non-monotonically.
%
% Inputs:
%   x               quantity of interest (vector)
%   t               time since step forcing was applied (first element should 
%                       be 0) (vector)
%   mark            fraction of difference between final and initial
%                       steady-state values of x (scalar; e.g., 0.95 for 
%                       reaching 95% of the total change in x)d
%   tstep           time at which step occurred in driver of x
%
% Output:
%   lagtime         time it takes to pass criterion.  Units are same as
%                       units of the input t vector.
%   x_atlagtime     value of x
%
% Underlying assumptions:
%   1. First element of x vector is the initial steady state value of x.
%   2. Final element of x vector is the final steady state value of x.


%% Compute lag.

% Assign the lag time to be the time at which the variable reaches an
% assigned fraction of the way from the old to the new steady state.

% Pre-step steady state value (initial conditions).
x_oldsteadystate = x(1);  % First point in x vector must be initial
    % steady state.    

% Post-step steady state value (final conditions).
x_newsteadystate = x(end);  % Last point in x vector must be final
    % steady state.

% One approach: Define mark for approaching final steady state as relative
% to the initial steady state.  First, compute absolute value between new
% and initial steady states.  Then compute difference between x and the
% final steady state value of x.  Starting from the end and working 
% backward in time, that difference starts at zero and grows.  When it
% passes the criterion for convergence, that's where to mark the lag time.
% E.g., if mark = 0.95, then find where x deviates from final steady state
% value by more than 5% of the difference between new and old steady 
% states.
x_absnewminusold = abs(x_newsteadystate - x_oldsteadystate);
x_absdiff2 = abs(x - x_newsteadystate);
x_absdiffratio2 = x_absdiff2/x_absnewminusold;

% To find the time at which the time series passes the mark, analyze the 
% time series in reverse.  Flip x from end to start, and find the first
% time (from the end) where the lag criterion is not met.  Things you need:
% 1) define the lag criterion; 2) compute whether criterion has been met at
% each point in the time series; 3) find the first point where the
% criterion goes from met to not met.

% Flip t vector.
if size(t,1)==1  % if t is a row vector, flip left to right
    flipped_t = fliplr(t);
else  % if it's a column vector, flip top to bottom
    flipped_t = flipud(t);
end

% Flip x vector.
if size(x_absdiffratio2,1)==1  % if x is a row vector, flip left to right
    flipped_x_absdiffratio = fliplr(x_absdiffratio2);
else  % if it's a column vector, flip top to bottom
    flipped_x_absdiffratio = flipud(x_absdiffratio2);
end

% Find time at which the criterion is first met.
criterionforlag = flipped_x_absdiffratio < (1 - mark); 
indexforcriterion = min(find(criterionforlag==0));  % index in vector
    % where the lag criterion is first met
if isempty(indexforcriterion)  % if criterion isn't met at any time, set
        % t_Wmark to the first value in the flipped_t vector, which is the
        % last element of t.
    t_Wmark = flipped_t(1);  % Set to maximum value of t
elseif indexforcriterion==1  % if criterion is met in first element of 
        % t_flipped
    t_Wmark = flipped_t(1);  % Set to maximum value of t
else
    t_Wmark = flipped_t(indexforcriterion-1);  % time at first moment the
        % lag criterion is met
end

% Compute time lag as difference between the time it hits the mark minus
% the time when the step happened.
lagtime = t_Wmark - tstep;

% Compute value of x the lag time.
x_atlagtime = x(min(find(t == tstep + lagtime)));


