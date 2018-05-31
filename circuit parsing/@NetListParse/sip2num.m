function [val,spl,sgf,tkn] = sip2num(obj,str,uni)
% Convert a metric-prefixed string into numeric values. (SI/engineering)
%
% (c) 2017 Stephen Cobeldick
%
% Convert a string containing numeric coefficients with metric prefixes into
% the equivalent numeric values, and also return the split string parts.
% For example the string '1k' is converted to the value 1000. The function
% identifies both symbol prefixes and full names, e.g. either 'k' or 'kilo'.
%
%%% Syntax:
%  val = sip2num(str)
%  val = sip2num(str,uni)
%  [val,spl,sgf] = sip2num(...)
%
% See also NUM2SIP BIP2NUM NUM2BIP WORDS2NUM STR2DOUBLE TEXTSCAN SSCANF NATSORT NATSORTFILES
%
%% Examples %%
%
% sip2num('10 k')  OR  sip2num('10.0 kilo')  OR  sip2num('10000')  OR  sip2num('1e4')
%   ans = 10000
%
% [val,spl] = sip2num('Power: 200 megawatt')
%   val = 200000000
%   spl = {'Power: ','watt'}
%
% [val,spl,sgf] = sip2num('from -3.6 MV to +1.24kV potential difference.')
%   val = [-3600000,1240]
%   spl = {'from ','V to ','V potential difference.'}
%   sgf = [2,3]
%
% [val,spl] = sip2num('100 meter','meter') % Try it without the second option.
%   val = 100
%   spl = {'','meter'}
%
% sip2num(num2sip(9*1000^4))
%   ans = 9000000000000  % == 9*1000^4
%
%% String Format %%
%
% * Any number of coefficients may occur in the string.
% * The coefficients may be any combination of digits, positive or negative,
%   integer or decimal, exponents may be included using E-notation (e/E).
% * An Inf or NaN value in the string will also be converted to a numeric.
% * The space-character between the coefficient and the prefix is optional.
% * The prefix is optional, either as the SI prefix symbol or full name.
% * By default checks first for prefix names, then symbols.
%
% Optional input <uni> controls the prefix/units recognition: if the units may
% contain the prefix characters, then this argument should be specified.
%
%% SI Prefix Strings %%
%
% Order  |1000^1 |1000^2 |1000^3 |1000^4 |1000^5 |1000^6 |1000^7 |1000^8 |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | kilo  | mega  | giga  | tera  | peta  | exa   | zetta | yotta |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol |   k   |   M   |   G   |   T   |   P   |   E   |   Z   |   Y   |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
%
% Order  |1000^-1|1000^-2|1000^-3|1000^-4|1000^-5|1000^-6|1000^-7|1000^-8|
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | milli | micro | nano  | pico  | femto | atto  | zepto | yocto |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol |   m   |   µ   |   n   |   p   |   f   |   a   |   z   |   y   |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
%
%% Input and Output Arguments %%
%
%%% Inputs (*=default):
%  str = String, with coefficients and prefixes to convert to numeric values.
%  uni = String, to specify the units that are given after the prefix.
%      = Logical Scalar, true/false -> match only the prefix name/symbol.
%      = *[], automagically check for prefix name or symbol, with any units.
%
%%% Outputs:
%  val = NumericVector, with values calculated from the coefficients and prefixes
%        given in <str>. The size is 1xN, N = number of detected coefficients.
%  spl = CellOfStrings, parts of <str> split by the detected coefficients(+prefixes).
%  sgf = NumericVector, same size as <val>, significant figures of each coefficient.
%
% [val,spl,sgf] = sip2num(str,*uni)

%% Disclaimer %%

% Copyright (c) 2017, Stephen Cobeldick
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%% Modifications From Copyright %%

% !@#$%^& Edited by Jared Baxter to use µ instead of u for micro *&^%$#@!%

%% Input Wrangling %%
%
pfC = {...
    'y',    'z',    'a',   'f',    'p',   'n',   'µ',    'm',    '','k',   'M',   'G',   'T',   'P',   'E',  'Z',    'Y';...
    'yocto','zepto','atto','femto','pico','nano','micro','milli','','kilo','mega','giga','tera','peta','exa','zetta','yotta'};
vPw = [ -24,    -21,   -18,    -15,   -12,    -9,     -6,     -3,+0,    +3,   +6,     +9,   +12,   +15,  +18,    +21,    +24];
pfB = 10;
%
pfx = 'yocto|zepto|atto|femto|pico|nano|micro|milli|kilo|mega|giga|tera|peta|exa|zetta|yotta';
sym = 'y|z|a|f|p|n|µ|m|k|M|G|T|P|E|Z|Y';
sep = '|';
sfx = '';
%
% Determine the prefix+unit combination:
if nargin<3||(isnumeric(uni)&&isempty(uni))
	% Name/symbol prefix, any units.
elseif ischar(uni)&&isrow(uni)
	% Units are the given string:
	sfx = ['(?=',regexptranslate('escape',uni),')'];
elseif islogical(uni)&&isscalar(uni)
	sep = '';
	if uni % Prefix names only.
		sym = '';
	else   % Prefix symbols only.
		pfx = '';
	end
else
	error('Second input <uni> must be a logical scalar, a string, or empty numeric.')
end
assert(ischar(str)&&isrow(str),'First input <str> must be a 1xN char.')
%
%% String Parsing %%
%
% Try to locate a coefficient, possibly with a prefix:
rgx = ['([+-]?(NaN|Inf|\d+\.?\d*)([eE][+-]?\d+)?)\s?(',pfx,sep,sym,')?',sfx];
[tkn,spl] = regexp(str,rgx,'tokens','split');
tkn = vertcat(tkn{:});
%
if isempty(tkn)
	% No coefficient found:
	val = [];
	sgf = [];
else
	% Calculate values from the coefficients:
	for k = size(tkn,1):-1:1
		if isempty(tkn{k,2})
			val(k) = sscanf(tkn{k,1},'%f');
		else
			[~,idx] = find(strcmp(tkn{k,2},pfC));
			val(k) = sscanf(tkn{k,1},'%f') * pfB^vPw(idx);
		end
	end
	if nargout>2
		rgx = {'^[+-]','Inf|NaN','[eE].*$','^0+\.?0*(?=[1-9])','^0+\.','\.'};
		sgf = cellfun('length',regexprep(tkn(:,1),rgx,'')).';
	end
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sip2num