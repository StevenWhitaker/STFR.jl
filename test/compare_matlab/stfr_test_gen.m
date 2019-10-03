% Generate .mat files for comparing STFR in Julia to STFR in MATLAB (they should
% be the same).

%% stfr1_test1
clear

M0 = 1;
T1 = 1000;
T2 = 80;
wf = 1 * 2*pi;
kappa = 1.1;
Tfree = 7.55;
Tg = 4.2;
alpha = 10 * 2*pi;
beta = 10 * 2*pi;
phi = 30 * 2*pi;
spoil = true;
% result uses the single-compartment STFR because ff = 0 (the second argument of
% stfr_2comp_forward_model(...))
% The transpose is necessary because the output of stfr_2comp_forward_model(...)
% is [N,P], whereas the output of the Julia version is [P,N]
result = stfr_2comp_forward_model(M0, 0, 1, 1, T1, T2, wf, 0, kappa, Tfree, ...
	                              Tg, alpha, beta, phi, 'spoil', spoil).';

save('stfr1_test1.mat', 'M0', 'T1', 'T2', 'wf', 'kappa', 'Tfree', 'Tg', ...
	 'alpha', 'beta', 'phi', 'spoil', 'result', '-v7.3');

%% stfr1_test2
clear

M0 = [1; 0.8; 1.1];
T1 = [1000; 900; 1500];
T2 = [80; 120; 300];
wf = [1; 2; 3] * 2*pi;
kappa = [1.1; 0.8; 1];
Tfree = 7.55;
Tg = 4.2;
alpha = 10 * 2*pi;
beta = 10 * 2*pi;
phi = 30 * 2*pi;
spoil = false;
result = stfr_2comp_forward_model(M0, [0; 0; 0], [1; 1; 1], [1; 1; 1], T1, ...
	                              T2, wf, [0; 0; 0], kappa, Tfree, ...
	                              Tg, alpha, beta, phi, 'spoil', spoil).';

save('stfr1_test2.mat', 'M0', 'T1', 'T2', 'wf', 'kappa', 'Tfree', 'Tg', ...
	 'alpha', 'beta', 'phi', 'spoil', 'result', '-v7.3');

%% stfr1_test3
clear

M0 = 1.05;
T1 = 2000;
T2 = 110;
wf = 3 * 2*pi;
kappa = 1.01;
Tfree = [7.45; 6];
Tg = [4.1; 3];
alpha = [11; 20] * 2*pi;
beta = [11; 15] * 2*pi;
phi = [31; 2] * 2*pi;
spoil = false;
result = stfr_2comp_forward_model(M0, 0, 1, 1, T1, T2, wf, 0, kappa, Tfree, ...
	                              Tg, alpha, beta, phi, 'spoil', spoil).';

save('stfr1_test3.mat', 'M0', 'T1', 'T2', 'wf', 'kappa', 'Tfree', 'Tg', ...
	 'alpha', 'beta', 'phi', 'spoil', 'result', '-v7.3');

%% stfr1_test4
clear

M0 = [0.95; 0.8; 1.1];
T1 = [1100; 900; 1500];
T2 = [90; 120; 300];
wf = [11; 2; 3] * 2*pi;
kappa = [1.3; 0.8; 1];
Tfree = [7.45; 6];
Tg = [4.1; 3];
alpha = [11; 20] * 2*pi;
beta = [11; 15] * 2*pi;
phi = [31; 2] * 2*pi;
spoil = true;
result = stfr_2comp_forward_model(M0, [0; 0; 0], [1; 1; 1], [1; 1; 1], T1, ...
	                              T2, wf, [0; 0; 0], kappa, Tfree, ...
	                              Tg, alpha, beta, phi, 'spoil', spoil).';

save('stfr1_test4.mat', 'M0', 'T1', 'T2', 'wf', 'kappa', 'Tfree', 'Tg', ...
	 'alpha', 'beta', 'phi', 'spoil', 'result', '-v7.3');

%% stfr2_test1
clear

M0 = 1;
ff = 0.2;
T1f = 400;
T1s = 1000;
T2f = 20;
T2s = 80;
wf = 1 * 2*pi;
wff = 10 * 2*pi;
kappa = 1.1;
Tfree = 7.55;
Tg = 4.2;
alpha = 10 * 2*pi;
beta = 10 * 2*pi;
phi = 30 * 2*pi;
spoil = true;
result = stfr_2comp_forward_model(M0, ff, T1f, T2f, T1s, T2s, wf, wff, ...
	                              kappa, Tfree, Tg, alpha, beta, phi, ...
								  'spoil', spoil).';

save('stfr2_test1.mat', 'M0', 'ff', 'T1f', 'T2f', 'T1s', 'T2s', 'wf', 'wff', ...
	 'kappa', 'Tfree', 'Tg', 'alpha', 'beta', 'phi', 'spoil', 'result', ...
	 '-v7.3');

%% stfr2_test2
clear

M0 = [1; 0.8; 1.1];
ff = [0.2; 0.12; 0.25];
T1f = [400; 300; 350];
T1s = [1000; 900; 1500];
T2f = [20; 25; 30];
T2s = [80; 120; 300];
wf = [1; 2; 3] * 2*pi;
wff = [10; 20; 30] * 2*pi;
kappa = [1.1; 0.8; 1];
Tfree = 7.55;
Tg = 4.2;
alpha = 10 * 2*pi;
beta = 10 * 2*pi;
phi = 30 * 2*pi;
spoil = false;
result = stfr_2comp_forward_model(M0, ff, T1f, T2f, T1s, T2s, wf, wff, ...
	                              kappa, Tfree, Tg, alpha, beta, phi, ...
								  'spoil', spoil).';

save('stfr2_test2.mat', 'M0', 'ff', 'T1f', 'T2f', 'T1s', 'T2s', 'wf', 'wff', ...
	 'kappa', 'Tfree', 'Tg', 'alpha', 'beta', 'phi', 'spoil', 'result', ...
	 '-v7.3');

%% stfr2_test3
clear

M0 = 1;
ff = 0.2;
T1f = 400;
T1s = 1000;
T2f = 20;
T2s = 80;
wf = 1 * 2*pi;
wff = 10 * 2*pi;
kappa = 1.1;
Tfree = [7.45; 6];
Tg = [4.1; 3];
alpha = [11; 20] * 2*pi;
beta = [11; 15] * 2*pi;
phi = [31; 2] * 2*pi;
spoil = false;
result = stfr_2comp_forward_model(M0, ff, T1f, T2f, T1s, T2s, wf, wff, ...
	                              kappa, Tfree, Tg, alpha, beta, phi, ...
								  'spoil', spoil).';

save('stfr2_test3.mat', 'M0', 'ff', 'T1f', 'T2f', 'T1s', 'T2s', 'wf', 'wff', ...
	 'kappa', 'Tfree', 'Tg', 'alpha', 'beta', 'phi', 'spoil', 'result', ...
	 '-v7.3');

%% stfr2_test4
clear

M0 = [1; 0.8; 1.1];
ff = [0.2; 0.12; 0.25];
T1f = [400; 300; 350];
T1s = [1000; 900; 1500];
T2f = [20; 25; 30];
T2s = [80; 120; 300];
wf = [1; 2; 3] * 2*pi;
wff = [10; 20; 30] * 2*pi;
kappa = [1.1; 0.8; 1];
Tfree = [7.45; 6];
Tg = [4.1; 3];
alpha = [11; 20] * 2*pi;
beta = [11; 15] * 2*pi;
phi = [31; 2] * 2*pi;
spoil = true;
result = stfr_2comp_forward_model(M0, ff, T1f, T2f, T1s, T2s, wf, wff, ...
	                              kappa, Tfree, Tg, alpha, beta, phi, ...
								  'spoil', spoil).';

save('stfr2_test4.mat', 'M0', 'ff', 'T1f', 'T2f', 'T1s', 'T2s', 'wf', 'wff', ...
	 'kappa', 'Tfree', 'Tg', 'alpha', 'beta', 'phi', 'spoil', 'result', ...
	 '-v7.3');


function [M] = stfr_2comp_forward_model(M0, ff, T1f, T2f, T1s, T2s, wf, dwf, ...
                                    	kappa, Tfree, Tg, alpha, beta, phi, ...
                                    	varargin)

	% Set default option values
	arg.spoil = true;
	arg.TE    = Tfree / 2;

	% Use provided option values
	arg = use_varargin(arg, varargin);

	% Make sure the provided option values are valid
	check = {{arg.spoil, 'logical'}};
	check_input(check);

	% Grab the dimensions of latent/known parameters and scan parameters
	N = length(M0);
	P = length(Tfree);

	% Preallocate space for output
	M = zeros(N,P);

	% Fill output array
	for n = 1:N
		tmp = wf(n) + dwf(n);
		for p = 1:P
		    M(n,p) = stfr_2comp(phi(p), tmp * Tfree(p)/1000, wf(n) * ...
		                        Tfree(p)/1000, alpha(p), beta(p), ff(n), ...
		                        T1f(n), T1s(n), T2f(n), T2s(n), Tfree(p), ...
		                        Tg(p), M0(n), kappa(n), arg.TE(p), arg.spoil);
		end
	end

end

function [M] = stfr_2comp(phi, thetaff, thetafs, alpha, beta, ff, T1f, T1s, ...
                          T2f, T2s, Tfree, Tsg, M0, kappa, TE, useRFSpoiling)

    M =      ff  * stfr(phi, thetaff, alpha, beta, T1f, T2f, Tfree, Tsg, ...
                        M0, kappa, useRFSpoiling) * exp(-TE/T2f) ...
                        * exp(-1j * thetaff * TE / Tfree) + ...
        (1 - ff) * stfr(phi, thetafs, alpha, beta, T1s, T2s, Tfree, Tsg, ...
                        M0, kappa, useRFSpoiling) * exp(-TE/T2s) ...
                        * exp(-1j * thetafs * TE / Tfree);

end

function [M] = stfr(phi, thetaf, alpha, beta, T1, T2, Tfree, Tsg, M0, ...
		            kappa, useRFSpoiling)

    % Update alpha and beta according to kappa
    alpha = kappa * alpha;
    beta  = kappa * beta;

	% Check if RF spoiling should be employed. This affects the signal equation
	% that will be used.
	if useRFSpoiling

		Ts = Tsg;

		num = exp(-Ts/T1) * (1 - exp(-Tfree/T1)) * sin(alpha) * cos(beta) + ...
			  (1 - exp(-Ts/T1)) * sin(alpha);

		denom = 1 - exp(-Ts/T1 - Tfree/T2) * sin(alpha) * sin(beta) * ...
                cos(thetaf - phi) - exp(-Ts/T1 - Tfree/T1) * cos(alpha) * ...
                cos(beta);

		M = M0 * num / denom;

	else

		Tg = Tsg;
		Ef1 = exp(-Tfree/T1);
		Ef2 = exp(-Tfree/T2);
		Eg1 = exp(-Tg/T1);
		Eg2 = exp(-Tg/T2);

		a = -1i * Eg2 * (Ef2 * (-1 + Eg1 + (-1 + Ef1) * Eg1 * cos(beta)) * ...
            cos(thetaf - phi) * sin(alpha) + (Ef1 * (-1 + Eg1) + (-1 + ...
            Ef1) * cos(alpha)) * sin(beta) + 1i * Ef2 * (-1 + Eg1 + (-1 + ...
            Ef1) * Eg1 * cos(beta)) * sin(alpha) * sin(thetaf - phi));

		b = Eg2 * (Ef2 * ((-1 + Ef1) * Eg1 + (-1 + Eg1) * cos(beta)) * ...
            cos(thetaf - phi) * sin(alpha) - (-1 + Ef1 + Ef1 * (-1 + Eg1) * ...
			cos(alpha)) * sin(beta) + 1i * Ef2 * ((-1 + Ef1) * Eg1 + (-1 + ...
            Eg1) * cos(beta)) * sin(alpha) * sin(thetaf - phi));

		c = 1i * ((-1 + Eg1 + (-1 + Ef1) * Eg1 * cos(beta)) * sin(alpha) + ...
            Ef2 * Eg2^2 * (Ef1 * (-1 + Eg1) + (-1 + Ef1) * cos(alpha)) * ...
            sin(beta) * (cos(thetaf - phi) + 1i * sin(thetaf - phi)));

		d = Eg2 * (-Ef2 * (-1 + Ef1 * Eg1) * (1 + cos(alpha) * cos(beta)) * ...
            cos(thetaf - phi) + (Ef1 - Ef2^2 * Eg1) * sin(alpha) * sin(beta));

		e = Ef2 * (-1 + Ef1 * Eg1) * Eg2 * (cos(alpha) + cos(beta)) * ...
            sin(thetaf - phi);

		f = -1 + Ef1 * Ef2^2 * Eg1 * Eg2^2 + (Ef1 * Eg1 - Ef2^2 * Eg2^2) * ...
            cos(alpha) * cos(beta) + Ef2 * (Eg1 - Ef1 * Eg2^2) * ...
            cos(thetaf - phi) * sin(alpha) * sin(beta);

        s = sqrt(f^2 - d^2 - e^2);

		first  = -c / s;
		second = (a * d + b * e) / (d^2 + e^2);
		third  = (-f - s) / s;

		M = M0 * (first - second * third);

	end

end

function [arg] = use_varargin(arg, vararg)

    if mod(length(vararg), 2) ~= 0
        error('Options must consist of option-value pairs.');
    end
    for n = 1:2:length(vararg)
        if isfield(arg, vararg{n})
            arg.(vararg{n}) = vararg{n+1};
        else
            warning('Unknown option ''%s'' specified. Will be ignored.', ...
                    vararg{n});
        end
    end

end

function check_input(check)

    % Make sure check is a cell array
    if ~iscell(check)
        error('check must be a cell array.');
    end

    % Vectorize the inputs
    check = check(:);

    % Iterate through each check
    for n = 1:length(check)

        % Make sure the check is a cell
        if ~iscell(check{n})
            error('Each check must be a cell. Error at index %d.', n);
        end

        % Make sure each check has at least 2 elements
        if length(check{n}) < 2
            error(['Each check must have at least an input variable paired ' ...
                   'with a check to make. Error at index %d.'], n);
        end

        % Grab the current input and check
        in = check{n}{1};
        t  = check{n}(2:end);

        % Make sure the first entry in t is a char
        if ~ischar(t{1})
            error(['The second entry of check{n} must be of type char. ' ...
                   'Error at index %d.'], n);
        end

        % Check the first entry in t to see what check to make
        if strcmp(t{1}, 'struct')
            % Make sure in is a struct
            if ~isstruct(in)
                error('Expected struct at index %d.', n);
            end
            % See if struct fields were specified
            if length(t) > 1
                % Make sure fields are in a cell array
                if ~iscell(t{2})
                    warning(['Expected cell array as extra struct ' ...
                             'argument. Ignoring extra argument at index ' ...
                             '%d.'], n);
                else
                    % Iterate through each struct field
                    fields = t{2};
                    for f = 1:numel(fields)
                        % Make sure all entries in fields are chars
                        if ~ischar(fields{f})
                            warning(['Expected char for struct field. ' ...
                                     'Ignoring field number %d at index ' ...
                                     '%d.'], f, n);
                        else
                            % Make sure struct has the specified field
                            if ~isfield(in, fields{f})
                                error(['Expected field ' fields{f} ' in ' ...
                                       'struct at index %d.'], n);
                            end
                        end
                    end
                end
            end
            % Show a warning if extra type arguments were provided
            if length(t) > 2
                warning(['Ignoring extra struct arguments at index %d. ' ...
                         'Expected 1 or 2, but got %d.'], n, length(t));
            end
        elseif strcmp(t{1}, 'size')
            % Make sure the size to compare the input to was provided
            if length(t) < 2
                error(['''size'' type option must be accompanied with a ' ...
                       'size to compare the input to. Error at index %d.'], n);
            end
            % Make sure the size is a numeric array
            if ~isnumeric(t{2})
                error('size must be a numeric array. Error at index %d.', n);
            end
            % Convert t{2} into a 1-D array (if not already)
            s = t{2};
            s = s(:);
            % Make sure the input is the correct size
            if length(size(in)) ~= length(s) || sum(size(in)' ~= s) > 0
                error('Size mismatch at index %d.', n);
            end
            % Show a warning if extra type arguments were provided
            if length(t) > 2
                warning(['Ignoring extra size arguments at index %d. ' ...
                         'Expected 2, but got %d.'], n, length(t));
            end
        elseif strcmp(t{1}, 'length')
            % Make sure the length to compare the input to was provided
            if length(t) < 2
                error(['''length'' type option must be accompanied with a ' ...
                       'length to compare the input to. Error at index ' ...
                       '%d.'], n);
            end
            % Make sure the length is a scalar
            len = t{2};
            if ~isscalar(len)
                error('len must be a scalar value. Error at index %d.', n);
            end
            % Make sure the input is the correct length
            if length(in) ~= len
                error('Length mismatch at index %d.', n);
            end
            % Show a warning if extra type arguments were provided
            if length(t) > 2
                warning(['Ignoring extra length arguments at index %d. ' ...
                         'Expected 2, but got %d.'], n, length(t));
            end
        elseif strcmp(t{1}, 'logical')
            % Make sure the input is a logical
            if ~islogical(in)
                error('Expected logical at index %d.', n);
            end
            % Show a warning if extra type arguments were provided
            if length(t) > 1
                warning(['Ignoring extra logical arguments at index %d. ' ...
                         'Expected 1, but got %d.'], n, length(t));
            end
        elseif strcmp(t{1}, 'char')
            % Make sure the input is a char
            if ~ischar(in)
                error('Expected char at index %d.', n);
            end
            % Show a warning if extra type arguments were provided
            if length(t) > 1
                warning(['Ignoring extra char arguments at index %d. ' ...
                         'Expected 1, but got %d.'], n, length(t));
            end
        elseif strcmp(t{1}, 'cell')
            % Make sure the input is a cell
            if ~iscell(in)
                error('Expected cell at index %d.', n);
            end
            % Show a warning if extra type arguments were provided
            if length(t) > 1
                warning(['Ignoring extra cell arguments at index %d. ' ...
                         'Expected 1, but got %d.'], n, length(t));
            end
        elseif strcmp(t{1}, 'scalar')
            % Make sure the input is a scalar
            if ~isscalar(in)
                error('Expected scalar at index %d.', n);
            end
            % Show a warning if extra type arguments were provided
            if length(t) > 1
                warning(['Ignoring extra scalar arguments at index %d. ' ...
                         'Expected 1, but got %d.'], n, length(t));
            end
        else
            error(['Invalid type ' t{1} ' at index %d'], n);
        end

    end

end
