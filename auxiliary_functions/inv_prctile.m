% PERCENTILE = inv_prctile(DATA,VALUE)
%
% Basically the opposite of MATLAB's defaul "prctile". "inv_prctile" will
% give the approximate percentile of a value you give it compared to the
% given data set:
%
% IN:   DATA  = Vector of data with an arbitrary distribution. Note that,
%               unlike the built-in function "prctile", this function has
%               no option to format "DATA" as a matrix. If "DATA" is a
%               matrix and you want to use all elements of it, you will
%               have to reshape it.
% IN:   VALUE = Some value whose percentile in "DATA" is sought. If "VALUE"
%               is greater than max(DATA), then "PERCENTILE" = 100. If it's
%               less than min(DATA), "PERCENTILE" = 0.
%
% OUT:  PERCENTILE = The exact percentile in "DATA" which "VALUE" has if it
%                    is an element of "DATA", or the percentile it would
%                    have if it were in data.
%
%
%
% NOTE: If there are multiple occurances of "VALUE" in "DATA", then the
% value of "PERCENTILE" given will be the "most middle" possible, e.g.:
%
% inv_prctile([1 2 3 3 3],3) gives 70. See the folowing example for
% clarification for why this isn't 60:
%
% [[10:10:100]; prctile([1 2 3 3 3],[10:10:100])]
%
% Basically (I think), it's the average of 60 (because 60% are below the
% middle "3") and 80 (because 100%-80% = 20% are above the middle "3")
% 
% For short "DATA" vectors, the output might not make a lot of sense; feel
% free to tweak this function as necessary.

function PERCENTILE = inv_prctile(DATA,VALUE)

if VALUE < min(DATA)
    PERCENTILE = 0;
elseif VALUE > max(DATA)
    PERCENTILE = 100;
else
    
    if isrow(DATA)
        DATA = DATA';
    else
    end
    
    DATA = sort(DATA);
    
    if sum(DATA == VALUE) == 1 % i.e. if VALUE is in DATA once
        PERCENTILE = 100*(dsearchn(DATA,VALUE)/length(DATA)) - 50/length(DATA);
    elseif sum(DATA == VALUE) > 1 % i.e. if VALUE is in DATA more than once
        % This is janky, but it should work:
        N_occurances = sum(DATA == VALUE);
        dDATA = abs(diff(DATA));
        dDATA(dDATA == 0) = nan;
        min_diff = min(dDATA); % to prevent overshoot
        min_prc = 100*(dsearchn(DATA,VALUE)/length(DATA));
        for ii = 1:(N_occurances-1) % slightly shrink all but one occurance of VALUE
            DATA(dsearchn(DATA,VALUE)) = VALUE - 0.5*min_diff;
        end
        max_prc = 100*(dsearchn(DATA,VALUE)/length(DATA));
        PERCENTILE = mean([min_prc max_prc]) - 50/length(DATA);
    else % i.e. if VALUE is not in DATA
        PERCENTILE = 100*(dsearchn(DATA,VALUE)/length(DATA)) - 50/length(DATA);
    end
    
end
end

%% Test for bias (it passes)

% for ii = 1:1000
%     foo(ii) = inv_prctile(randn(10000,1),0);
% end
% figure;histogram(foo)
% disp(['standard error of mean = ',num2str(std(foo)/sqrt(length(foo)))])
% disp(['difference between observed mean and true mean = ',num2str(abs(50 - mean(foo)))])
