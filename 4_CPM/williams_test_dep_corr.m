function [tval, pval] = williams_test_dep_corr(r12, r13, r23, n)
% Williams/Steiger test for comparing two dependent correlations
% that share one variable (Steiger, 1980; "Williams' t").
% r12 = corr(X1,Y), r13 = corr(X2,Y), r23 = corr(X1,X2), n = sample size.
% Returns t (df = n-3) and two-sided p-value.

    if n < 6
        tval = NaN; pval = NaN; return;
    end

    % Guard numerical issues
    r12 = max(min(r12, 0.999999), -0.999999);
    r13 = max(min(r13, 0.999999), -0.999999);
    r23 = max(min(r23, 0.999999), -0.999999);

    % Williams test statistic (dependent/overlapping case)
    num = (r12 - r13) * sqrt((n - 3) * (1 + r23));
    den_inner = 2 * (1 - r23) * (1 + r23 - r12^2 - r13^2 - 2*r12*r13*r23);
    den = sqrt(den_inner);

    if ~(den > 0)
        tval = NaN; pval = NaN; return;
    end

    tval = num / den;
    df = n - 3;
    pval = 2 * (1 - tcdf(abs(tval), df));
end
