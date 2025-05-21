function H = matrixEntropy(M)
% MATRIXENTROPY  Compute the Shannon entropy of the entries of a matrix.

    p = M(:);
    p = p ./ sum(p);          % normalise so sum(p)=1
    epsH = 1e-12;
    H = -sum( p .* log2(p + epsH) );
end
