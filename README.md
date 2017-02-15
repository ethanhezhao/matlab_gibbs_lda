# matlab_gibbs_lda

LDA with Gibbs samplings is implemented with MATLAB.

The hyper-parameter sampling of the Dirichlet distribution on document topic is based on Mallet, but re-implemented with pure MATLAB codes.

Input data format:

A N by 3 matrix X, for a row i, X(i,1) is the document id, X(i,2) is the token id, X(i,3) is the word counts for this token in this document.
