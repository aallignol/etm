#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

using namespace arma;

RcppExport SEXP cov_aj(SEXP __time,
		       SEXP __est,
		       SEXP __nrisk,
		       SEXP __nevent,
		       SEXP __dna)

{

    Rcpp::NumericVector _time(__time), _est(__est), _nevent(__nevent), _dna(__dna);
    Rcpp::IntegerVector dims = _est.att("dim");
    
    const int lt = _time.size();
    const int nstate = dims(0);

    cube est(_est.begin(), nstate, nstate, lt, false);
    cube nevent(_nevent.begin(), nstate, nstate, lt, false);
    cube dna(_dna.begin(), nstate, nstate, lt, false);
    vec time(_time.begin(), _time.size(), false);
    
}
