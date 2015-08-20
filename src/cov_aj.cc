#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

using namespace arma;

mat cov_dna(const mat & dna, const vec & nrisk, int d, int D);

RcppExport SEXP cov_aj(SEXP __time,
		       SEXP __est,
		       SEXP __nrisk,
		       SEXP __nevent,
		       SEXP __dna)

{

    Rcpp::NumericVector _time(__time), _est(__est), _nevent(__nevent), _dna(__dna);
    Rcpp::IntegerVector dims = _est.att("dim");
    Rcpp::NumericMatrix _nrisk(__nrisk);
    
    const int lt = _time.size();
    const int nstate = dims(0);
    const int D = pow(nstate, 2);

    cube est(_est.begin(), nstate, nstate, lt, false);
    cube nevent(_nevent.begin(), nstate, nstate, lt, false);
    cube dna(_dna.begin(), nstate, nstate, lt, false);
    mat nrisk(_nrisk.begin(), lt, nstate, false);
    vec time(_time.begin(), _time.size(), false);

    cube cov_etm(D, D, lt);
    cov_etm.zeros();

    mat I(nstate, nstate, fill::eye);
    mat II(D, D, fill::eye);
    mat cov_deltaNA(D, D, fill::zeros);

    // first iteration
    cov_deltaNA = cov_dna(dna.slice(0),
			  nrisk.row(0),
			  nstate,
			  D);
    cov_etm.slice(0) = II * cov_deltaNA * II;
    
    for (int t = 1; t < lt, ++t) {
	
	mat temp_dna(dna.slice(t).begin(), nstate, nstate, false);
	mat temp_est(est.slice(t).begin(), nstate, nstate, false);
	
	cov_deltaNA = cov_dna(temp_dna,
			      nrisk.row(t),
			      nstate,
			      D);
	
	cov_etm.slice(t) = kron((I + temp_dna).t(), II) * cov_etm.slice(t - 1) * kron((I+temp_dna),II) +
	    kron(I, temp_est) * cov_deltaNA * kron(I, temp_est.t());
	
    }

    return cov_etm;
}
