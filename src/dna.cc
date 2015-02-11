#include <RcppArmadillo.h>

// using namespace std;
using namespace arma;

RcppExport SEXP dna(SEXP _times,
		    SEXP _entry,
		    SEXP _exit,
		    SEXP _from,
		    SEXP _to,
		    SEXP _nstate)

{


    Rcpp::NumericVector entry(_entry), exit(_exit), times(_times);
    Rcpp::IntegerVector from(_from), to(_to);

    const int lt = times.size();
    const int n = entry.size();
    const int nstate = Rcpp::as<int>(_nstate);

    // define the matrices we need
    icube nrisk(nstate, nstate, lt); nrisk.zeros();
    icube nev(nstate, nstate, lt); nev.zeros();
    cube dna(nstate, nstate, lt); dna.zeros();

    for (int t=0; t < lt; ++t) {
	for (int i=0; i < n; ++i) {
	    if (entry[i] < times[t] && exit[i] >= times[t])
	    	nrisk.slice(t).row(from[i] - 1) += 1;
	    if (exit[i] == times[t] && to[i] != 0)
	    	nev.at(from[i] - 1, to[i] - 1, t) += 1;
	}
	    
	mat n = conv_to<mat>::from(nev.slice(t));
	mat y = conv_to<mat>::from(nrisk.slice(t));
	mat tmp(dna.slice(t).begin(), nstate, nstate, false);
	
	tmp = n / y;
	tmp.elem(find_nonfinite(tmp)).zeros();
	vec d = sum(tmp, 1);
	tmp.diag() = -d;
    }
    
    return Rcpp::List::create(Rcpp::Named("n.risk") = nrisk,
			      Rcpp::Named("n.event") = nev,
			      Rcpp::Named("dna") = dna);
	
	}
	
    
