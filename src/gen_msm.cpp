#include <RcppArmadillo.h>
//#include <gperftools/profiler.h>

using namespace arma;

cube prodint(cube & dna, int nstate, int ltimes);



RcppExport SEXP gen_msm(SEXP _times,
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

    //const int cova = Rcpp::as<int>(_covariance);

    // define the matrices we need
    mat nrisk(nstate, lt); nrisk.zeros();
    cube nev(nstate, nstate, lt); nev.zeros();
    cube dna(nstate, nstate, lt); dna.zeros();
//    ProfilerStart("/tmp/gen_msm.prof");

    for (int i=0; i < n; ++i) {
	for (int t=0; t < lt; ++t) {

	    if (entry[i] < times[t] && exit[i] >= times[t]) {
		nrisk.at(from[i] - 1, t) += 1;
	    }
	    if (exit[i] == times[t] && to[i] != 0) {
		nev.at(from[i] - 1, to[1] - 1, lt) += 1;
		break;
	    }
	}
    }
    
    for (int t=0; t < lt; ++t) {
	for (int i=0; i < nstate; ++i) {
	    for (int j=0; j<nstate; ++j) {
		if (nrisk.at(i, t) != 0) {
		    dna.at(i, j, t) = nev.at(i, j, t) / nrisk.at(i, t);
		}
	    }
	}
    }
			
     	// mat tmp(dna.slice(t).begin(), nstate, nstate, false);
	// mat dn(nev.slice(t).begin(), nstate, nstate, false);
	// mat y(nrisk.slice(t).begin(), nstate, nstate, false);

    	// tmp = dn / y;
	    
    	// tmp.elem(find_nonfinite(tmp)).zeros();
    	// vec d = sum(tmp, 1);
    	// tmp.diag() = -d;
    

    cube est = prodint(dna, nstate, lt);	
//    ProfilerStop();
    
    return Rcpp::List::create(Rcpp::Named("n.risk") = nrisk,
			      Rcpp::Named("n.event") = nev,
			      Rcpp::Named("dna") = dna,
			      Rcpp::Named("est") = est);
	
}
	

