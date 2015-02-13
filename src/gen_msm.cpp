#include <RcppArmadillo.h>
//#include <gperftools/profiler.h>

using namespace arma;

cube prodint(const cube & dna, int nstate, int ltimes);
cube deltaNA(const cube & nev, const mat & nrisk, int nstate, int ltimes);


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
    mat nrisk(lt, nstate); nrisk.zeros();
    cube nev(nstate, nstate, lt); nev.zeros();
//    ProfilerStart("/tmp/gen_msm.prof");

    for (int i=0; i < n; ++i) {
	for (int t=0; t < lt; ++t) {

	    if (entry[i] < times[t] && exit[i] >= times[t]) {
		nrisk.at(t, from[i] - 1) += 1;
	    }
	    if (exit[i] == times[t] && to[i] != 0) {
		nev.at(from[i] - 1, to[i] - 1, t) += 1;
		break;
	    }
	}
    }


    cube dna = deltaNA(nev, nrisk, nstate, lt);
    cube est = prodint(dna, nstate, lt);	
//    ProfilerStop();
    
    return Rcpp::List::create(Rcpp::Named("n.risk") = nrisk,
			      Rcpp::Named("n.event") = nev,
			      Rcpp::Named("dna") = dna,
			      Rcpp::Named("est") = est);
	
}
	

