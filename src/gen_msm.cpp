#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

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

    Rcpp::NumericVector __entry(_entry), __exit(_exit), times(_times);
    Rcpp::IntegerVector __from(_from), __to(_to);

    // ProfilerStart("/tmp/gen_msm.prof");
    
    vec Tentry(__entry.begin(), __entry.size(), false);
    vec Texit(__exit.begin(), __exit.size(), false);
    ivec Tfrom(__from.begin(), __from.size(), false);
    ivec Tto(__to.begin(), __to.size(), false);

    // vec entry(Tentry), exit(Texit);
    // ivec to(Tto), from_entry(Tfrom), from_exit(Tfrom);
    
    // do some sorting
    std::sort(times.begin(), times.end());
    uvec ind_entry = sort_index(Tentry);
    uvec ind_exit = sort_index(Texit);
    vec entry = Tentry.elem(ind_entry);
    ivec from_entry = Tfrom.elem(ind_entry);
    vec exit = Texit.elem(ind_exit);
    ivec to = Tto.elem(ind_exit);
    ivec from_exit = Tfrom.elem(ind_exit);
    
    const int lt = times.size();
    const int n = entry.size();
    const int nstate = Rcpp::as<int>(_nstate);

    //const int cova = Rcpp::as<int>(_covariance);
    
    // define the matrices we need
    mat nrisk(lt, nstate); nrisk.zeros();
    cube nev(nstate, nstate, lt); nev.zeros();


    if (to[0] != 0) nev.at(from_exit[0] - 1, to[0] - 1, 0) += 1;
    nrisk.at(1, from_exit[0] - 1) -= 1;
	
    // the events
    int t = 0;
    for (int i = 1; i<n; ++i) {
	
	if (exit[i] == exit[i - 1]) {
	    if (to[i] != 0) nev.at(from_exit[i] - 1, to[i] - 1, t) += 1;
	    if (t < lt - 1) nrisk.at(t + 1, from_exit[i] - 1) -= 1;
	} else {
	    ++t;
	    if (to[i] != 0) nev.at(from_exit[i] - 1, to[i] - 1, t) += 1;
	    if (t < lt - 1) nrisk.at(t + 1, from_exit[i] - 1) -= 1;
	}
    }

    // the entries
    int l = 0;
    for (int j = 0; j < n; ++j) {
    	if (entry[j] < times[l]) {
    	    nrisk.at(l, from_entry[j] - 1) += 1;
    	}
    	else {
	    do {
		++l;
		// nrisk(l, span::all) = nrisk(l - 1, span::all);
	    } while (entry[j] >= times[l]);
	    nrisk.at(l, from_entry[j] - 1) += 1;
	}
    }
    
    mat y = cumsum(nrisk);
	
    cube dna = deltaNA(nev, y, nstate, lt);
    cube est = prodint(dna, nstate, lt);	
    // ProfilerStop();
    
    return Rcpp::List::create(Rcpp::Named("n.risk") = y,
			      Rcpp::Named("n.event") = nev,
			      Rcpp::Named("dna") = dna,
			      Rcpp::Named("est") = est,
			      Rcpp::Named("time") = times);
	
}
	

