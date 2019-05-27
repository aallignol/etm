// #define ARMA_NO_DEBUG

// pbl in the new stuff for "1" event: Probably
// does something weird when last data point with ties

#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

using namespace arma;

cube prodint(const cube & dna, int nstate, int ltimes);
cube deltaNA(const cube & nev, const mat & nrisk, int nstate, int ltimes);
cube deltaNA_LY(const cube & nev, const mat & nrisk, const mat & which_compute, int nstate, int ltimes);


RcppExport SEXP gen_msm(SEXP _times,
			SEXP _entry,
			SEXP _exit,
			SEXP _from,
			SEXP _to,
			SEXP _nstate,
			SEXP _const_modif)

{

    Rcpp::NumericVector __entry(_entry), __exit(_exit), times(_times);
    Rcpp::IntegerVector __from(_from), __to(_to);
    Rcpp::IntegerMatrix __const_modif(_const_modif);

    // ProfilerStart("/tmp/gen_msm.prof");
    
    vec Tentry(__entry.begin(), __entry.size(), false);
    vec Texit(__exit.begin(), __exit.size(), false);
    ivec Tfrom(__from.begin(), __from.size(), false);
    ivec Tto(__to.begin(), __to.size(), false);
    imat const_modif(__const_modif.begin(), __const_modif.nrow(), __const_modif.ncol(), false);
    
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
    
    // define the matrices we need
    mat nrisk(lt, nstate, fill::zeros);
    
    cube nev(nstate, nstate, lt); nev.zeros();
    cube dna(nstate, nstate, lt); dna.zeros();


    // the entries
    int l = 0;
    for (int j = 0; j < n; ++j) {
    	if (entry[j] < times[l]) {
    	    nrisk(l, from_entry[j] - 1) += 1;
    	} else {
	    while (l < lt && entry[j] >= times[l] && entry[j] <= max(times)) {
		++l;
	    }
	    nrisk(l, from_entry[j] - 1) += 1;
	}
    }

    // the events
    int t = 0;

    // only one event time. And different code paths if 1 data point
    // or more
    if (lt == 1) {
	if (to[0] != 0) nev(from_exit[0] - 1, to[0] - 1, 0) += 1;
	if (n > 1) {
	    for (int i = 1; i<n; ++i) {
		if (exit[i] == exit[i - 1]) {
		    if (to[i] != 0) nev(from_exit[i] - 1, to[i] - 1, t) += 1;
		} else {
		    break;
		}
	    }
	}
	
    } else {

	// A code path for first event(s) censored
	int ii = 0;
	while (to[ii] == 0) {
	    nrisk(0, from_exit[ii] - 1) -= 1;
	    ++ii;
	}
	
	if (ii == 0) {
	    if (to[0] != 0) nev(from_exit[0] - 1, to[0] - 1, 0) += 1;
	    if (n > 1) {
		nrisk(1, from_exit[0] - 1) -= 1;
	    }
	    ii = 1;
	} else {
	    if (to[ii] != 0) nev(from_exit[ii] - 1, to[ii] - 1, 0) += 1;
	    nrisk(1, from_exit[ii] - 1) -= 1;
	    ++ii;
	}
	
	for (int i = ii; i<n; ++i) {
	
	    if (exit[i] == exit[i - 1]) {
		if (to[i] != 0) nev(from_exit[i] - 1, to[i] - 1, t) += 1;
		if (t < lt - 1) nrisk(t + 1, from_exit[i] - 1) -= 1;
	    } else {
		if (t < lt - 1) {
		    if (exit[i] == times[t+1]) {
			++t;
			if (to[i] != 0) nev(from_exit[i] - 1, to[i] - 1, t) += 1;
			if (t < lt - 1) nrisk(t + 1, from_exit[i] - 1) -= 1;
		    } else {
			nrisk(t + 1, from_exit[i] - 1) -= 1;
		    }
		}
		else {
		    if (exit[i] > times[t]) break;
		    if (to[i] != 0) nev(from_exit[i] - 1, to[i] - 1, t) += 1;
		}
	    }
	}
    }
    mat y = cumsum(nrisk);

    irowvec cc = const_modif.row(0);
    if (any(cc)) {
	umat tmp = (y >= const_modif);
	mat  which_compute = conv_to<mat>::from(tmp);
	// Nelson-Aalen (the increments) Lai and Ying
	dna = deltaNA_LY(nev, y, which_compute, nstate, lt);
    }
    else {
	// Nelson-Aalen (the increments) original
	dna = deltaNA(nev, y, nstate, lt);
    }
    
    cube est = prodint(dna, nstate, lt);	
    // ProfilerStop();
    
    return Rcpp::List::create(Rcpp::Named("n.risk") = y,
			      Rcpp::Named("n.event") = nev,
			      Rcpp::Named("dna") = dna,
			      Rcpp::Named("est") = est,
			      Rcpp::Named("time") = times);

}
	

