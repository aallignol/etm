// diff start at .begin() + 1

#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

using namespace arma;

RcppExport SEXP los_nocp(SEXP __times,
			 SEXP __tr_mat,
			 SEXP __surv,
			 SEXP __tau,
			 SEXP __aw,
			 SEXP __ratio) {

    Rcpp::NumericVector _times(__times), _surv(__surv), _tr_mat(__tr_mat);
    Rcpp::IntegerVector Dims = _tr_mat.attr("dim");
    const int aw = Rcpp::as<int>(__aw);
    const int ratio = Rcpp::as<int>(__ratio);
    const double tau = Rcpp::as<double>(__tau);
    const int lt = Dims[2];
    const int nstate = Dims[0];
    cube tr_mat(_tr_mat.begin(), nstate, nstate, lt, false);
    vec times(_times.begin(), _times.size(), false);
    vec surv(_surv.begin(), _surv.size(), false);
    
    // new stuffs we'll need
    mat::fixed<3, 3> I;
    I.eye();
    
    cube aj(nstate, nstate, lt); 
    
    for (int i = 0; i < lt; ++i) {
	tr_mat.slice(i) = tr_mat.slice(i) + I;
	aj.slice(i).eye();
    }
    
    vec los0(lt), los1(lt), a00(lt), a11(lt), a01(lt);;
    
    los0.fill(tau);
    los1.fill(tau);

    for (int t = (times.size() - 2); t >= 0; --t) {

	rowvec diff = times(span(t+2, lt)) - times(span(t+1, lt-1));

	for (int j = t; j < lt; ++j) {
	    aj.slice(j) = tr_mat.slice(t+1) * aj.slice(j);
	    a00[j] = aj.at(0, 0, j);
	    a01[j] = aj.at(0, 1, j);
	    a11[j] = aj.at(1, 1, j);
	}
	
	los0[t] = times[t+1] + as_scalar(diff * (a00(span(t, lt)) + a01(span(t, lt))));
	los1[t] = times[t+1] + as_scalar(diff * a11(span(t, lt)));
    }

    return Rcpp::List::create(Rcpp::Named("los0") = los0,
			      Rcpp::Named("los0") = los1);

}
    
    
