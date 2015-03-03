// diff start at .begin() + 1

#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

// using namespace arma;

// RcppExport SEXP los_nocp(SEXP __times,
// 			 SEXP __tr_mat,
// 			 SEXP __surv,
// 			 SEXP __tau,
// 			 SEXP __aw,
// 			 SEXP __ratio) {

//     Rcpp::NumericVector _times(__times), _surv(__surv), _tr_mat(__tr_mat);
//     Rcpp::IntegerVector Dims = _tr_mat.attr("dim");
//     const int aw = Rcpp::as<int>(__aw);
//     const int ratio = Rcpp::as<int>(__ratio);
//     const double = Rcpp::as<double>(__tau);
//     const int lt = Dims[2];
//     const int nstate = Dims[0];
//     cube tr_mat(_tr_mat.begin(), nstate, nstate, lt, false);
//     vec times(_times.begin(), _times.size(), false);
//     vec surv(_surv.begin(), _surv.size(), false);
    
//     // new stuffs we'll need
//     mat::fixed<3, 3> I;
//     I.eye();
    
//     cube aj(nstate, nstate, lt); 
    
//     for (int i = 0; i < lt; ++i) {
// 	tr_mat.slice(i) = tr_mat.slice(i) + I;
// 	aj.slice(i).eye();
//     }
    
//     vec los0(lt), los1(lt);
//     los0.fill(tau);
//     los1.fill(tau);

//     for (int t = (times.size() - 2); t >= 0; --t) {
	
	
	
//     }

    
    
    
