#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

//' Fast C++ implementation of the ASCAT distance grid calculation (BAF-only distance)
//' This avoids the memory explosion of creating large matrices in R and the overhead of forking.
//' @noRd
// [[Rcpp::export]]
NumericVector calculate_ascat_dist_matrix_cpp(
    NumericVector s_b, 
    NumericVector s_r, 
    NumericVector s_len, 
    NumericVector rho_vec, 
    NumericVector psi_vec,
    double gamma_param) {
    
    int n_seg = s_b.size();
    int n_grid = rho_vec.size();
    NumericVector results(n_grid);
    
    // PRE-CALCULATION: Pre-calculate segment-specific logR factors
    std::vector<double> logR_factors(n_seg);
    
    for (int i = 0; i < n_seg; ++i) {
        logR_factors[i] = std::pow(2.0, s_r[i] / gamma_param);
    }
    
    // Iterative calculation: One grid point at a time
    for (int g = 0; g < n_grid; ++g) {
        double rho = rho_vec[g];
        double psi = psi_vec[g];
        double two_one_minus_rho = 2.0 * (1.0 - rho);
        double one_minus_rho = 1.0 - rho;
        
        double total_weighted_dist = 0.0;
        
        for (int i = 0; i < n_seg; ++i) {
            
            double scale_factor = logR_factors[i] * (two_one_minus_rho + rho * psi);
            
            // Optimization: Find the best integer combination (k=1..4)
            double nMaj_raw = (rho - 1.0 + s_b[i] * scale_factor) / rho;
            // nMinor logic derived from nMajor + nMinor = total
            // But here we calculate independently based on BAF
            double nMin_raw = (rho - 1.0 + (1.0 - s_b[i]) * scale_factor) / rho;
            
            if (nMaj_raw < 0.01) nMaj_raw = 0.01;
            if (nMin_raw < 0.01) nMin_raw = 0.01;
            
            double J_f = std::floor(nMaj_raw);
            double J_c = std::ceil(nMaj_raw);
            double N_f = std::floor(nMin_raw);
            double N_c = std::ceil(nMin_raw);
            
            double best_d = 1e18; // Infinity
            double best_mu = 0.0;
            
            // 4 Integer combinations: (F, C), (C, C), (F, F), (C, F)
            double nMaj_opts[4] = {J_f, J_c, J_f, J_c};
            double nMin_opts[4] = {N_c, N_c, N_f, N_f};
            
            for (int k = 0; k < 4; ++k) {
                double nMaj = nMaj_opts[k];
                double nMin = nMin_opts[k];
                
                double denom = two_one_minus_rho + rho * (nMaj + nMin);
                if (denom < 1e-10) denom = 1e-10;
                
                double mu = (one_minus_rho + rho * nMaj) / denom;
                double dist = std::abs(mu - s_b[i]);
                
                if (dist < best_d) {
                    best_d = dist;
                    best_mu = mu;
                }
            }
            
            double diff = s_b[i] - best_mu;
            total_weighted_dist += (diff * diff) * s_len[i];
        }
        
        results[g] = total_weighted_dist;
        
        // Progress reporting every 500 grid points
        if ((g + 1) % 500 == 0 || g == n_grid - 1) {
            Rcpp::checkUserInterrupt(); // Allow user to cancel
        }
    }
    
    return results;
}
