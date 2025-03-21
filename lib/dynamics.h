//  last update: 21/03/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: General dynamics and Cherenkov specifcs kinematic helper fuynctions
//  TODO: -

//  Kinematics
//  --- Get Beta
inline float get_relativistic_beta(float gamma) { return TMath::Sqrt(1 - gamma * gamma); }
inline float get_relativistic_beta(float mass, float momentum) { return 1. / sqrt(1 + (mass * mass) / (momentum * momentum)); }
//  --- Get Gamma
inline float get_relativistic_gamma(float beta) { return 1. / (TMath::Sqrt(1 - beta * beta)); }
inline float get_relativistic_gamma(float mass, float momentum) { return get_relativistic_gamma(get_relativistic_beta(mass, momentum)); }

//  Cherenkov
//  --- Get Theta_c
inline float get_theta_cherenkov(float beta, float ref_index) { return TMath::ACos(1. / (beta * ref_index)); }
//  --- Get thresholds
inline float get_ch_thr_relativistic_beta(float ref_index) { return 1. / ref_index; }
inline float get_ch_thr_momentum(float ref_index, float mass) { return mass * get_ch_thr_relativistic_beta(ref_index) * get_relativistic_gamma(get_ch_thr_relativistic_beta(ref_index)); }
inline float get_ch_thr_mass(float ref_index, float momentum) { return momentum / (get_ch_thr_relativistic_beta(ref_index) * get_relativistic_gamma(get_ch_thr_relativistic_beta(ref_index))); }
//  --- Photons production

//  --- TBS
// inline float get_expected_radius(float mass, float momentum, float ref_index, float arm_length) { return TMath::Tan(get_theta_cherenkov(get_relativistic_beta(mass, momentum), ref_index)) * arm_length; }
// inline float get_expected_photons(float ref_index, float mass, float momentum, int z_charge = 1, float lambda_low = 350.e-9, float lambda_high = 500.e-9)