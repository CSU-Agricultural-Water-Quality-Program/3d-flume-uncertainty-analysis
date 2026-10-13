data{
     real sigma_H;
     int J;
     vector[30] H_obs;
     int N;
     vector[30] H_upper;
     vector[30] H_lower;
    array[30] int F;
     vector[30] Q_obs;
     real cv_Q_prior_log_sd;
     real cv_Q_prior_median;
     real sigma_Q_add_prior_median;
}
parameters{
     vector<lower=H_lower,upper=H_upper>[N] H_true;
     vector[N] Q_true;
     vector<lower=0>[J] a;
     vector<lower=0>[J] b;
     vector<lower=0>[J] sigma_proc;
     real<lower=0> cv_Q;
     vector<lower=0>[J] sigma_Q_add;
}
model{
     vector[30] mu;
    sigma_Q_add ~ lognormal( log(sigma_Q_add_prior_median) , 0.5 );
    cv_Q ~ lognormal( log(cv_Q_prior_median) , cv_Q_prior_log_sd );
    sigma_proc ~ exponential( 1 );
    b ~ lognormal( log(2.58) , 0.5 );
    a ~ lognormal( log(1.55) , 0.5 );
    Q_obs ~ normal( Q_true , sqrt(square(cv_Q * Q_true) + square(sigma_Q_add[F])) );
    for ( i in 1:30 ) {
        mu[i] = a[F[i]] * H_true[i]^b[F[i]];
    }
    Q_true ~ normal( mu , sigma_proc[F] );
    H_true ~ uniform( H_lower , H_upper );
}


