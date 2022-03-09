val binom_test :
  n_successes:int ->
  n_trials:int ->
  p:float ->
  float


type named_pval = {
  id : string ;
  pval : float ;
}

val false_discovery_rates : named_pval list -> named_pval list
