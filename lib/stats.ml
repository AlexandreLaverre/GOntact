let binom_test ~n_successes ~n_trials ~p =
  1. -. Gsl.Cdf.binomial_P ~k:(n_successes - 1) ~n:n_trials ~p

let%expect_test _ =
  print_endline {|
> binom.test(4, 5, 0.2,alternative="greater")

       Exact binomial test

data:  4 and 5
number of successes = 4, number of trials = 5, p-value = 0.00672
alternative hypothesis: true probability of success is greater than 0.2
95 percent confidence interval:
 0.3425917 1.0000000
sample estimates:
probability of success
                   0.8
|} ;
  Printf.printf "pvalue = %f\n" (binom_test ~n_successes:4 ~n_trials:5 ~p:0.2) ;
  [%expect {|
    > binom.test(4, 5, 0.2,alternative="greater")

           Exact binomial test

    data:  4 and 5
    number of successes = 4, number of trials = 5, p-value = 0.00672
    alternative hypothesis: true probability of success is greater than 0.2
    95 percent confidence interval:
     0.3425917 1.0000000
    sample estimates:
    probability of success
                       0.8

    pvalue = 0.006720 |}]
