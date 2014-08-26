#include <Rcpp.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string>

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;

// this file is adapted from Helge Roiders source code for TRAP
//
// this function slides the PWM over the sequence and computes for each
// position the number of expected TF binding events.
// the results are written to affinities which is an array of length
// seqlength - motiflength + 1
void affinity(Rcpp::NumericVector pwm, int motiflength, std::string sequence, int seqlength, double Rmax, double lambda, double* affinities, int bothstrands) {
  int n;

  /* some debugging here */
  /*
    printf("motiflength = %i\n", motiflength);
    printf("seqlength = %i\n", seqlength);
    printf("Rmax = %f, lambda = %f\n", Rmax, lambda);
    printf("sequence = '%s'\n", sequence);
    printf("PWM by columns = ");
    for (n=0; n<motiflength; n++)
    printf("%f", pwm[n]);
    printf("\n");
  */

  /* matrix for reverse strand */
  Rcpp::NumericVector complement(4 * motiflength);
  for (int m = 0; m < motiflength; m++) { //create complement to matrix
    complement[(motiflength-m-1) * 4 + A] = pwm[m * 4 + T];
    complement[(motiflength-m-1) * 4 + T] = pwm[m * 4 + A];
    complement[(motiflength-m-1) * 4 + G] = pwm[m * 4 + C];
    complement[(motiflength-m-1) * 4 + C] = pwm[m * 4 + G];
  }

  double P_combined = 0; //only palindrome correction
  //double P_uncorrected = 0; //uncorrected expected count

  int illegalBase = 0;

  /* LOOP OVER SEQUENCE */
  for(n = 0; n < seqlength - motiflength + 1; n++){
    double dE_forward = 0;
    double dE_compl = 0;

    /* LOOP OVER MOTIF */
    for(int m = 0; m < motiflength; m++) {
      int BASE = 0;
      illegalBase = 0;
      switch(sequence.at(n + m)) {
      case 'A':
        BASE = 0;break;
      case 'a':
        BASE = 0;break;
      case 'C':
        BASE = 1;break;
      case 'c':
        BASE = 1;break;
      case 'G':
        BASE = 2;break;
      case 'g':
        BASE = 2;break;
      case 'T':
        BASE = 3;break;
      case 't':
        BASE = 3;break;
      case 'N':
        illegalBase = 1;break;
      case 'n':
        illegalBase = 1;break;
      default:
        Rprintf("ILLEGAL CHARACTER FOUND!\n Only A/a, T/t, G/g, C/c, N/n are allowed!\n");
        illegalBase = 1;
      }

      dE_forward = dE_forward + pwm[m * 4 + BASE];
      dE_compl = dE_compl + complement[m * 4 + BASE];

    } /* loop over motif */

    /* CALCULATE P(BOUND) FOR CURRENT SITE */
    if(illegalBase == 0) {
      double product = Rmax * exp(-1*dE_forward);
      double P_bound_F = product/(1 + product);

      product = Rmax * exp(-1*dE_compl);
      double P_bound_C = product/(1 + product);

      /* correction for forward and reverse strand */
      if (bothstrands > 0) {
	affinities[n] = P_bound_F + (1 - P_bound_F) * P_bound_C;
      } else {
	affinities[n] = P_bound_F;
      }
      P_combined = P_combined + affinities[n];
      /* P_uncorrected = P_uncorrected + P_bound_F + P_bound_C;
         printf("%e\t%e\t%e\t%e\t%e\n", P_bound_F, P_bound_C, affinities[n], P_combined, P_uncorrected); */
    } else {
      affinities[n] = 0;
    }
  } /* loop over sequence */
}

//' Calculate the TRAP affinity for one or more sequences (one score per sequence)
//'
//' @param pwm The position weight matrix
//' @param motiflength The length of the PWM
//' @param sequence a list of one or more genomic sequences
//' @param seqlength the lengths of the sequences
//' @param Rmax rmax
//' @param lambda lambda
//' @param bothstrands should the affinity be calcuated on both strands (DNA) or just one (RNA)
//' @return a numeric vector of the scores for each sequence
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector R_affinity_sum_multi(Rcpp::NumericVector pwm, int motiflength, std::vector< std::string > sequence, Rcpp::IntegerVector seqlength, double Rmax, double lambda, int bothstrands) {
  int n;
  double sum;

  int len = sequence.size();
  Rcpp::NumericVector sums(len);

  for (int i = 0; i < len; i++) {
    sum = 0;
    Rcpp::NumericVector affinities(seqlength[i] - motiflength + 1);
    affinity(pwm, motiflength, sequence.at(i), seqlength[i], Rmax, lambda, affinities.begin(), bothstrands);
    for (n = 0; n < seqlength[i] - motiflength + 1; n++) {
      sum += affinities[n];
    }
    sums[i] = sum;
  }
  return(sums);
}

//' Calculate the TRAP affinity for each position along one or more sequences
//'
//' @inheritParams R_affinity_sum_multi
//' @return a list where each element is a numeric vector containing the
//' scores along the length of a single sequence.
//' @export
// [[Rcpp::export]]
Rcpp::List R_affinity_multi(Rcpp::NumericVector pwm, int motiflength, std::vector< std::string > sequence, Rcpp::IntegerVector seqlength, double Rmax, double lambda, int bothstrands) {
  int len = sequence.size();
  Rcpp::List affinity_list(len);

  for (int i = 0; i < len; i++) {
    affinity_list[i] = Rcpp::NumericVector(seqlength[i] - motiflength + 1);
    Rcpp::NumericVector affinities(seqlength[i] - motiflength + 1);
    affinity(pwm, motiflength, sequence.at(i), seqlength[i], Rmax, lambda, affinities.begin(), bothstrands);
    affinity_list[i] = affinities;
  }
  return(affinity_list);
}
