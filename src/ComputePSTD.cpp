//Includes/namespaces
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
//' @title
//' ComputeDistancesPSTD
//' @description
//' Function to compute the probabililistic distance between two gene tree probabilities for the Hellinger (dH, first), KL (dKL, second), and JS (dJS, third) distance measures of Garba et al 2018
//' 
//' IMPORTANT: This function assumes that the two gene tree probability vectors (vectorM1 and vectorM2) are sorted by to be the same. That is, the first gene tree in vectorM1 is also the first gene tree in vectorM2, and so on...
//' The gene tree probability output files of HYBRID-COAL are sorted automatically.
//' 
//' @param vectorM1 Vector of the gene tree probabilities for M1
//' 
//' @param vectorM2 Vector of the gene tree probabilities for M2
//' 
//' @examples 
//' 
//' vector.GeneTreeProbs.M1 <- Read.GeneTree.Probabilities(string.PathToGeneTreeProbs = '~/Desktop/EXAMPLE/SmallTree/M1.prob')
//' vector.GeneTreeProbs.M2 <- Read.GeneTree.Probabilities(string.PathToGeneTreeProbs = '~/Desktop/EXAMPLE/SmallTree/M2.prob')
//' 
//' ComputePSTD(vectorM1 = vector.GeneTreeProbs.M1, vectorM2 = vector.GeneTreeProbs.M2)
//' 
//' @details
//' \code{ComputeDistancesPSTD} Computes the probabilistic distances between two probability distributions
//' @export
// [[Rcpp::export]]
NumericVector ComputePSTD(NumericVector vectorM1, NumericVector vectorM2) {

		// init starting values
		int numGenes = vectorM1.size();
		double dHSummation = 0;
		double dKLSummation = 0;
		double dJS1Summation = 0;
		double dJS2Summation = 0;
		double Hdistance = 0;
		double KLdistance = 0;
		double JSdistance = 0;
		double probM1 = 0;
		double probM2 = 0;
		double probCombined = 0;
		Rcpp::NumericVector PSTD(3);

		
		// loop through gene trees
		for (int j = 0; j<numGenes;++j){
		
				// extract probability 
				probM1 = vectorM1[j];
				probM2 = vectorM2[j];
				probCombined = (probM1 + probM2)*0.50;

				// append to summations for dH and dKL
				dHSummation += pow(sqrt(probM1) - sqrt(probM2), 2);
				dKLSummation += (probM1 * log(probM1/probM2));
				
				// append to summations for dJS
				dJS1Summation += (probM1*log(probM1/probCombined));
				dJS2Summation += (probM2*log(probM2/probCombined));

			}
			
		Hdistance = 0.50*dHSummation;
		KLdistance = dKLSummation;
		JSdistance = (0.50*dJS1Summation) + (0.50*dJS2Summation);
		
		PSTD[0] = Hdistance;
		PSTD[1] = KLdistance;
		PSTD[2] = JSdistance;
		
		return PSTD;
		}