
// Enable C++11
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
#include <numeric>
#include <cstdlib>
//#include <iostream>

using namespace std;
using namespace Eigen;
using namespace Rcpp;

class Gamma
{
public:
    Gamma(MatrixXd x, VectorXd y, int n, int p, float h1, float h2, float c, float k, bool w);
    //Initializes gamma as a vector of all 0 (no included variables)
    ~Gamma(); //Destructor

    void updateAdd();
    void updateSub();
    //Updates gamma, matrix products and FAv given the flipped j

    void calculateH();

    void calculateFullCond();

    void calculateFlipRates();

    void updatePIP(double z);
    double calculateWeight();

    int whichFlip();

	int p, n;
	float h1, h2, c, k;
	bool w;
	//Outside information

	int *gamma;
//private:
    vector<int> included;
    int gamma_size;
    int j;
    //Information about which regressors are included

    MatrixXd F, A;
    RowVectorXd v;
    //FAv matrices

    double vFv, *aFa, *vFa;
    //Matrix products

    double *H;

    double *full_cond;

    double *flip_rates;

    double *PIP;
    double sumweight;

    double *diagXtX;
    vector<int> includedXtX;
    MatrixXd XtX;
    MatrixXd X;

    RowVectorXd ytX;
    double yty;
    //Needed global matrices (only computed once at the start)
};

//' auxiliary function that does weighted Tempered Gibbs Sampling from the posterior distribution of Gamma
//'
//' Performs weighted Tempered Gibbs Sampling 
//'
//' @param X_ a matrix of p regressors (independent variables) with dimension (nxp).
//' @param y_ a vector of n observations (dependent variable) with dimensions (nx1).
//' @param n_ number of observations n.
//' @param p_ number of regressors p.
//' @param n_iter a positive integer specifying the number of iterations for the Markov Chain.
//' @param burnin_ a positive integer specifying the number of burn in iterations for the Markov Chain.
//' @param h1_ first parameter of the beta distribution defining h. If value of h fixed, make it h1.
//' @param h2_ second parameter of the beta distribution defining h. If value of h fixed, make h2=0.
//' @param c_ constant of proportionality to prior covariace matrix.
//' @param k_ k_weighted parameter for weighted sampling.
//' @param weighted boolean indicating if weighted sampling.
//'
//' @return List with values for PIP and required elements to reproduce the samples. 
// [[Rcpp::export]]
List wTGS(SEXP X_, SEXP y_, SEXP n_, SEXP p_, SEXP n_iter, SEXP burnin_, SEXP h1_, SEXP h2_, SEXP c_, SEXP k_, SEXP weighted)
{
	typedef Map<MatrixXd> matrix;
	typedef Map<VectorXd> vector;
	
	const matrix X(as<matrix>(X_));
	const vector y(as<vector>(y_));
	const int n(as<int>(n_));
	const int p(as<int>(p_));
	const int burnin(as<int>(burnin_));
	const int n_it(as<int>(n_iter));
	const float h1(as<float>(h1_));
	const float h2(as<float>(h2_));
	const float c(as<float>(c_));
	const float k(as<float>(k_));
	const bool w(as<bool>(weighted));
	
	NumericVector PIP(p);
	NumericVector start(p);
	NumericVector sample_weights(n_it);
	NumericVector indices_sequence(n_it);
	
	Gamma main(X, y, n, p, h1, h2, c, k, w);
	
	//cout << "Start wTGS..." << endl;
	for (int i = 1; i <= n_it+burnin; i++)
	{
		//if (i < burnin)
			//cout << "\r\t(Burn in) Iteration " << i;
		//else 
		if (i == burnin)
		{
			//cout << "\r\t(Burn in) Iteration " << i << endl;
			for (int j = 0; j < p; j++)
				start[j] = main.gamma[j];
		}
		//else
			//cout << "\r\tIteration " << i-burnin;
		
		if (main.gamma[main.j] == 0)
			main.updateAdd();
		else
			main.updateSub();
		main.calculateH();
		main.calculateFullCond();
		main.calculateFlipRates();
		double z = main.calculateWeight();
		main.updatePIP(z);
		
		if (i > burnin)
		{
			sample_weights[i-burnin-1] = z;
			indices_sequence[i-burnin-1] = main.j+1;
		}
	}
	
	//cout << endl << "DONE" << endl;
	for (int i = 0; i < p; i++)
		PIP[i] = main.PIP[i]/main.sumweight;
	
	List output(4);
	output[0] = PIP;
	output[1] = start;
	output[2] = sample_weights/main.sumweight;
	output[3] = indices_sequence;
	return output;
}

//uses cmath
Gamma::Gamma(MatrixXd x, VectorXd y, int n, int p, float h1, float h2, float c, float k, bool w):
    p(p), n(n), h1(h1), h2(h2), c(c), k(k), w(w),//outside info
    gamma_size(0), //integer
    F(1,1), A(1,p), v(1), XtX(1,p), X(n,p), ytX(p) //matrices initialization
{
    X = x;
    gamma = new int[p];
    for (int i = 0; i < p; i++)
        gamma[i] = 0;
    //main vector (included stays empty)

    diagXtX = new double[p];
    for (int i = 0; i < p; i++)
        diagXtX[i] = X.col(i).transpose()*X.col(i);
    //XtX = X.transpose()*X;
    ytX = y.transpose()*X;
    yty = y.transpose()*y;

    //Matrices filled
    F(0,0) = 0;
    for (int i = 0; i < p; i++)
        A(0,i) = 0;
    v(0) = 0;

    //Matrix prod filled
    vFv = 0;
    aFa = new double[p];
    vFa = new double[p];
    for (int i = 0; i < p; i++)
    {
        aFa[i] = 0;
        vFa[i] = 0;
    }

    H = new double[2];
    if (h2 > 1e-6)
		H[0] = h1/(h2+p-1);
    else
        H[0] = H[1] = h1/(1-h1);

    full_cond = new double[p];
    flip_rates = new double[p];
    for (int i = 0; i < p; i++)
    {
        double S1 = yty - c/(1+c)*pow(ytX(i), 2)/diagXtX[i];
        double frac = H[0]*1/sqrt(1+c)*pow(yty/S1, n/2);
        full_cond[i] = 1/(1+frac);
        if (w)
            flip_rates[i] = (1-full_cond[i]+k/p)/full_cond[i];
        else
            flip_rates[i] = 1/full_cond[i];
    }

    PIP = new double[p];
    double z = calculateWeight();
    for (int i = 0; i < p; i++)
        PIP[i] = z*full_cond[i];

    includedXtX.push_back(j);
    XtX.row(0) = X.col(j).transpose()*X;
}

Gamma::~Gamma()
{
    delete[] gamma;
	delete[] aFa;
    delete[] vFa;
    delete[] H;
    delete[] full_cond;
    delete[] flip_rates;
    delete[] PIP;
    delete[] diagXtX;
}

//uses cmath
void Gamma::updateAdd()
{
    gamma[j] = 1;
    included.push_back(j);
    gamma_size++;
    //Gamma update

    unsigned int j_X = 0;
    while (includedXtX[j_X] != j)
    {
        j_X++;
        if (j_X == includedXtX.size())
        {
            XtX.conservativeResize(j_X+1, NoChange);
            XtX.row(j_X) = X.col(j).transpose()*X;
            includedXtX.push_back(j);
            break;
        }
    }
    //XtX update

    VectorXd F_hat = F*A.col(j);
    double d = diagXtX[j] - A.col(j).transpose()*F_hat;
    double vF = v*F_hat;

    vFv += pow(v*F_hat - ytX(j), 2)/d;

    for (int i = 0; i < p; i++)
    {
        double aF = A.col(i).transpose()*F_hat;
        aFa[i] += pow(aF - XtX(j_X,i), 2)/d;
        vFa[i] += (vF - ytX(j))*(aF - XtX(j_X,i))/d;
    }
    //Matrix products update

    if (gamma_size == 1)
    {
        F(0,0) = 1/d;
        A.row(0) = XtX.row(j_X);
        v(0) = ytX(j);
    }
    else
    {
        // VectorXd dF_hat = F_hat/d;
        // F += dF_hat*F_hat.transpose();
        // F.conservativeResize(NoChange, gamma_size);
        // F.col(gamma_size-1) = -dF_hat;
        // F.conservativeResize(gamma_size, NoChange);
        // F.row(gamma_size-1) << -dF_hat.transpose(), 1/d;

        A.conservativeResize(gamma_size, NoChange);
        A.row(gamma_size-1) = XtX.row(j_X);
		
		if (gamma_size < n)
		{
			VectorXd dF_hat = F_hat/d;
			F += dF_hat*F_hat.transpose();
			F.conservativeResize(NoChange, gamma_size);
			F.col(gamma_size-1) = -dF_hat;
			F.conservativeResize(gamma_size, NoChange);
			F.row(gamma_size-1) << -dF_hat.transpose(), 1/d;
		}
		else
		{
			F.resize(gamma_size, gamma_size);
			for (int i = 0; i < gamma_size; i++)
				F.col(i) = A.col(included[i]);
			// F = F.completeOrthogonalDecomposition().pseudoInverse();
			JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
			F = svd.matrixV()*MatrixXd( (svd.singularValues().diagonal().array().abs() > 1e-6*max(F.cols(), F.rows())*svd.singularValues().array().abs().maxCoeff()).select(svd.singularValues().diagonal().array().inverse(), 0) ).asDiagonal()*svd.matrixU().adjoint();
		}

        v.conservativeResize(NoChange, gamma_size);
        v(gamma_size-1) = ytX(j);
    }
    //FAv update
}

//uses cmath
void Gamma::updateSub()
{
    gamma[j] = 0;
    int j_new = 0;
    while (included[j_new] != j)
        j_new++;
    included.erase(included.begin()+j_new);
    gamma_size--;
    //Gamma update

    VectorXd F_tilde = F.col(j_new);
    double f = F(j_new, j_new);
    double vF = v*F_tilde;

    vFv -= pow(vF, 2)/f;

    for (int i = 0; i < p; i++)
    {
        double aF = A.col(i).transpose()*F_tilde;
        aFa[i] -= pow(aF, 2)/f;
        vFa[i] -= (vF*aF)/f;
    }
    //Matrix products update

    if (gamma_size == 0)
    {
        F(0,0) = 0;
        for (int i = 0; i < p; i++)
            A(0,i) = 0;
        v(0) = 0;
    }
    else
    {
        // F.block(j_new, 0, gamma_size-j_new, gamma_size+1) = F.block(j_new+1, 0, gamma_size-j_new, gamma_size+1);
        // F.conservativeResize(gamma_size, NoChange);
        // VectorXd F_col = F.col(j_new);
        // F.block(0, j_new, gamma_size, gamma_size-j_new) = F.block(0, j_new+1, gamma_size, gamma_size-j_new);
        // F.conservativeResize(NoChange, gamma_size);
        // F -= F_col*F_col.transpose()/f;

        A.block(j_new, 0, gamma_size-j_new, p) = A.block(j_new+1, 0, gamma_size-j_new, p);
        A.conservativeResize(gamma_size, p);
		
		if (gamma_size < n)
		{
			F.block(j_new, 0, gamma_size-j_new, gamma_size+1) = F.block(j_new+1, 0, gamma_size-j_new, gamma_size+1);
			F.conservativeResize(gamma_size, NoChange);
			VectorXd F_col = F.col(j_new);
			F.block(0, j_new, gamma_size, gamma_size-j_new) = F.block(0, j_new+1, gamma_size, gamma_size-j_new);
			F.conservativeResize(NoChange, gamma_size);
			F -= F_col*F_col.transpose()/f;
		}
		else
		{
			F.resize(gamma_size, gamma_size);
			for (int i = 0; i < gamma_size; i++)
				F.col(i) = A.col(included[i]);
			// F = F.completeOrthogonalDecomposition().pseudoInverse();
			JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
			F = svd.matrixV()*MatrixXd( (svd.singularValues().diagonal().array().abs() > 1e-6*max(F.cols(), F.rows())*svd.singularValues().array().abs().maxCoeff()).select(svd.singularValues().diagonal().array().inverse(), 0) ).asDiagonal()*svd.matrixU().adjoint();
		}

        v.segment(j_new, gamma_size-j_new) = v.segment(j_new+1, gamma_size-j_new);
        v.conservativeResize(NoChange, gamma_size);
    }
    //FAv update
}

void Gamma::calculateH()
{
    if (h2 > 1e-6)
    {
        H[0] = (h1+gamma_size)/(h2+p-gamma_size-1);
        H[1] = (h1+gamma_size-1)/(h1+p-gamma_size);
    }
}

void Gamma::calculateFullCond()
{
    double S = yty - c/(1+c)*vFv;

    for (int i = 0; i < p; i++)
    {
        if (gamma[i] == 0)
        {
            double S1 = S - c/(1+c)*pow(vFa[i] - ytX(i), 2)/(diagXtX[i] - aFa[i]);
            double frac = H[0]*1/sqrt(1+c)*pow(S/S1, n/2);
            full_cond[i] = 1/(1+frac);
        }
        else if (gamma_size == 1)
        {
            double frac = 1/H[1]*sqrt(1+c)*pow(S/yty, n/2);
            full_cond[i] = 1/(1+frac);
        }
        else
        {
            int j_new = 0;
            while (included[j_new] != i)
                j_new++;
            double S0 = S + c/(1+c)*pow(v*F.col(j_new), 2)/F(j_new, j_new);
            double frac = 1/H[1]*sqrt(1+c)*pow(S/S0, n/2);
            full_cond[i] = 1/(1+frac);
        }
    }
}

void Gamma::calculateFlipRates()
{
    for (int i = 0; i < p; i++)
    {
        if (w)
        {
            if (gamma[i] == 0)
                flip_rates[i] = (1-full_cond[i]+k/p)/full_cond[i];
            else
                flip_rates[i] = 1 + (k/p)/full_cond[i];
        }
        else
            flip_rates[i] = 1/full_cond[i];
    }
}

//uses numeric and cstdlib
double Gamma::calculateWeight()
{
    double z;
    double sum_flip = accumulate(flip_rates, flip_rates+p, 0.0);
    z = p/(2*sum_flip);
    sumweight += z;

    double runif = ((RAND_MAX - rand())/static_cast<double>(RAND_MAX))*sum_flip-1e-6;
    j = 0;
    while (runif > 0)
    {
        runif -= flip_rates[j];
        if (runif > 0)
			j++;
    }
    //which variable to flip

    return z;
}

void Gamma::updatePIP(double z)
{
    for (int i = 0; i < p; i++)
    {
        PIP[i] += z*(gamma[i]*full_cond[i] + (1-gamma[i])*(1-full_cond[i]));
    }
}