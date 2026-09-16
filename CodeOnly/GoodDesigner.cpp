#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cctype>
#include <sstream>
#include <vector>
#include <cmath>
#include <complex>
#include <set>
#include <chrono>
#include <mpi.h>
#include <algorithm>
#include "Eigen/Dense"

#define PI 3.14159265359

typedef struct {

	std::string name;
	long double min, ave, max, stddv;

}DESCRIPTOR;

std::vector<DESCRIPTOR> descriptor;

typedef struct {

	std::string atomname;
	int iatom;

}NEIGHBOR;

typedef struct {

	std::string atomname;
	std::vector<long double> v;
	std::vector<long double> zv;
	std::vector<NEIGHBOR> neighbor;
	bool dopeonly;

}ATOMDAT;

std::vector<ATOMDAT> atomdat;

typedef struct {

	std::string name;
	std::vector<std::string> sitename;
	std::vector<long double> v;
	std::vector<std::string> strv;

}SDESCRIPTOR;

typedef struct {

	std::string strname;
	std::vector<std::string> sitename;
	std::vector<SDESCRIPTOR> sdescriptor;

}STRDAT;

std::vector<STRDAT> strdat;

typedef struct {

	std::string name, ASB, postfix;
	int ides;
	int** idesB;
	int UType_numerator_S0A1E2N3;
	int UType_denominator_S0A1E2N3;
	std::vector<std::string> UType_numerator_AtomKinds;
	std::vector<std::string> UType_denominator_AtomKinds;
	long double UType_numerator_value;
	long double UType_denominator_value;
	int UType_m0d1;
	long double min;
	long double min_ave, ave, max_ave;
	long double max;
	long double min_stdev, min_skew, min_kurto;
	long double ave_stdev, ave_skew, ave_kurto;
	long double max_stdev, max_skew, max_kurto;
	int count;
	std::string latex;
	std::string mathematica;

}ANADESCRIPTOR;

std::vector<ANADESCRIPTOR> anadescriptor;

typedef struct {
	
	int whichreg;
	std::string whichregfilename;
	int whichfuncchain;
	int whichfunc;
	long double zcoeff;

}WHERETOFIND;

typedef struct {

	std::vector<std::string> name;
	int count;
	long double min_absz, ave_absz, max_absz;
	std::vector<WHERETOFIND> wheretofind;
	
}TERMANA;

typedef struct {

	std::vector<std::string> name;
	int metacount, count;
	long double min_absz, ave_absz, max_absz;

}TERMANAFORBAG;

typedef struct {

	std::vector<int> it_chain_list;

}CAVAILLES;

typedef struct {

	std::string AtomName;
	bool face;
	std::vector<long double> vFeature;

}FACE;

std::vector<FACE> face;

typedef struct {

	std::string name;
	std::vector<long double> coeff_0;
	std::vector<long double> power_0;
	std::vector<long double> coeff_1;
	std::vector<long double> power_1;
	std::string nonlinear_funcname;
	bool sum_only;
	std::vector<long double> coeff_2;
	std::vector<long double> power_2;
	std::vector<int> ianades_2;
	std::vector<int> imamssk_2;
	std::vector<int> ilogln_2;
	std::vector<int> ipow10exp_2;
	std::vector<int> ierf_2;
	bool added;
	std::vector<long double> coeff_3;
	std::vector<long double> power_3;
	std::vector<int> ianades_3;
	std::vector<int> imamssk_3;
	std::vector<int> ilogln_3;
	std::vector<int> ipow10exp_3;
	std::vector<int> ierf_3;

	std::vector<long double> v_list;
	long double v_ave, v_stddv;

}FUNCTYPE;

std::vector<FUNCTYPE> functype;

typedef struct {

	std::vector<std::string> userdefined_name;
	std::string FName;
	std::vector<long double> coeff;

}ORTHFUNCTYPE;

typedef struct {

	bool identified;
	std::string name;
	int type;
	std::string CationName;
	int i0anades;
	int i0mamssk;
	int i1anades[2];
	int i1mamssk[2];
	int i1m0d1;
	int i2whichlastfunc;
	int i3whichlastfunc[2];
	int i3m0d1;
	std::string latex;
	std::string mathematica;

}FEATURE;

typedef struct {

	bool identified, feature_identified;
	std::string name;
	long double coeff, zcoeff;
	int whichorthfunctype;
	int whichfeature;
	std::string latex;
	std::string mathematica;

}FUNC;

typedef struct {

	std::string filename;
	std::vector<FEATURE> feature;
	std::vector<FUNC> func;
	long double intercept;
	long double score_inverse_slope, score_inverse_intercept;
	long double revise_slope, revise_intercept;
	bool allow_beta;

}FUNCCHAIN;

typedef struct {

	std::string name;
	std::string beta_target;
	int ibeta_target;
	std::vector<ORTHFUNCTYPE> orthfunctype;
	std::vector<FUNCCHAIN> funcchain;
	std::vector<FUNC> postfunc;
	std::vector<long double> v_list, regressed_v_list;
	std::vector<long double> onev_acc_list;
	std::vector<TERMANA> termana;

}REGTYPE;

std::vector<REGTYPE> regtype;

typedef struct {

	std::vector<long double> coeff;
	long double coeff_sum;

}STACKINGENSEMBLE;

typedef struct {

	std::string name;
	std::vector<int> which_regtype;
	std::vector<long double> v_list, regressed_v_list, regressed_stddv_v_list, regressed_post_v_list;
	long double v_ave, v_stddv;
	bool user_defined_v_ave;
	bool user_defined_v_stddv;
	std::vector<TERMANAFORBAG> termanaforbag;
	long double obs_v_min, obs_v_max;
	std::vector<STACKINGENSEMBLE> stackingensemble;
	
}BAGTYPE;

std::vector<BAGTYPE> bagtype;

typedef struct {
	
	std::string name, composition;
	std::string strname;
	int istr;
	std::vector<std::string> atomname;
	std::vector<int> iatom;
	std::vector<std::string> sitename;
	std::vector<int> isite;
	std::vector<long double> content;
	long double anion_content;
	long double OM;
	std::vector<long double> CO;

}MAT;

typedef struct {

	std::vector<std::string> atomname;
	std::vector<long double> content;

}SMAT;

typedef struct {

	std::vector<SMAT> smat_min; 
	std::vector<SMAT> smat_max;

}SMATLIST;

std::vector<SMATLIST> smatlist;

typedef struct {

	std::string name;
	std::vector<long double> v;
	std::vector<int> yhit;

}DATASET;

typedef struct {

	Eigen::MatrixXd MatrixX;
	Eigen::MatrixXd MatrixY;
	Eigen::MatrixXd MatrixYexpected;
	Eigen::MatrixXd MatrixXTXi;
	Eigen::MatrixXd MatrixBeta;
	Eigen::MatrixXd MatrixBetaStandardized;
	std::vector<long double> fit_coeff, fit_coeff_standardized, fit_coeff_SE, fit_coeff_t, fit_coeff_p, fit_coeff_CI_start, fit_coeff_CI_end;
	long double SStot, SSres;
	long double MSE, RMSE, MAE;
	long double det_r2, det_radj2;
	int df;
	long double Fvalue, p_Fvalue;
	bool full_Frequentist;
	std::vector<int> fit_coeff_porder, fit_coeff_porder_on_arr;

}REGELE;

typedef struct {

	std::vector<long double> v;

}BACKUPV;

typedef struct {

	std::string anionname;
	long double anionvalence;

}ANIONS;

std::vector<ANIONS> Anions;

int pr_size, ierr, pr_id;
signed long long int* pr_job_size;
signed long long int** pr_job;
long double* pr_job_size_ld;
long double** pr_job_ld;
signed long long int jobsize;
long double jobsize_ld;
bool xpr_ld_get;
bool finite_job;
signed long long int* SigmaMinusXk;
long double* SigmaMinusXk_ld;
int Nfeatures;
std::vector<int> arrange;
std::vector<int> occup;
std::vector<signed long long int> jobsize_step;
std::vector<long double> jobsize_step_ld;

std::string AtomDatListFilename;
std::string FaceFilename;
std::string StructureListFilename;
int UseStrName;
std::string MaterialsFilename;
int HeadExists;
std::vector<std::string> PostTagNames;
std::vector<std::string> RegModelFilename;
std::string MapOutSpaceFilename;
std::vector<long double> TargetVReg;
std::vector<long double> TargetStdDvReg;
std::vector<long double> ZDistWeight;
std::vector<long double> UseAve;
std::vector<bool> bUseAve;
std::vector<long double> UseStdDv;
std::vector<bool> bUseStdDv;
std::vector<long double> ModulateContentFactor;
std::vector<long double> DopeFactor;
int UseOffset;
int MaxIteration;
std::vector<std::string> NeighborSelectSpace;
long double NeighborZCrit;
std::vector<long double> NeighborCrit;
std::vector<std::string> ExcludeNeighbors;
std::vector<std::string> DoNotFindNeighbors;
std::vector<bool> DopeOnly;
std::string AnionName;
long double modulation_limit = 0.1;
long double OM_min, OM_ave, OM_max, OM_count;
std::vector<std::string> CationNames;
std::vector<long double> CO_min, CO_ave, CO_max, CO_count;
long double AnionValence;
long double VInfinite;
bool verbose;
std::string OutputFilename;
std::string AnalysisOutputFilename;
int AnalysisMesh = 100;
std::vector<std::string> DoNotAnalyze;
int NVariableLimit[2];
int NCavailles;
int SkipMinMaxKurto;

bool FoldMode;
int NFolds, TakeFold, Train0Test1, Straight0Rotate1;
long double TestModelStDevLimit = -1.0;

bool mapout = false;
bool FurtherMaterialsDesign = true;
int mapout_preNRandom, mapout_NRandom, mapout_Verbose;
std::vector<std::string> sdescriptorname;
std::vector<std::string> sitename;
std::vector<std::string> FaceFeatureName;
signed long long int Nmat;
long double** xbox;
long double xbox_OM;
long double* xbox_CO;
long double* xbox_user;
std::vector<long double> zTargetVReg;
int CValence = -1;
int AValence = -1;
int tFTestMesh = 1000;
bool get_floating_v = false;
std::vector<long double> floating_v, floating_v_max, floating_v_min;

std::string sbuf;
long double ldbuf;
std::ifstream readpara, readatomdat, readstrdat, readmatdat, readreg, readregchain, readbpregressor, readfuncchain, readmapout, readface, readanalysis, readanalysisR, readnb, readfold, readoutputpr;
std::ofstream writeoutput, writeoutputpr, writeanalysis, writeanalysis_second, writeonevtest, writezc, writezcext, writelatex, writemathematica;

bool chk_eql(long double v, long double crit) {

	if (fabsl(v - crit) / fabsl(crit) < 0.001) {
		return true;
	}
	else {
		return false;
	}

}

std::string get_addstring(long double* addv) {

	std::string addstr;
	if (chk_eql(*addv, -0.05)) {
		addstr = "-\\frac{1}{20}";
	}
	else if (chk_eql(*addv, -0.1)) {
		addstr = "-\\frac{1}{10}";
	}
	else if (chk_eql(*addv, -0.5)) {
		addstr = "-\\frac{1}{2}";
	}
	else if (chk_eql(*addv, -0.875)) {
		addstr = "-\\frac{7}{8}";
	}
	else if (chk_eql(*addv, -1.0)) {
		addstr = "-1";
	}
	else if (chk_eql(*addv, -5.0)) {
		addstr = "-5";
	}
	else if (chk_eql(*addv, -10.0)) {
		addstr = "-10";
	}
	else if (chk_eql(*addv, -50.0)) {
		addstr = "-50";
	}
	else if (chk_eql(*addv, -100.0)) {
		addstr = "-100";
	}
	else if (chk_eql(*addv, -1000.0)) {
		addstr = "-1000";
	}
	else if (chk_eql(*addv, -5000.0)) {
		addstr = "-5000";
	}
	else if (chk_eql(*addv, -10000.0)) {
		addstr = "-10000";
	}
	else if (chk_eql(*addv, 0.05)) {
		addstr = "+\\frac{1}{20}";
	}
	else if (chk_eql(*addv, 0.1)) {
		addstr = "+\\frac{1}{10}";
	}
	else if (chk_eql(*addv, 0.5)) {
		addstr = "+\\frac{1}{2}";
	}
	else if (chk_eql(*addv, 0.875)) {
		addstr = "+\\frac{7}{8}";
	}
	else if (chk_eql(*addv, 1.0)) {
		addstr = "+1";
	}
	else if (chk_eql(*addv, 5.0)) {
		addstr = "+1";
	}
	else if (chk_eql(*addv, 10.0)) {
		addstr = "+10";
	}
	else if (chk_eql(*addv, 50.0)) {
		addstr = "+50";
	}
	else if (chk_eql(*addv, 100.0)) {
		addstr = "+100";
	}
	else if (chk_eql(*addv, 1000.0)) {
		addstr = "+1000";
	}
	else if (chk_eql(*addv, 5000.0)) {
		addstr = "+5000";
	}
	else if (chk_eql(*addv, 10000.0)) {
		addstr = "+10000";
	}
	else {
		if (*addv < 0) {
			addstr = "-" + std::to_string(fabsl(*addv));
		}
		else {
			addstr = "+" + std::to_string(fabsl(*addv));
		}

	}

	return addstr;

}

std::string get_addstringM(long double* addv) {

	std::string addstr;
	if (chk_eql(*addv, -0.05)) {
		addstr = "-1/20";
	}
	else if (chk_eql(*addv, -0.1)) {
		addstr = "-1/10";
	}
	else if (chk_eql(*addv, -0.5)) {
		addstr = "-1/2";
	}
	else if (chk_eql(*addv, -0.875)) {
		addstr = "-7/8";
	}
	else if (chk_eql(*addv, -1.0)) {
		addstr = "-1";
	}
	else if (chk_eql(*addv, -5.0)) {
		addstr = "-5";
	}
	else if (chk_eql(*addv, -10.0)) {
		addstr = "-10";
	}
	else if (chk_eql(*addv, -50.0)) {
		addstr = "-50";
	}
	else if (chk_eql(*addv, -100.0)) {
		addstr = "-100";
	}
	else if (chk_eql(*addv, -1000.0)) {
		addstr = "-1000";
	}
	else if (chk_eql(*addv, -5000.0)) {
		addstr = "-5000";
	}
	else if (chk_eql(*addv, -10000.0)) {
		addstr = "-10000";
	}
	else if (chk_eql(*addv, 0.05)) {
		addstr = "+1/20";
	}
	else if (chk_eql(*addv, 0.1)) {
		addstr = "+1/10";
	}
	else if (chk_eql(*addv, 0.5)) {
		addstr = "+1/2";
	}
	else if (chk_eql(*addv, 0.875)) {
		addstr = "+7/8";
	}
	else if (chk_eql(*addv, 1.0)) {
		addstr = "+1";
	}
	else if (chk_eql(*addv, 5.0)) {
		addstr = "+1";
	}
	else if (chk_eql(*addv, 10.0)) {
		addstr = "+10";
	}
	else if (chk_eql(*addv, 50.0)) {
		addstr = "+50";
	}
	else if (chk_eql(*addv, 100.0)) {
		addstr = "+100";
	}
	else if (chk_eql(*addv, 1000.0)) {
		addstr = "+1000";
	}
	else if (chk_eql(*addv, 5000.0)) {
		addstr = "+5000";
	}
	else if (chk_eql(*addv, 10000.0)) {
		addstr = "+10000";
	}
	else {
		if (*addv < 0) {
			addstr = "-" + std::to_string(fabsl(*addv));
		}
		else {
			addstr = "+" + std::to_string(fabsl(*addv));
		}

	}

	return addstr;

}

std::string get_powerstring(long double* powv) {

	std::string powvstr;
	if (chk_eql(*powv, -0.5)) {
		powvstr = "-\\frac{1}{2}";
	}
	else if (chk_eql(*powv, 0.5)) {
		powvstr = "\\frac{1}{2}";
	}
	else if (chk_eql(*powv, -0.333333)) {
		powvstr = "-\\frac{1}{3}";
	}
	else if (chk_eql(*powv, 0.333333)) {
		powvstr = "\\frac{1}{3}";
	}
	else if (chk_eql(*powv, -1.0)) {
		powvstr = "-1";
	}
	else if (chk_eql(*powv, -2.0)) {
		powvstr = "-2";
	}
	else if (chk_eql(*powv, 2.0)) {
		powvstr = "2";
	}
	else if (chk_eql(*powv, -3.0)) {
		powvstr = "-3";
	}
	else if (chk_eql(*powv, 3.0)) {
		powvstr = "3";
	}
	else if (chk_eql(*powv, -4.0)) {
		powvstr = "-4";
	}
	else if (chk_eql(*powv, 4.0)) {
		powvstr = "4";
	}
	else if (chk_eql(*powv, -5.0)) {
		powvstr = "-5";
	}
	else if (chk_eql(*powv, 5.0)) {
		powvstr = "5";
	}
	else {
		powvstr = std::to_string(*powv);
	}

	return powvstr;

}

std::string get_powerstringM(long double* powv) {

	std::string powvstr;
	if (chk_eql(*powv, -0.5)) {
		powvstr = "-1/2";
	}
	else if (chk_eql(*powv, 0.5)) {
		powvstr = "1/2";
	}
	else if (chk_eql(*powv, -0.333333)) {
		powvstr = "-1/3";
	}
	else if (chk_eql(*powv, 0.333333)) {
		powvstr = "1/3";
	}
	else if (chk_eql(*powv, -1.0)) {
		powvstr = "-1";
	}
	else if (chk_eql(*powv, -2.0)) {
		powvstr = "-2";
	}
	else if (chk_eql(*powv, 2.0)) {
		powvstr = "2";
	}
	else if (chk_eql(*powv, -3.0)) {
		powvstr = "-3";
	}
	else if (chk_eql(*powv, 3.0)) {
		powvstr = "3";
	}
	else if (chk_eql(*powv, -4.0)) {
		powvstr = "-4";
	}
	else if (chk_eql(*powv, 4.0)) {
		powvstr = "4";
	}
	else if (chk_eql(*powv, -5.0)) {
		powvstr = "-5";
	}
	else if (chk_eql(*powv, 5.0)) {
		powvstr = "5";
	}
	else {
		powvstr = std::to_string(*powv);
	}

	return powvstr;

}

std::string get_mulstring(long double* mulv) {

	std::string mulstr;
	if (chk_eql(*mulv, 1.57079632679)) {
		mulstr = "\\frac{\\pi}{2}";
	}
	else if (chk_eql(*mulv, 3.14159265359)) {
		mulstr = "\\pi";
	}
	else if (chk_eql(*mulv, 0.00157079632679)) {
		mulstr = "\\frac{\\pi}{2000}";
	}
	else if (chk_eql(*mulv, 0.00314159265359)) {
		mulstr = "\\frac{\\pi}{1000}";
	}
	else if (chk_eql(*mulv, 0.0157079632679)) {
		mulstr = "\\frac{\\pi}{200}";
	}
	else if (chk_eql(*mulv, 0.0314159265359)) {
		mulstr = "\\frac{\\pi}{100}";
	}
	else if (chk_eql(*mulv, 0.157079632679)) {
		mulstr = "\\frac{\\pi}{20}";
	}
	else if (chk_eql(*mulv, 0.314159265359)) {
		mulstr = "\\frac{\\pi}{10}";
	}
	else if (chk_eql(*mulv, 15.7079632679)) {
		mulstr = "5{\\pi}";
	}
	else if (chk_eql(*mulv, 31.4159265359)) {
		mulstr = "10{\\pi}";
	}
	else if (chk_eql(*mulv, 157.079632679)) {
		mulstr = "50{\\pi}";
	}
	else if (chk_eql(*mulv, 314.159265359)) {
		mulstr = "100{\\pi}";
	}
	else if (chk_eql(*mulv, 1570.79632679)) {
		mulstr = "500{\\pi}";
	}
	else if (chk_eql(*mulv, 3141.59265359)) {
		mulstr = "1000{\\pi}";
	}
	else if (chk_eql(*mulv, 0.001)) {
		mulstr = "\\frac{1}{1000}";
	}
	else if (chk_eql(*mulv, 0.01)) {
		mulstr = "\\frac{1}{100}";
	}
	else if (chk_eql(*mulv, 0.1)) {
		mulstr = "\\frac{1}{10}";
	}
	else if (chk_eql(*mulv, 10.0)) {
		mulstr = "10";
	}
	else if (chk_eql(*mulv, 100.0)) {
		mulstr = "100";
	}
	else if (chk_eql(*mulv, 1000.0)) {
		mulstr = "1000";
	}
	else if (chk_eql(*mulv, -0.001)) {
		mulstr = "-\\frac{1}{1000}";
	}
	else if (chk_eql(*mulv, -0.01)) {
		mulstr = "-\\frac{1}{100}";
	}
	else if (chk_eql(*mulv, -0.1)) {
		mulstr = "-\\frac{1}{10}";
	}
	else if (chk_eql(*mulv, -10.0)) {
		mulstr = "-10";
	}
	else if (chk_eql(*mulv, -100.0)) {
		mulstr = "-100";
	}
	else if (chk_eql(*mulv, -1000.0)) {
		mulstr = "-1000";
	}
	else {
		mulstr = std::to_string(*mulv);
	}

	return mulstr;

}

std::string get_mulstringM(long double* mulv) {

	std::string mulstr;
	if (chk_eql(*mulv, 1.57079632679)) {
		mulstr = "Pi/2";
	}
	else if (chk_eql(*mulv, 3.14159265359)) {
		mulstr = "Pi";
	}
	else if (chk_eql(*mulv, 0.00157079632679)) {
		mulstr = "Pi/2000";
	}
	else if (chk_eql(*mulv, 0.00314159265359)) {
		mulstr = "Pi/1000";
	}
	else if (chk_eql(*mulv, 0.0157079632679)) {
		mulstr = "Pi/200";
	}
	else if (chk_eql(*mulv, 0.0314159265359)) {
		mulstr = "Pi/100";
	}
	else if (chk_eql(*mulv, 0.157079632679)) {
		mulstr = "Pi/20";
	}
	else if (chk_eql(*mulv, 0.314159265359)) {
		mulstr = "Pi/10";
	}
	else if (chk_eql(*mulv, 15.7079632679)) {
		mulstr = "5Pi";
	}
	else if (chk_eql(*mulv, 31.4159265359)) {
		mulstr = "10Pi";
	}
	else if (chk_eql(*mulv, 157.079632679)) {
		mulstr = "50Pi";
	}
	else if (chk_eql(*mulv, 314.159265359)) {
		mulstr = "100Pi";
	}
	else if (chk_eql(*mulv, 1570.79632679)) {
		mulstr = "500Pi";
	}
	else if (chk_eql(*mulv, 3141.59265359)) {
		mulstr = "1000Pi";
	}
	else if (chk_eql(*mulv, 0.001)) {
		mulstr = "1/1000";
	}
	else if (chk_eql(*mulv, 0.01)) {
		mulstr = "1/100";
	}
	else if (chk_eql(*mulv, 0.1)) {
		mulstr = "1/10";
	}
	else if (chk_eql(*mulv, 10.0)) {
		mulstr = "10";
	}
	else if (chk_eql(*mulv, 100.0)) {
		mulstr = "100";
	}
	else if (chk_eql(*mulv, 1000.0)) {
		mulstr = "1000";
	}
	else if (chk_eql(*mulv, -0.001)) {
		mulstr = "-1/1000";
	}
	else if (chk_eql(*mulv, -0.01)) {
		mulstr = "-1/100";
	}
	else if (chk_eql(*mulv, -0.1)) {
		mulstr = "-1/10";
	}
	else if (chk_eql(*mulv, -10.0)) {
		mulstr = "-10";
	}
	else if (chk_eql(*mulv, -100.0)) {
		mulstr = "-100";
	}
	else if (chk_eql(*mulv, -1000.0)) {
		mulstr = "-1000";
	}
	else {
		mulstr = std::to_string(*mulv);
	}

	return mulstr;

}

std::string to_times10(long double x, int decimals) {

	std::ostringstream oss;
	oss << std::scientific << std::setprecision(decimals) << x;
	std::string s = oss.str();

	auto pos = s.find('e');
	std::string mant = s.substr(0, pos);
	std::string exp = s.substr(pos + 1);

	if (mant.find('.') != std::string::npos) {
		while (!mant.empty() && mant.back() == '0') mant.pop_back();
		if (!mant.empty() && mant.back() == '.') mant.pop_back();
	}

	int e = std::stoi(exp);
	return mant + "*10^" + std::to_string(e);

}

std::string to_times_e(long double x, int decimals) {

	std::ostringstream oss;
	oss << std::scientific << std::setprecision(decimals) << x;
	std::string s = oss.str();

	auto pos = s.find('e');
	std::string mant = s.substr(0, pos);
	std::string exp = s.substr(pos + 1);

	if (mant.find('.') != std::string::npos) {
		while (!mant.empty() && mant.back() == '0') mant.pop_back();
		if (!mant.empty() && mant.back() == '.') mant.pop_back();
	}

	int e = std::stoi(exp);
	if (e >= 0) {
		return mant + "e+" + std::to_string(e);
	}
	else {
		return mant + "e" + std::to_string(e);
	}
	

}

void initialize_regression(REGELE* RegEle, std::vector<DATASET>* dataset, int* Noccup) {

	RegEle->MatrixX = Eigen::MatrixXd::Random(dataset->size(), *Noccup + 1);
	RegEle->MatrixYexpected = Eigen::MatrixXd::Random(dataset->size(), 1);
	RegEle->MatrixXTXi = Eigen::MatrixXd::Random(*Noccup + 1, *Noccup + 1);
	RegEle->MatrixBeta = Eigen::MatrixXd::Random(*Noccup + 1, 1);
	RegEle->MatrixBetaStandardized = Eigen::MatrixXd::Random(*Noccup + 1, 1);

	if (RegEle->fit_coeff.size() > 0) {
		std::vector<long double>().swap(RegEle->fit_coeff);
	}
	if (RegEle->fit_coeff_standardized.size() > 0) {
		std::vector<long double>().swap(RegEle->fit_coeff_standardized);
	}
	if (RegEle->fit_coeff_SE.size() > 0) {
		std::vector<long double>().swap(RegEle->fit_coeff_SE);
	}
	if (RegEle->fit_coeff_t.size() > 0) {
		std::vector<long double>().swap(RegEle->fit_coeff_t);
	}
	if (RegEle->fit_coeff_p.size() > 0) {
		std::vector<long double>().swap(RegEle->fit_coeff_p);
	}
	if (RegEle->fit_coeff_CI_start.size() > 0) {
		std::vector<long double>().swap(RegEle->fit_coeff_CI_start);
	}
	if (RegEle->fit_coeff_CI_end.size() > 0) {
		std::vector<long double>().swap(RegEle->fit_coeff_CI_end);
	}
	if (RegEle->fit_coeff_porder.size() > 0) {
		std::vector<int>().swap(RegEle->fit_coeff_porder);
	}
	if (RegEle->fit_coeff_porder_on_arr.size() > 0) {
		std::vector<int>().swap(RegEle->fit_coeff_porder_on_arr);
	}

}

void set_data(REGELE* RegEle, std::vector<DATASET>* dataset, int* Noccup) {

	for (int i = 0; i < (signed)dataset->size(); i++) {
		RegEle->MatrixX(i, 0) = 1.0;
		for (int j = 0; j < *Noccup; j++) {
			RegEle->MatrixX(i, j + 1) = (*dataset)[i].v[j + 1];
		}
	}
	RegEle->MatrixY = Eigen::MatrixXd::Random(dataset->size(), 1);
	for (int i = 0; i < (signed)dataset->size(); i++) {
		(RegEle->MatrixY)(i, 0) = (*dataset)[i].v[0];
	}

}

void fit_any(REGELE* RegEle, std::vector<DATASET>* dataset, int* Noccup) {

	RegEle->MatrixXTXi = (RegEle->MatrixX.transpose() * RegEle->MatrixX).inverse();
	RegEle->MatrixBeta = RegEle->MatrixXTXi * RegEle->MatrixX.transpose() * RegEle->MatrixY;

	long double y_average = 0.0;
	long double y_stdev = 0.0;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		y_average += RegEle->MatrixY(i, 0);
	}
	y_average /= (long double)((signed)dataset->size());
	for (int i = 0; i < (signed)dataset->size(); i++) {
		y_stdev += powl(RegEle->MatrixY(i, 0) - y_average, 2.0);
	}
	y_stdev /= (long double)((signed)dataset->size());
	y_stdev = sqrtl(y_stdev);

	for (int j = 0; j < *Noccup + 1; j++) {
		long double x_average = 0.0;
		long double x_stdev = 0.0;
		for (int i = 0; i < (signed)dataset->size(); i++) {
			x_average += RegEle->MatrixX(i, j);
		}
		x_average /= (long double)((signed)dataset->size());
		for (int i = 0; i < (signed)dataset->size(); i++) {
			x_stdev += powl(RegEle->MatrixX(i, j) - x_average, 2.0);
		}
		x_stdev /= (long double)((signed)dataset->size());
		x_stdev = sqrtl(x_stdev);
		RegEle->MatrixBetaStandardized(j, 0) = RegEle->MatrixBeta(j, 0) * x_stdev / y_stdev;
	}

	for (int j = 0; j < *Noccup + 1; j++) {
		RegEle->fit_coeff.push_back(RegEle->MatrixBeta(j, 0));
		RegEle->fit_coeff_standardized.push_back(RegEle->MatrixBetaStandardized(j, 0));
	}

}

void measure_SStot_SSres(REGELE* RegEle, std::vector<DATASET>* dataset) {

	RegEle->MatrixYexpected = RegEle->MatrixX * RegEle->MatrixBeta;

	long double ybar = 0.0;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		ybar += RegEle->MatrixY(i, 0);
	}
	ybar /= (long double)(signed)dataset->size();

	RegEle->SStot = 0.0;
	RegEle->SSres = 0.0;
	RegEle->MSE = 0.0;
	RegEle->RMSE = 0.0;
	RegEle->MAE = 0.0;
	if (dataset->size() > 0) {
		for (int i = 0; i < (signed)dataset->size(); i++) {
			RegEle->SStot += powl(RegEle->MatrixY(i, 0) - ybar, 2.0);
			RegEle->SSres += powl(RegEle->MatrixY(i, 0) - RegEle->MatrixYexpected(i, 0), 2.0);
			RegEle->MAE += fabsl(RegEle->MatrixY(i, 0) - RegEle->MatrixYexpected(i, 0));
		}
		RegEle->MSE = RegEle->SSres / (long double)((signed)dataset->size());
		RegEle->RMSE = sqrtl(RegEle->MSE);
		RegEle->MAE /= (long double)((signed)dataset->size());
	}

}

void measuere_r2(REGELE* RegEle, std::vector<DATASET>* dataset) {

	RegEle->det_r2 = 1.0 - RegEle->SSres / RegEle->SStot;
	if (RegEle->df != 0) {
		RegEle->det_radj2 = 1.0 - (1.0 - RegEle->det_r2) * ((long double)((signed)dataset->size() - 1) / (long double)RegEle->df);
	}
	else {
		RegEle->det_radj2 = 0.0;
	}
	if (!(RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0)) {
		RegEle->det_r2 = -1;
		RegEle->det_radj2 = -1;
	}

}

void measure_t(REGELE* RegEle, int* Noccup) {

	for (int i = 0; i < *Noccup + 1; i++) {
		long double SE = sqrtl(RegEle->SSres * RegEle->MatrixXTXi(i, i) / (long double)RegEle->df);
		RegEle->fit_coeff_SE.push_back(SE);
		RegEle->fit_coeff_t.push_back(RegEle->MatrixBeta(i, 0) / SE);
	}

}

long double beta_func_integrand(long double* x, long double* a, long double* b) {

	return powl(*x, *a - 1.0) * powl(1.0 - *x, *b - 1.0);

}

long double cumulative_distribution_function_tCurve(REGELE* RegEle, long double* t) {

	long double I_subscript = (long double)RegEle->df / (powl(fabsl(*t), 2.0) + (long double)RegEle->df);
	long double dtTestMesh = 1.0 / (long double)tFTestMesh;
	long double incomplete_beta = 0.0;
	long double complete_beta = 0.0;
	for (long double x = 0.5 * dtTestMesh; x < 1.0; x += dtTestMesh) {
		long double a = 0.5 * (long double)RegEle->df;
		long double b = 0.5;
		complete_beta += beta_func_integrand(&x, &a, &b) * dtTestMesh;
		if (x < I_subscript) {
			incomplete_beta = complete_beta;
		}
	}

	long double CDF;
	if (*t == 0.0) {
		CDF = 0.5;
	}
	else {
		CDF = 1.0 - 0.5 * incomplete_beta / complete_beta;
	}

	if (*t >= 0.0) {
		return CDF;
	}
	else {
		return 1.0 - CDF;
	}

}

long double get_beta_PDF(long double* x, long double* a, long double* b) {

	long double dtTestMesh = 1.0 / (long double)tFTestMesh;
	long double incomplete_beta = 0.0;
	long double complete_beta = 0.0;
	for (long double xpos = 0.5 * dtTestMesh; xpos < 1.0; xpos += dtTestMesh) {
		complete_beta += beta_func_integrand(&xpos, a, b) * dtTestMesh;
	}

	return beta_func_integrand(x, a, b) / complete_beta;

}

long double probability_distribution_function_tCurve(REGELE* RegEle, long double* t) {

	long double dtTestMesh = 1.0 / (long double)tFTestMesh;
	long double complete_beta = 0.0;
	for (long double x = 0.5 * dtTestMesh; x < 1.0; x += dtTestMesh) {
		long double a = 0.5;
		long double b = 0.5 * (long double)RegEle->df;
		complete_beta += beta_func_integrand(&x, &a, &b) * dtTestMesh;
	}

	return (1.0 / (sqrtl((long double)RegEle->df) * complete_beta)) * powl(1.0 + powl(*t, 2.0) / (long double)RegEle->df, -((long double)RegEle->df + 1.0) / 2.0);

}

long double point_to_CIpoint(REGELE* RegEle, long double* CI_level) {

	long double dtTestMesh = 1.0 / (long double)tFTestMesh;
	std::vector<long double> incomplete_beta;
	long double complete_beta = 0.0;
	for (long double x = 0.5 * dtTestMesh; x < 1.0; x += dtTestMesh) {
		long double a = 0.5 * (long double)RegEle->df;
		long double b = 0.5;
		complete_beta += beta_func_integrand(&x, &a, &b) * dtTestMesh;
		incomplete_beta.push_back(complete_beta);
	}

	long double t;
	for (int i = (signed)incomplete_beta.size() - 1; i >= 0; i--) {
		long double CDF = 1.0 - 0.5 * incomplete_beta[i] / complete_beta;
		if (CDF - 0.5 > 0.5 * *CI_level) {
			long double target_x0 = (long double)(i - 1) * dtTestMesh + 0.5 * dtTestMesh;
			long double CDF_0 = CDF;
			long double target_x1 = (long double)i * dtTestMesh + 0.5 * dtTestMesh;
			long double CDF_1 = 1.0 - 0.5 * incomplete_beta[i + 1] / complete_beta;
			long double target_x = (0.5 * *CI_level + 0.5 - CDF_0) * (target_x1 - target_x0) / (CDF_1 - CDF_0) + target_x0;
			t = sqrtl((long double)RegEle->df / target_x - (long double)RegEle->df);
			break;
		}
	}

	return t;

}

void measure_p(REGELE* RegEle) {

	for (int i = 0; i < (signed)RegEle->fit_coeff_t.size(); i++) {
		long double t = -fabsl(RegEle->fit_coeff_t[i]);
		long double p_value = 2.0 * cumulative_distribution_function_tCurve(RegEle, &t);
		RegEle->fit_coeff_p.push_back(p_value);
	}

}

void measure_confidence_interval(REGELE* RegEle, int* Noccup, long double* CI_level) {

	long double test_t;
	if (RegEle->df > 0) {
		test_t = point_to_CIpoint(RegEle, CI_level);
	}
	else {
		test_t = 100000000.0;
	}

	for (int i = 0; i < *Noccup + 1; i++) {
		long double CI_half_delta = test_t * RegEle->fit_coeff_SE[i];
		RegEle->fit_coeff_CI_start.push_back(RegEle->MatrixBeta(i, 0) - CI_half_delta);
		RegEle->fit_coeff_CI_end.push_back(RegEle->MatrixBeta(i, 0) + CI_half_delta);
	}

}

void measure_F(REGELE* RegEle, int* Noccup) {

	long double d1 = (long double)*Noccup;
	long double d2 = (long double)RegEle->df;
	RegEle->Fvalue = (RegEle->SStot - RegEle->SSres) / d1;
	RegEle->Fvalue /= RegEle->SSres / d2;
	long double I_subscript = (d1 * RegEle->Fvalue) / (d1 * RegEle->Fvalue + d2);
	long double dtTestMesh = 1.0 / (long double)tFTestMesh;
	long double incomplete_beta = 0.0;
	long double complete_beta = 0.0;
	for (long double x = 0.5 * dtTestMesh; x < 1.0; x += dtTestMesh) {
		long double a = 0.5 * d1;
		long double b = 0.5 * d2;
		complete_beta += beta_func_integrand(&x, &a, &b) * dtTestMesh;
		if (x < I_subscript) {
			incomplete_beta = complete_beta;
		}
	}

	RegEle->p_Fvalue = 1.0 - incomplete_beta / complete_beta;

}

void check_Frequentist(REGELE* RegEle, long double* critic_p) {

	RegEle->full_Frequentist = true;

	if (RegEle->p_Fvalue > *critic_p) {
		RegEle->full_Frequentist = false;
	}

	if (RegEle->full_Frequentist) {
		for (int i = 0; i < (signed)RegEle->fit_coeff_p.size(); i++) {
			if (RegEle->fit_coeff_p[i] > *critic_p) {
				RegEle->full_Frequentist = false;
				break;
			}
		}
	}

}

void do_regression(REGELE* RegEle, std::vector<DATASET>* dataset, int* Noccup, long double* CI_level, long double* critic_p, bool* for_test) {

	initialize_regression(RegEle, dataset, Noccup);
	set_data(RegEle, dataset, Noccup);
	RegEle->df = (signed)dataset->size() - (*Noccup + 1);
	if (!*for_test) {
		fit_any(RegEle, dataset, Noccup);
	}
	measure_SStot_SSres(RegEle, dataset);
	measuere_r2(RegEle, dataset);
	if (!*for_test) {
		measure_t(RegEle, Noccup);
		measure_p(RegEle);
		measure_confidence_interval(RegEle, Noccup, CI_level);
		measure_F(RegEle, Noccup);
		check_Frequentist(RegEle, critic_p);
	}

}

long double get_stirling_approxiamtion(int* n) {

	return (long double)(*n) * logl((long double)(*n)) - (long double)(*n);

}

bool too_big_multinomial(int* n, std::vector<int>* m) {

	long double stirling = get_stirling_approxiamtion(n);
	for (int im = 0; im < m->size(); im++) {
		stirling -= get_stirling_approxiamtion(&(*m)[im]);
	}

	if (stirling > logl(powl(10.0, 16.0))) {
		return true;
	}
	else {
		return false;
	}

}

signed long long int get_multinomial(int* n, std::vector<int>* m) {

	int biggest_m = -1;
	int target_first_m = 1;
	for (int im = 0; im < m->size(); im++) {
		if ((*m)[im] > biggest_m) {
			biggest_m = (*m)[im];
			target_first_m = im;
		}
	}

	std::vector<int> up, down;
	std::vector<bool> down_chk;
	for (int i = biggest_m + 1; i <= *n; i++) {
		up.push_back(i);
	}
	for (int im = 0; im < m->size(); im++) {
		if (im != target_first_m) {
			for (int i = 2; i <= (*m)[im]; i++) {
				down.push_back(i);
				down_chk.push_back(false);
			}
		}
	}

	while (true) {
		int count = 0;
		for (int iu = 0; iu < (signed)up.size(); iu++) {
			for (int id = 0; id < (signed)down.size(); id++) {
				if (!down_chk[id] && up[iu] % down[id] == 0) {
					up[iu] /= down[id];
					down_chk[id] = true;
					count++;
				}
			}
		}
		if (count == 0) {
			break;
		}
	}

	signed long long int bino = 1;
	for (int iu = 0; iu < (signed)up.size(); iu++) {
		bino *= up[iu];
		for (int id = 0; id < (signed)down.size(); id++) {
			if (!down_chk[id] && bino % down[id] == 0) {
				bino /= down[id];
				down_chk[id] = true;
			}
		}
	}

	signed long long int down_ttl = 1;
	for (int id = 0; id < (signed)down.size(); id++) {
		if (!down_chk[id]) {
			down_ttl *= down[id];
		}
	}

	bino /= down_ttl;

	if (up.size() > 0) {
		std::vector<int>().swap(up);
	}
	if (down.size() > 0) {
		std::vector<int>().swap(down);
	}
	if (down_chk.size() > 0) {
		std::vector<bool>().swap(down_chk);
	}

	return bino;

}

long double get_multinomial_ld(int* n, std::vector<int>* m) {

	int biggest_m = -1;
	int target_first_m = 1;
	for (int im = 0; im < m->size(); im++) {
		if ((*m)[im] > biggest_m) {
			biggest_m = (*m)[im];
			target_first_m = im;
		}
	}

	std::vector<int> up, down;
	std::vector<bool> down_chk;
	for (int i = biggest_m + 1; i <= *n; i++) {
		up.push_back(i);
	}
	for (int im = 0; im < m->size(); im++) {
		if (im != target_first_m) {
			for (int i = 2; i <= (*m)[im]; i++) {
				down.push_back(i);
				down_chk.push_back(false);
			}
		}
	}

	while (true) {
		int count = 0;
		for (int iu = 0; iu < (signed)up.size(); iu++) {
			for (int id = 0; id < (signed)down.size(); id++) {
				if (!down_chk[id] && up[iu] % down[id] == 0) {
					up[iu] /= down[id];
					down_chk[id] = true;
					count++;
				}
			}
		}
		if (count == 0) {
			break;
		}
	}

	long double bino = 1.0;
	for (int iu = 0; iu < (signed)up.size(); iu++) {
		bino *= (long double)up[iu];
		for (int id = 0; id < (signed)down.size(); id++) {
			if (!down_chk[id] && remainderl(bino, (long double)down[id]) == 0) {
				bino /= (long double)down[id];
				down_chk[id] = true;
			}
		}
	}

	long double down_ttl = 1.0;
	for (int id = 0; id < (signed)down.size(); id++) {
		if (!down_chk[id]) {
			down_ttl *= (long double)down[id];
		}
	}

	bino /= down_ttl;

	if (up.size() > 0) {
		std::vector<int>().swap(up);
	}
	if (down.size() > 0) {
		std::vector<int>().swap(down);
	}
	if (down_chk.size() > 0) {
		std::vector<bool>().swap(down_chk);
	}

	return bino;

}

long double folder(long double v) {

	long double vv;
	if (v < 0.0) {
		while (true) {
			v += 1.0;
			if (v >= 0.0 && v < 1.0) {
				break;
			}
		}
		vv = v;
	}
	else if (v >= 0.0 && v < 1.0) {
		vv = v;
	}
	else {
		while (true) {
			v -= 1.0;
			if (v >= 0.0 && v < 1.0) {
				break;
			}
		}
		vv = v;
	}

	return vv;

}

long double erfinv(long double* x) {

	long double w, p;
	w = -logl((1.0 - *x) * (1.0 + *x));
	if (w < 5.000000) {
		w = w - 2.500000;
		p = 2.81022636e-08;
		p = 3.43273939e-07 + p * w;
		p = -3.5233877e-06 + p * w;
		p = -4.39150654e-06 + p * w;
		p = 0.00021858087 + p * w;
		p = -0.00125372503 + p * w;
		p = -0.00417768164 + p * w;
		p = 0.246640727 + p * w;
		p = 1.50140941 + p * w;
	}
	else {
		w = sqrtl(w) - 3.000000;
		p = -0.000200214257;
		p = 0.000100950558 + p * w;
		p = 0.00134934322 + p * w;
		p = -0.00367342844 + p * w;
		p = 0.00573950773 + p * w;
		p = -0.0076224613 + p * w;
		p = 0.00943887047 + p * w;
		p = 1.00167406 + p * w;
		p = 2.83297682 + p * w;
	}

	return p * *x;

}

long double get_v_special(std::string* ftype, long double* x) {

	if (*ftype == "logit" || *ftype == "isigmoid") {
		return(logl(*x / (1.0 - *x)));
	}
	else if (*ftype == "ilogit" || *ftype == "sigmoid") {
		return(1.0 / (1.0 + exp(-*x)));
	}
	else if (*ftype == "probit" || *ftype == "icdf") {
		long double newx = 2.0 * *x - 1.0;
		return(sqrtl(2.0) * erfinv(&newx));
	}
	else if (*ftype == "iprobit" || *ftype == "cdf") {
		return(0.5 * (1.0 + erfl(*x / sqrtl(2.0))));
	}
	else if (*ftype == "cloglog" || *ftype == "icexpexp") {
		return(logl(-logl(1.0 - *x)));
	}
	else if (*ftype == "icloglog" || *ftype == "cexpexp") {
		return(1.0 - expl(-expl(*x)));
	}
	else if (*ftype == "cauchit" || *ftype == "itangenttype") {
		return(tanl(PI * (*x - 0.5)));
	}
	else if (*ftype == "icauchit" || *ftype == "tangenttype") {
		return(atanl(*x) / PI + 0.5);
	}
	else if (*ftype == "nloglog" || *ftype == "loglog" || *ftype == "iexpexp") {
		return(-logl(-logl(*x)));
	}
	else if (*ftype == "inloglog" || *ftype == "iloglog" || *ftype == "expexp") {
		return(expl(-expl(-*x)));
	}
	else {
		return 0.0;
	}

}

long double get_v(ORTHFUNCTYPE* functype, long double* x) {

	if (functype->FName == "pow") {
		return powl(*x + functype->coeff[0], functype->coeff[1]);
	}
	else if (functype->FName == "log") {
		return powl(logl(*x + functype->coeff[0]) / logl(functype->coeff[1]), functype->coeff[2]);
	}
	else if (functype->FName == "sin") {
		return powl(sinl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "cos") {
		return powl(cosl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "tan") {
		return powl(tanl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "sinh") {
		return powl(sinhl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "cosh") {
		return powl(coshl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "tanh") {
		return powl(tanhl(powl(functype->coeff[1] * (*x + functype->coeff[0]), functype->coeff[2])), functype->coeff[3]);
	}
	else if (functype->FName == "exp") {
		return powl(powl(functype->coeff[3], functype->coeff[2] * powl(*x + functype->coeff[0], functype->coeff[1])), functype->coeff[4]);
	}
	else if (functype->FName == "erf") {
		return powl(erfl(functype->coeff[1] * (*x + functype->coeff[0])), functype->coeff[2]);
	}
	else if (functype->FName == "erfc") {
		return powl(erfcl(functype->coeff[1] * (*x + functype->coeff[0])), functype->coeff[2]);
	}
	else if (functype->FName == "abs") {
		return fabsl(powl(*x + functype->coeff[0], functype->coeff[1]));
	}
	else {
		return 0;
	}

}

bool isNumber(const std::string& s) {

	if (s.empty()) { 
		return false; 
	}

	std::istringstream iss(s);
	double d;
	char c;

	if (!(iss >> d)) {
		return false;
	}

	return !(iss >> c);

}

bool is_subset(std::vector<std::string>* superset, std::vector<std::string>* subset) {

	bool subset_check = true;

	for (int i = 0; i < (signed)subset->size(); i++) {
		bool exists = false;
		for (int j = 0; j < (signed)superset->size(); j++) {
			if ((*superset)[j] == (*subset)[i]) {
				exists = true;
				break;
			}
		}
		if (!exists) {
			subset_check = false;
			break;
		}
	}

	return subset_check;

}

void estimate_jobsize(bool* xpr_ld_get, signed long long int* jobsize, long double* jobsize_ld, int* Nfeatures, int* NVariables) {

	*xpr_ld_get = false;
	jobsize_step.clear();
	jobsize_step_ld.clear();

	if (NVariables[0] < 1) {
		NVariables[0] = 1;
	}
	if (NVariables[1] > *Nfeatures) {
		NVariables[1] = *Nfeatures;
	}

	*jobsize = 0;
	*jobsize_ld = 0.0;

	for (int taken = NVariables[0]; taken <= NVariables[1]; taken++) {
		std::vector<int> recipe_jobcontrol;
		recipe_jobcontrol.push_back(taken);
		recipe_jobcontrol.push_back(*Nfeatures - taken);
		signed long long int djobsize = get_multinomial(Nfeatures, &recipe_jobcontrol);
		*jobsize += djobsize;
		if (too_big_multinomial(Nfeatures, &recipe_jobcontrol) || *jobsize < 0) {
			*jobsize_ld += get_multinomial_ld(Nfeatures, &recipe_jobcontrol);
			*xpr_ld_get = true;
		}
		else {
			*jobsize_ld += (long double)djobsize;
		}
		jobsize_step.push_back(*jobsize);
		jobsize_step_ld.push_back(*jobsize_ld);
		if (recipe_jobcontrol.size() > 0) {
			std::vector<int>().swap(recipe_jobcontrol);
		}
	}

}

void set_pr_jobs(bool* xpr_ld_get, signed long long int* jobsize, long double* jobsize_ld) {

	pr_job = new signed long long int* [pr_size];
	pr_job_size = new signed long long int[pr_size];
	for (int i = 0; i < pr_size; i++) {
		pr_job[i] = new signed long long int[2];
		pr_job[i][0] = -1;
		pr_job[i][1] = -1;
		pr_job_size[i] = 0;
	}

	pr_job_ld = new long double* [pr_size];
	pr_job_size_ld = new long double[pr_size];
	for (int i = 0; i < pr_size; i++) {
		pr_job_ld[i] = new long double[2];
		pr_job_ld[i][0] = -1.0;
		pr_job_ld[i][1] = -1.0;
		pr_job_size_ld[i] = 0.0;
	}

	if (!*xpr_ld_get) {

		signed long long int prjobsize_per_prrank = *jobsize / (signed long long int)pr_size;
		signed long long int remain_prranks = *jobsize - prjobsize_per_prrank * (signed long long int)pr_size;
		for (int i = 0; i < pr_size; i++) {
			pr_job_size[i] = prjobsize_per_prrank;
			if ((signed long long int)i < remain_prranks) {
				pr_job_size[i]++;
			}
		}

		pr_job[0][0] = 0;
		for (int i = 1; i < pr_size; i++) {
			pr_job[i][0] = 0;
			for (int j = 0; j < i; j++) {
				pr_job[i][0] += pr_job_size[j];
			}
		}

		for (int i = 0; i < pr_size - 1; i++) {
			pr_job[i][1] = pr_job[i + 1][0] - 1;
		}
		pr_job[pr_size - 1][1] = *jobsize - 1;

		int pr_job_i1;
		for (int i = 0; i < pr_size; i++) {
			if (pr_job_size[i] > 0) {
				pr_job_i1 = pr_job[i][1];
			}
			else {
				pr_job[i][0] = pr_job_i1;
				pr_job[i][1] = pr_job_i1;
				pr_job_size[i] = 1;
			}
		}

		finite_job = true;
		if (pr_job[pr_id][0] == -1 || pr_job[pr_id][1] == -1) {
			finite_job = false;
		}

	}
	else {

		long double prjobsize_per_prrank_ld = *jobsize_ld / (long double)pr_size;
		long double remain_prranks_ld = *jobsize_ld - prjobsize_per_prrank_ld * (long double)pr_size;
		for (int i = 0; i < pr_size; i++) {
			pr_job_size_ld[i] = prjobsize_per_prrank_ld;
			if ((long double)i < remain_prranks_ld) {
				pr_job_size_ld[i]++;
			}
		}

		pr_job_ld[0][0] = 0;
		for (int i = 1; i < pr_size; i++) {
			pr_job_ld[i][0] = 0;
			for (int j = 0; j < i; j++) {
				pr_job_ld[i][0] += pr_job_size_ld[j];
			}
		}

		for (int i = 0; i < pr_size - 1; i++) {
			pr_job_ld[i][1] = pr_job_ld[i + 1][0] - 1;
		}
		pr_job_ld[pr_size - 1][1] = *jobsize_ld - 1;

		for (int i = 0; i < pr_size; i++) {
			if (pr_job_size_ld[i] <= 0) {
				pr_job_ld[i][0] = -1;
				pr_job_ld[i][1] = -1;
			}
		}

		finite_job = true;
		if (pr_job_ld[pr_id][0] == -1 || pr_job_ld[pr_id][1] == -1) {
			finite_job = false;
		}

	}

	
	SigmaMinusXk = new signed long long int[2];
	SigmaMinusXk_ld = new long double[2];

}

void initialize_occp(std::vector<int>* arrange, std::vector<int>* occup, int* Nfeatures) {

	arrange->clear();
	for (int i = 0; i < *Nfeatures; i++) {
		arrange->push_back(-1);
	}
	occup->clear();
	occup->push_back(0);
	occup->push_back(0);

}

void write_arrange(std::vector<REGTYPE>* regtype, std::vector<int>* arrange) {

	for (int ireg = 0; ireg < (signed)regtype->size(); ireg++) {
		TERMANA objTERMANA;
		objTERMANA.count = 0;
		objTERMANA.min_absz = VInfinite;
		objTERMANA.ave_absz = 0.0;
		objTERMANA.max_absz = -VInfinite;
		for (int iarr = 0; iarr < (signed)arrange->size(); iarr++) {
			if ((*arrange)[iarr] == 1) {
				objTERMANA.name.push_back((*regtype)[0].termana[iarr].name[0]);
			}
		}
		(*regtype)[ireg].termana.push_back(objTERMANA);
	}

}

void call_XPR(std::vector<int>* arrange, signed long long int* jobsize, std::vector<int>* occup, signed long long int* index) {

	signed long long int Ni[2];
	signed long long int Nbar;
	signed long long int index_dynamic = *index;

	for (int iX = 0; iX < (signed)occup->size(); iX++) {
		SigmaMinusXk[iX] = 0;
	}

	for (signed long long int i = 1; i <= arrange->size(); i++) {

		if (i == 1) {
			Ni[0] = *jobsize;
			Ni[1] = (signed long long int) arrange->size();
		}
		else {
			Ni[0] = Nbar;
			Ni[1] = (signed long long int) arrange->size() - i + 1;
		}

		signed long long int Nixk = 0;
		signed long long int Nixk_buf = 0;
		signed long long int DNiXm = 0;
		for (int m = 1; m <= (signed)occup->size(); m++) {
			signed long long int multiplier = (signed long long int)(*occup)[m - 1] - SigmaMinusXk[m - 1];
			if (multiplier < 0) {
				multiplier = 0;
			}
			DNiXm = Ni[0] * multiplier;
			if (DNiXm > 0) {
				DNiXm /= Ni[1];
			}
			else {
				long double DNiXm_ld = (long double)Ni[0] * (long double)multiplier / (long double)Ni[1];
				DNiXm = (signed long long int)DNiXm_ld;
			}
			Nixk_buf = Nixk;
			Nixk += DNiXm;
			if (index_dynamic > Nixk_buf && index_dynamic <= Nixk) {
				SigmaMinusXk[m - 1]++;
				(*arrange)[i - 1] = m - 1;
				index_dynamic -= Nixk_buf;
				Nbar = DNiXm;
				break;
			}
		}

	}

}

bool call_XPR_ld(std::vector<int>* arrange, long double* jobsize, std::vector<int>* occup, long double* index) {

	long double Ni[2];
	long double Nbar;
	long double index_dynamic = *index;

	for (int iX = 0; iX < (signed)occup->size(); iX++) {
		SigmaMinusXk_ld[iX] = 0;
	}

	bool stable = true;

	for (signed long long int i = 1; i <= arrange->size(); i++) {

		if (i == 1) {
			Ni[0] = *jobsize;
			Ni[1] = (signed long long int) arrange->size();
		}
		else {
			Ni[0] = Nbar;
			Ni[1] = (signed long long int) arrange->size() - i + 1;
		}

		long double Nixk = 0;
		long double Nixk_buf = 0;
		long double DNiXm = 0;
		bool m_taken = false;
		for (int m = 1; m <= (signed)occup->size(); m++) {
			signed long long int multiplier = (signed long long int)(*occup)[m - 1] - SigmaMinusXk_ld[m - 1];
			if (multiplier < 0) {
				multiplier = 0;
			}
			DNiXm = Ni[0] * multiplier;
			if (DNiXm > 0) {
				DNiXm /= Ni[1];
			}
			else {
				long double DNiXm_ld = (long double)Ni[0] * (long double)multiplier / (long double)Ni[1];
				DNiXm = (signed long long int)DNiXm_ld;
			}
			Nixk_buf = Nixk;
			Nixk += DNiXm;
			if (index_dynamic > Nixk_buf && index_dynamic <= Nixk) {
				SigmaMinusXk_ld[m - 1]++;
				(*arrange)[i - 1] = m - 1;
				index_dynamic -= Nixk_buf;
				Nbar = DNiXm;
				m_taken = true;
				break;
			}
		}
		if (!m_taken) {
			int unoccup_count = 0;
			for (int iarr = 0; iarr < (signed)arrange->size(); iarr++) {
				if (unoccup_count < (*occup)[0]) {
					(*arrange)[iarr] = 0;
					unoccup_count++;
				}
				else {
					(*arrange)[iarr] = 1;
				}
			}
			break;
			stable = false;
		}

	}

	return stable;

}

void run_through_models(std::vector<REGTYPE>* regtype, bool* xpr_ld_get, int* Nfeatures, int* NVariables) {

	signed long long int RunningIndexDelta = 1;
	long double RunningIndexDelta_ld = 1.0;

	if (!*xpr_ld_get) {

		for (signed long long int running_index = pr_job[pr_id][0]; running_index <= pr_job[pr_id][1]; running_index += RunningIndexDelta) {
			int whichjobsize = -1;
			for (int ijobsize = 0; ijobsize < (signed)jobsize_step.size(); ijobsize++) {
				int sjobsize;
				if (ijobsize == 0) {
					sjobsize = 0;
				}
				else {
					sjobsize = jobsize_step[ijobsize - 1];
				}
				if (running_index >= sjobsize && running_index < jobsize_step[ijobsize]) {
					whichjobsize = ijobsize;
					break;
				}
			}
			if (whichjobsize == -1) {
				whichjobsize = (signed)jobsize_step.size() - 1;
			}
			occup[0] = *Nfeatures - (whichjobsize + NVariables[0]);
			occup[1] = whichjobsize + NVariables[0];
			signed long long int jobsize_step_taken;
			if (whichjobsize == 0) {
				jobsize_step_taken = 0;
			}
			else {
				jobsize_step_taken = jobsize_step[whichjobsize - 1];
			}
			signed long long int sub_running_index = running_index - jobsize_step_taken + 1;
			signed long long int sub_jobsize = jobsize_step[whichjobsize] - jobsize_step_taken;
			call_XPR(&arrange, &sub_jobsize, &occup, &sub_running_index);
			write_arrange(regtype, &arrange);
		}

	}
	else {
	
		for (long double running_index_ld = pr_job_ld[pr_id][1]; running_index_ld >= pr_job_ld[pr_id][0]; running_index_ld -= RunningIndexDelta_ld) {
			int whichjobsize = -1;
			for (int ijobsize = 0; ijobsize < (signed)jobsize_step_ld.size(); ijobsize++) {
				long double sjobsize;
				if (ijobsize == 0) {
					sjobsize = 0.0;
				}
				else {
					sjobsize = jobsize_step_ld[ijobsize - 1];
				}
				if (running_index_ld >= sjobsize && running_index_ld < jobsize_step_ld[ijobsize]) {
					whichjobsize = ijobsize;
					break;
				}
			}
			if (whichjobsize == -1) {
				whichjobsize = (signed)jobsize_step_ld.size() - 1;
			}
			occup[0] = *Nfeatures - (whichjobsize + NVariables[0]);
			occup[1] = whichjobsize + NVariables[0];
			long double jobsize_step_taken_ld;
			if (whichjobsize == 0) {
				jobsize_step_taken_ld = 0;
			}
			else {
				jobsize_step_taken_ld = jobsize_step_ld[whichjobsize - 1];
			}
			long double sub_running_index_ld = running_index_ld - jobsize_step_taken_ld + 1;
			long double sub_jobsize_ld = jobsize_step_ld[whichjobsize] - jobsize_step_taken_ld;
			if (!call_XPR_ld(&arrange, &sub_jobsize_ld, &occup, &sub_running_index_ld)) {
				if (pr_id == 0) {
					std::cout << pr_id << " !call_XPR_ld\n";
					std::cout.flush();
				}
				long double temporary_RunningIndexDelta_ld = 0.0001 * (pr_job_ld[pr_id][1] - pr_job_ld[pr_id][0]);
				int error_count = 0;
				while (true) {
					sub_running_index_ld -= temporary_RunningIndexDelta_ld;
					if (call_XPR_ld(&arrange, &sub_jobsize_ld, &occup, &sub_running_index_ld)) {
						break;
					}
					else {
						error_count++;
					}
					if (error_count >= 100) {
						break;
					}
				}
				if (error_count >= 100) {
					sub_running_index_ld = pr_job_ld[pr_id][0];
					call_XPR_ld(&arrange, &sub_jobsize_ld, &occup, &sub_running_index_ld);
					write_arrange(regtype, &arrange);
				}
			}
		}
	
	}

}

void get_mamssk(std::vector<long double>* mamssk, std::vector<long double>* dv, long double* VInfinite, MAT* mat, std::vector<std::string>* DoNotAnalyze, bool* AveOnly) {

	std::vector<long double>().swap(*mamssk);

	std::vector<long double> finite_dv;
	std::vector<long double> finite_content;
	for (int i = 0; i < (signed)dv->size(); i++) {
		bool analyze = true;
		for (int j = 0; j < (signed)DoNotAnalyze->size(); j++) {
			if (mat->atomname[i] == (*DoNotAnalyze)[j]) {
				analyze = false;
				break;
			}
		}
		if (!*AveOnly) {
			if (analyze && (*dv)[i] != 0.0 && mat->content[i] != 0.0) {
				finite_dv.push_back((*dv)[i]);
				finite_content.push_back(mat->content[i]);
			}
		}
		else {
			if (analyze && mat->content[i] != 0.0) {
				finite_dv.push_back((*dv)[i]);
				finite_content.push_back(mat->content[i]);
			}
		}
	}

	long double min = *VInfinite;
	long double ave = 0.0;
	long double max = -*VInfinite;

	long double content_sum = 0.0;
	for (int i = 0; i < (signed)finite_content.size(); i++) {
		content_sum += finite_content[i];
	}

	for (int i = 0; i < (signed)finite_dv.size(); i++) {
		if (finite_dv[i] < min) {
			min = finite_dv[i];
		}
		ave += finite_content[i] * finite_dv[i];
		if (finite_dv[i] > max) {
			max = finite_dv[i];
		}
	}
	ave /= content_sum;

	long double stddv = 0.0;
	for (int i = 0; i < (signed)finite_dv.size(); i++) {
		stddv += finite_content[i] * powl(finite_dv[i] - ave, 2.0);
	}
	stddv = sqrtl(stddv / content_sum);

	long double skew = 0.0;
	long double kurto = 0.0;
	if (fabsl(stddv) > 1 / *VInfinite) {
		for (int i = 0; i < (signed)finite_dv.size(); i++) {
			skew += finite_content[i] * powl((finite_dv[i] - ave) / stddv, 3.0);
			kurto += finite_content[i] * powl((finite_dv[i] - ave) / stddv, 4.0);
		}
		skew /= content_sum;
		kurto /= content_sum;
	}
	else {
		skew = 0.0;
		kurto = *VInfinite;
	}

	mamssk->push_back(min);
	mamssk->push_back(ave);
	mamssk->push_back(max);
	mamssk->push_back(stddv);
	mamssk->push_back(skew);
	mamssk->push_back(kurto);

	std::vector<long double>().swap(finite_dv);
	std::vector<long double>().swap(finite_content);

}

int get_WhichAnaDescriptor(std::string DesName, std::vector<ANADESCRIPTOR>* anadescriptor) {

	int WAD = -1;
	for (int k = 0; k < (signed)anadescriptor->size(); k++) {
		if ((*anadescriptor)[k].name == DesName) {
			WAD = k;
			break;
		}
	}

	return WAD;

}

void get_dv(std::vector<long double>* dv, MAT* mat, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, int* WhichAnaDescriptor, bool* AveOnly) {

	if ((*anadescriptor)[*WhichAnaDescriptor].ASB == "A") {
		if ((*anadescriptor)[*WhichAnaDescriptor].ides >= 0) {
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				dv->push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[*WhichAnaDescriptor].ides]);
			}
		}
		/* USER-DEFINED SUPERVISED FUNCS */
		else if ((*anadescriptor)[*WhichAnaDescriptor].ides == -1) {
			
			std::vector<long double> v1, mamssk1;
			
			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("OuterFillingnonquenched", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AveOuterFillingnonquenched = mamssk1[1];
			
			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("PoissonRatio", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AvePoissonRatio = mamssk1[1];

			long double AveOuterFillingnonquenchedAvePoissonRatio = AveOuterFillingnonquenched * AvePoissonRatio;
			long double AveOuterFillingnonquenchedAvePoissonRatioX = -195.7 * powl(AveOuterFillingnonquenchedAvePoissonRatio, 2.0) + 28.516 * AveOuterFillingnonquenchedAvePoissonRatio + 0.8503;

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("Filling", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AveFilling = mamssk1[1];
			
			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("Radius06", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double StdDvRadius06 = mamssk1[3];

			long double SmoothCategory = 0.163981508328047 + 0.530154930618852 * AveOuterFillingnonquenchedAvePoissonRatioX - 3.66811654584888 * powl(erfl(AveFilling - 0.1), 2.0) - 1.48394565038536 * powl(erfl(StdDvRadius06 - 0.5), 2.0);
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				dv->push_back(SmoothCategory);
			}

			std::vector<long double>().swap(v1);
			std::vector<long double>().swap(mamssk1);

		}
		/* USER-DEFINED SUPERVISED FUNCS */
		else if ((*anadescriptor)[*WhichAnaDescriptor].ides == -2) {

			std::vector<long double> v1, mamssk1;

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("OuterFillingnonquenched", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AveOuterFillingnonquenched = mamssk1[1];

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("PoissonRatio", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AvePoissonRatio = mamssk1[1];

			long double AveOuterFillingnonquenchedAvePoissonRatio = AveOuterFillingnonquenched * AvePoissonRatio;

			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				dv->push_back(AveOuterFillingnonquenchedAvePoissonRatio);
			}

			std::vector<long double>().swap(v1);
			std::vector<long double>().swap(mamssk1);

		}
		/* USER-DEFINED SUPERVISED FUNCS */
		else if ((*anadescriptor)[*WhichAnaDescriptor].ides == -3) {

			std::vector<long double> v1, mamssk1;

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("OuterFillingnonquenched", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AveOuterFillingnonquenched = mamssk1[1];

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("PoissonRatio", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AvePoissonRatio = mamssk1[1];

			long double AveOuterFillingnonquenchedAvePoissonRatio = AveOuterFillingnonquenched * AvePoissonRatio;
			long double AveOuterFillingnonquenchedAvePoissonRatioX = -195.7 * powl(AveOuterFillingnonquenchedAvePoissonRatio, 2.0) + 28.516 * AveOuterFillingnonquenchedAvePoissonRatio + 0.8503;

			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				dv->push_back(AveOuterFillingnonquenchedAvePoissonRatioX);
			}

			std::vector<long double>().swap(v1);
			std::vector<long double>().swap(mamssk1);

		}
		/* USER-DEFINED SUPERVISED FUNCS */
		else if ((*anadescriptor)[*WhichAnaDescriptor].ides == -4) {

			std::vector<long double> v1, mamssk1;

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("OuterFillingnonquenched", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AveOuterFillingnonquenched = mamssk1[1];

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("PoissonRatio", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AvePoissonRatio = mamssk1[1];

			long double AveOuterFillingnonquenchedAvePoissonRatio = AveOuterFillingnonquenched * AvePoissonRatio;
			long double AveOuterFillingnonquenchedAvePoissonRatioX = -195.7 * powl(AveOuterFillingnonquenchedAvePoissonRatio, 2.0) + 28.516 * AveOuterFillingnonquenchedAvePoissonRatio + 0.8503;

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("Mass", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AveMass = mamssk1[1];

			long double AveOuterFillingnonquenchedAvePoissonRatioXAveMass = AveOuterFillingnonquenchedAvePoissonRatioX * AveMass;
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				dv->push_back(AveOuterFillingnonquenchedAvePoissonRatioXAveMass);
			}

			std::vector<long double>().swap(v1);
			std::vector<long double>().swap(mamssk1);

		}
		/* USER-DEFINED SUPERVISED FUNCS */
		else if ((*anadescriptor)[*WhichAnaDescriptor].ides == -5) {

			std::vector<long double> v1, mamssk1;

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("OuterFillingnonquenched", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AveOuterFillingnonquenched = mamssk1[1];

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("PoissonRatio", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double AvePoissonRatio = mamssk1[1];

			long double AveOuterFillingnonquenchedAvePoissonRatio = AveOuterFillingnonquenched * AvePoissonRatio;
			long double AveOuterFillingnonquenchedAvePoissonRatioX = -195.7 * powl(AveOuterFillingnonquenchedAvePoissonRatio, 2.0) + 28.516 * AveOuterFillingnonquenchedAvePoissonRatio + 0.8503;

			v1.clear();
			mamssk1.clear();
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				v1.push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[get_WhichAnaDescriptor("Mass", anadescriptor)].ides]);
			}
			get_mamssk(&mamssk1, &v1, &VInfinite, mat, &DoNotAnalyze, AveOnly);
			long double StdDvMass = mamssk1[3];

			long double AveOuterFillingnonquenchedAvePoissonRatioXStdDvMass = AveOuterFillingnonquenchedAvePoissonRatioX * StdDvMass;
			for (int k = 0; k < (signed)mat->iatom.size(); k++) {
				dv->push_back(AveOuterFillingnonquenchedAvePoissonRatioXStdDvMass);
			}

			std::vector<long double>().swap(v1);
			std::vector<long double>().swap(mamssk1);

		}
	}
	else if ((*anadescriptor)[*WhichAnaDescriptor].ASB == "S") {
		for (int k = 0; k < (signed)mat->isite.size(); k++) {
			dv->push_back((*strdat)[mat->istr].sdescriptor[(*anadescriptor)[*WhichAnaDescriptor].ides].v[mat->isite[k]]);
		}
	}
	else if ((*anadescriptor)[*WhichAnaDescriptor].ASB == "B") {
		for (int k = 0; k < (signed)mat->iatom.size(); k++) {
			dv->push_back((*atomdat)[mat->iatom[k]].v[(*anadescriptor)[*WhichAnaDescriptor].idesB[mat->istr][mat->isite[k]]]);
		}
	}

}

long double calc_reg(FUNCTYPE* functype) {

	/*xbox_CO!!!*/

	long double value = 0.0;
	for (int i = 0; i < (signed)functype->coeff_2.size(); i++) {
		long double sub_value = functype->coeff_2[i];
		if (functype->power_2[i] != 0.0) {
			if (functype->ianades_2[i] == -1) {
				if (functype->ilogln_2[i] == 1) {
					sub_value = sub_value * powl(log10l(xbox_OM), functype->power_2[i]);
				}
				else if (functype->ilogln_2[i] == 2) {
					sub_value = sub_value * powl(logl(xbox_OM), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 1) {
					sub_value = sub_value * powl(powl(10.0, xbox_OM), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 2) {
					sub_value = sub_value * powl(powl(2.71828, xbox_OM), functype->power_2[i]);
				}
				else if (functype->ierf_2[i] != -1000) {
					sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_2[i]) * xbox_OM), functype->power_2[i]);
				}
				else {
					sub_value = sub_value * powl(xbox_OM, functype->power_2[i]);
				}
			}
			else { 
				if (functype->ilogln_2[i] == 1) {
					sub_value = sub_value * powl(log10l(xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ilogln_2[i] == 2) {
					sub_value = sub_value * powl(logl(xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 1) {
					sub_value = sub_value * powl(powl(10.0, xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ipow10exp_2[i] == 2) {
					sub_value = sub_value * powl(powl(2.71828, xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else if (functype->ierf_2[i] != -1000) {
					sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_2[i]) * xbox[functype->ianades_2[i]][functype->imamssk_2[i]]), functype->power_2[i]);
				}
				else {
					sub_value = sub_value * powl(xbox[functype->ianades_2[i]][functype->imamssk_2[i]], functype->power_2[i]);
				}
			}
		}
		value += sub_value;
	}

	long double value_sofar = value;
	if (functype->nonlinear_funcname != "") {
		long double super_value = 0.0;
		for (int i = 0; i < (signed)functype->coeff_1.size(); i++) {
			long double sub_value = functype->coeff_1[i];
			if (functype->power_1[i] != 0.0) {
				sub_value = sub_value * powl(get_v_special(&functype->nonlinear_funcname, &value), functype->power_1[i]);
			}
			super_value += sub_value;
		}
		long double ssuper_value = 0.0;
		for (int i = 0; i < (signed)functype->coeff_0.size(); i++) {
			long double sub_value = functype->coeff_0[i];
			if (functype->power_0[i] != 0.0) {
				sub_value = sub_value * powl(super_value, functype->power_0[i]);
			}
			ssuper_value += sub_value;
		}
		value_sofar = ssuper_value;
	}

	if (functype->added) {
		long double addedv = 0.0;
		for (int i = 0; i < (signed)functype->coeff_3.size(); i++) {
			long double sub_value = functype->coeff_3[i];
			if (functype->power_3[i] != 0.0) {
				if (functype->ianades_3[i] == -1) {
					if (functype->ilogln_3[i] == 1) {
						sub_value = sub_value * powl(log10l(xbox_OM), functype->power_3[i]);
					}
					else if (functype->ilogln_3[i] == 2) {
						sub_value = sub_value * powl(logl(xbox_OM), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 1) {
						sub_value = sub_value * powl(powl(10.0, xbox_OM), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 2) {
						sub_value = sub_value * powl(powl(2.71828, xbox_OM), functype->power_3[i]);
					}
					else if (functype->ierf_3[i] != -1000) {
						sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_3[i]) * xbox_OM), functype->power_3[i]);
					}
					else {
						sub_value = sub_value * powl(xbox_OM, functype->power_3[i]);
					}
				}
				else {
					if (functype->ilogln_3[i] == 1) {
						sub_value = sub_value * powl(log10l(xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ilogln_3[i] == 2) {
						sub_value = sub_value * powl(logl(xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 1) {
						sub_value = sub_value * powl(powl(10.0, xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ipow10exp_3[i] == 2) {
						sub_value = sub_value * powl(powl(2.71828, xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else if (functype->ierf_3[i] != -1000) {
						sub_value = sub_value * powl(erfl(powl(10.0, (long double)functype->ierf_3[i]) * xbox[functype->ianades_3[i]][functype->imamssk_3[i]]), functype->power_3[i]);
					}
					else {
						sub_value = sub_value * powl(xbox[functype->ianades_3[i]][functype->imamssk_3[i]], functype->power_3[i]);
					}
				}
			}
			addedv += sub_value;
		}
		value_sofar += addedv;
	}

	return value_sofar;

}

long double get_v_in_the_set(std::vector<int>* anaset, std::vector<int>* mamsskset, int* whichana, int* whichmamssk, std::vector<ANADESCRIPTOR>* anadescriptor, std::vector<long double>* test_anades_degree) {

	long double value = 0.0;
	bool exists = false;
	
	int target_ana = -1;
	for (int iana = 0; iana < (signed)anaset->size(); iana++) {
		if ((*anaset)[iana] == *whichana) {
			exists = true;
			target_ana = iana;
			break;
		}
	}

	if (exists) {
		if (*whichmamssk != (*mamsskset)[target_ana]) {
			exists = false;
		}
	}

	if (exists) {
		if (*whichmamssk == 0) {
			value = (*anadescriptor)[*whichana].min;
			if (get_floating_v) {
				floating_v_min[target_ana] = (*anadescriptor)[*whichana].min;
				floating_v_max[target_ana] = (*anadescriptor)[*whichana].min;
			}
		}
		else if (*whichmamssk == 1) {
			value = (*anadescriptor)[*whichana].min_ave + (*test_anades_degree)[target_ana] * ((*anadescriptor)[*whichana].max_ave - (*anadescriptor)[*whichana].min_ave);
			if (get_floating_v) {
				floating_v_min[target_ana] = (*anadescriptor)[*whichana].min_ave;
				floating_v_max[target_ana] = (*anadescriptor)[*whichana].max_ave;
			}
		}
		else if (*whichmamssk == 2) {
			value = (*anadescriptor)[*whichana].max;
			if (get_floating_v) {
				floating_v_min[target_ana] = (*anadescriptor)[*whichana].max;
				floating_v_max[target_ana] = (*anadescriptor)[*whichana].max;
			}
		}
		else if (*whichmamssk == 3) {
			value = (*anadescriptor)[*whichana].min_stdev + (*test_anades_degree)[target_ana] * ((*anadescriptor)[*whichana].max_stdev - (*anadescriptor)[*whichana].min_stdev);
			if (get_floating_v) {
				floating_v_min[target_ana] = (*anadescriptor)[*whichana].min_stdev;
				floating_v_max[target_ana] = (*anadescriptor)[*whichana].max_stdev;
			}
		}
		else if (*whichmamssk == 4) {
			value = (*anadescriptor)[*whichana].min_skew + (*test_anades_degree)[target_ana] * ((*anadescriptor)[*whichana].max_skew - (*anadescriptor)[*whichana].min_skew);
			if (get_floating_v) {
				floating_v_min[target_ana] = (*anadescriptor)[*whichana].min_skew;
				floating_v_max[target_ana] = (*anadescriptor)[*whichana].max_skew;
			}
		}
		else if (*whichmamssk == 5) {
			value = (*anadescriptor)[*whichana].min_kurto + (*test_anades_degree)[target_ana] * ((*anadescriptor)[*whichana].max_kurto - (*anadescriptor)[*whichana].min_kurto);
			if (get_floating_v) {
				floating_v_min[target_ana] = (*anadescriptor)[*whichana].min_kurto;
				floating_v_max[target_ana] = (*anadescriptor)[*whichana].max_kurto;
			}
		}
		if (get_floating_v) {
			floating_v[target_ana] = value;	
		}
	}
	else {
		if (*whichmamssk == 0) {
			value = (*anadescriptor)[*whichana].min;
		}
		else if (*whichmamssk == 1) {
			value = (*anadescriptor)[*whichana].ave;
		}
		else if (*whichmamssk == 2) {
			value = (*anadescriptor)[*whichana].max;
		}
		else if (*whichmamssk == 3) {
			value = (*anadescriptor)[*whichana].ave_stdev;
		}
		else if (*whichmamssk == 4) {
			value = (*anadescriptor)[*whichana].ave_skew;
		}
		else if (*whichmamssk == 5) {
			value = (*anadescriptor)[*whichana].ave_kurto;
		}
	}

	
	return value;

}

long double calc_reg_chain(REGTYPE* regtype, MAT* objMAT, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::vector<std::string>* CationNames, std::vector<int>* test_anades, std::vector<int>* test_mamssk, bool *test_OM, int *test_CO, std::vector<long double>* test_anades_degree, long double* test_OM_degree, long double* test_CO_degree, bool *verbose, bool* AveOnly) {

	std::vector<long double> funcele, funcele_last;
	for (int i = 0; i < (signed)regtype->funcchain.size(); i++) {
		funcele.clear();
		for (int j = 0; j < regtype->funcchain[i].func.size(); j++) {
			long double v = 0.0;
			if (regtype->funcchain[i].func[j].identified && regtype->funcchain[i].func[j].feature_identified) {
				int ift = regtype->funcchain[i].func[j].whichfeature;
				if (regtype->funcchain[i].feature[ift].type == 0) {
					std::vector<long double> dv, mamssk;
					if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
						get_dv(&dv, objMAT, atomdat, strdat, anadescriptor, &regtype->funcchain[i].feature[ift].i0anades, AveOnly);
						get_mamssk(&mamssk, &dv, VInfinite, objMAT, DoNotAnalyze, AveOnly);
						v = mamssk[regtype->funcchain[i].feature[ift].i0mamssk];
					}
					else {
						v = get_v_in_the_set(test_anades, test_mamssk, &regtype->funcchain[i].feature[ift].i0anades, &regtype->funcchain[i].feature[ift].i0mamssk, anadescriptor, test_anades_degree);
					}
					if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_A]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}
				else if (regtype->funcchain[i].feature[ift].type == 1) {
					long double v0 = 0.0;
					long double v1 = 0.0;
					if (regtype->funcchain[i].feature[ift].i1mamssk[0] >= 0) {
						std::vector<long double> dv, mamssk;
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							get_dv(&dv, objMAT, atomdat, strdat, anadescriptor, &regtype->funcchain[i].feature[ift].i1anades[0], AveOnly);
							get_mamssk(&mamssk, &dv, VInfinite, objMAT, DoNotAnalyze, AveOnly);
							v0 = mamssk[regtype->funcchain[i].feature[ift].i1mamssk[0]];
						}
						else {
							v0 = get_v_in_the_set(test_anades, test_mamssk, &regtype->funcchain[i].feature[ift].i1anades[0], &regtype->funcchain[i].feature[ift].i1mamssk[0], anadescriptor, test_anades_degree);
						}
						if (*verbose && (std::isnan(v0) || std::isinf(v0) || std::isinf(-v0))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_B]v0 is " << v0 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
						std::vector<long double>().swap(dv);
						std::vector<long double>().swap(mamssk);
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[0] == -1) {
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							long double cation_sum = 0.0;
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								cation_sum += objMAT->content[ia];
							}
							v0 = objMAT->anion_content / cation_sum;
						}
						else {
							if (!*test_OM) {
								v0 = OM_ave;
							}
							else {
								v0 = OM_min + *test_OM_degree * (OM_max - OM_min);
							}
						}
						if (*verbose && (std::isnan(v0) || std::isinf(v0) || std::isinf(-v0))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_C]v0 is " << v0 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[0] == -2) {
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							long double cation_sum = 0.0;
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								if (objMAT->atomname[ia] == (*CationNames)[regtype->funcchain[i].feature[ift].i1anades[0]]) {
									cation_sum += objMAT->content[ia];
								}
							}
							v0 = cation_sum / objMAT->anion_content;
						}
						else {
							if (*test_CO == -100) {
								v0 = CO_ave[regtype->funcchain[i].feature[ift].i1anades[0]];
							}
							else {
								v0 = CO_min[regtype->funcchain[i].feature[ift].i1anades[0]] + *test_CO_degree * (CO_max[regtype->funcchain[i].feature[ift].i1anades[0]] - CO_min[regtype->funcchain[i].feature[ift].i1anades[0]]);
							}
						}
						if (*verbose && (std::isnan(v0) || std::isinf(v0) || std::isinf(-v0))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_D]v0 is " << v0 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[0] == -3) {
						int ianades = regtype->funcchain[i].feature[ift].i1anades[0];
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							long double vv0 = 0.0;
							if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
											vv0 += objMAT->content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									vv0 += objMAT->content[ia];
								}
							}
							else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									bool found = false;
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
											found = true;
											break;
										}
									}
									if (!found) {
										vv0 += objMAT->content[ia];
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 3) {
								vv0 = (*anadescriptor)[ianades].UType_numerator_value;
							}
							long double vv1 = 0.0;
							if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
											vv1 += objMAT->content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									vv1 += objMAT->content[ia];
								}
							}
							else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									bool found = false;
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
											found = true;
											break;
										}
									}
									if (!found) {
										vv1 += objMAT->content[ia];
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 3) {
								vv1 = (*anadescriptor)[ianades].UType_denominator_value;
							}
							if ((*anadescriptor)[ianades].UType_m0d1 == 0) {
								v0 = vv0 * vv1;
							}
							else {
								v0 = vv0 / vv1;
							}
						}
						else {
							int null_mamssk = 1;
							v0 = get_v_in_the_set(test_anades, test_mamssk, &ianades, &null_mamssk, anadescriptor, test_anades_degree);
						}
						if (*verbose && (std::isnan(v0) || std::isinf(v0) || std::isinf(-v0))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_E]v0 is " << v0 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
					}
					if (regtype->funcchain[i].feature[ift].i1mamssk[1] >= 0) {
						std::vector<long double> dv, mamssk;
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							get_dv(&dv, objMAT, atomdat, strdat, anadescriptor, &regtype->funcchain[i].feature[ift].i1anades[1], AveOnly);
							get_mamssk(&mamssk, &dv, VInfinite, objMAT, DoNotAnalyze, AveOnly);
							v1 = mamssk[regtype->funcchain[i].feature[ift].i1mamssk[1]];
						}
						else {
							v1 = get_v_in_the_set(test_anades, test_mamssk, &regtype->funcchain[i].feature[ift].i1anades[1], &regtype->funcchain[i].feature[ift].i1mamssk[1], anadescriptor, test_anades_degree);
						}
						if (*verbose && (std::isnan(v1) || std::isinf(v1) || std::isinf(-v1))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_F]v1 is " << v1 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
						std::vector<long double>().swap(dv);
						std::vector<long double>().swap(mamssk);
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[1] == -1) {
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							long double cation_sum = 0.0;
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								cation_sum += objMAT->content[ia];
							}
							v1 = objMAT->anion_content / cation_sum;
						}
						else {
							if (!*test_OM) {
								v1 = OM_ave;
							}
							else {
								v1 = OM_min + *test_OM_degree * (OM_max - OM_min);
							}
						}
						if (*verbose && (std::isnan(v1) || std::isinf(v1) || std::isinf(-v1))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_G]v1 is " << v1 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[1] == -2) {
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							long double cation_sum = 0.0;
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								if (objMAT->atomname[ia] == (*CationNames)[regtype->funcchain[i].feature[ift].i1anades[1]]) {
									cation_sum += objMAT->content[ia];
								}
							}
							v1 = cation_sum / objMAT->anion_content;
						}
						else {
							if (*test_CO == -100) {
								v1 = CO_ave[regtype->funcchain[i].feature[ift].i1anades[1]];
							}
							else {
								v1 = CO_min[regtype->funcchain[i].feature[ift].i1anades[1]] + *test_CO_degree * (CO_max[regtype->funcchain[i].feature[ift].i1anades[1]] - CO_min[regtype->funcchain[i].feature[ift].i1anades[1]]);
							}
						}
						if (*verbose && (std::isnan(v1) || std::isinf(v1) || std::isinf(-v1))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_H]v1 is " << v1 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
					}
					else if (regtype->funcchain[i].feature[ift].i1mamssk[1] == -3) {
						int ianades = regtype->funcchain[i].feature[ift].i1anades[1];
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							long double vv0 = 0.0;
							if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
											vv0 += objMAT->content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									vv0 += objMAT->content[ia];
								}
							}
							else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									bool found = false;
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
											found = true;
											break;
										}
									}
									if (!found) {
										vv0 += objMAT->content[ia];
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 3) {
								vv0 = (*anadescriptor)[ianades].UType_numerator_value;
							}
							long double vv1 = 0.0;
							if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
											vv1 += objMAT->content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									vv1 += objMAT->content[ia];
								}
							}
							else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
									bool found = false;
									for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
											found = true;
											break;
										}
									}
									if (!found) {
										vv1 += objMAT->content[ia];
									}
								}
							}
							else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 3) {
								vv1 = (*anadescriptor)[ianades].UType_denominator_value;
							}
							if ((*anadescriptor)[ianades].UType_m0d1 == 0) {
								v1 = vv0 * vv1;
							}
							else {
								v1 = vv0 / vv1;
							}
						}
						else {
							int null_mamssk = 1;
							v1 = get_v_in_the_set(test_anades, test_mamssk, &ianades, &null_mamssk, anadescriptor, test_anades_degree);
						}
						if (*verbose && (std::isnan(v1) || std::isinf(v1) || std::isinf(-v1))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_I]v1 is " << v1 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
					}
					if (regtype->funcchain[i].feature[ift].i1m0d1 == 0) {
						v = v0 * v1;
					}
					else {
						v = v0 / v1;
					}
					if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_J]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
				}
				else if (regtype->funcchain[i].feature[ift].type == 2) {
					v = funcele_last[regtype->funcchain[i].feature[ift].i2whichlastfunc];
					if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_K]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
				}
				else if (regtype->funcchain[i].feature[ift].type == 3) {
					long double v0 = funcele_last[regtype->funcchain[i].feature[ift].i3whichlastfunc[0]];
					if (*verbose && (std::isnan(v0) || std::isinf(v0) || std::isinf(-v0))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_L]v0 is " << v0 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
					long double v1 = funcele_last[regtype->funcchain[i].feature[ift].i3whichlastfunc[1]];
					if (*verbose && (std::isnan(v1) || std::isinf(v1) || std::isinf(-v1))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_M]v1 is " << v1 << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
					if (regtype->funcchain[i].feature[ift].i3m0d1 == 0) {
						v = v0 * v1;
					}
					else {
						v = v0 / v1;
					}
					if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_N]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
				}
				else if (regtype->funcchain[i].feature[ift].type == -1) {
					if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
						long double cation_sum = 0.0;
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							cation_sum += objMAT->content[ia];
						}
						v = objMAT->anion_content / cation_sum;
					}
					else {
						if (!*test_OM) {
							v = OM_ave;
						}
						else {
							v = OM_min + *test_OM_degree * (OM_max - OM_min);
						}
					}
					if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_O]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
				}
				else if (regtype->funcchain[i].feature[ift].type == -2) {
					if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
						long double cation_sum = 0.0;
						for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
							if (objMAT->atomname[ia] == regtype->funcchain[i].feature[ift].CationName) {
								cation_sum += objMAT->content[ia];
							}
						}
						v = cation_sum / objMAT->anion_content;
					}
					else {
						int found_target_iC = -1;
						for (int ia = 0; ia < (signed)CationNames->size(); ia++) {
							if (regtype->funcchain[i].feature[ift].CationName == (*CationNames)[ia]) {
								found_target_iC = ia;
								break;
							}
						}
						if (*test_CO == -100) {
							v = CO_ave[found_target_iC];
						}
						else {
							v = CO_min[found_target_iC] + *test_CO_degree * (CO_max[found_target_iC] - CO_min[found_target_iC]);
						}
					}
					if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_P]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
				}
				else if (regtype->funcchain[i].feature[ift].type == -3) {
					int ianades = regtype->funcchain[i].feature[ift].i0anades;
					if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
						long double v0 = 0.0;
						if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 0) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
										v0 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 1) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								v0 += objMAT->content[ia];
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 2) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								bool found = false;
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_numerator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_numerator_AtomKinds[ja]) {
										found = true;
										break;
									}
								}
								if (!found) {
									v0 += objMAT->content[ia];
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_numerator_S0A1E2N3 == 3) {
							v0 = (*anadescriptor)[ianades].UType_numerator_value;
						}
						long double v1 = 0.0;
						if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 0) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
										v1 += objMAT->content[ia];
										break;
									}
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 1) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								v1 += objMAT->content[ia];
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 2) {
							for (int ia = 0; ia < (signed)objMAT->content.size(); ia++) {
								bool found = false;
								for (int ja = 0; ja < (signed)(*anadescriptor)[ianades].UType_denominator_AtomKinds.size(); ja++) {
									if (objMAT->atomname[ia] == (*anadescriptor)[ianades].UType_denominator_AtomKinds[ja]) {
										found = true;
										break;
									}
								}
								if (!found) {
									v1 += objMAT->content[ia];
								}
							}
						}
						else if ((*anadescriptor)[ianades].UType_denominator_S0A1E2N3 == 3) {
							v1 = (*anadescriptor)[ianades].UType_denominator_value;
						}
						if ((*anadescriptor)[ianades].UType_m0d1 == 0) {
							v = v0 * v1;
						}
						else {
							v = v0 / v1;
						}
						if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
							if (pr_id == 0) {
								std::cout << "[ERROR_Q]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
								std::cout.flush();
							}
						}
					}
					else {
						int null_mamssk = 1;
						v = get_v_in_the_set(test_anades, test_mamssk, &ianades, &null_mamssk, anadescriptor, test_anades_degree);
					
					}
					if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
						if (pr_id == 0) {
							std::cout << "[ERROR_R]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
							std::cout.flush();
						}
					}
				}
				if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
					if (pr_id == 0) {
						std::cout << "[ERROR_S1]v is " << v << "... at funcchain " << i << " feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
						std::cout.flush();
					}
				}
				long double vbuf = v;
				v = get_v(&(regtype->orthfunctype)[regtype->funcchain[i].func[j].whichorthfunctype], &v);
				if (*verbose && (std::isnan(v) || std::isinf(v) || std::isinf(-v))) {
					if (pr_id == 0) {
						std::cout << "[ERROR_S2]vbuf is " << vbuf << " , v is " << v << "... at funcchain " << i << " : " << (regtype->orthfunctype)[regtype->funcchain[i].func[j].whichorthfunctype].FName << " : feature : " << ift << " : " << regtype->funcchain[i].feature[ift].name << "\n";
						std::cout.flush();
					}
				}
			}
			funcele.push_back(v);
		}
		funcele_last.clear();
		funcele_last = funcele;
	}

	long double value = 0.0;
	for (int j = 0; j < regtype->funcchain[regtype->funcchain.size() - 1].func.size(); j++) {
		if (*verbose && (std::isnan(funcele[j]) || std::isinf(funcele[j]) || std::isinf(-funcele[j]))) {
			if (pr_id == 0) {
				std::cout << "[ERROR_T]v is " << funcele[j] << " ... at funcchain " << regtype->funcchain.size() - 1 << " at func : " << j << "\n";
				std::cout.flush();
			}
		}
		value += regtype->funcchain[regtype->funcchain.size() - 1].func[j].coeff * funcele[j];
		if (*verbose && (std::isnan(value) || std::isinf(value) || std::isinf(-value))) {
			if (pr_id == 0) {
				std::cout << "[ERROR_U]v is " << value << " ... at funcchain " << regtype->funcchain.size() - 1 << " at func : " << j << "\n";
				std::cout.flush();
			}
		}
	}
	value += regtype->funcchain[regtype->funcchain.size() - 1].intercept;
	
	std::vector<long double>().swap(funcele);
	std::vector<long double>().swap(funcele_last);

	if (regtype->funcchain[regtype->funcchain.size() - 1].allow_beta) {
		std::string ftype = "i" + regtype->beta_target;
		value = regtype->funcchain[regtype->funcchain.size() - 1].score_inverse_slope * get_v_special(&ftype, &value) + regtype->funcchain[regtype->funcchain.size() - 1].score_inverse_intercept;
		value = regtype->funcchain[regtype->funcchain.size() - 1].revise_slope * value + regtype->funcchain[regtype->funcchain.size() - 1].revise_intercept;
	}

	return value;

}

long double calc_bag_chain(std::vector<long double>* v_series, std::vector<long double>* v_series_acc, long double* v_stdev, std::vector<REGTYPE>* regtype, BAGTYPE* bagtype, MAT* objMAT, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::vector<std::string>* CationNames, std::vector<int>* test_anades, std::vector<int>* test_mamssk, bool* test_OM, int* test_CO, std::vector<long double>* test_anades_degree, long double* test_OM_degree, long double* test_CO_degree, bool* verbose, bool* post_if_exists, bool* AveOnly, long double *post_v) {

	std::vector<long double>().swap(*v_series);
	std::vector<long double>().swap(*v_series_acc);
	v_series->clear();
	v_series_acc->clear();

	*post_v = 0.0;
	if (*post_if_exists) {
		for (int pk = 0; pk < (signed)bagtype->which_regtype.size(); pk++) {
			int k = bagtype->which_regtype[pk];
			for (int pkk = 0; pkk < (signed)(*regtype)[k].postfunc.size(); pkk++) {
				if ((*regtype)[k].postfunc[pkk].identified && (*regtype)[k].postfunc[pkk].feature_identified) {
					int ift = (*regtype)[k].postfunc[pkk].whichfeature;
					if ((*regtype)[k].funcchain[0].feature[ift].type == 0) {
						std::vector<long double> dv, mamssk;
						long double post_v_element = 0.0;
						if (test_anades->size() == 0 && !*test_OM && *test_CO == -100) {
							get_dv(&dv, objMAT, atomdat, strdat, anadescriptor, &(*regtype)[k].funcchain[0].feature[ift].i0anades, AveOnly);
							get_mamssk(&mamssk, &dv, VInfinite, objMAT, DoNotAnalyze, AveOnly);
							post_v_element = mamssk[(*regtype)[k].funcchain[0].feature[ift].i0mamssk];
						}
						else {
							post_v_element = get_v_in_the_set(test_anades, test_mamssk, &(*regtype)[k].funcchain[0].feature[ift].i0anades, &(*regtype)[k].funcchain[0].feature[ift].i0mamssk, anadescriptor, test_anades_degree);
						}
						long double post_v_element_funcized = get_v(&((*regtype)[k].orthfunctype)[(*regtype)[k].postfunc[pkk].whichorthfunctype], &post_v_element);
						*post_v += post_v_element_funcized;
					}
				}
			}
		}
	}

	for (int pk = 0; pk < (signed)bagtype->which_regtype.size(); pk++) {
		int k = bagtype->which_regtype[pk];
		v_series->push_back(calc_reg_chain(&(*regtype)[k], objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, test_anades, test_mamssk, test_OM, test_CO, test_anades_degree, test_OM_degree, test_CO_degree, verbose, AveOnly)/* + post_v*/);
	}
	
	int count = 0;
	long double v = 0.0;
	if (bagtype->stackingensemble.size() == 0) {
		for (int i = 0; i < (signed)v_series->size(); i++) {
			if (!std::isnan((*v_series)[i]) && !std::isinf((*v_series)[i]) && !std::isinf(-(*v_series)[i])) {
				count++;
				v += (*v_series)[i];
			}
			long double acc = 0.0;
			int acc_count = 0;
			for (int ii = 0; ii <= i; ii++) {
				if (!std::isnan((*v_series)[ii]) && !std::isinf((*v_series)[ii]) && !std::isinf(-(*v_series)[ii])) {
					acc_count++;
					acc += (*v_series)[ii];
				}
			}
			acc /= (long double)acc_count;
			v_series_acc->push_back(acc);
		}
		v /= (long double)count;
		*v_stdev = 0.0;
		for (int i = 0; i < (signed)v_series->size(); i++) {
			if (!std::isnan((*v_series)[i]) && !std::isinf((*v_series)[i]) && !std::isinf(-(*v_series)[i])) {
				*v_stdev += powl((*v_series)[i] - v, 2.0);
			}
		}
		*v_stdev /= (long double)count;
		*v_stdev = sqrtl(*v_stdev);
	}
	else {
		v = bagtype->stackingensemble[bagtype->stackingensemble.size() - 1].coeff[0];
		long double simple_ave = 0.0;
		bool nan_detected = false;
		for (int i = 0; i < (signed)v_series->size(); i++) {
			v += bagtype->stackingensemble[bagtype->stackingensemble.size() - 1].coeff[i + 1] * (*v_series)[i];
			simple_ave += (*v_series)[i];
			if (std::isnan((*v_series)[i]) || std::isinf((*v_series)[i]) || std::isinf(-(*v_series)[i])) {
				nan_detected = true;
				break;
			}
		}
		if (!nan_detected) {
			simple_ave /= (long double)(signed)v_series->size();
			for (int i = 0; i < (signed)v_series->size(); i++) {
				long double acc = bagtype->stackingensemble[i].coeff[0];
				for (int j = 0; j <= i; j++) {
					acc += bagtype->stackingensemble[i].coeff[j + 1] * (*v_series)[j];
				}
				v_series_acc->push_back(acc);
			}
			*v_stdev = 0.0;
			for (int i = 0; i < (signed)v_series->size(); i++) {
				*v_stdev += fabsl(bagtype->stackingensemble[bagtype->stackingensemble.size() - 1].coeff[i + 1]) * powl((*v_series)[i] - simple_ave, 2.0) / bagtype->stackingensemble[bagtype->stackingensemble.size() - 1].coeff_sum;
			}
			*v_stdev = sqrtl(*v_stdev);
			if (std::isnan(*v_stdev) || std::isinf(*v_stdev) || std::isinf(-*v_stdev)) {
				*v_stdev = *VInfinite;
			}
		}
		else {
			v = 0.0;
			count = 0;
			for (int i = 0; i < (signed)v_series->size(); i++) {
				if (!std::isnan((*v_series)[i]) && !std::isinf((*v_series)[i]) && !std::isinf(-(*v_series)[i]) && (*v_series)[i] >= bagtype->obs_v_min && (*v_series)[i] <= bagtype->obs_v_max) {
					count++;
					v += (*v_series)[i];
				}
				long double acc = 0.0;
				int acc_count = 0;
				for (int ii = 0; ii <= i; ii++) {
					if (!std::isnan((*v_series)[ii]) && !std::isinf((*v_series)[ii]) && !std::isinf(-(*v_series)[ii]) && (*v_series)[ii] >= bagtype->obs_v_min && (*v_series)[ii] <= bagtype->obs_v_max) {
						acc_count++;
						acc += (*v_series)[ii];
					}
				}
				acc /= (long double)acc_count;
				v_series_acc->push_back(acc);
			}
			v /= (long double)count;
			*v_stdev = 0.0;
			for (int i = 0; i < (signed)v_series->size(); i++) {
				if (!std::isnan((*v_series)[i]) && !std::isinf((*v_series)[i]) && !std::isinf(-(*v_series)[i]) && (*v_series)[i] >= bagtype->obs_v_min && (*v_series)[i] <= bagtype->obs_v_max) {
					*v_stdev += powl((*v_series)[i] - v, 2.0);
				}
			}
			*v_stdev /= (long double)count;
			*v_stdev = sqrtl(*v_stdev);
		}
	}
	
	return v + *post_v;

}

void calc_bag_chain_partial_test(std::vector<long double>* v_series, std::vector<long double>* v_series_acc, long double* v_stdev, std::vector<REGTYPE>* regtype, BAGTYPE* bagtype, MAT* objMAT, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::string* AnionName, std::vector<std::string>* CationNames, std::vector<std::string>* test_name, std::vector<long double>* test_anades_degree, long double* test_OM_degree, long double* test_CO_degree, bool* verbose, bool* post_if_exists, bool* AveOnly) {

	std::vector<int> test_anades; 
	std::vector<int> test_mamssk; 
	std::vector<std::string> find_mamssk;
	if (!*AveOnly) {
		find_mamssk.push_back("Min");
		find_mamssk.push_back("Ave");
		find_mamssk.push_back("Max");
		find_mamssk.push_back("StdDv");
		find_mamssk.push_back("Kurto");
	}
	else {
		find_mamssk.push_back("");
	}
	
	bool test_OM = false;
	int test_CO = -100;

	for (int itn = 0; itn < (signed)test_name->size(); itn++) {
		bool found = false;
		for (int iana = 0; iana < (signed)anadescriptor->size(); iana++) {
			if ((*anadescriptor)[iana].ASB != "U") {
				for (int imamssk = 0; imamssk < (signed)find_mamssk.size(); imamssk++) {
					if ((*test_name)[itn] == find_mamssk[imamssk] + (*anadescriptor)[iana].name) {
						test_anades.push_back(iana);
						if (!*AveOnly) {
							test_mamssk.push_back(imamssk);
						}
						else {
							test_mamssk.push_back(1);
						}
						found = true;
					}
				}
			}
			else {
				if ((*test_name)[itn] == (*anadescriptor)[iana].name) {
					test_anades.push_back(iana);
					test_mamssk.push_back(1);
					found = true;
				}
			}
			if (found) {
				break;
			}
		}
		if (!found && (*test_name)[itn] == *AnionName + "M") {
			test_OM = true;
			found = true;
		}
		if (!found) {
			for (int ic = 0; ic < (signed)CationNames->size(); ic++) {
				if ((*test_name)[itn] == (*CationNames)[ic] + *AnionName) {
					test_CO = ic;
					found = true;
					break;
				}
			}
		}
	}

	if (test_anades.size() > 0 || test_OM == true || test_CO != -100) {
		long double post_v = 0.0;
		calc_bag_chain(v_series, v_series_acc, v_stdev, regtype, bagtype, objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, test_anades_degree, test_OM_degree, test_CO_degree, verbose, post_if_exists, AveOnly, &post_v);
	}
	else {
		for (int ireg = 0; ireg < (signed)bagtype->which_regtype.size(); ireg++) {
			v_series->push_back(0.0);
			v_series_acc->push_back(0.0);
		}
	}
	
	std::vector<int>().swap(test_anades);
	std::vector<int>().swap(test_mamssk);
	std::vector<std::string>().swap(find_mamssk);

}

void test_robustness(long double* SStot, long double* SSres, long double* MSE, long double* RMSE, long double* MAE, long double* xbar, long double* ybar, long double* R2, long double* PearsonR, int* invalid_count, std::vector<long double>* y_regressed, std::vector<long double>* y_raw, bool* use_limit, long double* minv, long double* maxv, int* starting_index, int* jumping_jack_flash_index) {

	*SStot = 0.0;
	*SSres = 0.0;
	*MSE = 0.0;
	*RMSE = 0.0;
	*MAE = 0.0;
	*xbar = 0.0;
	*ybar = 0.0;
	*R2 = 0.0;
	*PearsonR = 0.0;
	*invalid_count = 0;

	int icount = 0;
	for (int i = *starting_index; i < (signed)y_raw->size(); i += *jumping_jack_flash_index) {
		if (!std::isnan((*y_raw)[i]) && !std::isinf((*y_raw)[i]) && !std::isinf(-(*y_raw)[i]) && !std::isnan((*y_regressed)[i]) && !std::isinf((*y_regressed)[i]) && !std::isinf(-(*y_regressed)[i])) {
			if (!*use_limit) {
				*ybar += (*y_raw)[i];
				*xbar += (*y_regressed)[i];
				icount++;
			}
			else {
				if ((*y_raw)[i] > *minv && (*y_raw)[i] < *maxv && (*y_regressed)[i] > *minv && (*y_regressed)[i] < *maxv) {
					*ybar += (*y_raw)[i];
					*xbar += (*y_regressed)[i];
					icount++;
				}
				else {
					(*invalid_count)++;
				}
			}
		}
		else {
			(*invalid_count)++;
		}
	}

	if (icount > 0) {

		*ybar /= (long double)icount;
		*xbar /= (long double)icount;
		long double PearsonR_numerator = 0.0;
		long double PearsonR_denominator_x = 0.0;
		long double PearsonR_denominator_y = 0.0;

		for (int i = *starting_index; i < (signed)y_raw->size(); i += *jumping_jack_flash_index) {
			if (!std::isnan((*y_raw)[i]) && !std::isinf((*y_raw)[i]) && !std::isinf(-(*y_raw)[i]) && !std::isnan((*y_regressed)[i]) && !std::isinf((*y_regressed)[i]) && !std::isinf(-(*y_regressed)[i])) {
				if (!*use_limit) {
					*SStot += powl((*y_raw)[i] - *ybar, 2.0);
					*SSres += powl((*y_raw)[i] - (*y_regressed)[i], 2.0);
					*MAE += fabsl((*y_raw)[i] - (*y_regressed)[i]);
					PearsonR_numerator += ((*y_regressed)[i] - *xbar) * ((*y_raw)[i] - *ybar);
					PearsonR_denominator_x += powl((*y_regressed)[i] - *xbar, 2.0);
					PearsonR_denominator_y += powl((*y_raw)[i] - *ybar, 2.0);
				}
				else {
					if ((*y_raw)[i] > *minv && (*y_raw)[i] < *maxv && (*y_regressed)[i] > *minv && (*y_regressed)[i] < *maxv) {
						*SStot += powl((*y_raw)[i] - *ybar, 2.0);
						*SSres += powl((*y_raw)[i] - (*y_regressed)[i], 2.0);
						*MAE += fabsl((*y_raw)[i] - (*y_regressed)[i]);
						PearsonR_numerator += ((*y_regressed)[i] - *xbar) * ((*y_raw)[i] - *ybar);
						PearsonR_denominator_x += powl((*y_regressed)[i] - *xbar, 2.0);
						PearsonR_denominator_y += powl((*y_raw)[i] - *ybar, 2.0);
					}
				}
			}
		}

		*MSE = *SSres / (long double)icount;
		*RMSE = sqrtl(*MSE);
		*MAE /= (long double)icount;
		*R2 = 1.0 - *SSres / *SStot;
		*PearsonR = PearsonR_numerator / sqrtl(PearsonR_denominator_x * PearsonR_denominator_y);

	}

}

void read_para() {

	readpara.open("basic_properties_GoodDesigner.txt");

	readpara >> sbuf >> AtomDatListFilename;
	readpara >> sbuf >> FaceFilename;
	readpara >> sbuf >> StructureListFilename >> UseStrName;
	readpara >> sbuf >> MaterialsFilename >> HeadExists;
	std::getline(readpara, sbuf);
	std::stringstream lcstr;
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		if (sbuf != "") {
			PostTagNames.push_back(sbuf);
		}
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		if (sbuf == "-") {
			break;
		}
		if (sbuf != "!") {
			RegModelFilename.push_back(sbuf);
		}
		else {
			FurtherMaterialsDesign = false;
		}
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf >> MapOutSpaceFilename;
	if (MapOutSpaceFilename != "-") {
		mapout = true;
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		lcstr >> ldbuf;
		TargetVReg.push_back(ldbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		lcstr >> ldbuf;
		TargetStdDvReg.push_back(ldbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		lcstr >> ldbuf;
		ZDistWeight.push_back(ldbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		lcstr >> sbuf;
		if (sbuf == "-") {
			bUseAve.push_back(false);
			UseAve.push_back(0.0);
		}
		else {
			bUseAve.push_back(true);
			std::stringstream llcstr;
			llcstr.str("");
			llcstr.clear();
			llcstr << sbuf;
			llcstr >> ldbuf;
			UseAve.push_back(ldbuf);
		}
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		lcstr >> sbuf;
		if (sbuf == "-") {
			bUseStdDv.push_back(false);
			UseStdDv.push_back(0.0);
		}
		else {
			bUseStdDv.push_back(true);
			std::stringstream llcstr;
			llcstr.str("");
			llcstr.clear();
			llcstr << sbuf;
			llcstr >> ldbuf;
			UseStdDv.push_back(ldbuf);
		}
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> ldbuf;
		if (ldbuf > 0) {
			ModulateContentFactor.push_back(ldbuf);
		}
		else {
			DopeFactor.push_back(-ldbuf);
		}
		if (lcstr.eof()) {
			break;
		}
	}
	
	readpara >> sbuf >> UseOffset;
	readpara >> sbuf >> MaxIteration;

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		NeighborSelectSpace.push_back(sbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readpara >> sbuf >> NeighborZCrit;
	if (NeighborZCrit <= 0.0) {
		NeighborCrit.push_back(NeighborZCrit);
		std::getline(readpara, sbuf);
		lcstr.str("");
		lcstr.clear();
		lcstr << sbuf;
		while (true) {
			lcstr >> ldbuf;
			NeighborCrit.push_back(ldbuf);
			if (lcstr.eof()) {
				break;
			}
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	if (sbuf.find("-") == std::string::npos) {
		while (true) {
			lcstr >> sbuf;
			if (sbuf == "Toxic" || sbuf == "toxic") {
				ExcludeNeighbors.push_back("Hg");
				ExcludeNeighbors.push_back("Pb");
				ExcludeNeighbors.push_back("Cd");
				ExcludeNeighbors.push_back("As");
				ExcludeNeighbors.push_back("Tl");
				ExcludeNeighbors.push_back("Be");
				ExcludeNeighbors.push_back("Cr");
				ExcludeNeighbors.push_back("Cs");
				ExcludeNeighbors.push_back("Po");
				ExcludeNeighbors.push_back("At");
				ExcludeNeighbors.push_back("Rn");
				ExcludeNeighbors.push_back("Fr");
				ExcludeNeighbors.push_back("Ra");
				ExcludeNeighbors.push_back("Ac");
				ExcludeNeighbors.push_back("Th");
				ExcludeNeighbors.push_back("Pa");
				ExcludeNeighbors.push_back("U");
				ExcludeNeighbors.push_back("Np");
				ExcludeNeighbors.push_back("Pu");
				ExcludeNeighbors.push_back("Am");
				ExcludeNeighbors.push_back("Cm");
				ExcludeNeighbors.push_back("Bk");
				ExcludeNeighbors.push_back("Cf");
				ExcludeNeighbors.push_back("Es");
				ExcludeNeighbors.push_back("Fm");
				ExcludeNeighbors.push_back("Md");
				ExcludeNeighbors.push_back("No");
				ExcludeNeighbors.push_back("Lr");
				ExcludeNeighbors.push_back("Rf");
				ExcludeNeighbors.push_back("Db");
				ExcludeNeighbors.push_back("Sg");
				ExcludeNeighbors.push_back("Bh");
				ExcludeNeighbors.push_back("Hs");
				ExcludeNeighbors.push_back("Mt");
				ExcludeNeighbors.push_back("Ds");
				ExcludeNeighbors.push_back("Rg");
				ExcludeNeighbors.push_back("Cn");
				ExcludeNeighbors.push_back("Nh");
				ExcludeNeighbors.push_back("Fl");
				ExcludeNeighbors.push_back("Mc");
				ExcludeNeighbors.push_back("Lv");
				ExcludeNeighbors.push_back("Ts");
				ExcludeNeighbors.push_back("Og");
			}
			else if (sbuf.find("rare") != std::string::npos || sbuf.find("Rare") != std::string::npos) {
				ExcludeNeighbors.push_back("Ce");
				ExcludeNeighbors.push_back("Pr");
				ExcludeNeighbors.push_back("Nd");
				ExcludeNeighbors.push_back("Pm");
				ExcludeNeighbors.push_back("Sm");
				ExcludeNeighbors.push_back("Eu");
				ExcludeNeighbors.push_back("Gd");
				ExcludeNeighbors.push_back("Tb");
				ExcludeNeighbors.push_back("Dy");
				ExcludeNeighbors.push_back("Ho");
				ExcludeNeighbors.push_back("Er");
				ExcludeNeighbors.push_back("Tm");
				ExcludeNeighbors.push_back("Yb");
				ExcludeNeighbors.push_back("Lu");
			}
			else if (sbuf.find("Non-quenched") != std::string::npos || sbuf.find("NQ") != std::string::npos) {
				ExcludeNeighbors.push_back("Mn");
				ExcludeNeighbors.push_back("Fe");
				ExcludeNeighbors.push_back("Co");
				ExcludeNeighbors.push_back("Ni");
				ExcludeNeighbors.push_back("Cu");
				ExcludeNeighbors.push_back("Ru");
				ExcludeNeighbors.push_back("Rh");
				ExcludeNeighbors.push_back("Pd");
				ExcludeNeighbors.push_back("Ce");
				ExcludeNeighbors.push_back("Pr");
				ExcludeNeighbors.push_back("Nd");
				ExcludeNeighbors.push_back("Pm");
				ExcludeNeighbors.push_back("Sm");
				ExcludeNeighbors.push_back("Eu");
				ExcludeNeighbors.push_back("Gd");
				ExcludeNeighbors.push_back("Tb");
				ExcludeNeighbors.push_back("Dy");
				ExcludeNeighbors.push_back("Ho");
				ExcludeNeighbors.push_back("Er");
				ExcludeNeighbors.push_back("Tm");
				ExcludeNeighbors.push_back("Yb");
				ExcludeNeighbors.push_back("Os");
				ExcludeNeighbors.push_back("Ir");
				ExcludeNeighbors.push_back("Pt");
				ExcludeNeighbors.push_back("Au");
			}
			else {
				ExcludeNeighbors.push_back(sbuf);
			}
			if (lcstr.eof()) {
				break;
			}
		}
	}

	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf;
	while (true) {
		lcstr >> sbuf;
		if (sbuf != "-") {
			if (sbuf.find("-") == std::string::npos) {
				DoNotFindNeighbors.push_back(sbuf);
				DopeOnly.push_back(false);
			}
			else {
				DoNotFindNeighbors.push_back(sbuf.substr(1));
				DopeOnly.push_back(true);
			}	
		}
		if (lcstr.eof()) {
			break;
		}
	}
	
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf >> AnionName;
	OM_count = 0;
	if (AnionName != "*") {
		if (!lcstr.eof()) {
			while (true) {
				lcstr >> sbuf;
				CationNames.push_back(sbuf);
				CO_count.push_back(0);
				if (lcstr.eof()) {
					break;
				}
			}
		}
	}
	else {
		AnionName = "-";
		if (!lcstr.eof()) {
			while (true) {
				ANIONS objANIONS;
				lcstr >> objANIONS.anionname >> objANIONS.anionvalence;
				if (objANIONS.anionname.find("limit") != std::string::npos || objANIONS.anionname.find("Limit") != std::string::npos || objANIONS.anionname.find("LIMIT") != std::string::npos) {
					modulation_limit = objANIONS.anionvalence;
				}
				else {
					Anions.push_back(objANIONS);
				}
				if (lcstr.eof()) {
					break;
				}
			}
		}
	}
	
	readpara >> sbuf >> AnionValence;
	readpara >> sbuf >> VInfinite;
	OM_min = VInfinite;
	OM_ave = 0.0;
	OM_max = -VInfinite;
	for (int i = 0; i < (signed)CationNames.size(); i++) {
		CO_min.push_back(VInfinite);
		CO_ave.push_back(0.0);
		CO_max.push_back(-VInfinite);
	}
	int iverbose;
	readpara >> sbuf >> iverbose;
	if (iverbose == 0) {
		verbose = false;
	}
	else {
		verbose = true;
	}

	readpara >> sbuf >> OutputFilename;
	std::getline(readpara, sbuf);
	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf;
	lcstr >> AnalysisOutputFilename;
	if (!lcstr.eof()) {
		lcstr >> AnalysisMesh;
	}

	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf;
	if (!lcstr.eof()) {
		while (true) {
			lcstr >> sbuf;
			DoNotAnalyze.push_back(sbuf);
			if (lcstr.eof()) {
				break;
			}
		}
	}

	std::getline(readpara, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	lcstr >> sbuf;
	lcstr >> NVariableLimit[1];
	if (NVariableLimit[1] < 1) {
		NVariableLimit[1] = 1;
	}
	if (NVariableLimit[1] == 1) {
		NVariableLimit[0] = 1;
	}
	else if (NVariableLimit[1] > 1) {
		NVariableLimit[0] = 2;
	}
	if (!lcstr.eof()) {
		lcstr >> NCavailles;
	}
	else {
		NCavailles = 1;
	}

	readpara >> sbuf >> SkipMinMaxKurto;
	
	readpara >> sbuf;
	int Ucount = 0;
	while (true) {

		ANADESCRIPTOR objANADESCRIPTOR;
		readpara >> objANADESCRIPTOR.name;
		if (objANADESCRIPTOR.name == "STOP") {
			break;
		}

		objANADESCRIPTOR.min = VInfinite;
		objANADESCRIPTOR.min_ave = VInfinite;
		objANADESCRIPTOR.ave = 0.0;
		objANADESCRIPTOR.max_ave = -VInfinite;
		objANADESCRIPTOR.max = -VInfinite;
		objANADESCRIPTOR.min_stdev = VInfinite;
		objANADESCRIPTOR.min_skew = VInfinite;
		objANADESCRIPTOR.min_kurto = VInfinite;
		objANADESCRIPTOR.ave_stdev = 0.0;
		objANADESCRIPTOR.ave_skew = 0.0;
		objANADESCRIPTOR.ave_kurto = 0.0;
		objANADESCRIPTOR.max_stdev = -VInfinite;
		objANADESCRIPTOR.max_skew = -VInfinite;
		objANADESCRIPTOR.max_kurto = -VInfinite;
		objANADESCRIPTOR.count = 0;
		objANADESCRIPTOR.latex = "";
		objANADESCRIPTOR.mathematica = "";
		
		bool found_latex_tag = false;
		bool found_mathematica_tag = false;
		readpara >> objANADESCRIPTOR.ASB;
		/* USER-DEFINED SUPERVISED FUNCS */
		if (objANADESCRIPTOR.ASB == "P") {
			readpara >> objANADESCRIPTOR.ides;
		}

		if (objANADESCRIPTOR.ASB == "B") {
			readpara >> objANADESCRIPTOR.postfix;
		}
		else if (objANADESCRIPTOR.ASB == "U") {
			Ucount++;
			objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 0;
			objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 0;
			while (true) {
				readpara >> sbuf;
				if (sbuf == "*" || sbuf == "/") {
					if (sbuf == "*") {
						objANADESCRIPTOR.UType_m0d1 = 0;
					}
					else {
						objANADESCRIPTOR.UType_m0d1 = 1;
					}
					break;
				}
				else {
					if (sbuf == "All") {
						objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 1;
					}
					else if (sbuf == "!") {
						objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 2;
					}
					else {
						if (isNumber(sbuf)) {
							objANADESCRIPTOR.UType_numerator_S0A1E2N3 = 3;
							std::stringstream ldbuf;
							ldbuf.str("");
							ldbuf.clear();
							ldbuf << sbuf;
							ldbuf >> objANADESCRIPTOR.UType_numerator_value;
						}
						else {
							objANADESCRIPTOR.UType_numerator_AtomKinds.push_back(sbuf);
						}
					}
				}
			}
			std::getline(readpara, sbuf);
			std::stringstream decompbuf;
			decompbuf.str("");
			decompbuf.clear();
			decompbuf << sbuf;
			while (true) {
				decompbuf >> sbuf;
				if (sbuf == "L=") {
					found_latex_tag = true;
					bool initial_latex = true;
					while (true) {
						decompbuf >> sbuf;
						if (sbuf == "M=") {
							found_mathematica_tag = true;
							bool initial_mathematica = true;
							while (true) {
								decompbuf >> sbuf;
								if (!initial_mathematica) {
									objANADESCRIPTOR.mathematica += " ";
								}
								objANADESCRIPTOR.mathematica += sbuf;
								initial_mathematica = false;
								if (decompbuf.eof()) {
									break;
								}
							}
						}
						if (found_mathematica_tag && decompbuf.eof()) {
							break;
						}
						if (!initial_latex) {
							objANADESCRIPTOR.latex += " ";
						}
						objANADESCRIPTOR.latex += sbuf;
						initial_latex = false;
						if (decompbuf.eof()) {
							break;
						}
					}
					break;
				}
				if (sbuf == "All") {
					objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 1;
				}
				else if (sbuf == "!") {
					objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 2;
				}
				else {
					if (isNumber(sbuf)) {
						objANADESCRIPTOR.UType_denominator_S0A1E2N3 = 3;
						std::stringstream ldbuf;
						ldbuf.str("");
						ldbuf.clear();
						ldbuf << sbuf;
						ldbuf >> objANADESCRIPTOR.UType_denominator_value;
					}
					else {
						objANADESCRIPTOR.UType_denominator_AtomKinds.push_back(sbuf);
					}
				}
				if (decompbuf.eof()) {
					break;
				}
			}
		}

		if (!found_latex_tag) {
			std::getline(readpara, sbuf);
			std::stringstream latexbuf;
			latexbuf.str("");
			latexbuf.clear();
			latexbuf << sbuf;
			latexbuf >> sbuf;
			if (sbuf == "L=") {
				found_latex_tag = true;
				bool initial_latex = true;
				while (true) {
					latexbuf >> sbuf;
					if (sbuf == "M=") {
						found_mathematica_tag = true;
						break;
					}
					if (!initial_latex) {
						objANADESCRIPTOR.latex += " ";
					}
					objANADESCRIPTOR.latex += sbuf;
					initial_latex = false;
					if (latexbuf.eof()) {
						break;
					}
				}
			}
			if (found_latex_tag && found_mathematica_tag) {
				bool initial_mathematica = true;
				while (true) {
					latexbuf >> sbuf;
					if (!initial_mathematica) {
						objANADESCRIPTOR.mathematica += " ";
					}
					objANADESCRIPTOR.mathematica += sbuf;
					initial_mathematica = false;
					if (latexbuf.eof()) {
						break;
					}
				}
			}
		}

		if (!found_latex_tag) {
			objANADESCRIPTOR.latex = "a" + std::to_string(anadescriptor.size() + 1);
		}
		if (!found_mathematica_tag) {
			objANADESCRIPTOR.mathematica = "a" + std::to_string( anadescriptor.size() + 1);
		}

		if (!mapout) {
			anadescriptor.push_back(objANADESCRIPTOR);
		}
		else if (objANADESCRIPTOR.ASB == "A" || objANADESCRIPTOR.ASB == "U") {
			anadescriptor.push_back(objANADESCRIPTOR);
		}

		if (readpara.eof()) {
			break;
		}

	}

	xbox = new long double *[anadescriptor.size()];
	for (int i = 0; i < (signed)anadescriptor.size(); i++) {
		xbox[i] = new long double[6];
		for (int j = 0; j < 6; j++) {
			xbox[i][j] = 0.0;
		}
	}
	xbox_OM = 0.0;
	xbox_CO = new long double[CationNames.size()];
	for (int i = 0; i < (signed)CationNames.size(); i++) {
		xbox_CO[i] = 0.0;
	}
	xbox_user = new long double[Ucount];
	for (int i = 0; i < Ucount; i++) {
		xbox_user[i] = 0.0;
	}

	readpara.close();

}

void read_atomdat(std::vector<ATOMDAT>* atomdat, std::vector<DESCRIPTOR> * descriptor, std::vector<ANADESCRIPTOR> * anadescriptor, std::string * AtomDatListFilename, std::vector<std::string>* DoNotAnalyze) {

	readatomdat.open(*AtomDatListFilename);

	std::getline(readatomdat, sbuf);
	std::stringstream lcstr;
	lcstr << sbuf;
	lcstr >> sbuf;
	while (true) {
		DESCRIPTOR objDESCRIPTOR;
		lcstr >> objDESCRIPTOR.name;
		if (objDESCRIPTOR.name == "Valence") {
			CValence = (signed)descriptor->size();
		}
		else if (objDESCRIPTOR.name.find("Valence") != std::string::npos) {
			AValence = (signed)descriptor->size();
		}
		descriptor->push_back(objDESCRIPTOR);
		if (lcstr.eof()) {
			break;
		}
	}

	for (int ii = 0; ii < (signed)anadescriptor->size(); ii++) {
		if ((*anadescriptor)[ii].ASB == "A") {
			for (int i = 0; i < (signed)descriptor->size(); i++) {
				if ((*descriptor)[i].name == (*anadescriptor)[ii].name) {
					(*anadescriptor)[ii].ides = i;
					break;
				}
			}
		}
		if ((*anadescriptor)[ii].ASB == "P") {
			(*anadescriptor)[ii].ASB = "A";
		}
	}
	
	while (true) {
		ATOMDAT objATOMDAT;
		readatomdat >> objATOMDAT.atomname;
		if (objATOMDAT.atomname == "STOP") {
			break;
		}
		for(int i = 0; i < (signed)descriptor->size(); i++){
			readatomdat >> ldbuf;
			objATOMDAT.v.push_back(ldbuf);
		}
		atomdat->push_back(objATOMDAT);
		if (readatomdat.eof()) {
			break;
		}
	}

	readatomdat.close();

	for (int j = 0; j < (signed)descriptor->size(); j++) {
		(*descriptor)[j].min = 10000000000000000000000.0;
		(*descriptor)[j].ave = 0.0;
		(*descriptor)[j].max = -10000000000000000000000.0;
		for (int i = 0; i < (signed)atomdat->size(); i++) {
			if ((*atomdat)[i].v[j] < (*descriptor)[j].min) {
				(*descriptor)[j].min = (*atomdat)[i].v[j];
			}
			(*descriptor)[j].ave += (*atomdat)[i].v[j];
			if ((*atomdat)[i].v[j] > (*descriptor)[j].max) {
				(*descriptor)[j].max = (*atomdat)[i].v[j];
			}
		}
		(*descriptor)[j].ave /= (long double)(signed)atomdat->size();
		(*descriptor)[j].stddv = 0.0;
		for (int i = 0; i < (signed)atomdat->size(); i++) {
			(*descriptor)[j].stddv += powl((*atomdat)[i].v[j] - (*descriptor)[j].ave, 2.0);
		}
		(*descriptor)[j].stddv /= (long double)(signed)atomdat->size();
		(*descriptor)[j].stddv = sqrtl((*descriptor)[j].stddv);
		for (int i = 0; i < (signed)atomdat->size(); i++) {
			(*atomdat)[i].zv.push_back(((*atomdat)[i].v[j] - (*descriptor)[j].ave) / (*descriptor)[j].stddv);
		}
	}

	if (pr_id == 0) {
		writeoutput << "Atom\t";
		writeoutput.flush();
		for (int j = 0; j < (signed)descriptor->size(); j++) {
			writeoutput << "z" << (*descriptor)[j].name;
			writeoutput.flush();
			if (j != (signed)descriptor->size() - 1) {
				writeoutput << "\t";
				writeoutput.flush();
			}
			else {
				writeoutput << "\n";
				writeoutput.flush();
			}
		}
		for (int i = 0; i < (signed)atomdat->size(); i++) {
			writeoutput << (*atomdat)[i].atomname << "\t";
			writeoutput.flush();
			for (int j = 0; j < (signed)descriptor->size(); j++) {
				writeoutput << (*atomdat)[i].zv[j];
				writeoutput.flush();
				if (j != (signed)descriptor->size() - 1) {
					writeoutput << "\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\n";
					writeoutput.flush();
				}
			}
		}

		writeoutput << "=================================================================================\n";
		writeoutput.flush();
	}
	
}

void read_face(std::vector<FACE>* face, std::vector<std::string>* FaceFeatureName, std::string* FaceFilename) {

	readface.open(*FaceFilename);

	std::getline(readface, sbuf);
	std::stringstream lcstream;
	lcstream << sbuf;
	lcstream >> sbuf >> sbuf;
	while (true) {
		lcstream >> sbuf;
		FaceFeatureName->push_back(sbuf);
		if (lcstream.eof()) {
			break;
		}
	}

	while (true) {
		std::getline(readface, sbuf);
		if (readface.eof() || sbuf.find("STOP") != std::string::npos) {
			break;
		}
		else {
			lcstream.str("");
			lcstream.clear();
			lcstream << sbuf;
			FACE objFACE;
			lcstream >> objFACE.AtomName;
			lcstream >> sbuf;
			if (sbuf == "-") {
				objFACE.face = false;
			}
			else {
				objFACE.face = true;
			}
			long double ldbuf;
			for (int i = 0; i < (signed)FaceFeatureName->size(); i++) {
				lcstream >> ldbuf;
				objFACE.vFeature.push_back(ldbuf);
			}
			face->push_back(objFACE);
		}
	}

	readface.close();

}

void pick_neighbor(std::vector<ATOMDAT>* atomdat, std::vector<std::string>* NeighborSelectSpace, long double* NeighborZCrit, std::vector<long double>* NeighborCrit, std::vector<std::string>* ExcludeNeighbors, std::vector<std::string>* DoNotFindNeighbors, std::vector<bool>* DopeOnly, std::vector<DESCRIPTOR>* descriptor) {

	std::vector<int> ides;
	for (int i = 0; i < (signed)NeighborSelectSpace->size(); i++) {
		for (int j = 0; j < (signed)descriptor->size(); j++) {
			if ((*NeighborSelectSpace)[i] == (*descriptor)[j].name) {
				ides.push_back(j);
				break;
			}
		}
	}

	for (int ii = 0; ii < (signed)atomdat->size(); ii++) {
		readnb.open((*atomdat)[ii].atomname + ".nb.txt");
		if (readnb.is_open()) {
			while (true) {
				readnb >> sbuf;
				for (int jj = 0; jj < (signed)atomdat->size(); jj++) {
					if ((*atomdat)[jj].atomname != (*atomdat)[ii].atomname && (*atomdat)[jj].atomname == sbuf) {
						NEIGHBOR objNEIGHBOR;
						objNEIGHBOR.atomname = (*atomdat)[jj].atomname;
						objNEIGHBOR.iatom = jj;
						(*atomdat)[ii].neighbor.push_back(objNEIGHBOR);
						if (pr_id == 0) {
							std::cout << "Overwrite NEIGHBOR: " << (*atomdat)[ii].atomname << "\t...NEIGHBOR...\t" << (*atomdat)[jj].atomname << "\n";
							std::cout.flush();
						}
						break;
					}
				}
				if (readnb.eof()) {
					break;
				}
			}
			readnb.close();
			for (int jj = 0; jj < (signed)DoNotFindNeighbors->size(); jj++) {
				if ((*atomdat)[ii].atomname == (*DoNotFindNeighbors)[jj]) {
					if ((*DopeOnly)[jj]) {
						(*atomdat)[ii].dopeonly = true;
					}
					break;
				}
			}
		}
		else {
			readnb.close();
			bool needtofind = true;
			for (int jj = 0; jj < (signed)DoNotFindNeighbors->size(); jj++) {
				if ((*atomdat)[ii].atomname == (*DoNotFindNeighbors)[jj]) {
					needtofind = false;
					if ((*DopeOnly)[jj]) {
						needtofind = true;
						(*atomdat)[ii].dopeonly = true;
					}
					break;
				}
			}
			if (needtofind) {
				for (int jj = 0; jj < (signed)atomdat->size(); jj++) {
					bool exclude = false;
					for (int kk = 0; kk < (signed)ExcludeNeighbors->size(); kk++) {
						if ((*atomdat)[jj].atomname == (*ExcludeNeighbors)[kk]) {
							exclude = true;
							break;
						}
					}
					if (ii != jj && !exclude) {
						if (*NeighborZCrit > 0.0) {
							long double zdist = 0.0;
							for (int i = 0; i < (signed)ides.size(); i++) {
								zdist += powl((*atomdat)[jj].zv[ides[i]] - (*atomdat)[ii].zv[ides[i]], 2.0);
							}
							zdist = sqrtl(zdist);
							if (zdist <= *NeighborZCrit) {
								NEIGHBOR objNEIGHBOR;
								objNEIGHBOR.atomname = (*atomdat)[jj].atomname;
								objNEIGHBOR.iatom = jj;
								(*atomdat)[ii].neighbor.push_back(objNEIGHBOR);
							}
						}
						else {
							bool include_neighbor = true;
							for (int i = 0; i < (signed)ides.size(); i++) {
								if (fabsl((*atomdat)[jj].v[ides[i]] - (*atomdat)[ii].v[ides[i]]) > fabsl((*NeighborCrit)[i])) {
									include_neighbor = false;
									break;
								}
							}
							if (include_neighbor) {
								NEIGHBOR objNEIGHBOR;
								objNEIGHBOR.atomname = (*atomdat)[jj].atomname;
								objNEIGHBOR.iatom = jj;
								(*atomdat)[ii].neighbor.push_back(objNEIGHBOR);
							}
						}
					}
				}
			}
		}

	}
	
	if (pr_id == 0) {
		for (int ii = 0; ii < (signed)atomdat->size(); ii++) {
			writeoutput << (*atomdat)[ii].atomname << "\t...NEIGHBOR...";
			writeoutput.flush();
			if ((*atomdat)[ii].neighbor.size() > 0) {
				writeoutput << "\t";
				writeoutput.flush();
				for (int jj = 0; jj < (signed)(*atomdat)[ii].neighbor.size(); jj++) {
					writeoutput << (*atomdat)[ii].neighbor[jj].atomname;
					writeoutput.flush();
					if (jj != (signed)(*atomdat)[ii].neighbor.size() - 1) {
						writeoutput << "\t";
						writeoutput.flush();
					}
					else {
						writeoutput << "\n";
						writeoutput.flush();
					}
				}
			}
			else {
				writeoutput << "\n";
				writeoutput.flush();
			}
		}

		writeoutput << "=================================================================================\n";
		writeoutput.flush();
	}
	
	std::vector<int>().swap(ides);
	
}

void read_strdat(std::vector<STRDAT> * strdat, std::vector<ANADESCRIPTOR>* anadescriptor, std::vector<std::string> * sdescriptorname, std::vector<std::string>* sitename, std::string * StructureListFilename, std::vector<DESCRIPTOR>* descriptor, std::vector<ATOMDAT>* atomdat, bool *usestrname) {

	readstrdat.open(*StructureListFilename);

	readstrdat >> sbuf;
	std::getline(readstrdat, sbuf);
	std::stringstream lcstr;
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		sdescriptorname->push_back(sbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	readstrdat >> sbuf;
	std::getline(readstrdat, sbuf);
	lcstr.str("");
	lcstr.clear();
	lcstr << sbuf;
	while (true) {
		lcstr >> sbuf;
		sitename->push_back(sbuf);
		if (lcstr.eof()) {
			break;
		}
	}

	std::getline(readstrdat, sbuf);

	while (true) {
		readstrdat >> sbuf;
		if (sbuf == "STOP" || readstrdat.eof()) {
			break;
		}
		else {
			STRDAT objSTRDAT;
			if (*usestrname) {
				objSTRDAT.strname = sbuf;
			}
			else {
				objSTRDAT.strname = "None";
			}
			for (int i = 0; i < (signed)sdescriptorname->size(); i++) {
				SDESCRIPTOR objSDESCRIPTOR;
				objSDESCRIPTOR.name = (*sdescriptorname)[i];
				for (int j = 0; j < (signed)sitename->size(); j++) {
					readstrdat >> sbuf;
					if (sbuf != "X") {
						objSDESCRIPTOR.sitename.push_back((*sitename)[j]);
						objSDESCRIPTOR.strv.push_back(sbuf);
						std::stringstream lcstr;
						lcstr << sbuf;
						lcstr >> ldbuf;
						objSDESCRIPTOR.v.push_back(ldbuf);
					}
				}
				objSTRDAT.sdescriptor.push_back(objSDESCRIPTOR);
			}
			strdat->push_back(objSTRDAT);
		}
		if (readstrdat.eof()) {
			break;
		}
	}

	for (int i = 0; i < (signed)strdat->size(); i++) {
		for (int j = 0; j < (signed)(*strdat)[i].sdescriptor[0].sitename.size(); j++) {
			(*strdat)[i].sitename.push_back((*strdat)[i].sdescriptor[0].sitename[j]);
		}
	}

	readstrdat.close();

	for (int ii = 0; ii < (signed)anadescriptor->size(); ii++) {
		if ((*anadescriptor)[ii].ASB == "S") {
			for (int i = 0; i < (signed)sdescriptorname->size(); i++) {
				if ((*sdescriptorname)[i] == (*anadescriptor)[ii].name) {
					(*anadescriptor)[ii].ides = i;
					break;
				}
			}
		}
		else if ((*anadescriptor)[ii].ASB == "B") {
			(*anadescriptor)[ii].idesB = new int* [strdat->size()];
			for (int i = 0; i < (signed)strdat->size(); i++) {
				(*anadescriptor)[ii].idesB[i] = new int[sitename->size()];
				for (int j = 0; j < (signed)sitename->size(); j++) {
					(*anadescriptor)[ii].idesB[i][j] = -1;
					std::string target_name = (*anadescriptor)[ii].name;
					int target_ipostfix = -1;
					for (int k = 0; k < (signed)sdescriptorname->size(); k++) {
						if ((*sdescriptorname)[k] == (*anadescriptor)[ii].postfix) {
							target_ipostfix = k;
							break;
						}
					}
					if (target_ipostfix != -1 && j < (*strdat)[i].sdescriptor[target_ipostfix].strv.size()) {
						target_name += (*strdat)[i].sdescriptor[target_ipostfix].strv[j];
					}
					for (int k = 0; k < (signed)descriptor->size(); k++) {
						if ((*descriptor)[k].name == target_name) {
							(*anadescriptor)[ii].idesB[i][j] = k;
							break;
						}
					}
				}
			}
		}
	}

}

void find_nonlinear(std::string * functype, std::string * termline) {

	*functype = "";

	if ((termline->find("logit") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("sigmoid") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "logit";
	}
	else if ((termline->find("sigmoid") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("logit") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "sigmoid";
	}
	else if ((termline->find("probit") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cdf") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "probit";
	}
	else if ((termline->find("cdf") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("probit") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cdf";
	}
	else if ((termline->find("cloglog") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cexpexp") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cloglog";
	}
	else if ((termline->find("cexpexp") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cloglog") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cexpexp";
	}
	else if ((termline->find("cauchit") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("tangenttype") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "cauchit";
	}
	else if ((termline->find("tangenttype") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("cauchit") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "tangenttype";
	}
	else if ((termline->find("nloglog") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("expexp") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "nloglog";
	}
	else if ((termline->find("expexp") != std::string::npos && termline->find("inverse") == std::string::npos) || (termline->find("nloglog") != std::string::npos && termline->find("inverse") != std::string::npos)) {
		*functype = "expexp";
	}

}

void analyze_term(long double* coeff, long double* power, std::vector<std::string>* term) {

	if (term->size() == 0) {
		*coeff = 0.0;
		*power = 0.0;
	}
	else if (term->size() == 1) {
		std::stringstream lcstr;
		lcstr << (*term)[0];
		lcstr >> ldbuf;
		*coeff = ldbuf;
		*power = 0.0;
	}
	else {
		if ((*term)[0] == (*term)[1]) {
			(*term)[0] = "1.0";
		}
		std::stringstream lcstr;
		lcstr << (*term)[0];
		lcstr >> ldbuf;
		*coeff = ldbuf;
		if ((*term)[term->size() - 1].find("^") == std::string::npos) {
			*power = 1.0;
		}
		else {
			sbuf = (*term)[term->size() - 1].substr((*term)[term->size() - 1].find("^") + 1);
			lcstr.str("");
			lcstr.clear();
			lcstr << sbuf;
			lcstr >> ldbuf;
			*power = ldbuf;
		}
	}

}

void identify_anades(int* ianades, int* imamssk, int* ilogln, int* ipow10exp, int* ierf, std::string * term, std::vector<ANADESCRIPTOR> * anadescriptor) {

	*ianades = -1; // OM
	size_t max_length = 0;
	for (int i = 0; i < (signed)anadescriptor->size(); i++) {
		if((term->find("^") != std::string::npos && term->find((*anadescriptor)[i].name + "^") != std::string::npos) || (term->find("^") == std::string::npos && term->find((*anadescriptor)[i].name) != std::string::npos)) {
			if ((*anadescriptor)[i].name.length() > max_length) {
				*ianades = i;
				max_length = (*anadescriptor)[i].name.length();
			}
		}
	}

	*ipow10exp = -1;
	if (term->find("Min") != std::string::npos) {
		*imamssk = 0;
		if (term->find("10Min") != std::string::npos && term->find("log10Min") == std::string::npos && term->find("Log10Min") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expMin") != std::string::npos || term->find("ExpMin") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Ave") != std::string::npos) {
		*imamssk = 1;
		if (term->find("10Ave") != std::string::npos && term->find("log10Ave") == std::string::npos && term->find("Log10Ave") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expAve") != std::string::npos || term->find("ExpAve") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Max") != std::string::npos) {
		*imamssk = 2;
		if (term->find("10Max") != std::string::npos && term->find("log10Max") == std::string::npos && term->find("Log10Max") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expMax") != std::string::npos || term->find("ExpMax") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("StdDv") != std::string::npos) {
		*imamssk = 3;
		if (term->find("10StdDv") != std::string::npos && term->find("log10StdDv") == std::string::npos && term->find("Log10StdDv") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expStdDv") != std::string::npos || term->find("ExpStdDv") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Skew") != std::string::npos) {
		*imamssk = 4;
		if (term->find("10Skew") != std::string::npos && term->find("log10Skew") == std::string::npos && term->find("Log10Skew") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expSkew") != std::string::npos || term->find("ExpSkew") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else if (term->find("Kurto") != std::string::npos) {
		*imamssk = 5;
		if (term->find("10Kurto") != std::string::npos && term->find("log10Kurto") == std::string::npos && term->find("Log10Kurto") == std::string::npos) {
			*ipow10exp = 1;
		}
		else if (term->find("expKurto") != std::string::npos || term->find("ExpKurto") != std::string::npos) {
			*ipow10exp = 2;
		}
	}
	else {
		*imamssk = -1; //intercept
	}

	if (term->find("Log10") != std::string::npos || term->find("log10") != std::string::npos) {
		*ilogln = 1;
	}
	else if (term->find("Log") != std::string::npos || term->find("log") != std::string::npos || term->find("Ln") != std::string::npos || term->find("ln") != std::string::npos) {
		*ilogln = 2;
	}
	else {
		*ilogln = -1;
	}

	*ierf = -1000;
	if (term->find("erf") != std::string::npos || term->find("Erf") != std::string::npos) {
		std::string terminology = "erf";
		if (term->find("Erf") != std::string::npos) {
			terminology = "Erf";
		}
		*ierf = 0;
		for (int ierf_find = -9; ierf_find <= 9; ierf_find++) {
			if (term->find(terminology + std::to_string(ierf_find)) != std::string::npos) {
				*ierf = ierf_find;
				break;
			}
		}
	}

}

bool identify_anades_chain(int* itype, int* ianades, int* imamssk, std::string* latex, std::string* mathematica, std::string* CationName, std::string* term, std::vector<ANADESCRIPTOR>* anadescriptor, std::string* AnionName, std::vector<std::string>* CationNames, bool* AveOnly) {

	bool found = false;

	if (!*AveOnly) {

		if (*term == *AnionName + "M") {
			*imamssk = -1;
			*ianades = -1;
			*itype = -1;
			*latex = *AnionName + "_{M}";
			*mathematica = "Subscript[" + *AnionName + ", M]";
			found = true;
		}
		else {
			for (int i = 0; i < (signed)anadescriptor->size(); i++) {
				if (*term == "Min" + (*anadescriptor)[i].name) {
					*imamssk = 0;
					*ianades = i;
					*itype = 0;
					*latex = " m \\left( " + (*anadescriptor)[i].latex + " \\right) ";
					*mathematica = "Subscript[" + (*anadescriptor)[i].mathematica + ", min]";
					found = true;
				}
				else if (*term == "Ave" + (*anadescriptor)[i].name) {
					*imamssk = 1;
					*ianades = i;
					*itype = 0;
					*latex = " \\langle " + (*anadescriptor)[i].latex + " \\rangle ";
					*mathematica = "Subscript[" + (*anadescriptor)[i].mathematica + ", ave]";
					found = true;
				}
				else if (*term == "Max" + (*anadescriptor)[i].name) {
					*imamssk = 2;
					*ianades = i;
					*itype = 0;
					*latex = " M \\left( " + (*anadescriptor)[i].latex + " \\right) ";
					*mathematica = "Subscript[" + (*anadescriptor)[i].mathematica + ", max]";
					found = true;
				}
				else if (*term == "StdDv" + (*anadescriptor)[i].name) {
					*imamssk = 3;
					*ianades = i;
					*itype = 0;
					*latex = " \\sigma \\left( " + (*anadescriptor)[i].latex + " \\right) ";
					*mathematica = "Subscript[" + (*anadescriptor)[i].mathematica + ", stdev]";
					found = true;
				}
				else if (*term == "Skew" + (*anadescriptor)[i].name) {
					*imamssk = 4;
					*ianades = i;
					*itype = 0;
					*latex = " r \\left( " + (*anadescriptor)[i].latex + " \\right) ";
					*mathematica = "Subscript[" + (*anadescriptor)[i].mathematica + ", skew]";
					found = true;
				}
				else if (*term == "Kurto" + (*anadescriptor)[i].name) {
					*imamssk = 5;
					*ianades = i;
					*itype = 0;
					*latex = " K \\left( " + (*anadescriptor)[i].latex + " \\right) ";
					*mathematica = "Subscript[" + (*anadescriptor)[i].mathematica + ", kurto]";
					found = true;
				}
				else if (*term == (*anadescriptor)[i].name && (*anadescriptor)[i].ASB == "U") {
					*imamssk = -3;
					*ianades = i;
					*itype = -3;
					*latex = (*anadescriptor)[i].latex;
					*mathematica = (*anadescriptor)[i].mathematica;
					found = true;
				}
				if (found) {
					break;
				}
			}
		}

		if (!found) {
			for (int i = 0; i < CationNames->size(); i++) {
				if (*term == (*CationNames)[i] + *AnionName) {
					*imamssk = -2;
					*ianades = i;
					*itype = -2;
					*CationName = (*CationNames)[i];
					*latex = *CationName + "_{" + *AnionName + "}";
					*mathematica = "Subscript[" + *CationName + ", " + *AnionName + "]";
					found = true;
					break;
				}
			}
		}
	
	}
	else {

		for (int i = 0; i < (signed)anadescriptor->size(); i++) {
			if (*term == (*anadescriptor)[i].name) {
				*imamssk = 1;
				*ianades = i;
				*itype = 0;
				*latex = (*anadescriptor)[i].latex;
				*mathematica = (*anadescriptor)[i].mathematica;
				found = true;
			}
		}
	
	}
	
	return found;

}

void read_reg(std::vector<FUNCTYPE> * functype, std::string * RegFilename, std::vector<ANADESCRIPTOR>* anadescriptor) {

	readreg.open(*RegFilename);

	if (readreg.is_open()) {
		FUNCTYPE objFUNCTYPE;
		objFUNCTYPE.added = false;
		readreg >> sbuf >> objFUNCTYPE.name;
		std::getline(readreg, sbuf);
		std::getline(readreg, sbuf);
		std::vector<std::string> slinebuf;
		while (true) {
			std::getline(readreg, sbuf);
			if (sbuf.find("=====") == std::string::npos) {
				slinebuf.push_back(sbuf);
			}
			else {
				break;
			}
		}
		if (slinebuf.size() >= 2) {
			objFUNCTYPE.sum_only = false;
			find_nonlinear(&objFUNCTYPE.nonlinear_funcname, &slinebuf[1]);
			for (int iline = 0; iline < 2; iline++) {
				std::stringstream lcstr;
				lcstr << slinebuf[iline];
				while (true) {
					lcstr >> sbuf;
					if (sbuf == "=") {
						break;
					}
				}
				while (true) {
					std::vector<std::string> term;
					while (true) {
						lcstr >> sbuf;
						if (sbuf == "+") {
							break;
						}
						else {
							if (sbuf != "*") {
								term.push_back(sbuf);
							}
						}
						if (lcstr.eof()) {
							break;
						}
					}
					long double coeff, power;
					analyze_term(&coeff, &power, &term);
					if (iline == 0) {
						objFUNCTYPE.coeff_0.push_back(coeff);
						objFUNCTYPE.power_0.push_back(power);
					}
					else {
						objFUNCTYPE.coeff_1.push_back(coeff);
						objFUNCTYPE.power_1.push_back(power);
					}
					std::vector<std::string>().swap(term);
					if (lcstr.eof()) {
						break;
					}
				}
			}
		}
		else {
			objFUNCTYPE.sum_only = true;
		}
		std::vector<std::string>().swap(slinebuf);
		std::getline(readreg, sbuf);
		bool addfeatures_exist = false;
		while (true) {
			std::vector<std::string> term;
			std::getline(readreg, sbuf);
			if (sbuf.find("STOP") != std::string::npos || sbuf.find("*****") != std::string::npos || readreg.eof()) {
				if (sbuf.find("*****") != std::string::npos) {
					addfeatures_exist = true;
				}
				break;
			}
			std::stringstream lcstr;
			lcstr << sbuf;
			lcstr >> sbuf;
			if (sbuf != "intercept") {
				term.push_back(sbuf);
				term.push_back(sbuf);
			}
			lcstr >> sbuf;
			if (term.size() != 0) {
				term[0] = sbuf;
			}
			else {
				term.push_back(sbuf);
			}
			long double coeff, power;
			analyze_term(&coeff, &power, &term);
			objFUNCTYPE.coeff_2.push_back(coeff);
			objFUNCTYPE.power_2.push_back(power);
			if (term.size() == 1) {
				objFUNCTYPE.ianades_2.push_back(-1);
				objFUNCTYPE.imamssk_2.push_back(-1);
				objFUNCTYPE.ilogln_2.push_back(-1);
				objFUNCTYPE.ipow10exp_2.push_back(-1);
				objFUNCTYPE.ierf_2.push_back(-1000);
			}
			else {
				int ianades, imamssk, ilogln, ipow10exp, ierf;
				identify_anades(&ianades, &imamssk, &ilogln, &ipow10exp, &ierf, &term[1], anadescriptor);
				objFUNCTYPE.ianades_2.push_back(ianades);
				objFUNCTYPE.imamssk_2.push_back(imamssk);
				objFUNCTYPE.ilogln_2.push_back(ilogln);
				objFUNCTYPE.ipow10exp_2.push_back(ipow10exp);
				objFUNCTYPE.ierf_2.push_back(ierf);
			}
			std::vector<std::string>().swap(term);
		}
		if (addfeatures_exist) {
			objFUNCTYPE.added = true;
			std::getline(readreg, sbuf);
			while (true) {
				std::vector<std::string> term;
				std::getline(readreg, sbuf);
				if (sbuf.find("STOP") != std::string::npos || readreg.eof()) {
					break;
				}
				std::stringstream lcstr;
				lcstr << sbuf;
				lcstr >> sbuf;
				if (sbuf != "intercept") {
					term.push_back(sbuf);
					term.push_back(sbuf);
				}
				lcstr >> sbuf;
				if (term.size() != 0) {
					term[0] = sbuf;
				}
				else {
					term.push_back(sbuf);
				}
				long double coeff, power;
				analyze_term(&coeff, &power, &term);
				objFUNCTYPE.coeff_3.push_back(coeff);
				objFUNCTYPE.power_3.push_back(power);
				if (term.size() == 1) {
					objFUNCTYPE.ianades_3.push_back(-1);
					objFUNCTYPE.imamssk_3.push_back(-1);
					objFUNCTYPE.ilogln_3.push_back(-1);
					objFUNCTYPE.ipow10exp_3.push_back(-1);
					objFUNCTYPE.ierf_3.push_back(-1000);
				}
				else {
					int ianades, imamssk, ilogln, ipow10exp, ierf;
					identify_anades(&ianades, &imamssk, &ilogln, &ipow10exp, &ierf, &term[1], anadescriptor);
					objFUNCTYPE.ianades_3.push_back(ianades);
					objFUNCTYPE.imamssk_3.push_back(imamssk);
					objFUNCTYPE.ilogln_3.push_back(ilogln);
					objFUNCTYPE.ipow10exp_3.push_back(ipow10exp);
					objFUNCTYPE.ierf_3.push_back(ierf);
				}
				std::vector<std::string>().swap(term);
			}
		}
		functype->push_back(objFUNCTYPE);
	}

	readreg.close();

}

REGTYPE read_reg_chain(std::string* RegFilename, std::vector<ANADESCRIPTOR>* anadescriptor, std::string* AnionName, std::vector<std::string>* CationNames, bool* verbose, bool* AveOnly) {

	readregchain.open(*RegFilename);
	REGTYPE objREGTYPE;
	std::stringstream localstr;

	if (readregchain.is_open()) {

		readregchain >> sbuf >> objREGTYPE.name;
		readregchain >> sbuf >> objREGTYPE.beta_target;
		
		std::string bpregressor_filename;
		readregchain >> sbuf >> bpregressor_filename;
		readbpregressor.open(bpregressor_filename);
		if (readbpregressor.is_open()) {
			while (true) {
				readbpregressor >> sbuf;
				if (sbuf == "Functypes") {
					break;
				}
			}
			while (true) {
				std::getline(readbpregressor, sbuf);
				if (sbuf.find("BetaRegression") == std::string::npos) {
					ORTHFUNCTYPE objORTHFUNCTYPE;
					localstr.str("");
					localstr.clear();
					localstr << sbuf;
					localstr >> sbuf >> sbuf;
					std::stringstream getname(sbuf);
					std::string seg;
					while (std::getline(getname, seg, '^')) {
						objORTHFUNCTYPE.userdefined_name.push_back(seg);
					}
					localstr >> objORTHFUNCTYPE.FName;
					int NCoeff = 0;
					if (objORTHFUNCTYPE.FName == "pow" || objORTHFUNCTYPE.FName == "abs") {
						NCoeff = 2;
					}
					else if (objORTHFUNCTYPE.FName == "log" || objORTHFUNCTYPE.FName == "erf" || objORTHFUNCTYPE.FName == "erfc") {
						NCoeff = 3;
					}
					else if (objORTHFUNCTYPE.FName == "sin" || objORTHFUNCTYPE.FName == "cos" || objORTHFUNCTYPE.FName == "tan" || objORTHFUNCTYPE.FName == "sinh" || objORTHFUNCTYPE.FName == "cosh" || objORTHFUNCTYPE.FName == "tanh") {
						NCoeff = 4;
					}
					else if (objORTHFUNCTYPE.FName == "exp") {
						NCoeff = 5;
					}
					for (int j = 0; j < NCoeff; j++) {
						long double ldbuf;
						localstr >> sbuf >> ldbuf;
						objORTHFUNCTYPE.coeff.push_back(ldbuf);
					}
					objREGTYPE.orthfunctype.push_back(objORTHFUNCTYPE);
				}
				else {
					break;
				}
			}
			localstr.str("");
			localstr.clear();
			localstr << sbuf;
			localstr >> sbuf;
			objREGTYPE.ibeta_target = 0;
			if (objREGTYPE.beta_target != "none") {
				objREGTYPE.ibeta_target = 1;
				while (true) {
					localstr >> sbuf;
					if (sbuf == objREGTYPE.beta_target || localstr.eof()) {
						break;
					}
					else {
						objREGTYPE.ibeta_target++;
					}
				}
			}
		}
		readbpregressor.close();

		readregchain >> sbuf;
		std::getline(readregchain, sbuf);
		while (true) {
			std::getline(readregchain, sbuf);
			if (sbuf == "STOP" || readregchain.eof()) {
				break;
			}
			else {
				FUNCCHAIN objFUNCCHAIN;
				localstr.str("");
				localstr.clear();
				localstr << sbuf;
				localstr >> objFUNCCHAIN.filename;
				objREGTYPE.funcchain.push_back(objFUNCCHAIN);
			}
		}

		if (!readregchain.eof()) {
			while (true) {
				std::getline(readregchain, sbuf);
				if (sbuf.find("POST") != std::string::npos || sbuf.find("post") != std::string::npos || sbuf.find("Post") != std::string::npos) {
					FUNC objPFUNC;
					localstr.str("");
					localstr.clear();
					localstr << sbuf;
					objPFUNC.feature_identified = false;
					objPFUNC.identified = false;
					localstr >> sbuf >> objPFUNC.name >> objPFUNC.coeff;
					objREGTYPE.postfunc.push_back(objPFUNC);
				}
				if (sbuf.find("STOP") != std::string::npos || readregchain.eof()) {
					break;
				}
			}
		}

	}

	for (int ichain = 0; ichain < (signed)objREGTYPE.funcchain.size(); ichain++) {

		std::vector<std::string> dataset_firstcolumn;
		std::vector<std::string> beta_score_inverse;
		std::vector<std::string> beta_revised;
		readfuncchain.open(objREGTYPE.funcchain[ichain].filename);
		beta_revised.push_back("");
		if (readfuncchain.is_open()) {
			while (true) {
				readfuncchain >> sbuf;
				if (readfuncchain.eof()) {
					break;
				}
				else if (sbuf == "******************************DATASET******************************") {
					std::getline(readfuncchain, sbuf);
					readfuncchain >> sbuf;
					if (sbuf.find("Name") != std::string::npos) {
						readfuncchain >> sbuf;
					}
					std::getline(readfuncchain, sbuf);
					localstr.str("");
					localstr.clear();
					localstr << sbuf;
					std::string dataset_firstcolumn_element = "";
					while (true) {
						localstr >> sbuf;
						dataset_firstcolumn_element += sbuf;
						if (localstr.eof()) {
							break;
						}
						else {
							dataset_firstcolumn_element += " ";
						}
					}
					dataset_firstcolumn.push_back(dataset_firstcolumn_element);
				}
				else if (sbuf == "ModelDF") {
					readfuncchain >> sbuf;
					std::getline(readfuncchain, sbuf);
					beta_score_inverse.push_back(sbuf);
				}
				else if (sbuf == "det_r2_revised") {
					readfuncchain >> sbuf >> sbuf;
					if (sbuf == ":") {
						std::getline(readfuncchain, sbuf);
						beta_revised.push_back(sbuf);
					}
				}
			}
			
		}
		readfuncchain.close();
		
		std::string target_dataset_firstcolumn = "";
		std::string target_beta_score_inverse = "";
		std::string target_beta_revised = "";
		int ibeta_t = objREGTYPE.ibeta_target;
		if (objREGTYPE.ibeta_target >= (signed)dataset_firstcolumn.size()) {
			ibeta_t = 0;
		}
		target_dataset_firstcolumn = dataset_firstcolumn[ibeta_t];
		ibeta_t = objREGTYPE.ibeta_target;
		if (objREGTYPE.ibeta_target >= (signed)beta_score_inverse.size()) {
			ibeta_t = 0;
		}
		target_beta_score_inverse = beta_score_inverse[ibeta_t];
		ibeta_t = objREGTYPE.ibeta_target;
		if (objREGTYPE.ibeta_target >= (signed)beta_revised.size()) {
			ibeta_t = 0;
		}
		target_beta_revised = beta_revised[ibeta_t];

		objREGTYPE.funcchain[ichain].allow_beta = false;
		if (target_beta_score_inverse != "" && target_beta_revised != "") {
			objREGTYPE.funcchain[ichain].allow_beta = true;
		}

		if (target_beta_score_inverse != "") {
			localstr.str("");
			localstr.clear();
			localstr << target_beta_score_inverse;
			localstr >> sbuf >> sbuf >> objREGTYPE.funcchain[ichain].score_inverse_slope >> sbuf >> sbuf >> sbuf >> objREGTYPE.funcchain[ichain].score_inverse_intercept;
		}
		if (target_beta_revised != "") {
			localstr.str("");
			localstr.clear();
			localstr << target_beta_revised;
			localstr >> sbuf >> sbuf >> objREGTYPE.funcchain[ichain].revise_intercept >> sbuf >> objREGTYPE.funcchain[ichain].revise_slope >> sbuf;
		}

		std::vector<std::string>().swap(dataset_firstcolumn);
		std::vector<std::string>().swap(beta_score_inverse);
		std::vector<std::string>().swap(beta_revised);
		localstr.str("");
		localstr.clear();
		localstr << target_dataset_firstcolumn;
		while (true) {
			FEATURE objFEATURE;
			localstr >> objFEATURE.name;
			objFEATURE.identified = identify_anades_chain(&objFEATURE.type, &objFEATURE.i0anades, &objFEATURE.i0mamssk, &objFEATURE.latex, &objFEATURE.mathematica, & objFEATURE.CationName, &objFEATURE.name, anadescriptor, AnionName, CationNames, AveOnly);
			objREGTYPE.funcchain[ichain].feature.push_back(objFEATURE);
			if (localstr.eof()) {
				break;
			}
		}

		for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
			if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
				if (*verbose) {
					if (pr_id == 0) {
						std::cout << *RegFilename << " -- identifying funcchain " << ichain << " : " << ift << "\n";
						std::cout.flush();
					}
				}
				int starting_jft = 0;
				int starting_kft = 1;
				for (int jft = starting_jft; jft < (signed)objREGTYPE.funcchain[ichain].feature.size(); jft++) {
					if (objREGTYPE.funcchain[ichain].feature[jft].identified && objREGTYPE.funcchain[ichain].feature[ift].name.find(objREGTYPE.funcchain[ichain].feature[jft].name) != std::string::npos) {
						for (int kft = starting_kft; kft < (signed)objREGTYPE.funcchain[ichain].feature.size(); kft++) {
							if (objREGTYPE.funcchain[ichain].feature[kft].identified && objREGTYPE.funcchain[ichain].feature[ift].name.find(objREGTYPE.funcchain[ichain].feature[kft].name) != std::string::npos) {
								if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain].feature[jft].name + "*" + objREGTYPE.funcchain[ichain].feature[kft].name + ")") {
									objREGTYPE.funcchain[ichain].feature[ift].identified = true;
									objREGTYPE.funcchain[ichain].feature[ift].type = 1;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[0] = objREGTYPE.funcchain[ichain].feature[jft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[0] = objREGTYPE.funcchain[ichain].feature[jft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[1] = objREGTYPE.funcchain[ichain].feature[kft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[1] = objREGTYPE.funcchain[ichain].feature[kft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1m0d1 = 0;
									objREGTYPE.funcchain[ichain].feature[ift].latex = objREGTYPE.funcchain[ichain].feature[jft].latex + objREGTYPE.funcchain[ichain].feature[kft].latex;
									objREGTYPE.funcchain[ichain].feature[ift].mathematica = "((" + objREGTYPE.funcchain[ichain].feature[jft].mathematica + ")*(" + objREGTYPE.funcchain[ichain].feature[kft].mathematica + "))";
									starting_jft = jft;
									starting_kft = kft;
									break;
								}
								else if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain].feature[jft].name + "/" + objREGTYPE.funcchain[ichain].feature[kft].name + ")") {
									objREGTYPE.funcchain[ichain].feature[ift].identified = true;
									objREGTYPE.funcchain[ichain].feature[ift].type = 1;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[0] = objREGTYPE.funcchain[ichain].feature[jft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[0] = objREGTYPE.funcchain[ichain].feature[jft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1mamssk[1] = objREGTYPE.funcchain[ichain].feature[kft].i0mamssk;
									objREGTYPE.funcchain[ichain].feature[ift].i1anades[1] = objREGTYPE.funcchain[ichain].feature[kft].i0anades;
									objREGTYPE.funcchain[ichain].feature[ift].i1m0d1 = 1;
									objREGTYPE.funcchain[ichain].feature[ift].latex = "\\frac{" + objREGTYPE.funcchain[ichain].feature[jft].latex + "}{" + objREGTYPE.funcchain[ichain].feature[kft].latex + "}";
									objREGTYPE.funcchain[ichain].feature[ift].mathematica = "((" + objREGTYPE.funcchain[ichain].feature[jft].mathematica + ")/(" + objREGTYPE.funcchain[ichain].feature[kft].mathematica + "))";
									starting_jft = jft;
									starting_kft = kft;
									break;
								}
							}
						}
					}
					if (objREGTYPE.funcchain[ichain].feature[ift].identified) {
						break;
					}
				}
			}
		}

		if (ichain > 0) {

			for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
				if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
					if (*verbose) {
						if (pr_id == 0) {
							std::cout << *RegFilename << " -- identifying funcchain " << ichain << " : " << ift << "\n";
							std::cout.flush();
						}
					}
					for (int ifc = 0; ifc < (signed)objREGTYPE.funcchain[ichain - 1].func.size(); ifc++) {
						if (objREGTYPE.funcchain[ichain - 1].func[ifc].identified && objREGTYPE.funcchain[ichain - 1].func[ifc].feature_identified) {
							if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain - 1].func[ifc].name + ")") {
								objREGTYPE.funcchain[ichain].feature[ift].identified = true;
								objREGTYPE.funcchain[ichain].feature[ift].type = 2;
								objREGTYPE.funcchain[ichain].feature[ift].i2whichlastfunc = ifc;
								objREGTYPE.funcchain[ichain].feature[ift].latex = objREGTYPE.funcchain[ichain - 1].func[ifc].latex;
								objREGTYPE.funcchain[ichain].feature[ift].mathematica = "(" + objREGTYPE.funcchain[ichain - 1].func[ifc].mathematica + ")";
								break;
							}
						}
					}
				}
			}

			for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
				if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
					if (*verbose) {
						if (pr_id == 0) {
							std::cout << *RegFilename << " -- identifying funcchain " << ichain << " : " << ift << "\n";
							std::cout.flush();
						}
					}
					int starting_ifc = 0;
					int starting_jfc = 1;
					for (int ifc = starting_ifc; ifc < (signed)objREGTYPE.funcchain[ichain - 1].func.size(); ifc++) {
						if (objREGTYPE.funcchain[ichain - 1].func[ifc].identified && objREGTYPE.funcchain[ichain - 1].func[ifc].feature_identified && objREGTYPE.funcchain[ichain].feature[ift].name.find(objREGTYPE.funcchain[ichain - 1].func[ifc].name) != std::string::npos) {
							for (int jfc = starting_jfc; jfc < (signed)objREGTYPE.funcchain[ichain - 1].func.size(); jfc++) {
								if (objREGTYPE.funcchain[ichain - 1].func[jfc].identified && objREGTYPE.funcchain[ichain - 1].func[jfc].feature_identified && objREGTYPE.funcchain[ichain].feature[ift].name.find(objREGTYPE.funcchain[ichain - 1].func[jfc].name) != std::string::npos) {
									if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain - 1].func[ifc].name + "*" + objREGTYPE.funcchain[ichain - 1].func[jfc].name + ")") {
										objREGTYPE.funcchain[ichain].feature[ift].identified = true;
										objREGTYPE.funcchain[ichain].feature[ift].type = 3;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[0] = ifc;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[1] = jfc;
										objREGTYPE.funcchain[ichain].feature[ift].i3m0d1 = 0;
										objREGTYPE.funcchain[ichain].feature[ift].latex = objREGTYPE.funcchain[ichain - 1].func[ifc].latex + objREGTYPE.funcchain[ichain - 1].func[jfc].latex;
										objREGTYPE.funcchain[ichain].feature[ift].mathematica = "((" + objREGTYPE.funcchain[ichain - 1].func[ifc].mathematica + ")*(" + objREGTYPE.funcchain[ichain - 1].func[jfc].mathematica + "))";
										starting_ifc = ifc;
										starting_jfc = jfc;
										break;
									}
									else if (objREGTYPE.funcchain[ichain].feature[ift].name == "(" + objREGTYPE.funcchain[ichain - 1].func[ifc].name + "/" + objREGTYPE.funcchain[ichain - 1].func[jfc].name + ")") {
										objREGTYPE.funcchain[ichain].feature[ift].identified = true;
										objREGTYPE.funcchain[ichain].feature[ift].type = 3;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[0] = ifc;
										objREGTYPE.funcchain[ichain].feature[ift].i3whichlastfunc[1] = jfc;
										objREGTYPE.funcchain[ichain].feature[ift].i3m0d1 = 1;
										objREGTYPE.funcchain[ichain].feature[ift].latex = "\\frac{" + objREGTYPE.funcchain[ichain - 1].func[ifc].latex + "}{" + objREGTYPE.funcchain[ichain - 1].func[jfc].latex + "}";
										objREGTYPE.funcchain[ichain].feature[ift].mathematica = "((" + objREGTYPE.funcchain[ichain - 1].func[ifc].mathematica + ")/(" + objREGTYPE.funcchain[ichain - 1].func[jfc].mathematica + "))";
										starting_ifc = ifc;
										starting_jfc = jfc;
										break;
									}
								}
							}
						}
						if (objREGTYPE.funcchain[ichain].feature[ift].identified) {
							break;
						}
					}
				}
			}

		}

		for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
			if (!objREGTYPE.funcchain[ichain].feature[ift].identified) {
				if (pr_id == 0) {
					std::cout << "[ERROR] " << *RegFilename << " ... " << objREGTYPE.funcchain[ichain].filename << "\n";
					std::cout.flush();
					std::cout << "[ERROR] " << *RegFilename << " ... " << objREGTYPE.funcchain[ichain].feature[ift].name << " cannot be identified!!\n";
					std::cout.flush();
				}
			}
		}

		readfuncchain.open(objREGTYPE.funcchain[ichain].filename);
		if (readfuncchain.is_open()) {

			int beta_count = -1;
			while (true) {
				readfuncchain >> sbuf;
				if (readfuncchain.eof()) {
					break;
				}
				else if (sbuf == "fit_variable") {
					beta_count++;
					if (beta_count == objREGTYPE.ibeta_target) {
						break;
					}
				}
			}
			std::getline(readfuncchain, sbuf);
			while (true) {
				std::getline(readfuncchain, sbuf);
				if (sbuf.find("dataset.size()") != std::string::npos) {
					break;
				}
				else if (sbuf.find("intercept") == std::string::npos) {
					std::stringstream lcstr;
					lcstr.str("");
					lcstr.clear();
					lcstr << sbuf;
					FUNC objFUNC;
					objFUNC.latex = "UND";
					objFUNC.mathematica = "UND";
					lcstr >> objFUNC.name >> objFUNC.coeff >> objFUNC.zcoeff;
					if (pr_id == 0) {
						std::cout << *RegFilename << " -- identifying funcchain " << ichain << " : " << objFUNC.name << "\n";
						std::cout.flush();
					}
					objFUNC.identified = false;
					for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
						for (int jfc = 0; jfc < (signed)objREGTYPE.orthfunctype.size(); jfc++) {
							if (objREGTYPE.orthfunctype[jfc].userdefined_name.size() > 0) {
								std::string referred_name = objREGTYPE.orthfunctype[jfc].userdefined_name[0] + objREGTYPE.funcchain[ichain].feature[ift].name;
								if (objREGTYPE.orthfunctype[jfc].userdefined_name.size() > 1) {
									for (int iname = 1; iname < (signed)objREGTYPE.orthfunctype[jfc].userdefined_name.size(); iname++) {
										referred_name += "^" + objREGTYPE.orthfunctype[jfc].userdefined_name[iname];
									}
								}
								if (objFUNC.name == referred_name) {
									objFUNC.identified = true;
									objFUNC.whichfeature = ift;
									objFUNC.whichorthfunctype = jfc;
									objFUNC.feature_identified = objREGTYPE.funcchain[ichain].feature[ift].identified;
									if (objREGTYPE.orthfunctype[jfc].FName == "pow" || objREGTYPE.orthfunctype[jfc].FName == "abs") {
										objFUNC.latex = objREGTYPE.funcchain[ichain].feature[ift].latex;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.latex += get_addstring(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											if (objREGTYPE.orthfunctype[jfc].FName == "pow") {
												objFUNC.latex = "\\left(" + objFUNC.latex + "\\right)";
											}
											else {
												objFUNC.latex = "|" + objFUNC.latex + "|";
											}
											objFUNC.latex = objFUNC.latex + "^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "}";
										}
										objFUNC.mathematica = objREGTYPE.funcchain[ichain].feature[ift].mathematica;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.mathematica += get_addstringM(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											if (objREGTYPE.orthfunctype[jfc].FName == "abs") {
												objFUNC.mathematica = "Abs[" + objFUNC.mathematica + "]";
											}
											objFUNC.mathematica = "Power[" + objFUNC.mathematica + ", " + get_powerstringM(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "]";
										}
									}
									else if (objREGTYPE.orthfunctype[jfc].FName == "log") {
										objFUNC.latex = objREGTYPE.funcchain[ichain].feature[ift].latex;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.latex += get_addstring(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										std::string coeff1 = std::to_string(objREGTYPE.orthfunctype[jfc].coeff[1]);
										if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[1], 2.71828)) {
											coeff1 = "e";
										}
										else if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[1], 10.0)) {
											coeff1 = "10";
										}
										objFUNC.latex = "\\log_{" + coeff1 + "}(" + objFUNC.latex + ")";
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											objFUNC.latex = "(" + objFUNC.latex + ")^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "}";
										}
										objFUNC.mathematica = objREGTYPE.funcchain[ichain].feature[ift].mathematica;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.mathematica += get_addstringM(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										coeff1 = std::to_string(objREGTYPE.orthfunctype[jfc].coeff[1]);
										if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[1], 2.71828)) {
											coeff1 = "E";
										}
										else if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[1], 10.0)) {
											coeff1 = "10";
										}
										objFUNC.mathematica = "Log[" + coeff1 + ", " + objFUNC.mathematica + "]";
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											objFUNC.mathematica = "Power[" + objFUNC.mathematica + ", " + get_powerstringM(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "]";
										}
									}
									else if (objREGTYPE.orthfunctype[jfc].FName == "erf" || objREGTYPE.orthfunctype[jfc].FName == "erfc") {
										objFUNC.latex = objREGTYPE.funcchain[ichain].feature[ift].latex;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.latex += get_addstring(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
												objFUNC.latex = get_mulstring(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "\\left(" + objFUNC.latex + "\\right)";
											}
											else {
												objFUNC.latex = get_mulstring(&objREGTYPE.orthfunctype[jfc].coeff[1]) + objFUNC.latex;
											}
										}
										objFUNC.latex = "\\mathrm{" + objREGTYPE.orthfunctype[jfc].FName + "}\\left(" + objFUNC.latex + "\\right)";
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											objFUNC.latex = "\\left(" + objFUNC.latex + "\\right)^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "}";
										}
										objFUNC.mathematica = objREGTYPE.funcchain[ichain].feature[ift].mathematica;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.mathematica += get_addstringM(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											objFUNC.mathematica = get_mulstringM(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "(" + objFUNC.mathematica + ")";
										}
										if (objREGTYPE.orthfunctype[jfc].FName == "erf") {
											objFUNC.mathematica = "Erf[" + objFUNC.mathematica + "]";
										}
										else {
											objFUNC.mathematica = "Erfc[" + objFUNC.mathematica + "]";
										}
										
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											objFUNC.mathematica = "Power[" + objFUNC.mathematica + ", " + get_powerstringM(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "]";
										}
									}
									else if (objREGTYPE.orthfunctype[jfc].FName == "sin" || objREGTYPE.orthfunctype[jfc].FName == "cos" || objREGTYPE.orthfunctype[jfc].FName == "tan" || objREGTYPE.orthfunctype[jfc].FName == "sinh" || objREGTYPE.orthfunctype[jfc].FName == "cosh" || objREGTYPE.orthfunctype[jfc].FName == "tanh") {
										objFUNC.latex = objREGTYPE.funcchain[ichain].feature[ift].latex;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.latex += get_addstring(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
												objFUNC.latex = get_mulstring(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "\\left(" + objFUNC.latex + "\\right)";
											}
											else {
												objFUNC.latex = get_mulstring(&objREGTYPE.orthfunctype[jfc].coeff[1]) + objFUNC.latex;
											}
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											if (!(objREGTYPE.orthfunctype[jfc].coeff[0] == 0.0 && objREGTYPE.orthfunctype[jfc].coeff[1] == 1.0)) {
												objFUNC.latex = "\\left(" + objFUNC.latex + "\\right)^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "}";
											}
											else {
												objFUNC.latex = objFUNC.latex + "^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "}";
											}
										}
										objFUNC.latex = "\\" + objREGTYPE.orthfunctype[jfc].FName + "\\left(" + objFUNC.latex + "\\right)";
										if (objREGTYPE.orthfunctype[jfc].coeff[3] != 1.0) {
											objFUNC.latex = "\\left(" + objFUNC.latex + "\\right)^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[3]) + "}";
										}
										objFUNC.mathematica = objREGTYPE.funcchain[ichain].feature[ift].mathematica;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.mathematica += get_addstringM(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											objFUNC.mathematica = get_mulstringM(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "(" + objFUNC.mathematica + ")";
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											objFUNC.mathematica = "Power[" + objFUNC.mathematica + ", " + get_powerstringM(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "]";
										}
										if (objREGTYPE.orthfunctype[jfc].FName == "sin") {
											objFUNC.mathematica = "Sin[" + objFUNC.mathematica + "]";
										}
										else if (objREGTYPE.orthfunctype[jfc].FName == "cos") {
											objFUNC.mathematica = "Cos[" + objFUNC.mathematica + "]";
										}
										else if (objREGTYPE.orthfunctype[jfc].FName == "tan") {
											objFUNC.mathematica = "Tan[" + objFUNC.mathematica + "]";
										}
										else if (objREGTYPE.orthfunctype[jfc].FName == "sinh") {
											objFUNC.mathematica = "Sinh[" + objFUNC.mathematica + "]";
										}
										else if (objREGTYPE.orthfunctype[jfc].FName == "cosh") {
											objFUNC.mathematica = "Cosh[" + objFUNC.mathematica + "]";
										}
										else if (objREGTYPE.orthfunctype[jfc].FName == "tanh") {
											objFUNC.mathematica = "Tanh[" + objFUNC.mathematica + "]";
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[3] != 1.0) {
											objFUNC.mathematica = "Power[" + objFUNC.mathematica + ", " + get_powerstringM(&objREGTYPE.orthfunctype[jfc].coeff[3]) + "]";
										}
									}
									else if (objREGTYPE.orthfunctype[jfc].FName == "exp") {
										objFUNC.latex = objREGTYPE.funcchain[ichain].feature[ift].latex;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.latex += get_addstring(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
												objFUNC.latex = "\\left(" + objFUNC.latex + "\\right)^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "}";
											}
											else {
												objFUNC.latex = objFUNC.latex + "^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "}";
											}
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											if (!(objREGTYPE.orthfunctype[jfc].coeff[0] == 0.0 && objREGTYPE.orthfunctype[jfc].coeff[1] == 0.0)) {
												objFUNC.latex = get_mulstring(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "\\left(" + objFUNC.latex + "\\right)";
											}
											else {
												objFUNC.latex = get_mulstring(&objREGTYPE.orthfunctype[jfc].coeff[2]) + objFUNC.latex;
											}	
										}
										if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[3], 2.71828)) {
											objFUNC.latex = "\\exp\\left(" + objFUNC.latex + "\\right)";
										}
										else if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[3], 10.0)) {
											objFUNC.latex = "10^{" + objFUNC.latex + "}";
										}
										else {
											objFUNC.latex = std::to_string(objREGTYPE.orthfunctype[jfc].coeff[3]) + "^{" + objFUNC.latex + "}";
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[4] != 1.0) {
											objFUNC.latex = "\\left(" + objFUNC.latex + "\\right)^{" + get_powerstring(&objREGTYPE.orthfunctype[jfc].coeff[4]) + "}";
										}
										objFUNC.mathematica = objREGTYPE.funcchain[ichain].feature[ift].mathematica;
										if (objREGTYPE.orthfunctype[jfc].coeff[0] != 0.0) {
											objFUNC.mathematica += get_addstringM(&objREGTYPE.orthfunctype[jfc].coeff[0]);
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[1] != 1.0) {
											objFUNC.mathematica = "Power[" + objFUNC.mathematica + ", " + get_powerstringM(&objREGTYPE.orthfunctype[jfc].coeff[1]) + "]";
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[2] != 1.0) {
											objFUNC.mathematica = get_mulstringM(&objREGTYPE.orthfunctype[jfc].coeff[2]) + "(" + objFUNC.mathematica + ")";
										}
										if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[3], 2.71828)) {
											objFUNC.mathematica = "Exp[" + objFUNC.mathematica + "]";
										}
										else if (chk_eql(objREGTYPE.orthfunctype[jfc].coeff[3], 10.0)) {
											objFUNC.mathematica = "Power[10, " + objFUNC.mathematica + "]";
										}
										else {
											objFUNC.mathematica = "Power[" + std::to_string(objREGTYPE.orthfunctype[jfc].coeff[3]) + ", " + objFUNC.mathematica + "]";
										}
										if (objREGTYPE.orthfunctype[jfc].coeff[4] != 1.0) {
											objFUNC.mathematica = "Power[" + objFUNC.mathematica + ", " + get_powerstringM(&objREGTYPE.orthfunctype[jfc].coeff[4]) + "]";
										}
									}
									break;
								}
							}
						}
						if (objFUNC.identified) {
							break;
						}
					}
					objREGTYPE.funcchain[ichain].func.push_back(objFUNC);
				}
				else {
					std::stringstream lcstr;
					lcstr.str("");
					lcstr.clear();
					lcstr << sbuf;
					lcstr >> sbuf >> objREGTYPE.funcchain[ichain].intercept;
				}
			}

			std::vector<std::string> zerocoeff;
			std::vector<int> izerocoeff;
			for (int icoeff = 0; icoeff < (signed)objREGTYPE.funcchain[ichain].func.size(); icoeff++) {
				if (objREGTYPE.funcchain[ichain].func[icoeff].coeff == 0.0) {
					zerocoeff.push_back(objREGTYPE.funcchain[ichain].func[icoeff].name);
					izerocoeff.push_back(icoeff);
				}
			}

			if (zerocoeff.size() > 0) {

				std::vector<int> targetcolumn;
				std::getline(readfuncchain, sbuf);
				std::stringstream lcstr;
				lcstr.str("");
				lcstr.clear();
				lcstr << sbuf;
				int columncount = 0;
				while (true) {
					lcstr >> sbuf;
					if (sbuf == objREGTYPE.name) {
						targetcolumn.push_back(columncount);
					}
					else {
						int columnfound = -1;
						for (int icoeff = 0; icoeff < (signed)zerocoeff.size(); icoeff++) {
							if (sbuf == zerocoeff[icoeff]) {
								columnfound = icoeff;
								break;
							}
						}
						if (columnfound != -1) {
							targetcolumn.push_back(columncount);
						}
					}
					columncount++;
					if (lcstr.eof()) {
						break;
					}
				}

				std::vector<BACKUPV> backupv;
				while (true) {
					std::getline(readfuncchain, sbuf);
					if (sbuf.find("TEST") != std::string::npos) {
						break;
					}
					BACKUPV objBACKUPV;
					lcstr.str("");
					lcstr.clear();
					lcstr << sbuf;
					int ccount = 0;
					while (true) {
						lcstr >> sbuf;
						bool take_memo = false;
						bool touch_end = false;
						for (int ic = 0; ic < (signed)targetcolumn.size(); ic++) {
							if (ccount == targetcolumn[ic]) {
								take_memo = true;
								if (ic == targetcolumn[targetcolumn.size() - 1]) {
									touch_end = true;
								}
								break;
							}
						}
						if (take_memo) {
							std::stringstream lcstr2;
							lcstr2.str("");
							lcstr2.clear();
							lcstr2 << sbuf;
							long double ldbuf2;
							lcstr2 >> ldbuf2;
							objBACKUPV.v.push_back(ldbuf2);
						}
						ccount++;
						if (lcstr.eof() || touch_end) {
							break;
						}
					}
					backupv.push_back(objBACKUPV);
				}

				std::vector<long double> m_series, s_series;
				for (int itc = 0; itc < (signed)targetcolumn.size(); itc++) {
					m_series.push_back(0.0);
					s_series.push_back(0.0);
				}
				for (int ib = 0; ib < (signed)backupv.size(); ib++) {
					for (int itc = 0; itc < (signed)backupv[ib].v.size(); itc++) {
						m_series[itc] += backupv[ib].v[itc];
					}
				}
				for (int itc = 0; itc < (signed)targetcolumn.size(); itc++) {
					m_series[itc] /= (long double)backupv.size();
				}
				for (int ib = 0; ib < (signed)backupv.size(); ib++) {
					for (int itc = 0; itc < (signed)backupv[ib].v.size(); itc++) {
						s_series[itc] += powl(backupv[ib].v[itc] - m_series[itc], 2.0);
					}
				}
				for (int itc = 0; itc < (signed)targetcolumn.size(); itc++) {
					s_series[itc] /= (long double)backupv.size();
					s_series[itc] = sqrtl(s_series[itc]);
				}

				for (int icoeff = 0; icoeff < (signed)izerocoeff.size(); icoeff++) {
					objREGTYPE.funcchain[ichain].func[izerocoeff[icoeff]].coeff = objREGTYPE.funcchain[ichain].func[izerocoeff[icoeff]].zcoeff * s_series[0] / s_series[icoeff + 1];
					if (pr_id == 0) {
						std::cout << "[RECOVERED coeff] " << objREGTYPE.funcchain[ichain].func[izerocoeff[icoeff]].name << " " << objREGTYPE.funcchain[ichain].func[izerocoeff[icoeff]].coeff << "\n";
						std::cout.flush();
					}
				}

				std::vector<BACKUPV>().swap(backupv);
				std::vector<long double>().swap(m_series);
				std::vector<long double>().swap(s_series);

			}

			std::vector<std::string>().swap(zerocoeff);

			if (ichain == 0 && objREGTYPE.postfunc.size() > 0) {
				for (int ip = 0; ip < (signed)objREGTYPE.postfunc.size(); ip++) {
					if (pr_id == 0) {
						std::cout << *RegFilename << " -- identifying postfunc : " << objREGTYPE.postfunc[ip].name << "\n";
						std::cout.flush();
					}
					for (int ift = 0; ift < (signed)objREGTYPE.funcchain[ichain].feature.size(); ift++) {
						for (int jfc = 0; jfc < (signed)objREGTYPE.orthfunctype.size(); jfc++) {
							if (objREGTYPE.orthfunctype[jfc].userdefined_name.size() > 0) {
								std::string referred_name = objREGTYPE.orthfunctype[jfc].userdefined_name[0] + objREGTYPE.funcchain[ichain].feature[ift].name;
								if (objREGTYPE.orthfunctype[jfc].userdefined_name.size() > 1) {
									for (int iname = 1; iname < (signed)objREGTYPE.orthfunctype[jfc].userdefined_name.size(); iname++) {
										referred_name += "^" + objREGTYPE.orthfunctype[jfc].userdefined_name[iname];
									}
								}
								if (objREGTYPE.postfunc[ip].name == referred_name) {
									objREGTYPE.postfunc[ip].identified = true;
									objREGTYPE.postfunc[ip].whichfeature = ift;
									objREGTYPE.postfunc[ip].whichorthfunctype = jfc;
									objREGTYPE.postfunc[ip].feature_identified = objREGTYPE.funcchain[ichain].feature[ift].identified;
									break;
								}
							}
						}
						if (objREGTYPE.postfunc[ip].identified) {
							break;
						}
					}
				}
			}

		}
		readfuncchain.close();

	}

	readregchain.close();

	return objREGTYPE;

}

void generate_tuples(std::vector<int>* occs, int pos, int* imin, int* imax, std::vector<std::string>* test_name, BAGTYPE* bagtype, int* ibag, std::vector<REGTYPE>* regtype, std::vector<std::string>* RegModelFilename, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::string* AnionName, std::vector<std::string>* CationNames, bool* verbose, std::string* AnalysisOutputFilename, std::vector<int>* tcolumns, int* tcolumnR, bool* AveOnly) {

	long double increment = (long double)(*imax - *imin);

	for (int i = *imin; i <= *imax; i++) {
		(*occs)[pos] = i;
		if (pos != occs->size() - 1) {
			generate_tuples(occs, pos + 1, imin, imax, test_name, bagtype, ibag, regtype, RegModelFilename, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, AnionName, CationNames, verbose, AnalysisOutputFilename, tcolumns, tcolumnR, AveOnly);
		}
		else {
			std::vector<long double> voccs;
			std::vector<long double> test_anades_degree;
			long double test_OM_degree = 0.0;
			long double test_CO_degree = 0.0;
			for (int j = 0; j < (signed)occs->size(); j++) {
				voccs.push_back((long double)(*occs)[j] / increment);
			}
			for (int j = 0; j < (signed)voccs.size(); j++) {
				if (j < test_name->size()) {
					bool found = false;
					if ((*test_name)[j] == *AnionName + "M") {
						test_OM_degree = voccs[j];
						found = true;
					}
					else {
						for (int ic = 0; ic < (signed)CationNames->size(); ic++) {
							if ((*test_name)[j] == (*CationNames)[ic] + *AnionName) {
								test_CO_degree = voccs[j];
								found = true;
								break;
							}
						}
					}
					if (!found) {
						test_anades_degree.push_back(voccs[j]);
					}
				}
			}
			std::vector<long double> v_series;
			std::vector<long double> v_series_acc;
			v_series.clear();
			v_series_acc.clear();
			long double v_stdev;
			bool post_if_exists = true;
			MAT nullMAT;
			
			std::vector<long double>().swap(floating_v);
			std::vector<long double>().swap(floating_v_min);
			std::vector<long double>().swap(floating_v_max);
			floating_v.clear();
			floating_v_min.clear();
			floating_v_max.clear();
			for (int j = 0; j < voccs.size(); j++) {
				floating_v.push_back(0.0);
				floating_v_min.push_back(0.0);
				floating_v_max.push_back(0.0);
			}
			calc_bag_chain_partial_test(&v_series, &v_series_acc, &v_stdev, regtype, bagtype, &nullMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, AnionName, CationNames, test_name, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly);
			if (pr_id == 0) {
				for (int j = 0; j < voccs.size(); j++) {
					writezcext << voccs[j] << "\t";
					writezcext.flush();
				}
				for (int j = 0; j < floating_v.size(); j++) {
					writezcext << floating_v[j] << "\t";
					writezcext.flush();
				}
			}
			readanalysis.open(*AnalysisOutputFilename);
			readanalysisR.open("__R__." + * AnalysisOutputFilename);
			std::string sbufR;
			std::getline(readanalysis, sbuf);
			std::getline(readanalysisR, sbufR);
			std::vector<long double> real_values;
			std::vector<long double> reg_values;
			while (true) {
				std::getline(readanalysis, sbuf);
				std::getline(readanalysisR, sbufR);
				if (sbuf.find("===") != std::string::npos) {
					break;
				}
				std::stringstream linebuf, linebufR;
				linebuf.str("");
				linebuf.clear();
				linebuf << sbuf;
				linebufR.str("");
				linebufR.clear();
				linebufR << sbufR;
				int icol = -1;
				int icolR = -1;
				std::vector<long double> collect_v;
				long double collect_vR = 0.0;
				for (int ic = 0; ic < tcolumns->size(); ic++) {
					collect_v.push_back(0.0);
				}
				while (true) {
					linebuf >> sbuf;
					icol++;
					for (int ic = 0; ic < tcolumns->size(); ic++) {
						if ((*tcolumns)[ic] == icol) {
							std::stringstream cbuf;
							cbuf.str("");
							cbuf.clear();
							cbuf << sbuf;
							cbuf >> ldbuf;
							collect_v[ic] = ldbuf;
							break;
						}
					}
					if (linebuf.eof()) {
						break;
					}
				}
				while (true) {
					linebufR >> sbufR;
					icolR++;
					if (*tcolumnR == icolR) {
						std::stringstream cbuf;
						cbuf.str("");
						cbuf.clear();
						cbuf << sbufR;
						cbuf >> ldbuf;
						collect_vR = ldbuf;
						break;
					}
				}
				bool take = true;
				for (int iblock = 0; iblock < (signed)collect_v.size() - 1; iblock++) {
					long double delta = fabsl(0.5 * (floating_v_max[iblock] - floating_v_min[iblock]) / increment);
					if (!(collect_v[iblock + 1] >= floating_v[iblock] - delta && collect_v[iblock + 1] < floating_v[iblock] + delta)) {
						take = false;
						break;
					}
				}
				if (take) {
					real_values.push_back(collect_v[0]);
					reg_values.push_back(collect_vR);
				}
				std::vector<long double>().swap(collect_v);
				if (readanalysis.eof() || readanalysisR.eof()) {
					break;
				}
			}
			if (pr_id == 0) {
				writezcext << v_series_acc[v_series_acc.size() - 1] << "\t" << v_stdev << "\t";
				writezcext.flush();
			}
			if (real_values.size() == 0 || reg_values.size() == 0) {
				if (pr_id == 0) {
					writezcext << "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
					writezcext.flush();
				}
			}
			else {
				long double min = *VInfinite;
				long double max = -*VInfinite;
				long double ave = 0.0;
				long double stdv = 0.0;
				for (int ireal = 0; ireal < (signed)real_values.size(); ireal++) {
					if (real_values[ireal] > max) {
						max = real_values[ireal];
					}
					ave += real_values[ireal];
					if (real_values[ireal] < min) {
						min = real_values[ireal];
					}
				}
				ave /= (long double)real_values.size();
				for (int ireal = 0; ireal < (signed)real_values.size(); ireal++) {
					stdv += powl(real_values[ireal] - ave, 2.0);
				}
				stdv /= (long double)real_values.size();
				stdv = sqrtl(stdv);
				if (pr_id == 0) {
					writezcext << real_values.size() << "\t" << min << "\t" << ave << "\t" << stdv << "\t" << max << "\t";
					writezcext.flush();
				}
				min = *VInfinite;
				max = -*VInfinite;
				ave = 0.0;
				stdv = 0.0;
				for (int ireg = 0; ireg < (signed)reg_values.size(); ireg++) {
					if (reg_values[ireg] > max) {
						max = reg_values[ireg];
					}
					ave += reg_values[ireg];
					if (reg_values[ireg] < min) {
						min = reg_values[ireg];
					}
				}
				ave /= (long double)reg_values.size();
				for (int ireg = 0; ireg < (signed)reg_values.size(); ireg++) {
					stdv += powl(reg_values[ireg] - ave, 2.0);
				}
				stdv /= (long double)reg_values.size();
				stdv = sqrtl(stdv);
				if (pr_id == 0) {
					writezcext << reg_values.size() << "\t" << min << "\t" << ave << "\t" << stdv << "\t" << max;
					writezcext.flush();
				}
				if (real_values.size() == reg_values.size()) {
					long double average_absolute_error = 0.0;
					long double max_absolute_error = -1000000000000.0;
					for (int ireg = 0; ireg < (signed)reg_values.size(); ireg++) {
						long double absolute_error = fabsl(reg_values[ireg] - real_values[ireg]);
						average_absolute_error += absolute_error;
						if (absolute_error > max_absolute_error) {
							max_absolute_error = absolute_error;
						}
					}
					average_absolute_error /= (long double)real_values.size();
					if (pr_id == 0) {
						writezcext << "\t" << average_absolute_error << "\t" << max_absolute_error << "\n";
						writezcext.flush();
					}
				}
				else {
					if (pr_id == 0) {
						writezcext << "\t-\t-\n";
						writezcext.flush();
					}
				}
			}
			readanalysis.close();
			readanalysisR.close();
			std::vector<long double>().swap(voccs);
			std::vector<long double>().swap(test_anades_degree);
			std::vector<long double>().swap(real_values);
			std::vector<long double>().swap(reg_values);
		}
	}

}

void analyze_terms_zc(BAGTYPE* bagtype, int* ibag, std::vector<REGTYPE>* regtype, std::vector<std::string>* RegModelFilename, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::string* AnionName, std::vector<std::string>* CationNames, bool* verbose, int* NCavailles, std::string* AnalysisOutputFilename, int* AnalysisMesh, bool* AveOnly) {
		
	for (int ir = 0; ir < (signed)bagtype->which_regtype.size(); ir++) {
		for (int ifc = (signed)(*regtype)[bagtype->which_regtype[ir]].funcchain.size() - 1; ifc < (signed)(*regtype)[bagtype->which_regtype[ir]].funcchain.size(); ifc++) {
			for (int iff = 0; iff < (signed)(*regtype)[bagtype->which_regtype[ir]].funcchain[ifc].func.size(); iff++) {
	
				for (int it = 0; it < (signed)(*regtype)[bagtype->which_regtype[ir]].termana.size(); it++) {
					bool all_exist = true;
					for (int itn = 0; itn < (signed)(*regtype)[bagtype->which_regtype[ir]].termana[it].name.size(); itn++) {
						if ((*regtype)[bagtype->which_regtype[ir]].funcchain[ifc].func[iff].name.find((*regtype)[bagtype->which_regtype[ir]].termana[it].name[itn]) == std::string::npos) {
							all_exist = false;
							break;
						}
					}
					if (all_exist) {
						(*regtype)[bagtype->which_regtype[ir]].termana[it].count++;
						WHERETOFIND objWHERETOFIND;
						objWHERETOFIND.whichreg = bagtype->which_regtype[ir];
						objWHERETOFIND.whichregfilename = (*RegModelFilename)[bagtype->which_regtype[ir]];
						objWHERETOFIND.whichfuncchain = ifc;
						objWHERETOFIND.whichfunc = iff;
						objWHERETOFIND.zcoeff = (*regtype)[bagtype->which_regtype[ir]].funcchain[ifc].func[iff].zcoeff;
						(*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind.push_back(objWHERETOFIND);
					}
				}

			}
		}
	}

	for (int ir = 0; ir < (signed)bagtype->which_regtype.size(); ir++) {
		for (int it = 0; it < (signed)(*regtype)[bagtype->which_regtype[ir]].termana.size(); it++) {
			for (int izc = 0; izc < (signed)(*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind.size(); izc++) {
				if (fabsl((*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[izc].zcoeff) < (*regtype)[bagtype->which_regtype[ir]].termana[it].min_absz) {
					(*regtype)[bagtype->which_regtype[ir]].termana[it].min_absz = fabsl((*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[izc].zcoeff);
				}
				(*regtype)[bagtype->which_regtype[ir]].termana[it].ave_absz += bagtype->stackingensemble[bagtype->stackingensemble.size() - 1].coeff[ir + 1] * fabsl((*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[izc].zcoeff);
				if (fabsl((*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[izc].zcoeff) > (*regtype)[bagtype->which_regtype[ir]].termana[it].max_absz) {
					(*regtype)[bagtype->which_regtype[ir]].termana[it].max_absz = fabsl((*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[izc].zcoeff);
				}
			}
			if ((*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind.size() == 0) {
				(*regtype)[bagtype->which_regtype[ir]].termana[it].ave_absz = 0.0;
			}
		}
	}

	for (int ib = 0; ib < (signed)bagtype->termanaforbag.size(); ib++) {
		for (int ir = 0; ir < (signed)bagtype->which_regtype.size(); ir++) {
			bagtype->termanaforbag[ib].count += (*regtype)[bagtype->which_regtype[ir]].termana[ib].count;
			if ((*regtype)[bagtype->which_regtype[ir]].termana[ib].count > 0) {
				bagtype->termanaforbag[ib].metacount++;
			}
			if ((*regtype)[bagtype->which_regtype[ir]].termana[ib].min_absz < bagtype->termanaforbag[ib].min_absz) {
				bagtype->termanaforbag[ib].min_absz = (*regtype)[bagtype->which_regtype[ir]].termana[ib].min_absz;
			}
			bagtype->termanaforbag[ib].ave_absz += (*regtype)[bagtype->which_regtype[ir]].termana[ib].ave_absz * (long double)(*regtype)[bagtype->which_regtype[ir]].termana[ib].count;
			if ((*regtype)[bagtype->which_regtype[ir]].termana[ib].max_absz > bagtype->termanaforbag[ib].max_absz) {
				bagtype->termanaforbag[ib].max_absz = (*regtype)[bagtype->which_regtype[ir]].termana[ib].max_absz;
			}
		}
		bagtype->termanaforbag[ib].ave_absz /= (long double)bagtype->termanaforbag[ib].count;
	}

	if (pr_id == 0) {
		writezc.open("bag" + std::to_string(*ibag) + "_" + bagtype->name + ".zc.txt", std::ios::binary);
		writezc << std::fixed << std::setprecision(20);
		writezc << "TermName\tInteractionLevel\tMin(|Z|)\tAve(|Z|)\tMax(|Z|)\tInHowManyModelsAppear\tHowManyAppear\tWhereToFind(z@address)\n";
		writezc.flush();
	}

	std::vector<int> it_list;
	for (int it = 0; it < (signed)bagtype->termanaforbag.size(); it++) {
		if (bagtype->termanaforbag[it].count > 0) {
			it_list.push_back(it);
			if (pr_id == 0) {
				for (int itn = 0; itn < (signed)bagtype->termanaforbag[it].name.size(); itn++) {
					writezc << bagtype->termanaforbag[it].name[itn];
					writezc.flush();
					if (itn != (signed)bagtype->termanaforbag[it].name.size() - 1) {
						writezc << "-";
						writezc.flush();
					}
				}
				writezc << "\t" << bagtype->termanaforbag[it].name.size() << "\t" << bagtype->termanaforbag[it].min_absz << "\t" << bagtype->termanaforbag[it].ave_absz << "\t" << bagtype->termanaforbag[it].max_absz << "\t" << bagtype->termanaforbag[it].metacount << "\t" << bagtype->termanaforbag[it].count << "\t";
				writezc.flush();
				for (int ir = 0; ir < (signed)bagtype->which_regtype.size(); ir++) {
					for (int iwhere = 0; iwhere < (signed)(*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind.size(); iwhere++) {
						writezc << (*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[iwhere].zcoeff << "\t@\t" << (*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[iwhere].whichregfilename << "-" << (*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[iwhere].whichfuncchain << "-" << (*regtype)[bagtype->which_regtype[ir]].termana[it].wheretofind[iwhere].whichfunc << "\t";
						writezc.flush();
					}
				}
				writezc << "\n";
				writezc.flush();
			}
		}
	}

	std::vector<CAVAILLES> cavailles;
	for (int ic = 0; ic < *NCavailles; ic++) {
		CAVAILLES objCAVAILLES;
		cavailles.push_back(objCAVAILLES);
	}

	for (int ic = 0; ic < (signed)cavailles.size(); ic++) {
		int metacount_max = -1;
		long double aveabsz_max = 0.0;
		int target_it = -1;
		int it_nextlevelstart = -1;
		int it_current_level_search = 1;
		for (int it = 0; it < (signed)bagtype->termanaforbag.size(); it++) {
			if ((signed)bagtype->termanaforbag[it].name.size() == it_current_level_search) {
				if (bagtype->termanaforbag[it].metacount > metacount_max) {
					metacount_max = bagtype->termanaforbag[it].metacount;
				}
			}
			else {
				break;
			}
		}
		if (metacount_max != -1) {
			while (true) {
				for (int it = 0; it < (signed)bagtype->termanaforbag.size(); it++) {
					bool passover = false;
					if (ic > 0) {
						for (int jc = 0; jc <= ic - 1; jc++) {
							if (cavailles[jc].it_chain_list[0] == it) {
								passover = true;
								break;
							}
						}
					}
					if (!passover && (signed)bagtype->termanaforbag[it].name.size() == it_current_level_search) {
						if (bagtype->termanaforbag[it].metacount == metacount_max && bagtype->termanaforbag[it].ave_absz > aveabsz_max) {
							aveabsz_max = bagtype->termanaforbag[it].ave_absz;
							target_it = it;
						}
					}
					if ((signed)bagtype->termanaforbag[it].name.size() != it_current_level_search) {
						it_nextlevelstart = it;
						break;
					}
					if (it == (signed)bagtype->termanaforbag.size() - 1) {
						it_nextlevelstart = bagtype->termanaforbag.size();
					}
				}
				if (target_it != -1) {
					break;
				}
				else {
					metacount_max--;
					if (metacount_max == 0) {
						break;
					}
				}
			}
		}

		if (target_it != -1) {
			cavailles[ic].it_chain_list.push_back(target_it);
			if (it_nextlevelstart < (signed)bagtype->termanaforbag.size()) {
				while (true) {
					it_current_level_search++;
					int aveabsz_max = 0.0;
					int target_it = -1;
					for (int it = it_nextlevelstart; it < (signed)bagtype->termanaforbag.size(); it++) {
						if ((signed)bagtype->termanaforbag[it].name.size() == it_current_level_search) {
							if (bagtype->termanaforbag[it].metacount >= metacount_max - 2 && bagtype->termanaforbag[it].ave_absz > aveabsz_max && is_subset(&bagtype->termanaforbag[it].name, &bagtype->termanaforbag[cavailles[ic].it_chain_list[cavailles[ic].it_chain_list.size() - 1]].name)) {
								aveabsz_max = bagtype->termanaforbag[it].ave_absz;
								target_it = it;
							}
						}
						else {
							it_nextlevelstart = it;
							break;
						}
						if (it == (signed)bagtype->termanaforbag.size() - 1) {
							it_nextlevelstart = bagtype->termanaforbag.size();
						}
					}
					if (target_it != -1) {
						cavailles[ic].it_chain_list.push_back(target_it);
						metacount_max = bagtype->termanaforbag[target_it].metacount;
					}
					else {
						break;
					}
					if (it_nextlevelstart >= (signed)bagtype->termanaforbag.size()) {
						break;
					}
				}
			}
		}
	}
	
	if (pr_id == 0) {
		writezc << "===================================================================================\n";
		writezc.flush();
		writezc << "TermName\tInteractionLevel\tMin(|Z|)\tAve(|Z|)\tMax(|Z|)\tInHowManyModelsAppear\tHowManyAppear\tWhereToFind(z@address)\n";
		writezc.flush();
		for (int ic = 0; ic < (signed)cavailles.size(); ic++) {
			for (int it = 0; it < (signed)cavailles[ic].it_chain_list.size(); it++) {
				for (int itn = 0; itn < (signed)bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].name.size(); itn++) {
					writezc << bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].name[itn];
					writezc.flush();
					if (itn != (signed)bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].name.size() - 1) {
						writezc << "-";
						writezc.flush();
					}
				}
				writezc << "\t" << bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].name.size() << "\t" << bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].min_absz << "\t" << bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].ave_absz << "\t" << bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].max_absz << "\t" << bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].metacount << "\t" << bagtype->termanaforbag[cavailles[ic].it_chain_list[it]].count << "\t";
				writezc.flush();
				for (int ir = 0; ir < (signed)bagtype->which_regtype.size(); ir++) {
					for (int iwhere = 0; iwhere < (signed)(*regtype)[bagtype->which_regtype[ir]].termana[cavailles[ic].it_chain_list[it]].wheretofind.size(); iwhere++) {
						writezc << (*regtype)[bagtype->which_regtype[ir]].termana[cavailles[ic].it_chain_list[it]].wheretofind[iwhere].zcoeff << "\t@\t" << (*regtype)[bagtype->which_regtype[ir]].termana[cavailles[ic].it_chain_list[it]].wheretofind[iwhere].whichregfilename << "-" << (*regtype)[bagtype->which_regtype[ir]].termana[cavailles[ic].it_chain_list[it]].wheretofind[iwhere].whichfuncchain << "-" << (*regtype)[bagtype->which_regtype[ir]].termana[cavailles[ic].it_chain_list[it]].wheretofind[iwhere].whichfunc << "\t";
						writezc.flush();
					}
				}
				writezc << "\n";
				writezc.flush();
			}
		}

		writezc << "===================================================================================\n";
		writezc.flush();
		writezc << "v\t";
		writezc.flush();
		for (int it = 0; it < (signed)it_list.size(); it++) {
			for (int itn = 0; itn < (signed)bagtype->termanaforbag[it_list[it]].name.size(); itn++) {
				writezc << bagtype->termanaforbag[it_list[it]].name[itn];
				writezc.flush();
				if (itn != (signed)bagtype->termanaforbag[it_list[it]].name.size() - 1) {
					writezc << "-";
					writezc.flush();
				}
			}
			if (it != (signed)it_list.size() - 1) {
				writezc << "\t";
				writezc.flush();
			}
		}
		writezc << "\n";
		writezc.flush();
	}

	for (int itest = 0; itest <= 100; itest++) {
		std::vector<long double> test_anades_degree;
		long double test_OM_degree = (long double)itest / 100.0;
		long double test_CO_degree = (long double)itest / 100.0;
		if (pr_id == 0) {
			writezc << (long double)itest / 100.0 << "\t";
			writezc.flush();
		}
		for (int it = 0; it < (signed)it_list.size(); it++) {
			std::vector<long double> v_series;
			std::vector<long double> v_series_acc;
			v_series.clear();
			v_series_acc.clear();
			long double v_stdev;
			bool post_if_exists = true;
			MAT nullMAT;
			test_anades_degree.clear();
			int degree_count = (signed)bagtype->termanaforbag[it_list[it]].name.size();
			for (int in = 0; in < (signed)bagtype->termanaforbag[it_list[it]].name.size(); in++) {
				if (bagtype->termanaforbag[it_list[it]].name[in] == *AnionName + "M") {
					degree_count--;
				}
				else {
					bool found_CO = false;
					for (int ic = 0; ic < (signed)CationNames->size(); ic++) {
						if (bagtype->termanaforbag[it_list[it]].name[in] == (*CationNames)[ic] + *AnionName) {
							found_CO = true;
							break;
						}
					}
					if (found_CO) {
						degree_count--;
					}
				}
			}
			for (int idc = 0; idc < degree_count; idc++) {
				test_anades_degree.push_back((long double)itest / 100.0);
			}
			calc_bag_chain_partial_test(&v_series, &v_series_acc, &v_stdev, regtype, bagtype, &nullMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, AnionName, CationNames, &bagtype->termanaforbag[it_list[it]].name, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly);
			if (pr_id == 0) {
				writezc << v_series_acc[v_series_acc.size() - 1];
				writezc.flush();
				if (it != (signed)it_list.size() - 1) {
					writezc << "\t";
					writezc.flush();
				}
			}
			std::vector<long double>().swap(v_series);
			std::vector<long double>().swap(v_series_acc);
		}
		if (itest != 100) {
			if (pr_id == 0) {
				writezc << "\n";
				writezc.flush();
			}
		}
		std::vector<long double>().swap(test_anades_degree);
	}

	for (int ic = 0; ic < (signed)cavailles.size(); ic++) {
		for (int it = 0; it < (signed)it_list.size(); it++) {

			bool writeable = true;
			if (*verbose) {
				writeable = true;
			}
			else {
				writeable = false;
				for (int jt = 0; jt < (signed)cavailles[ic].it_chain_list.size(); jt++) {
					if (it_list[it] == cavailles[ic].it_chain_list[jt]) {
						writeable = true;
						break;
					}
				}
			}

			if (writeable) {

				std::string subfilename = "";
				for (int in = 0; in < (signed)bagtype->termanaforbag[it_list[it]].name.size(); in++) {
					subfilename += bagtype->termanaforbag[it_list[it]].name[in];
					if (in != (signed)bagtype->termanaforbag[it_list[it]].name.size() - 1) {
						subfilename += "-";
					}
				}
				if (pr_id == 0) {
					writezcext.open("_bag" + std::to_string(*ibag) + "_" + bagtype->name + "." + subfilename + ".txt", std::ios::binary);
					writezcext << std::fixed << std::setprecision(20);
				}
				get_floating_v = true;
				std::vector<int> point_anades;
				if (pr_id == 0) {
					for (int in = 0; in < (signed)bagtype->termanaforbag[it_list[it]].name.size(); in++) {
						writezcext << "v" << bagtype->termanaforbag[it_list[it]].name[in] << "\t";
						writezcext.flush();
					}
					for (int in = 0; in < (signed)bagtype->termanaforbag[it_list[it]].name.size(); in++) {
						writezcext << bagtype->termanaforbag[it_list[it]].name[in] << "\t";
						writezcext.flush();
					}
					writezcext << bagtype->name << "\tModelStdDv\tNsample\tMin\tAve\t+-\tMax\tNsample(R)\tMin(R)\tAve(R)\t+-(R)\tMax(R)\tMeanAE\tMaxAE\n";
					writezcext.flush();
				}
				std::vector<int> tcolumns;
				tcolumns.clear();
				tcolumns.push_back(-1);
				for (int in = 0; in < (signed)bagtype->termanaforbag[it_list[it]].name.size(); in++) {
					tcolumns.push_back(-1);
				}
				readanalysis.open(*AnalysisOutputFilename);
				std::getline(readanalysis, sbuf);
				std::stringstream linebuf;
				linebuf.str("");
				linebuf.clear();
				linebuf << sbuf;
				int icol = -1;
				while (true) {
					linebuf >> sbuf;
					icol++;
					if (sbuf == bagtype->name) {
						tcolumns[0] = icol;
					}
					for (int in = 0; in < (signed)bagtype->termanaforbag[it_list[it]].name.size(); in++) {
						if ((!*AveOnly && sbuf == bagtype->termanaforbag[it_list[it]].name[in]) || (*AveOnly && sbuf == "Ave" + bagtype->termanaforbag[it_list[it]].name[in])) {
							tcolumns[in + 1] = icol;
							break;
						}
					}
					if (linebuf.eof()) {
						break;
					}
				}
				readanalysis.close();

				int tcolumnR = -1;
				readanalysisR.open("__R__." + * AnalysisOutputFilename);
				std::getline(readanalysisR, sbuf);
				linebuf.str("");
				linebuf.clear();
				linebuf << sbuf;
				icol = -1;
				while (true) {
					linebuf >> sbuf;
					icol++;
					if (sbuf == "bag" + std::to_string(*ibag) + "_" + bagtype->name) {
						tcolumnR = icol;
						break;
					}
					if (linebuf.eof()) {
						break;
					}
				}
				readanalysisR.close();

				int imin = 0;
				int imax = *AnalysisMesh;
				std::vector<int> occs;
				for (int iocc = 0; iocc < (signed)bagtype->termanaforbag[it_list[it]].name.size(); iocc++) {
					occs.push_back(-1);
				}
				int init_pos = 0;
				generate_tuples(&occs, init_pos, &imin, &imax, &bagtype->termanaforbag[it_list[it]].name, bagtype, ibag, regtype, RegModelFilename, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, AnionName, CationNames, verbose, AnalysisOutputFilename, &tcolumns, &tcolumnR, AveOnly);
				get_floating_v = false;
				if (pr_id == 0) {
					writezcext.close();
				}
				std::vector<int>().swap(tcolumns);

			}

		}
	}

	std::vector<int>().swap(it_list);
	std::vector<CAVAILLES>().swap(cavailles);
	if (pr_id == 0) {
		writezc.close();
	}

}

void set_par(std::string* MaterialsFilename, bool *headexists) {

	readmatdat.open(*MaterialsFilename);

	if (readmatdat.is_open()) {

		if (*headexists) {
			std::getline(readmatdat, sbuf);
		}

		Nmat = 0;
		while (true) {
			std::getline(readmatdat, sbuf);
			if (sbuf.find("STOP") != std::string::npos) {
				break;
			}
			else {
				Nmat++;
			}
			if (readmatdat.eof()) {
				break;
			}
		}

		readmatdat.close();

		if (sbuf == "") {
			Nmat--;
		}

		if (pr_id == 0) {
			if (Nmat <= 1) {
				writeoutput << Nmat << " material is found...\n";
				writeoutput.flush();
			}
			else {
				writeoutput << Nmat << " materials are found...\n";
				writeoutput.flush();
			}
			writeoutput << "=================================================================================\n";
			writeoutput.flush();
		}

		pr_job_size = new signed long long int[pr_size];
		pr_job = new signed long long int* [pr_size];
		for (int i = 0; i < pr_size; i++) {
			pr_job_size[i] = 0;
			pr_job[i] = new signed long long int[2];
			for (int j = 0; j < 2; j++) {
				pr_job[i][j] = -1;
			}
		}

		int pr_pointer = 0;
		for (signed long long int i = 0; i < Nmat; i++) {
			pr_job_size[pr_pointer]++;
			pr_pointer++;
			if (pr_pointer == pr_size) {
				pr_pointer = 0;
			}
		}

		for (int i = 0; i < pr_size; i++) {
			if (pr_job_size[i] > 0) {
				if (i == 0) {
					pr_job[i][0] = 0;
				}
				else {
					pr_job[i][0] = pr_job[i - 1][1] + 1;
				}
				pr_job[i][1] = pr_job[i][0] + pr_job_size[i] - 1;
			}
		}

		if (pr_id == 0) {
			writeoutput << "Materials Design Job Distribution:\n";
			writeoutput.flush();
			writeoutput << "pr_id\tpr_job[i][0]\tpr_job[i][1]\n";
			writeoutput.flush();
			for (int i = 0; i < pr_size; i++) {
				if (pr_job[i][0] == -1 || pr_job[i][1] == -1) {
					writeoutput << i << "\t" << pr_job[i][0] << "\t" << pr_job[i][1] << "\n";
					writeoutput.flush();
				}
			}
			writeoutput << "=================================================================================\n";
			writeoutput.flush();
		}

	}

}

void designmate_iatom(MAT * mat, std::vector<ATOMDAT> * atomdat) {

	std::vector<int>().swap(mat->iatom);

	for (int i = 0; i < mat->atomname.size(); i++) {
		int iatom = -1;
		for (int l = 0; l < atomdat->size(); l++) {
			if ((*atomdat)[l].atomname == mat->atomname[i]) {
				iatom = l;
				break;
			}
		}
		mat->iatom.push_back(iatom);
	}
	
}

void designmate_isite(MAT * mat, std::vector<std::string> * sitename) {

	for (int i = 0; i < mat->sitename.size(); i++) {
		int isite = -1;
		for (int l = 0; l < sitename->size(); l++) {
			if ((*sitename)[l] == mat->sitename[i]) {
				isite = l;
				break;
			}
		}
		mat->isite.push_back(isite);
	}

}

void write_zdist(MAT * mat, long double * zdist, std::string* AnionName) {

	writeoutputpr << mat->composition << "\t";
	writeoutputpr.flush();
	for (int i = 0; i < (signed)mat->atomname.size(); i++) {
		writeoutputpr << mat->atomname[i] << " " << mat->content[i];
		writeoutputpr.flush();
		if (i != (signed)mat->atomname.size() - 1) {
			writeoutputpr << " ";
			writeoutputpr.flush();
		}
	}
	if (*AnionName != "-") {
		writeoutputpr << " " << * AnionName << " " << mat->anion_content;
		writeoutputpr.flush();
	}
	writeoutputpr << "\t";
	writeoutputpr.flush();
	writeoutputpr << "zDist=\t" << *zdist;
	writeoutputpr.flush();

}

void dope_mat(MAT* objMAT, long double* dope) {

	std::vector<std::string> r_atomname;
	std::vector<int> r_iatom;
	std::vector<std::string> r_sitename;
	std::vector<int> r_isite;
	std::vector<long double> r_content;

	for (int iiatom = 0; iiatom < (signed)objMAT->atomname.size(); iiatom++) {
		r_atomname.push_back(objMAT->atomname[iiatom]);
		r_iatom.push_back(objMAT->iatom[iiatom]);
		r_sitename.push_back(objMAT->sitename[iiatom]);
		r_isite.push_back(objMAT->isite[iiatom]);
		r_content.push_back(objMAT->content[iiatom]);
	}

	std::vector<std::string>().swap(objMAT->atomname);
	objMAT->atomname.clear();
	std::vector<int>().swap(objMAT->iatom);
	objMAT->iatom.clear();
	std::vector<std::string>().swap(objMAT->sitename);
	objMAT->sitename.clear();
	std::vector<int>().swap(objMAT->isite);
	objMAT->isite.clear();
	std::vector<long double>().swap(objMAT->content);
	objMAT->content.clear();

	for (int iiatom = 0; iiatom < (signed)r_atomname.size(); iiatom++) {
		objMAT->atomname.push_back(r_atomname[iiatom]);
		objMAT->atomname.push_back(r_atomname[iiatom]);
		objMAT->iatom.push_back(r_iatom[iiatom]);
		objMAT->iatom.push_back(r_iatom[iiatom]);
		objMAT->sitename.push_back(r_sitename[iiatom]);
		objMAT->sitename.push_back(r_sitename[iiatom]);
		objMAT->isite.push_back(r_isite[iiatom]);
		objMAT->isite.push_back(r_isite[iiatom]);
		objMAT->content.push_back((1.0 - *dope) * r_content[iiatom]);
		objMAT->content.push_back(*dope * r_content[iiatom]);
	}

}

void process_matdat(std::vector<long double>* zTargetVReg, std::vector<long double>* ZDistWeight, std::vector<long double>* TargetVReg, std::vector<long double>* TargetStdDvReg, std::vector<long double>* ModulateContentFactor, std::vector<long double>* DopeFactor, std::string* MaterialsFilename, std::vector<REGTYPE>* regtype, std::vector<BAGTYPE>* bagtype, std::vector<ATOMDAT>* atomdat, std::string* FaceFilename, std::vector<std::string>* FaceFeatureName, std::vector<FACE>* face, std::vector<STRDAT>* strdat, std::vector<std::string>* sitename, std::vector<ANADESCRIPTOR>* anadescriptor, std::string* AnionName, std::vector<ANIONS>* Anions, long double* modulation_limit, std::vector<std::string>* CationNames, long double* AnionValence, long double* VInfinite, bool* collect_only, bool* UseOffset, int* MaxIteration, bool* headexists, bool* usestrname, std::vector<std::string>* PostTagNames, std::vector<std::string>* DoNotAnalyze, bool* verbose, bool* AveOnly, bool* FoldMode, int* NFolds, int* TakeFold, int* Train0Test1, int* Straight0Rotate1, long double* TestModelStDevLimit) {

	readmatdat.open(*MaterialsFilename);
	
	if (readmatdat.is_open()) {

		if (*headexists) {
			std::getline(readmatdat, sbuf);
		}

		MAT gMAT;
		std::vector<long double> gmonitor_value, gmonitor_stddv;
		if (!*collect_only) {
			for (int k = 0; k < (signed)bagtype->size(); k++) {
				gmonitor_value.push_back(0.0);
				gmonitor_stddv.push_back(0.0);
			}
		}
		long double global_min_zdist = *VInfinite;
		std::string gstr, pgstr;

		signed long long int Nstart = 0;
		signed long long int Nend = Nmat - 1;
		bool run_matdat = true;
		if (!*collect_only) {
			Nstart = pr_job[pr_id][0];
			Nend = pr_job[pr_id][1];
			if (Nstart == -1 || Nend == -1) {
				run_matdat = false;
			}
			if (run_matdat) {
				writeoutputpr.open("MD." + std::to_string(pr_id) + "." + OutputFilename, std::ios::binary);
				writeoutputpr << std::fixed << std::setprecision(20);
			}
		}

		if (Nstart != 0) {
			for (signed long long int i = 0; i <= Nstart - 1; i++) {
				std::getline(readmatdat, sbuf);
			}
		}

		for (signed long long int i = Nstart; i <= Nend && run_matdat; i++) {

			std::getline(readmatdat, sbuf);
			std::string linebuf = sbuf;
			int WhichFold = -1;
			bool ReadDatumInFoldMode = false;
			if (*FoldMode) {
				if (*Straight0Rotate1 == 0) {
					WhichFold = (int)((long double)(*NFolds * i) / (long double)Nmat) + 1;
				}
				else {
					WhichFold = i % *NFolds + 1;
				}
				if (*Train0Test1 == 0 && WhichFold != *TakeFold) {
					ReadDatumInFoldMode = true;
				}
				else if (*Train0Test1 == 1 && WhichFold == *TakeFold) {
					ReadDatumInFoldMode = true;
				}
			}

			if ((!*FoldMode && i >= 0 && i <= Nmat - 1) || (*FoldMode && ReadDatumInFoldMode)) {

				std::vector<std::string> PostTags;
				if(PostTagNames->size() > 0){ 
					for (int ipost = 0; ipost < (signed)PostTagNames->size(); ipost++) {
						PostTags.push_back("");
					}
					std::stringstream linestream;
					linestream << linebuf;
					std::string subbuf;
					while (true) {
						linestream >> subbuf;
						if (PostTags.size() > 1) {
							for (unsigned int ii = 1; ii < PostTags.size(); ii++) {
								PostTags[ii - 1] = PostTags[ii];
							}
						}
						PostTags[PostTags.size() - 1] = subbuf;
						if (linestream.eof()) {
							break;
						}
					}
				}

				MAT objMAT;
				std::stringstream lcstr;
				lcstr << sbuf;
				lcstr >> objMAT.name;
				lcstr >> sbuf;
				lcstr >> objMAT.strname;
				if (!*usestrname) {
					objMAT.strname = "None";
				}
				objMAT.istr = -1;
				for (int j = 0; j < (signed)strdat->size(); j++) {
					if (objMAT.strname == (*strdat)[j].strname) {
						objMAT.istr = j;
						break;
					}
				}
				while (true) {
					lcstr >> sbuf;
					if (sbuf == "||") {
						break;
					}
				}
				objMAT.composition = "";
				long double cation_sum = 0.0;
				long double* scation_sum;
				scation_sum = new long double[CationNames->size()];
				for (int j = 0; j < (signed)CationNames->size(); j++) {
					scation_sum[j] = 0.0;
				}
				for (int j = 0; j < (signed)sitename->size(); j++) {
					while (true) {
						lcstr >> sbuf;
						if (sbuf != "||") {
							objMAT.atomname.push_back(sbuf);
							objMAT.sitename.push_back((*sitename)[j]);
							lcstr >> ldbuf;
							objMAT.content.push_back(ldbuf);
							cation_sum += ldbuf;
							for (int k = 0; k < (signed)CationNames->size(); k++) {
								if (objMAT.atomname[objMAT.atomname.size() - 1] == (*CationNames)[k]) {
									scation_sum[k] += ldbuf;
									break;
								}
							}
							objMAT.composition += sbuf;
							if (ldbuf != 1.0) {
								objMAT.composition += std::to_string(ldbuf).substr(0, 5);
							}
						}
						else {
							break;
						}
					}
				}
				designmate_iatom(&objMAT, atomdat);
				designmate_isite(&objMAT, sitename);
				if (*AnionName != "-") {
					lcstr >> sbuf >> objMAT.anion_content;
					objMAT.composition += sbuf;
					if (objMAT.anion_content != 1.0) {
						objMAT.composition += std::to_string(objMAT.anion_content).substr(0, 5);
					}
					objMAT.OM = objMAT.anion_content / cation_sum;
					for (int k = 0; k < (signed)CationNames->size(); k++) {
						objMAT.CO.push_back(scation_sum[k] / objMAT.anion_content);
					}
				}
				else {
					objMAT.OM = 0.0;
				}
				delete[] scation_sum;
				scation_sum = nullptr;

				if (*collect_only) {
					if (pr_id == 0) {
						writeoutput << i + 1 << "\t...\t" << objMAT.name << "\t...\t" << objMAT.composition << "\t...\t" << objMAT.strname << "\n";
						writeoutput.flush();
						writeanalysis << i + 1 << "\t" << objMAT.name << "\t" << objMAT.composition << "\t" << objMAT.strname << "\t";
						writeanalysis.flush();
						if (*AnionName != "-") {
							writeanalysis << objMAT.OM << "\t";
							writeanalysis.flush();
							for (int k = 0; k < (signed)CationNames->size(); k++) {
								writeanalysis << objMAT.CO[k] << "\t";
								writeanalysis.flush();
							}
						}
					}
				}
				else {
					writeoutputpr << i + 1 << "\t...\t" << objMAT.name << "\t...\t" << objMAT.composition << "\t...\t" << objMAT.strname << "\n";
					writeoutputpr.flush();
				}
				
				if (*collect_only) {	

					if (pr_id == 0) {
						if (*AnionName != "-") {
							writeoutput << *AnionName << "M\t" << objMAT.OM << "\n";
							writeoutput.flush();
							for (int k = 0; k < (signed)CationNames->size(); k++) {
								writeoutput << (*CationNames)[k] << *AnionName << "\t" << objMAT.CO[k];
								writeoutput.flush();
								if (k != (signed)CationNames->size() - 1) {
									writeoutput << "\t";
									writeoutput.flush();
								}
								else {
									writeoutput << "\n";
									writeoutput.flush();
								}
							}
						}
					}

					if (!std::isnan(objMAT.OM)) {
						if (objMAT.OM < OM_min) {
							OM_min = objMAT.OM;
						}
						OM_ave += objMAT.OM;
						if (objMAT.OM > OM_max) {
							OM_max = objMAT.OM;
						}
						OM_count++;
					}
					for (int k = 0; k < (signed)CationNames->size(); k++) {
						if (!std::isnan(objMAT.CO[k])) {
							if (objMAT.CO[k] < CO_min[k]) {
								CO_min[k] = objMAT.CO[k];
							}
							CO_ave[k] += objMAT.CO[k];
							if (objMAT.CO[k] > CO_max[k]) {
								CO_max[k] = objMAT.CO[k];
							}
							CO_count[k]++;
						}
					}

					if (pr_id == 0) {
						writeoutput << "Atomname\t";
						writeoutput.flush();
						for (int j = 0; j < (signed)objMAT.atomname.size(); j++) {
							writeoutput << objMAT.atomname[j];
							writeoutput.flush();
							if (j != (signed)objMAT.atomname.size() - 1) {
								writeoutput << "\t";
								writeoutput.flush();
							}
							else {
								writeoutput << "\n";
								writeoutput.flush();
							}
						}
						writeoutput << "Content\t";
						writeoutput.flush();
						for (int j = 0; j < (signed)objMAT.content.size(); j++) {
							writeoutput << objMAT.content[j];
							writeoutput.flush();
							if (j != (signed)objMAT.content.size() - 1) {
								writeoutput << "\t";
								writeoutput.flush();
							}
							else {
								writeoutput << "\n";
								writeoutput.flush();
							}
						}
					}

				}

				if (*FaceFilename != "") {
					
					std::string found_face = "";
					std::vector<long double> vface;
					for (int j = 0; j < (signed)FaceFeatureName->size(); j++) {
						vface.push_back(0.0);
					}

					long double find_max_content = -1000000000000000.0;
					for (int jj = 0; jj < (signed)objMAT.atomname.size(); jj++) {
						for (int j = 0; j < (signed)face->size(); j++) {
							if (objMAT.atomname[jj] == (*face)[j].AtomName && objMAT.content[jj] > find_max_content && (*face)[j].face) {
								found_face = (*face)[j].AtomName;
								find_max_content = objMAT.content[jj];
								for (int iv = 0; iv < (signed)FaceFeatureName->size(); iv++) {
									vface[iv] = (*face)[j].vFeature[iv];
								}
								break;
							}
						}
					}
					
					if (*collect_only && pr_id == 0) {
						if (found_face != "") {

							writeoutput << "FACE_at\t" << found_face << "\n";
							writeoutput.flush();
							for (int iv = 0; iv < (signed)FaceFeatureName->size(); iv++) {
								writeoutput << (*FaceFeatureName)[iv];
								writeoutput.flush();
								if (iv != (signed)FaceFeatureName->size() - 1) {
									writeoutput << "\t";
									writeoutput.flush();
								}
								else {
									writeoutput << "\n";
									writeoutput.flush();
								}
							}
							for (int iv = 0; iv < (signed)vface.size(); iv++) {
								writeoutput << vface[iv] << "\t";
								writeoutput.flush();
							}

						}

						for (int iv = 0; iv < (signed)vface.size(); iv++) {
							writeanalysis << vface[iv] << "\t";
							writeanalysis.flush();
						}
					}

					std::vector<long double>().swap(vface);

				}

				if (objMAT.istr != -1) {
					int notUlastcolumn = (signed)anadescriptor->size() - 1;
					for (int j = 0; j < (signed)anadescriptor->size(); j++) {
						if ((*anadescriptor)[j].ASB != "U") {
							notUlastcolumn = j;
						}
					}
					for (int j = 0; j < (signed)anadescriptor->size(); j++) {
						
						if ((*anadescriptor)[j].ASB != "U") {
							std::vector<long double> dv, mamssk;
							get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j, AveOnly);
							get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze, AveOnly);
							if (*collect_only) {
								if (pr_id == 0) {
									writeoutput << (*anadescriptor)[j].name << "\t";
									writeoutput.flush();
								}
								for (int k = 0; k < (signed)dv.size(); k++) {
									if (pr_id == 0) {
										writeoutput << dv[k];
										writeoutput.flush();
									}
									if (k != (signed)dv.size() - 1) {
										if (pr_id == 0) {
											writeoutput << "\t";
											writeoutput.flush();
										}
									}
									else {
										if (pr_id == 0) {
											writeoutput << "\t***\tMin=\t" << mamssk[0] << "\tAve=\t" << mamssk[1] << "\tMax=\t" << mamssk[2] << "\tStdDv=\t" << mamssk[3] << "\tSkew=\t" << mamssk[4] << "\tKurto=\t" << mamssk[5];
											writeoutput.flush();
											writeanalysis << mamssk[0] << "\t" << mamssk[1] << "\t" << mamssk[2] << "\t" << mamssk[3] << "\t" << mamssk[4] << "\t" << mamssk[5];
											writeanalysis.flush();
											if (j != notUlastcolumn) {
												writeanalysis << "\t";
												writeanalysis.flush();
											}
										}
										
										if (!std::isnan(mamssk[0]) && !std::isnan(mamssk[1]) && !std::isnan(mamssk[2]) && !std::isnan(mamssk[3]) && !std::isnan(mamssk[4]) && !std::isnan(mamssk[5])) {
											if (mamssk[0] < (*anadescriptor)[j].min) {
												(*anadescriptor)[j].min = mamssk[0];
											}
											if (mamssk[1] < (*anadescriptor)[j].min_ave) {
												(*anadescriptor)[j].min_ave = mamssk[1];
											}
											(*anadescriptor)[j].ave += mamssk[1];
											if (mamssk[1] > (*anadescriptor)[j].max_ave) {
												(*anadescriptor)[j].max_ave = mamssk[1];
											}
											if (mamssk[2] > (*anadescriptor)[j].max) {
												(*anadescriptor)[j].max = mamssk[2];
											}
											if (mamssk[3] < (*anadescriptor)[j].min_stdev) {
												(*anadescriptor)[j].min_stdev = mamssk[3];
											}
											(*anadescriptor)[j].ave_stdev += mamssk[3];
											if (mamssk[3] > (*anadescriptor)[j].max_stdev) {
												(*anadescriptor)[j].max_stdev = mamssk[3];
											}
											if (mamssk[4] < (*anadescriptor)[j].min_skew) {
												(*anadescriptor)[j].min_skew = mamssk[4];
											}
											(*anadescriptor)[j].ave_skew += mamssk[4];
											if (mamssk[4] > (*anadescriptor)[j].max_skew) {
												(*anadescriptor)[j].max_skew = mamssk[4];
											}
											if (mamssk[5] < (*anadescriptor)[j].min_kurto) {
												(*anadescriptor)[j].min_kurto = mamssk[5];
											}
											(*anadescriptor)[j].ave_kurto += mamssk[5];
											if (mamssk[5] > (*anadescriptor)[j].max_kurto) {
												(*anadescriptor)[j].max_kurto = mamssk[5];
											}
											(*anadescriptor)[j].count++;
										}
									}
								}
							}
							for (int k = 0; k < 6; k++) {
								xbox[j][k] = mamssk[k];
							}
							std::vector<long double>().swap(dv);
							std::vector<long double>().swap(mamssk);
						}
						else {
							long double v = 0.0;
							long double v0 = 0.0;
							if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] == (*anadescriptor)[j].UType_numerator_AtomKinds[ja]) {
											v0 += objMAT.content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									v0 += objMAT.content[ia];
								}
							}
							else if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									bool found = false;
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_numerator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] == (*anadescriptor)[j].UType_numerator_AtomKinds[ja]) {
											found = true;
											break;
										}
									}
									if (!found) {
										v0 += objMAT.content[ia];
									}
								}
							}
							else if ((*anadescriptor)[j].UType_numerator_S0A1E2N3 == 3) {
								v0 = (*anadescriptor)[j].UType_numerator_value;
							}
							long double v1 = 0.0;
							if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 0) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] == (*anadescriptor)[j].UType_denominator_AtomKinds[ja]) {
											v1 += objMAT.content[ia];
											break;
										}
									}
								}
							}
							else if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 1) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									v1 += objMAT.content[ia];
								}
							}
							else if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 2) {
								for (int ia = 0; ia < (signed)objMAT.content.size(); ia++) {
									bool found = false;
									for (int ja = 0; ja < (signed)(*anadescriptor)[j].UType_denominator_AtomKinds.size(); ja++) {
										if (objMAT.atomname[ia] == (*anadescriptor)[j].UType_denominator_AtomKinds[ja]) {
											found = true;
											break;
										}
									}
									if (!found) {
										v1 += objMAT.content[ia];
									}
								}
							}
							else if ((*anadescriptor)[j].UType_denominator_S0A1E2N3 == 3) {
								v1 = (*anadescriptor)[j].UType_denominator_value;
							}
							if ((*anadescriptor)[j].UType_m0d1 == 0) {
								v = v0 * v1;
							}
							else {
								v = v0 / v1;
							}
							if (*collect_only) {
								if (pr_id == 0) {
									writeoutput << (*anadescriptor)[j].name << "\t" << v;
									writeoutput.flush();
									writeanalysis << "\t" << v;
									writeanalysis.flush();
								}
								if(!std::isnan(v)) {
									if (v < (*anadescriptor)[j].min) {
										(*anadescriptor)[j].min = v;
										
									}
									if (v < (*anadescriptor)[j].min_ave) {
										(*anadescriptor)[j].min_ave = v;
									}
									(*anadescriptor)[j].ave += v;
									if (v > (*anadescriptor)[j].max_ave) {
										(*anadescriptor)[j].max_ave = v;
									}
									if (v > (*anadescriptor)[j].max) {
										(*anadescriptor)[j].max = v;
									}
									(*anadescriptor)[j].min_stdev = 0.0;
									(*anadescriptor)[j].ave_stdev += 0.0;
									(*anadescriptor)[j].max_stdev = 0.0;
									(*anadescriptor)[j].min_skew = 0.0;
									(*anadescriptor)[j].ave_skew += 0.0;
									(*anadescriptor)[j].max_skew = 0.0;
									(*anadescriptor)[j].min_kurto = *VInfinite;
									(*anadescriptor)[j].ave_kurto += *VInfinite;
									(*anadescriptor)[j].max_kurto = *VInfinite;
									(*anadescriptor)[j].count++;
								}
							}
						}
						if (*collect_only) {
							if (pr_id == 0) {
								writeoutput << "\n";
								writeoutput.flush();
							}
						}
					}

					std::vector<long double> regv, regv2, bagv, bagstdev;
					regv.clear();
					regv2.clear();
					bagv.clear();
					bagstdev.clear();
					if (*collect_only) {

						if (regtype->size() > 0 && bagtype->size()) {
							for (int k = 0; k < (signed)bagtype->size(); k++) {
								std::vector<int> test_anades;
								std::vector<int> test_mamssk;
								std::vector<long double> test_anades_degree;
								long double test_OM_degree = 0.0;
								long double test_CO_degree = 0.0;
								bool test_OM = false;
								int test_CO = -100;
								bool post_if_exists = true;
								std::vector<long double> vacc;
								vacc.clear();
								long double v_stdev;
								long double post_v = 0.0;
								bagv.push_back(calc_bag_chain(&regv2, &vacc, &v_stdev, regtype, &(*bagtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v));
								bagstdev.push_back(v_stdev);
								if (!(*FoldMode && *Train0Test1 == 1 && *TestModelStDevLimit > 0.0 && v_stdev > *TestModelStDevLimit)) {
									(*bagtype)[k].regressed_v_list.push_back(bagv[bagv.size() - 1]);
								}
								else {
									(*bagtype)[k].regressed_v_list.push_back(std::nanl(""));
								}
								(*bagtype)[k].regressed_stddv_v_list.push_back(0.0);
								(*bagtype)[k].regressed_post_v_list.push_back(post_v);
								for (int kk = 0; kk < (signed)regv2.size(); kk++) {
									regv.push_back(regv2[kk]);
									(*regtype)[(*bagtype)[k].which_regtype[kk]].regressed_v_list.push_back(regv2[kk]);
								}
								std::vector<long double>().swap(vacc);
								std::vector<int>().swap(test_anades);
								std::vector<int>().swap(test_mamssk);
								std::vector<long double>().swap(test_anades_degree);
							}
						}
						if (PostTags.size() > 0) {
							for (int ipost = 0; ipost < (signed)PostTags.size(); ipost++) {
								if (pr_id == 0) {
									writeanalysis << "\t" << PostTags[ipost];
									writeanalysis.flush();
								}
							}
						}
						if (pr_id == 0) {
							writeanalysis << "\n";
							writeanalysis.flush();
						}
					}

					xbox_OM = 0.0;
					if (*AnionName != "-") {
						xbox_OM = objMAT.OM;
						for (int k = 0; k < (signed)CationNames->size(); k++) {
							xbox_CO[k] = objMAT.CO[k];
						}
					}
					if ((signed)bagtype->size() > 0) {
						if (*collect_only) {

							if (pr_id == 0) {
								writeoutput << "!!!\t";
								writeoutput.flush();
							}
							if (*UseOffset) {
								lcstr >> sbuf;
							}
							for (int k = 0; k < (signed)bagtype->size(); k++) {
								long double value = 0.0; 
								if (!*UseOffset) {
									value = bagv[k];
								}
								else {
									lcstr >> value;
									std::string target_ref_v_str = "";
									bool target_ref_v_found = false;
									for (int ipost = 0; ipost < (signed)PostTagNames->size(); ipost++) {
										if ((*PostTagNames)[ipost] == (*bagtype)[k].name) {
											target_ref_v_str = PostTags[ipost];
											target_ref_v_found = true;
											break;
										}
									}
									if (target_ref_v_found) {
										std::stringstream glcstr;
										glcstr.str("");
										glcstr.clear();
										glcstr << target_ref_v_str;
										glcstr >> value;
									}
								}
								(*bagtype)[k].v_list.push_back(value);
								for (int kk = 0; kk < (signed)(*bagtype)[k].which_regtype.size(); kk++) {
									(*regtype)[(*bagtype)[k].which_regtype[kk]].v_list.push_back(value);
								}
								if (pr_id == 0) {
									writeoutput << (*bagtype)[k].name << "=\t" << value;
									writeoutput.flush();
									if (!*UseOffset) {
										writeoutput << "\t+-\t" << bagstdev[k];
										writeoutput.flush();
									}
								}
								if (value < (*bagtype)[k].obs_v_min) {
									(*bagtype)[k].obs_v_min = value;
								}
								if (value > (*bagtype)[k].obs_v_max) {
									(*bagtype)[k].obs_v_max = value;
								}
								if (pr_id == 0) {
									if (k != (signed)bagtype->size() - 1) {
										writeoutput << "\t";
										writeoutput.flush();
									}
									else {
										writeoutput << "\n";
										writeoutput.flush();
									}
								}
							}

							std::vector<long double>().swap(regv);
							std::vector<long double>().swap(regv2);
							std::vector<long double>().swap(bagv);
							std::vector<long double>().swap(bagstdev);

						}
						else {

							MAT bufMAT = objMAT;

							for (int idope = 0; idope < (signed)DopeFactor->size() + 1; idope++) {

								if (idope > 0) {
									objMAT = bufMAT;
									dope_mat(&objMAT, &(*DopeFactor)[idope - 1]);
								}
							
								std::vector<long double> init_zvalue;
								std::vector<long double> monitor_value;
								std::vector<long double> monitor_stddv;
								std::vector<long double> target_monitor_value;
								std::vector<long double> target_monitor_stddv;
								std::vector<long double> init_monitor_value;
								std::vector<long double> init_monitor_stddv;
								std::vector<long double> offset;
								init_zvalue.clear();
								monitor_value.clear();
								monitor_stddv.clear();
								target_monitor_value.clear();
								target_monitor_stddv.clear();
								init_monitor_value.clear();
								init_monitor_stddv.clear();
								offset.clear();
								long double init_zdist = 0.0;
								writeoutputpr << "Init:\t";
								writeoutputpr.flush();
								pgstr = "Init:\t";
								if (*UseOffset) {
									lcstr >> sbuf;
								}
								for (int kk = 0; kk < (signed)bagtype->size(); kk++) {
									std::vector<int> test_anades;
									std::vector<int> test_mamssk;
									std::vector<long double> test_anades_degree;
									long double test_OM_degree = 0.0;
									long double test_CO_degree = 0.0;
									bool test_OM = false;
									int test_CO = -100;
									bool post_if_exists = true;
									std::vector<long double> rv, rvacc;
									long double v_stdev;
									long double post_v = 0.0;
									long double init_value = calc_bag_chain(&rv, &rvacc, &v_stdev, regtype, &(*bagtype)[kk], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v);
									std::vector<long double>().swap(rv);
									std::vector<long double>().swap(rvacc);
									std::vector<int>().swap(test_anades);
									std::vector<int>().swap(test_mamssk);
									std::vector<long double>().swap(test_anades_degree);
									long double offset_value = 0.0;
									if (*UseOffset) {
										lcstr >> ldbuf;
										std::string target_ref_v_str = "";
										bool target_ref_v_found = false;
										for (int ipost = 0; ipost < (signed)PostTagNames->size(); ipost++) {
											if ((*PostTagNames)[ipost] == (*bagtype)[kk].name) {
												target_ref_v_str = PostTags[ipost];
												target_ref_v_found = true;
												break;
											}
										}
										if (target_ref_v_found) {
											std::stringstream glcstr;
											glcstr.str("");
											glcstr.clear();
											glcstr << target_ref_v_str;
											glcstr >> ldbuf;
										}
										offset_value = ldbuf - init_value;
									}
									long double inz = (init_value + offset_value - (*bagtype)[kk].v_ave) / (*bagtype)[kk].v_stddv;
									init_zvalue.push_back(inz);
									monitor_value.push_back(init_value + offset_value);
									target_monitor_value.push_back(init_value + offset_value);
									init_monitor_value.push_back(init_value + offset_value);
									monitor_stddv.push_back(v_stdev);
									target_monitor_stddv.push_back(v_stdev);
									init_monitor_stddv.push_back(v_stdev);
									offset.push_back(offset_value);
									init_zdist += (*ZDistWeight)[kk] * powl(inz - (*zTargetVReg)[kk], 2.0);
								}
								init_zdist = sqrtl(init_zdist);
								writeoutputpr << "zDist=\t" << init_zdist << "\t";
								writeoutputpr.flush();
								pgstr += "zDist=\t" + std::to_string(init_zdist) + "\t";
								for (int k = 0; k < (signed)init_zvalue.size(); k++) {
									writeoutputpr << (*bagtype)[k].name << "=\t" << monitor_value[k] << "\t+-\t" << monitor_stddv[k];
									writeoutputpr.flush();
									pgstr += (*bagtype)[k].name + "=\t" + std::to_string(monitor_value[k]) + "\t+-\t" + std::to_string(monitor_stddv[k]);
									if (k != (signed)init_zvalue.size() - 1) {
										writeoutputpr << "\t";
										writeoutputpr.flush();
										pgstr += "\t";
									}
									else {
										writeoutputpr << "\n";
										writeoutputpr.flush();
									}
								}

								int super_monitor = 0;
								MAT updMAT, startMAT;
								bool update = false;
								startMAT = objMAT;

								while (true) {

									for (int itest = 0; itest < 2; itest++) {
										
										int monitor = 0;

										while (true) {

											int target_t = -1;
											int target_tt = -1;
											long double min_upd_zdist = *VInfinite;
											bool success = false;
											MAT targetMAT = objMAT;

											for (int t = 0; t < (signed)objMAT.iatom.size(); t++) {

												int tt_size = 0;
												if (itest == 0) {
													tt_size = (signed)(*atomdat)[startMAT.iatom[t]].neighbor.size();
													if ((*atomdat)[startMAT.iatom[t]].dopeonly) {
														if (idope == 0) {
															tt_size = 0;
														}
														else {
															if (t % 2 == 0) {
																tt_size = 0;
															}
														}													}
												}
												else {
													tt_size = (signed)ModulateContentFactor->size();
												}

												for (int tt = 0; tt < tt_size; tt++) {
													
													if (itest == 0 || itest == 1 && (*ModulateContentFactor)[tt] != 1.0) {

														updMAT = objMAT;
														if (itest == 0) {
															updMAT.atomname[t] = (*atomdat)[startMAT.iatom[t]].neighbor[tt].atomname;
															updMAT.iatom[t] = (*atomdat)[startMAT.iatom[t]].neighbor[tt].iatom;
														}
														else {
															updMAT.content[t] = (*ModulateContentFactor)[tt] * startMAT.content[t];
														}

														designmate_iatom(&objMAT, atomdat);
														long double upd_CationSum = 0.0;
														long double upd_AnionSum = 0.0;
														long double* upd_scation_sum;
														upd_scation_sum = new long double[CationNames->size()];
														for (int upd_j = 0; upd_j < (signed)CationNames->size(); upd_j++) {
															upd_scation_sum[upd_j] = 0.0;
														}
														bool impossible = false;

														if (*AnionName != "-") {
															long double CationValenceTotal = 0.0;
															for (int ttt = 0; ttt < (signed)updMAT.iatom.size(); ttt++) {
																CationValenceTotal += updMAT.content[ttt] * (*atomdat)[updMAT.iatom[ttt]].v[CValence];
																upd_CationSum += updMAT.content[ttt];
															}
															updMAT.anion_content = -CationValenceTotal / *AnionValence;
															for (int ttt = 0; ttt < (signed)updMAT.atomname.size(); ttt++) {
																for (int upd_k = 0; upd_k < (signed)CationNames->size(); upd_k++) {
																	if (updMAT.atomname[ttt] == (*CationNames)[upd_k]) {
																		upd_scation_sum[upd_k] += updMAT.content[ttt];
																		break;
																	}
																}
															}
														}
														else if (Anions->size() > 0) {

															updMAT.OM = 0.0;
															long double CationValenceTotal = 0.0;
															long double AnionValenceTotal = 0.0;
															bool any_anion_exists = false;
															for (int ttt = 0; ttt < (signed)updMAT.iatom.size(); ttt++) {
																bool is_anion = false;
																long double anion_valence = 0.0;
																for (int tttt = 0; tttt < (signed)Anions->size(); tttt++) {
																	if ((*Anions)[tttt].anionname == updMAT.atomname[ttt]) {
																		is_anion = true;
																		any_anion_exists = true;
																		anion_valence = (*Anions)[tttt].anionvalence;
																		break;
																	}
																}
																if (!is_anion) {
																	if ((*atomdat)[updMAT.iatom[ttt]].v[CValence] != 0.0) {
																		CationValenceTotal += updMAT.content[ttt] * (*atomdat)[updMAT.iatom[ttt]].v[CValence];
																	}
																	else {
																		CationValenceTotal += updMAT.content[ttt] * (*atomdat)[updMAT.iatom[ttt]].v[AValence];
																	}
																	upd_CationSum += updMAT.content[ttt];
																}
																else {
																	AnionValenceTotal += updMAT.content[ttt] * anion_valence;
																	upd_AnionSum += updMAT.content[ttt];
																}
															}
															long double modulate_factor = -CationValenceTotal / AnionValenceTotal;
															if (any_anion_exists && AnionValenceTotal != 0.0 && modulate_factor > 0.0 && !(std::isnan(AnionValenceTotal) || std::isinf(AnionValenceTotal) || std::isinf(AnionValenceTotal)) && !(std::isnan(modulate_factor) || std::isinf(modulate_factor) || std::isinf(modulate_factor))) {
																for (int ttt = 0; ttt < (signed)updMAT.iatom.size(); ttt++) {
																	bool is_anion = false;
																	for (int tttt = 0; tttt < (signed)Anions->size(); tttt++) {
																		if ((*Anions)[tttt].anionname == updMAT.atomname[ttt]) {
																			is_anion = true;
																			break;
																		}
																	}
																	if (is_anion) {
																		long double candidate_content = updMAT.content[ttt] * modulate_factor;
																		if (candidate_content / startMAT.content[ttt] > 1.0 + *modulation_limit || candidate_content / startMAT.content[ttt] < 1.0 - *modulation_limit) {
																			impossible = true;
																		}
																		updMAT.content[ttt] = candidate_content;
																	}
																}
															}
															else {
																impossible = true;
															}
														}
														else {
															updMAT.OM = 0.0;
														}

														if (!impossible) {
															updMAT.composition = "";
															for (int ttt = 0; ttt < (signed)updMAT.atomname.size(); ttt++) {
																updMAT.composition += updMAT.atomname[ttt];
																if (updMAT.content[ttt] != 1.0) {
																	updMAT.composition += std::to_string(updMAT.content[ttt]).substr(0, 5);
																}
															}
															if (*AnionName != "-") {
																updMAT.composition += *AnionName;
																if (updMAT.anion_content != 1.0) {
																	updMAT.composition += std::to_string(updMAT.anion_content).substr(0, 5);
																}
																updMAT.OM = updMAT.anion_content / upd_CationSum;
																updMAT.CO.clear();
																for (int upd_k = 0; upd_k < (signed)CationNames->size(); upd_k++) {
																	updMAT.CO.push_back(upd_scation_sum[upd_k] / updMAT.anion_content);
																}
															}
															else {
																updMAT.OM = 0.0;
															}

															std::vector<long double> dv, mamssk;
															for (int j = 0; j < (signed)anadescriptor->size(); j++) {
																get_dv(&dv, &updMAT, atomdat, strdat, anadescriptor, &j, AveOnly);
																get_mamssk(&mamssk, &dv, VInfinite, &updMAT, DoNotAnalyze, AveOnly);
																for (int k = 0; k < 6; k++) {
																	xbox[j][k] = mamssk[k];
																}
																std::vector<long double>().swap(dv);
																std::vector<long double>().swap(mamssk);
															}
															xbox_OM = 0.0;
															if (*AnionName != "-") {
																xbox_OM = updMAT.OM;
																for (int k = 0; k < (signed)CationNames->size(); k++) {
																	xbox_CO[k] = updMAT.CO[k];
																}
															}
															std::vector<long double> upd_zvalue;
															long double upd_zdist = 0.0;
															bool stdev_test = true;
															for (int kk = 0; kk < (signed)bagtype->size(); kk++) {
																std::vector<int> test_anades;
																std::vector<int> test_mamssk;
																std::vector<long double> test_anades_degree;
																long double test_OM_degree = 0.0;
																long double test_CO_degree = 0.0;
																bool test_OM = false;
																int test_CO = -100;
																bool post_if_exists = true;
																std::vector<long double> rv, rvacc;
																long double v_stdev;
																long double post_v = 0.0;
																long double upd_value = calc_bag_chain(&rv, &rvacc, &v_stdev, regtype, &(*bagtype)[kk], &updMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v);
																std::vector<long double>().swap(rv);
																std::vector<long double>().swap(rvacc);
																std::vector<int>().swap(test_anades);
																std::vector<int>().swap(test_mamssk);
																std::vector<long double>().swap(test_anades_degree);
																monitor_value[kk] = upd_value + offset[kk];
																monitor_stddv[kk] = v_stdev;
																long double upz = (upd_value + offset[kk] - (*bagtype)[kk].v_ave) / (*bagtype)[kk].v_stddv;
																upd_zvalue.push_back(upz);
																upd_zdist += (*ZDistWeight)[kk] * powl(upz - (*zTargetVReg)[kk], 2.0);
																if (v_stdev > (*TargetStdDvReg)[kk]) {
																	stdev_test = false;
																}
															}

															std::vector<long double>().swap(upd_zvalue);
															upd_zdist = sqrtl(upd_zdist);
															if (upd_zdist < min_upd_zdist && upd_zdist < init_zdist && stdev_test) {
																min_upd_zdist = upd_zdist;
																target_t = t;
																target_tt = tt;
																targetMAT = updMAT;
																success = true;
																update = true;
																for (int kk = 0; kk < (signed)bagtype->size(); kk++) {
																	target_monitor_value[kk] = monitor_value[kk];
																	target_monitor_stddv[kk] = monitor_stddv[kk];
																}
																if (*verbose) {
																	write_zdist(&updMAT, &min_upd_zdist, AnionName);
																	writeoutputpr << "\t";
																	writeoutputpr.flush();
																	for (int k = 0; k < (signed)bagtype->size(); k++) {
																		writeoutputpr << (*bagtype)[k].name << "=\t" << monitor_value[k] << "\t+-\t" << monitor_stddv[k] << "\t";
																		writeoutputpr.flush();
																		if (k != (signed)bagtype->size() - 1) {
																			writeoutputpr << "\t";
																			writeoutputpr.flush();
																		}
																		else {
																			writeoutputpr << "\n";
																			writeoutputpr.flush();
																		}
																	}
																}

															}
														}
													
													}

												}
											}

											if (success) {
												objMAT = targetMAT;
												init_zdist = min_upd_zdist;
												if (min_upd_zdist < global_min_zdist) {
													global_min_zdist = min_upd_zdist;
													gMAT = targetMAT;
													for (int k = 0; k < (signed)bagtype->size(); k++) {
														gmonitor_value[k] = target_monitor_value[k];
														gmonitor_stddv[k] = target_monitor_stddv[k];
													}
													gstr = std::to_string(i + 1) + "\t...\t" + startMAT.name + "\t...\t" + startMAT.composition + "\t...\t" + startMAT.strname;
													gstr += "\n" + pgstr;
												}
												write_zdist(&objMAT, &min_upd_zdist, AnionName);
												writeoutputpr << "\t";
												writeoutputpr.flush();
												for (int k = 0; k < (signed)bagtype->size(); k++) {
													writeoutputpr << (*bagtype)[k].name << "=\t" << target_monitor_value[k] << "\t+-\t" << target_monitor_stddv[k] << "\t";
													writeoutputpr.flush();
													if (k != (signed)bagtype->size() - 1) {
														writeoutputpr << "\t";
														writeoutputpr.flush();
													}
													else {
														writeoutputpr << "\n";
														writeoutputpr.flush();
													}
												}
											}

											monitor++;
											if (monitor == *MaxIteration || !success) {
												break;
											}

										}

									}

									super_monitor++;
									if (super_monitor == *MaxIteration || !update) {
										break;
									}

								}

								if (update) {
									writeoutputpr << "@@@\t";
									writeoutputpr.flush();
									write_zdist(&objMAT, &init_zdist, AnionName);
									writeoutputpr << "\t";
									writeoutputpr.flush();
									for (int k = 0; k < (signed)bagtype->size(); k++) {
										writeoutputpr << (*bagtype)[k].name << "=\t" << target_monitor_value[k] << "\t+-\t" << target_monitor_stddv[k] << "\t";
										writeoutputpr.flush();
										if (k != (signed)bagtype->size() - 1) {
											writeoutputpr << "\t";
											writeoutputpr.flush();
										}
										else {
											writeoutputpr << "\n";
											writeoutputpr.flush();
										}
									}
								}
								else {
									if (init_zdist < global_min_zdist) {
										global_min_zdist = init_zdist;
										gMAT = startMAT;
										for (int k = 0; k < (signed)bagtype->size(); k++) {
											gmonitor_value[k] = init_monitor_value[k];
											gmonitor_stddv[k] = init_monitor_stddv[k];
										}
										gstr = std::to_string(i + 1) + "\t...\t" + startMAT.name + "\t...\t" + startMAT.composition + "\t...\t" + startMAT.strname;
										gstr += "\n" + pgstr;
									}
								}

								std::vector<long double>().swap(init_zvalue);
								std::vector<long double>().swap(monitor_value);
								std::vector<long double>().swap(monitor_stddv);
								std::vector<long double>().swap(target_monitor_value);
								std::vector<long double>().swap(target_monitor_stddv);
								std::vector<long double>().swap(init_monitor_value);
								std::vector<long double>().swap(init_monitor_stddv);
								std::vector<long double>().swap(offset);

							}
							
						}
					}
					if (*collect_only) {
						writeoutput << "=================================================================================\n";
						writeoutput.flush();
					}
					else {
						writeoutputpr << "=================================================================================\n";
						writeoutputpr.flush();
					}
				}
				else {
					if (pr_id == 0) {
						std::cout << "UseStrName was not identified ... Check StructureListFilename or UseStrName...\n";
						std::cout.flush();
					}
				}

				std::vector<std::string>().swap(PostTags);

			}

			if (i == Nmat - 1) {
				break;
			}

		}

		if (!*collect_only) {
			if (run_matdat) {
				writeoutputpr << gstr << "\n@@@\t";
				writeoutputpr.flush();
				write_zdist(&gMAT, &global_min_zdist, AnionName);
				writeoutputpr << "\t";
				writeoutputpr.flush();
				for (int k = 0; k < (signed)bagtype->size(); k++) {
					writeoutputpr << (*bagtype)[k].name << "=\t" << gmonitor_value[k] << "\t+-\t" << gmonitor_stddv[k] << "\t";
					writeoutputpr.flush();
					if (k != (signed)bagtype->size() - 1) {
						writeoutputpr << "\t";
						writeoutputpr.flush();
					}
					else {
						writeoutputpr << "\n";
						writeoutputpr.flush();
					}
				}
			}
		}
		else {
			
			for (int j = 0; j < (signed)anadescriptor->size(); j++) {
				(*anadescriptor)[j].ave /= (long double)(*anadescriptor)[j].count;
				(*anadescriptor)[j].ave_stdev /= (long double)(*anadescriptor)[j].count;
				(*anadescriptor)[j].ave_skew /= (long double)(*anadescriptor)[j].count;
				(*anadescriptor)[j].ave_kurto /= (long double)(*anadescriptor)[j].count;
				if ((*anadescriptor)[j].min_ave < 0.0 && (*anadescriptor)[j].max_ave < 0.0) {
					long double minvbuf = (*anadescriptor)[j].min_ave;
					(*anadescriptor)[j].min_ave = (*anadescriptor)[j].max_ave;
					(*anadescriptor)[j].max_ave = minvbuf;
				}
			}
			OM_ave /= (signed)OM_count;
			for (int k = 0; k < (signed)CationNames->size(); k++) {
				CO_ave[k] /= (signed)CO_count[k];
			}

		}

		readmatdat.close();

		if (!*collect_only && run_matdat) {
			writeoutputpr.close();
		}
	
	}

	if (*collect_only && bagtype->size() > 0 && regtype->size() > 0) {

		for (int k = 0; k < (signed)bagtype->size(); k++) {
			int Nsample = 0;
			if (!(*bagtype)[k].user_defined_v_ave) {
				(*bagtype)[k].v_ave = 0.0;
				int nancount = 0;
				for (int kk = 0; kk < (signed)(*bagtype)[k].v_list.size(); kk++) {
					if (!std::isnan((*bagtype)[k].v_list[kk])) {
						(*bagtype)[k].v_ave += (*bagtype)[k].v_list[kk];
					}
					else {
						nancount++;
					}
				}
				Nsample = (signed)(*bagtype)[k].v_list.size() - nancount;
				(*bagtype)[k].v_ave /= (long double)Nsample;
			}
			if (!(*bagtype)[k].user_defined_v_stddv) {
				(*bagtype)[k].v_stddv = 0.0;
				for (int kk = 0; kk < (signed)(*bagtype)[k].v_list.size(); kk++) {
					if (!std::isnan((*bagtype)[k].v_list[kk])) {
						(*bagtype)[k].v_stddv += powl((*bagtype)[k].v_list[kk] - (*bagtype)[k].v_ave, 2.0);
					}
				}
				(*bagtype)[k].v_stddv = sqrtl((*bagtype)[k].v_stddv / (long double)Nsample);
			}
			zTargetVReg->push_back(((*TargetVReg)[k] - (*bagtype)[k].v_ave) / (*bagtype)[k].v_stddv);
			if (pr_id == 0) {
				writeoutput << "___\t" << (*bagtype)[k].name << "\t:\tAve=\t" << (*bagtype)[k].v_ave << "\tStdDv=\t" << (*bagtype)[k].v_stddv << "\tzTarget=\t" << (*zTargetVReg)[k] << "\n";
				writeoutput.flush();
			}
		}
		if (pr_id == 0) {
			writeoutput << "=================================================================================\n";
			writeoutput.flush();

			writeanalysis << "=================================================================================\n";
			writeanalysis.flush();
			writeanalysis << "ModelName\tSStot\tSSres\tMSE\tRMSE\tMAE\txbar\tybar\tR2\tPearsonR\tPearsonR2\tinvalid_count\n";
			writeanalysis.flush();
		}
		for (int k = 0; k < (signed)regtype->size(); k++) {
			long double SStot, SSres, MSE, RMSE, MAE, xbar, ybar, R2, PearsonR;
			int invalid_count;
			bool use_limit = false;
			long double mVInfinite = -*VInfinite;
			int si = 0;
			int ji = 1;
			test_robustness(&SStot, &SSres, &MSE, &RMSE, &MAE, &xbar, &ybar, &R2, &PearsonR, &invalid_count, &(*regtype)[k].regressed_v_list, &(*regtype)[k].v_list, &use_limit, &mVInfinite, VInfinite, &si, &ji);
			if (pr_id == 0) {
				writeanalysis << "reg" << k << "_" << (*regtype)[k].name << "\t" << SStot << "\t" << SSres << "\t" << MSE << "\t" << RMSE << "\t" << MAE << "\t" << xbar << "\t" << ybar << "\t" << R2 << "\t" << PearsonR << "\t" << powl(PearsonR, 2.0) << "\t" << invalid_count << "\n";
				writeanalysis.flush();
			}
		}
		for (int k = 0; k < (signed)bagtype->size(); k++) {
			std::string secoeff_filename = (*bagtype)[k].name + ".se.txt";
			std::ifstream readsecoeff;
			readsecoeff.open(secoeff_filename);
			if (!readsecoeff.is_open()) {
				for (int kk = 0; kk < (signed)(*bagtype)[k].which_regtype.size(); kk++) {
					std::vector<DATASET> stackingemsemble_datset;
					stackingemsemble_datset.clear();
					for (int iele = 0; iele < (signed)(*regtype)[(*bagtype)[k].which_regtype[kk]].regressed_v_list.size(); iele++) {
						DATASET objDATASET;
						objDATASET.v.push_back((*bagtype)[k].v_list[iele]);
						bool appendable = true;
						if (std::isnan((*bagtype)[k].v_list[kk]) || std::isinf((*bagtype)[k].v_list[kk]) || std::isinf((*bagtype)[k].v_list[kk])) {
							appendable = false;
						}
						if (appendable) {
							for (int kkk = 0; kkk <= kk; kkk++) {
								objDATASET.v.push_back((*regtype)[(*bagtype)[k].which_regtype[kkk]].regressed_v_list[iele]);
								if (std::isnan((*regtype)[(*bagtype)[k].which_regtype[kkk]].regressed_v_list[iele]) || std::isinf((*regtype)[(*bagtype)[k].which_regtype[kkk]].regressed_v_list[iele]) || std::isinf((*regtype)[(*bagtype)[k].which_regtype[kkk]].regressed_v_list[iele])) {
									appendable = false;
									break;
								}
							}
						}
						if (appendable) {
							stackingemsemble_datset.push_back(objDATASET);
						}
					}
					REGELE objREGELE;
					int Noccup = kk + 1;
					long double CI_level = 0.95;
					long double critic_p = 0.05;
					bool for_test = false;
					do_regression(&objREGELE, &stackingemsemble_datset, &Noccup, &CI_level, &critic_p, &for_test);
					STACKINGENSEMBLE objSTACKINGENSEMBLE;
					objSTACKINGENSEMBLE.coeff_sum = 0.0;
					for (int icoeff = 0; icoeff < (signed)objREGELE.fit_coeff.size(); icoeff++) {
						objSTACKINGENSEMBLE.coeff.push_back(objREGELE.fit_coeff[icoeff]);
						if (icoeff > 0) {
							objSTACKINGENSEMBLE.coeff_sum += fabsl(objREGELE.fit_coeff[icoeff]);
						}
					}
					(*bagtype)[k].stackingensemble.push_back(objSTACKINGENSEMBLE);
					std::vector<DATASET>().swap(stackingemsemble_datset);
				}
			}
			else {
				(*bagtype)[k].stackingensemble.clear();
				if (pr_id == 0) {
					std::cout << "OVERWRITING Coefficients for stacking ensembles... " << (*bagtype)[k].name << " from " << secoeff_filename << "\n";
					std::cout.flush();
					writeanalysis << "OVERWRITING Coefficients for stacking ensembles... " << (*bagtype)[k].name << " from " << secoeff_filename << "\n";
					writeanalysis.flush();
				}
				std::getline(readsecoeff, sbuf);
				for (int kk = 0; kk < (signed)(*bagtype)[k].which_regtype.size(); kk++) {
					STACKINGENSEMBLE objSTACKINGENSEMBLE;
					objSTACKINGENSEMBLE.coeff_sum = 0.0;
					for (int kkk = 0; kkk < kk + 2; kkk++) {
						if (kkk == 0 || kkk == 1) {
							readsecoeff >> sbuf;
						}
						else {
							readsecoeff >> sbuf >> sbuf;
						}
						long double ldbuf;
						readsecoeff >> ldbuf;
						objSTACKINGENSEMBLE.coeff.push_back(ldbuf);
						if (objSTACKINGENSEMBLE.coeff.size() > 1) {
							objSTACKINGENSEMBLE.coeff_sum += fabsl(ldbuf);
						}
					}
					std::getline(readsecoeff, sbuf);
					(*bagtype)[k].stackingensemble.push_back(objSTACKINGENSEMBLE);
				}
			}
			readsecoeff.close();
		}
		for (int k = 0; k < (signed)bagtype->size(); k++) {
			for (int iele = 0; iele < (signed)(*bagtype)[k].regressed_v_list.size(); iele++) {
				long double ldbuf = (*bagtype)[k].regressed_v_list[iele];
				(*bagtype)[k].regressed_v_list[iele] = (*bagtype)[k].stackingensemble[(*bagtype)[k].stackingensemble.size() - 1].coeff[0];
				long double simple_ave = 0.0;
				for (int kk = 0; kk < (signed)(*bagtype)[k].which_regtype.size(); kk++) {
					(*bagtype)[k].regressed_v_list[iele] += (*bagtype)[k].stackingensemble[(*bagtype)[k].stackingensemble.size() - 1].coeff[kk + 1] * (*regtype)[(*bagtype)[k].which_regtype[kk]].regressed_v_list[iele];
					simple_ave += (*regtype)[(*bagtype)[k].which_regtype[kk]].regressed_v_list[iele];
				}
				simple_ave /= (long double)(signed)(*bagtype)[k].which_regtype.size();
				for (int kk = 0; kk < (signed)(*bagtype)[k].which_regtype.size(); kk++) {
					(*bagtype)[k].regressed_stddv_v_list[iele] += fabsl((*bagtype)[k].stackingensemble[(*bagtype)[k].stackingensemble.size() - 1].coeff[kk + 1]) * powl((*regtype)[(*bagtype)[k].which_regtype[kk]].regressed_v_list[iele] - simple_ave, 2.0) / (*bagtype)[k].stackingensemble[(*bagtype)[k].stackingensemble.size() - 1].coeff_sum;
				}
				(*bagtype)[k].regressed_stddv_v_list[iele] = sqrtl((*bagtype)[k].regressed_stddv_v_list[iele]);
				if (std::isnan((*bagtype)[k].regressed_stddv_v_list[iele]) || std::isinf((*bagtype)[k].regressed_stddv_v_list[iele]) || std::isinf(-(*bagtype)[k].regressed_stddv_v_list[iele])) {
					(*bagtype)[k].regressed_stddv_v_list[iele] = *VInfinite;
				}
				if (std::isnan(ldbuf)) {
					(*bagtype)[k].regressed_v_list[iele] = ldbuf;
				}
			}
		}
		for (int k = 0; k < (signed)bagtype->size(); k++) {
			long double SStot, SSres, MSE, RMSE, MAE, xbar, ybar, R2, PearsonR;
			int invalid_count;
			bool use_limit = false;
			long double mVInfinite = -*VInfinite;
			int si = 0;
			int ji = 1;
			test_robustness(&SStot, &SSres, &MSE, &RMSE, &MAE, &xbar, &ybar, &R2, &PearsonR, &invalid_count, &(*bagtype)[k].regressed_v_list, &(*bagtype)[k].v_list, &use_limit, &mVInfinite, VInfinite, &si, &ji);
			if (pr_id == 0) {
				writeanalysis << "bag" << k << "_" << (*bagtype)[k].name << "\t" << SStot << "\t" << SSres << "\t" << MSE << "\t" << RMSE << "\t" << MAE << "\t" << xbar << "\t" << ybar << "\t" << R2 << "\t" << PearsonR << "\t" << powl(PearsonR, 2.0) << "\t" << invalid_count << "\n";
				writeanalysis.flush();
			}
			for (int kk = 0; kk < (signed)(*bagtype)[k].which_regtype.size(); kk++) {
				std::vector<long double> acc_regressed_v_list;
				acc_regressed_v_list.clear();
				for (int iele = 0; iele < (signed)(*regtype)[(*bagtype)[k].which_regtype[kk]].regressed_v_list.size(); iele++) {
					long double ele_v = (*bagtype)[k].stackingensemble[kk].coeff[0];
					for (int kkk = 0; kkk <= kk; kkk++) {
						long double test_v = (*regtype)[(*bagtype)[k].which_regtype[kkk]].regressed_v_list[iele];
						ele_v += (*bagtype)[k].stackingensemble[kk].coeff[kkk + 1] * test_v;
					}
					if (!std::isnan((*bagtype)[k].regressed_v_list[iele])) {
						acc_regressed_v_list.push_back(ele_v);
					}
					else {
						acc_regressed_v_list.push_back(std::nanl(""));
					}
				}
				test_robustness(&SStot, &SSres, &MSE, &RMSE, &MAE, &xbar, &ybar, &R2, &PearsonR, &invalid_count, &acc_regressed_v_list, &(*bagtype)[k].v_list, &use_limit, &mVInfinite, VInfinite, &si, &ji);
				if (pr_id == 0) {
					writeanalysis << "...bag" << k << "_" << (*bagtype)[k].name << "_acclevel_" << kk << "\t" << SStot << "\t" << SSres << "\t" << MSE << "\t" << RMSE << "\t" << MAE << "\t" << xbar << "\t" << ybar << "\t" << R2 << "\t" << PearsonR << "\t" << powl(PearsonR, 2.0) << "\t" << invalid_count << "\n";
					writeanalysis.flush();
				}
				std::vector<long double>().swap(acc_regressed_v_list);
			}
		}
		if (pr_id == 0) {
			writeanalysis << "=================================================================================\n";
			writeanalysis.flush();
			writeanalysis << "StackingEnsemble\n";
			writeanalysis.flush();
			for (int k = 0; k < (signed)bagtype->size(); k++) {
				for (int kk = 0; kk < (signed)(*bagtype)[k].stackingensemble.size(); kk++) {
					writeanalysis << "...bag" << k << "_" << (*bagtype)[k].name << "_acclevel_" << kk << "=";
					writeanalysis.flush();
					for (int kkk = 0; kkk < (signed)(*bagtype)[k].stackingensemble[kk].coeff.size(); kkk++) {
						writeanalysis << "\t" << (*bagtype)[k].stackingensemble[kk].coeff[kkk];
						writeanalysis.flush();
						if (kkk > 0) {
							writeanalysis << "\t*" << "reg" << (*bagtype)[k].which_regtype[kkk - 1] << "_" << (*regtype)[(*bagtype)[k].which_regtype[kkk - 1]].name;
							writeanalysis.flush();
						}
						if (kkk != (signed)(*bagtype)[k].stackingensemble[kk].coeff.size() - 1) {
							writeanalysis << "\t+";
							writeanalysis.flush();
						}
					}
					writeanalysis << "\n";
					writeanalysis.flush();
				}
			}
			for (int iele = 0; iele < (signed)(*regtype)[(*bagtype)[0].which_regtype[0]].regressed_v_list.size(); iele++) {
				for (int k = 0; k < (signed)bagtype->size(); k++) {
					for (int kk = 0; kk < (signed)(*bagtype)[k].which_regtype.size(); kk++) {
						writeanalysis_second << (*regtype)[(*bagtype)[k].which_regtype[kk]].regressed_v_list[iele] << "\t";
						writeanalysis_second.flush();
					}
				}
				for (int k = 0; k < (signed)bagtype->size(); k++) {
					writeanalysis_second << (*bagtype)[k].regressed_v_list[iele] + (*bagtype)[k].regressed_post_v_list[iele] << "\t" << (*bagtype)[k].regressed_stddv_v_list[iele] << "\t";
					writeanalysis_second.flush();
				}
				for (int k = 0; k < (signed)bagtype->size(); k++) {
					writeanalysis_second << (*bagtype)[k].v_list[iele];
					writeanalysis_second.flush();
					if (k != (signed)bagtype->size() - 1) {
						writeanalysis_second << "\t";
						writeanalysis_second.flush();
					}
					else {
						writeanalysis_second << "\n";
						writeanalysis_second.flush();
					}
				}
			}
		}

		if (pr_id == 0) {

			for (int k = 0; k < (signed)bagtype->size(); k++) {

				std::string grand_meth = "";
				writeonevtest.open("bag" + std::to_string(k) + "_" + (*bagtype)[k].name + ".onevtest.txt", std::ios::binary);
				writeonevtest << std::fixed << std::setprecision(20);
				writeonevtest << "v\t";
				writeonevtest.flush();
				std::vector<std::string> feature_name;
				feature_name.clear();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].name;
					writeonevtest.flush();
					feature_name.push_back((*anadescriptor)[j].name);
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << *AnionName << "M";
					writeonevtest.flush();
					feature_name.push_back(*AnionName + "M");
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << (*CationNames)[j] << "O";
						writeonevtest.flush();
						feature_name.push_back((*CationNames)[j] + "O");
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "min\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].min;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << OM_min;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << CO_min[j];
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "min_ave\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].min_ave;
					writeonevtest.flush();
					if (!*AveOnly) {
						grand_meth += "Subscript[" + (*anadescriptor)[j].mathematica + ", min] = " + to_times10((*anadescriptor)[j].min_ave, 20) + ";\n";
					}
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << OM_ave;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << CO_ave[j];
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest << "ave\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].ave;
					writeonevtest.flush();
					if (!*AveOnly) {
						grand_meth += "Subscript[" + (*anadescriptor)[j].mathematica + ", ave] = " + to_times10((*anadescriptor)[j].ave, 20) + ";\n";
					}
					else {
						grand_meth += (*anadescriptor)[j].mathematica + " = " + to_times10((*anadescriptor)[j].ave, 20) + ";\n";
					}
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << OM_ave;
					writeonevtest.flush();
					if (!*AveOnly) {
						grand_meth += "Subscript[" + *AnionName + ", M] = " + to_times10(OM_ave, 20) + ";\n";
					}
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << CO_ave[j];
						writeonevtest.flush();
						if (!*AveOnly) {
							grand_meth += "Subscript[" + (*CationNames)[j] + ", " + *AnionName + "]= " + to_times10(CO_ave[j], 20) + ";\n";
						}
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "max_ave\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].max_ave;
					writeonevtest.flush();
					if (!*AveOnly) {
						grand_meth += "Subscript[" + (*anadescriptor)[j].mathematica + ", max] = " + to_times10((*anadescriptor)[j].max_ave, 20) + ";\n";
					}
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << OM_ave;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << CO_ave[j];
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest << "max\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].max;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << OM_max;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << CO_max[j];
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();

				long double zld = 0.0;

				writeonevtest << "min_stdev\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].min_stdev;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << zld;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << zld;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "ave_stdev\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].ave_stdev;
					writeonevtest.flush();
					if (!*AveOnly) {
						grand_meth += "Subscript[" + (*anadescriptor)[j].mathematica + ", stdev] = " + to_times10((*anadescriptor)[j].ave_stdev, 20) + ";\n";
					}
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << zld;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << zld;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "max_stdev\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].max_stdev;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << zld;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << zld;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "min_skew\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].min_skew;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << zld;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << zld;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "ave_skew\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].ave_skew;
					writeonevtest.flush();
					if (!*AveOnly) {
						grand_meth += "Subscript[" + (*anadescriptor)[j].mathematica + ", skew] = " + to_times10((*anadescriptor)[j].ave_skew, 20) + ";\n";
					}
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << zld;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << zld;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "max_skew\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].max_skew;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << zld;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << zld;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "min_kurto\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].min_kurto;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << *VInfinite;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << *VInfinite;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "ave_kurto\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].ave_kurto;
					writeonevtest.flush();
					if (!*AveOnly) {
						grand_meth += "Subscript[" + (*anadescriptor)[j].mathematica + ", kurto] = " + to_times10((*anadescriptor)[j].ave_kurto, 20) + ";\n";
					}
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << *VInfinite;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << *VInfinite;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "max_kurto\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					writeonevtest << (*anadescriptor)[j].max_kurto;
					writeonevtest.flush();
					if (j != (signed)anadescriptor->size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t" << *VInfinite;
					writeonevtest.flush();
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						writeonevtest << *VInfinite;
						writeonevtest.flush();
						if (j != (signed)CationNames->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				writeonevtest << "v\t";
				writeonevtest.flush();
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
						writeonevtest << (*anadescriptor)[j].name << "_" << jj << "\t";
						writeonevtest.flush();
					}
					for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
						writeonevtest << (*anadescriptor)[j].name << "_acc_" << jj;
						writeonevtest.flush();
						if (!(j == (signed)anadescriptor->size() - 1 && jj == (signed)(*bagtype)[k].which_regtype.size() - 1)) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				if (*AnionName != "-") {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
						writeonevtest << *AnionName << "M_" << jj << "\t";
						writeonevtest.flush();
					}
					for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
						writeonevtest << *AnionName << "M_acc_" << jj;
						writeonevtest.flush();
						if (jj != (signed)(*bagtype)[k].which_regtype.size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
				}
				if (CationNames->size() > 0) {
					writeonevtest << "\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)CationNames->size(); j++) {
						for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
							writeonevtest << (*CationNames)[j] << "O_" << jj << "\t";
							writeonevtest.flush();
						}
						for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
							writeonevtest << (*CationNames)[j] << "O_acc_" << jj;
							writeonevtest.flush();
							if (!(j == (signed)CationNames->size() - 1 && jj == (signed)(*bagtype)[k].which_regtype.size() - 1)) {
								writeonevtest << "\t";
								writeonevtest.flush();
							}
						}
					}
				}
				writeonevtest << "\t";
				writeonevtest.flush();
				for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
					writeonevtest << "v\t";
					writeonevtest.flush();
					for (int j = 0; j < (signed)anadescriptor->size(); j++) {
						writeonevtest << (*anadescriptor)[j].name << "_acc_" << jj;
						writeonevtest.flush();
						if (j != (signed)anadescriptor->size() - 1) {
							writeonevtest << "\t";
							writeonevtest.flush();
						}
					}
					if (*AnionName != "-") {
						writeonevtest << "\t";
						writeonevtest.flush();
						writeonevtest << *AnionName << "M_acc_" << jj;
						writeonevtest.flush();
					}
					if (CationNames->size() > 0) {
						writeonevtest << "\t";
						writeonevtest.flush();
						for (int j = 0; j < (signed)CationNames->size(); j++) {
							writeonevtest << (*CationNames)[j] << "O_acc_" << jj;
							writeonevtest.flush();
							if (j != (signed)CationNames->size() - 1) {
								writeonevtest << "\t";
								writeonevtest.flush();
							}
						}
					}
					if (jj != (signed)(*bagtype)[k].which_regtype.size() - 1) {
						writeonevtest << "\t";
						writeonevtest.flush();
					}
				}
				writeonevtest << "\n";
				writeonevtest.flush();
				for (int itest = 0; itest <= 100; itest++) {
					std::vector<long double> test_anades_degree;
					test_anades_degree.push_back((long double)itest / 100.0);
					writeonevtest << test_anades_degree[0] << "\t";
					writeonevtest.flush();
					std::vector<long double> v_series, v_series_acc, v_series_acc_grande;
					int count_feature = 0;
					for (int j = 0; j < (signed)anadescriptor->size(); j++) {
						count_feature++;
						v_series.clear();
						v_series_acc.clear();
						long double v_stdev;
						std::vector<int> test_anades;
						test_anades.push_back(j);
						std::vector<int> test_mamssk;
						test_mamssk.push_back(1);
						bool test_OM = false;
						int test_CO = -100;
						bool post_if_exists = true;
						long double test_OM_degree = 0.0;
						long double test_CO_degree = 0.0;
						MAT nullMAT;
						long double post_v = 0.0;
						calc_bag_chain(&v_series, &v_series_acc, &v_stdev, regtype, &(*bagtype)[k], &nullMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v);
						for (int jj = 0; jj < (signed)v_series.size(); jj++) {
							writeonevtest << v_series[jj] << "\t";
							writeonevtest.flush();
						}
						for (int jj = 0; jj < (signed)v_series_acc.size(); jj++) {
							writeonevtest << v_series_acc[jj];
							writeonevtest.flush();
							v_series_acc_grande.push_back(v_series_acc[jj]);
							(*regtype)[(*bagtype)[k].which_regtype[jj]].onev_acc_list.push_back(v_series_acc[jj]);
							if (!(j == (signed)anadescriptor->size() - 1 && jj == (signed)v_series_acc.size() - 1)) {
								writeonevtest << "\t";
								writeonevtest.flush();
							}
						}
						std::vector<long double>().swap(v_series);
						std::vector<long double>().swap(v_series_acc);
						std::vector<int>().swap(test_anades);
						std::vector<int>().swap(test_mamssk);
					}
					test_anades_degree.clear();
					if (*AnionName != "-") {
						writeonevtest << "\t";
						writeonevtest.flush();
						count_feature++;
						v_series.clear();
						v_series_acc.clear();
						long double v_stdev;
						std::vector<int> test_anades;
						std::vector<int> test_mamssk;
						bool test_OM = true;
						int test_CO = -100;
						long double test_OM_degree = (long double)itest / 100.0;
						long double test_CO_degree = 0.0;
						bool post_if_exists = true;
						MAT nullMAT;
						long double post_v = 0.0;
						calc_bag_chain(&v_series, &v_series_acc, &v_stdev, regtype, &(*bagtype)[k], &nullMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v);
						for (int jj = 0; jj < (signed)v_series.size(); jj++) {
							writeonevtest << v_series[jj] << "\t";
							writeonevtest.flush();
						}
						for (int jj = 0; jj < (signed)v_series_acc.size(); jj++) {
							writeonevtest << v_series_acc[jj];
							writeonevtest.flush();
							v_series_acc_grande.push_back(v_series_acc[jj]);
							(*regtype)[(*bagtype)[k].which_regtype[jj]].onev_acc_list.push_back(v_series_acc[jj]);
							if (jj != (signed)v_series_acc.size() - 1) {
								writeonevtest << "\t";
								writeonevtest.flush();
							}
						}
						std::vector<long double>().swap(v_series);
						std::vector<long double>().swap(v_series_acc);
						std::vector<int>().swap(test_anades);
						std::vector<int>().swap(test_mamssk);
					}
					if (CationNames->size() > 0) {
						writeonevtest << "\t";
						writeonevtest.flush();
						for (int j = 0; j < (signed)CationNames->size(); j++) {
							count_feature++;
							v_series.clear();
							v_series_acc.clear();
							long double v_stdev;
							std::vector<int> test_anades;
							std::vector<int> test_mamssk;
							bool test_OM = false;
							int test_CO = j;
							long double test_OM_degree = 0.0;
							long double test_CO_degree = (long double)itest / 100.0;
							bool post_if_exists = true;
							MAT nullMAT;
							long double post_v = 0.0;
							calc_bag_chain(&v_series, &v_series_acc, &v_stdev, regtype, &(*bagtype)[k], &nullMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v);
							for (int jj = 0; jj < (signed)v_series.size(); jj++) {
								writeonevtest << v_series[jj] << "\t";
								writeonevtest.flush();
							}
							for (int jj = 0; jj < (signed)v_series_acc.size(); jj++) {
								writeonevtest << v_series_acc[jj];
								writeonevtest.flush();
								v_series_acc_grande.push_back(v_series_acc[jj]);
								(*regtype)[(*bagtype)[k].which_regtype[jj]].onev_acc_list.push_back(v_series_acc[jj]);
								if (!(j == (signed)CationNames->size() - 1 && jj == (signed)v_series_acc.size() - 1)) {
									writeonevtest << "\t";
									writeonevtest.flush();
								}
							}
							std::vector<long double>().swap(v_series);
							std::vector<long double>().swap(v_series_acc);
							std::vector<int>().swap(test_anades);
							std::vector<int>().swap(test_mamssk);
						}
					}
					std::vector<long double>().swap(test_anades_degree);
					if (v_series_acc_grande.size() > 0) {
						writeonevtest << "\t";
						writeonevtest.flush();
						for (int jj = 0; jj < (signed)(*bagtype)[k].which_regtype.size(); jj++) {
							writeonevtest << (long double)itest / 100.0 << "\t";
							writeonevtest.flush();
							for (int jjj = 0; jjj < count_feature; jjj++) {
								writeonevtest << v_series_acc_grande[jj + jjj * (signed)(*bagtype)[k].which_regtype.size()];
								writeonevtest.flush();
								if (!(jj == (signed)(*bagtype)[k].which_regtype.size() - 1 && jjj == count_feature - 1)) {
									writeonevtest << "\t";
									writeonevtest.flush();
								}
							}
						}
					}
					std::vector<long double>().swap(v_series);
					std::vector<long double>().swap(v_series_acc);
					std::vector<long double>().swap(v_series_acc_grande);
					writeonevtest << "\n";
					writeonevtest.flush();
				}

				writeonevtest << "=================================================================================\n";
				writeonevtest.flush();
				writeonevtest << "ModelName\tSStot\tSSres\tMSE\tRMSE\tMAE\txbar\tybar\tR2\tPearsonR\tPearsonR2\tinvalid_count\n";
				writeonevtest.flush();
				if ((*bagtype)[k].which_regtype.size() >= 2) {
					for (int ireg = 0; ireg < (signed)(*bagtype)[k].which_regtype.size() - 1; ireg++) {
						long double SStot, SSres, MSE, RMSE, MAE, xbar, ybar, R2, PearsonR;
						int invalid_count;
						bool use_limit = false;
						long double mVInfinite = -*VInfinite;
						int si = 0;
						int ji = 1;
						test_robustness(&SStot, &SSres, &MSE, &RMSE, &MAE, &xbar, &ybar, &R2, &PearsonR, &invalid_count, &(*regtype)[(*bagtype)[k].which_regtype[ireg]].onev_acc_list, &(*regtype)[(*bagtype)[k].which_regtype[ireg + 1]].onev_acc_list, &use_limit, &mVInfinite, VInfinite, &si, &ji);
						writeonevtest << "acc" << ireg << "-" << "acc" << ireg + 1 << "\t" << SStot << "\t" << SSres << "\t" << MSE << "\t" << RMSE << "\t" << MAE << "\t" << xbar << "\t" << ybar << "\t" << R2 << "\t" << PearsonR << "\t" << powl(PearsonR, 2.0) << "\t" << invalid_count << "\n";
						writeonevtest.flush();
					}
				}

				writeonevtest << "=================================================================================\n";
				writeonevtest.flush();
				writeonevtest << "y:[\t" << (*bagtype)[k].obs_v_min << "\t,\t" << (*bagtype)[k].obs_v_max << "\t]\n";
				writeonevtest.flush();
				writeonevtest << "ModelName\tRestrained_SStot\tRestrained_SSres\tRestrained_MSE\tRestrained_RMSE\tRestrained_MAE\tRestrained_xbar\tRestrained_ybar\tRestrained_R2\tRestrained_PearsonR\tRestrained_PearsonR2\tRestrained_invalid_count\n";
				writeonevtest.flush();
				if ((*bagtype)[k].which_regtype.size() >= 2) {
					for (int ireg = 0; ireg < (signed)(*bagtype)[k].which_regtype.size() - 1; ireg++) {
						long double SStot, SSres, MSE, RMSE, MAE, xbar, ybar, R2, PearsonR;
						int invalid_count;
						bool use_limit = true;
						int si = 0;
						int ji = 1;
						test_robustness(&SStot, &SSres, &MSE, &RMSE, &MAE, &xbar, &ybar, &R2, &PearsonR, &invalid_count, &(*regtype)[(*bagtype)[k].which_regtype[ireg]].onev_acc_list, &(*regtype)[(*bagtype)[k].which_regtype[ireg + 1]].onev_acc_list, &use_limit, &(*bagtype)[k].obs_v_min, &(*bagtype)[k].obs_v_max, &si, &ji);
						writeonevtest << "acc" << ireg << "-" << "acc" << ireg + 1 << "\t" << SStot << "\t" << SSres << "\t" << MSE << "\t" << RMSE << "\t" << MAE << "\t" << xbar << "\t" << ybar << "\t" << R2 << "\t" << PearsonR << "\t" << powl(PearsonR, 2.0) << "\t" << invalid_count << "\n";
						writeonevtest.flush();
					}
				}

				for (int ifeat = 0; ifeat < (signed)feature_name.size(); ifeat++) {
					writeonevtest << "=================================================================================\n";
					writeonevtest.flush();
					writeonevtest << "TARGET:\t" << feature_name[ifeat] << "\n";
					writeonevtest.flush();
					writeonevtest << "y:[\t" << (*bagtype)[k].obs_v_min << "\t,\t" << (*bagtype)[k].obs_v_max << "\t]\n";
					writeonevtest.flush();
					writeonevtest << "ModelName\tRestrained_SStot\tRestrained_SSres\tRestrained_MSE\tRestrained_RMSE\tRestrained_MAE\tRestrained_xbar\tRestrained_ybar\tRestrained_R2\tRestrained_PearsonR\tRestrained_PearsonR2\tRestrained_invalid_count\n";
					writeonevtest.flush();
					if ((*bagtype)[k].which_regtype.size() >= 2) {
						for (int ireg = 0; ireg < (signed)(*bagtype)[k].which_regtype.size() - 1; ireg++) {
							long double SStot, SSres, MSE, RMSE, MAE, xbar, ybar, R2, PearsonR;
							int invalid_count;
							bool use_limit = true;
							int si = ifeat;
							int ji = (signed)feature_name.size();
							test_robustness(&SStot, &SSres, &MSE, &RMSE, &MAE, &xbar, &ybar, &R2, &PearsonR, &invalid_count, &(*regtype)[(*bagtype)[k].which_regtype[ireg]].onev_acc_list, &(*regtype)[(*bagtype)[k].which_regtype[ireg + 1]].onev_acc_list, &use_limit, &(*bagtype)[k].obs_v_min, &(*bagtype)[k].obs_v_max, &si, &ji);
							writeonevtest << "acc" << ireg << "-" << "acc" << ireg + 1 << "\t" << SStot << "\t" << SSres << "\t" << MSE << "\t" << RMSE << "\t" << MAE << "\t" << xbar << "\t" << ybar << "\t" << R2 << "\t" << PearsonR << "\t" << powl(PearsonR, 2.0) << "\t" << invalid_count << "\n";
							writeonevtest.flush();
						}
					}
				}

				writeonevtest << "=================================================================================\n";
				writeonevtest << grand_meth;

				writeonevtest.close();

				std::vector<std::string>().swap(feature_name);

			}
		
		}
		
	}

}

void process_mapout(std::vector<long double>* zTargetVReg, std::vector<long double>* ZDistWeight, std::vector<long double>* TargetVReg, std::string* MapoutFilename, std::vector<REGTYPE>* regtype, std::vector<BAGTYPE>* bagtype, std::vector<ATOMDAT>* atomdat, std::vector<STRDAT>* strdat, std::vector<ANADESCRIPTOR>* anadescriptor, long double* VInfinite, std::vector<std::string>* DoNotAnalyze, std::vector<std::string>* CationNames, bool* verbose, bool* AveOnly) {

	readmapout.open(*MapoutFilename);

	if (readmapout.is_open()) {

		readmapout >> sbuf >> mapout_preNRandom >> mapout_NRandom;
		readmapout >> sbuf >> mapout_Verbose;
		std::getline(readmapout, sbuf);
		std::getline(readmapout, sbuf);
		int system_count = 0;
		long double global_min;
		global_min = (*VInfinite);
		std::string global_min_str;
		SMATLIST objSMATLIST;
		std::vector<long double> ave, stdev;
		for (int ifunc = 0; ifunc < (signed)bagtype->size(); ifunc++) {
			ave.push_back(0.0);
			stdev.push_back(0.0);
		}
		int max_comp = -1;

		int sample_count = 0;
		if (pr_id == 0) {
			writeoutput << "Gathering ave vector...\n";
			writeoutput.flush();
		}
		while (true) {

			std::getline(readmapout, sbuf);
			if (sbuf == "STOP" || readmapout.eof()) {
				break;
			}

			MAT objMAT;
			std::stringstream lcstr;
			lcstr << sbuf;
			if (sbuf.find("*") == std::string::npos) {
				while (true) {
					lcstr >> sbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(0.0);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			else {
				lcstr >> sbuf;
				while (true) {
					long double ldbuf;
					lcstr >> sbuf >> ldbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(ldbuf);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			
			system_count++;
			designmate_iatom(&objMAT, atomdat);
			if ((signed)objMAT.content.size() > max_comp) {
				max_comp = (signed)objMAT.content.size();
			}

			for (int i = 0; i < (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom + mapout_preNRandom; i++) {

				if (i < objMAT.content.size()) {
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii == i) {
							objMAT.content[ii] = 1.0;
						}
						else {
							objMAT.content[ii] = 0.0;
						}
					}
				}
				else if (objMAT.content.size() > 2 && ((signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom)) {
					int isetzero = (i - (signed)objMAT.content.size()) / mapout_preNRandom;
					if (isetzero < 0) {
						isetzero = 0;
					}
					else if (isetzero >= (signed)objMAT.content.size()) {
						isetzero = (signed)objMAT.content.size() - 1;
					}
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii != isetzero) {
							objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						}
						else {
							objMAT.content[ii] = 0.0;
						}
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}
				else {
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}

				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					std::vector<long double> dv, mamssk;
					get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j, AveOnly);
					get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze, AveOnly);
					for (int k = 0; k < 6; k++) {
						xbox[j][k] = mamssk[k];
					}
					xbox_OM = 0.0;
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}

				sample_count++;
				for (int k = 0; k < (signed)bagtype->size(); k++) {
					std::vector<int> test_anades;
					std::vector<int> test_mamssk;
					std::vector<long double> test_anades_degree;
					bool test_OM = false;
					int test_CO = -100;
					long double test_OM_degree = 0.0;
					long double test_CO_degree = 0.0;
					bool post_if_exists = true;
					std::vector<long double> v, va;
					long double v_stdev;
					long double post_v = 0.0;
					ave[k] += calc_bag_chain(&v, &va, &v_stdev, regtype, &(*bagtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v);
					std::vector<long double>().swap(v);
					std::vector<long double>().swap(va);
					std::vector<int>().swap(test_anades);
					std::vector<int>().swap(test_mamssk);
					std::vector<long double>().swap(test_anades_degree);
				}

			}

		}

		for (int ifunc = 0; ifunc < (signed)bagtype->size(); ifunc++) {
			ave[ifunc] /= (long double)sample_count;
			if (pr_id == 0) {
				writeoutput << (*bagtype)[ifunc].name << "\t[AVE]\t" << ave[ifunc] << "\n";
				writeoutput.flush();
			}
		}

		readmapout.close();
		readmapout.open(*MapoutFilename);
		for (int iline = 0; iline < 3; iline++) {
			std::getline(readmapout, sbuf);
		}

		if (pr_id == 0) {
			writeoutput << "Gathering stdev vector...\n";
			writeoutput.flush();
		}
		for (int isystem = 0; isystem < system_count; isystem++) {

			std::getline(readmapout, sbuf);
			MAT objMAT;
			std::stringstream lcstr;
			lcstr << sbuf;
			if (sbuf.find("*") == std::string::npos) {
				while (true) {
					lcstr >> sbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(0.0);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			else {
				lcstr >> sbuf;
				while (true) {
					long double ldbuf;
					lcstr >> sbuf >> ldbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(ldbuf);
					if (lcstr.eof()) {
						break;
					}
				}
			}

			designmate_iatom(&objMAT, atomdat);

			for (int i = 0; i < (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom + mapout_preNRandom; i++) {

				if (i < objMAT.content.size()) {
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii == i) {
							objMAT.content[ii] = 1.0;
						}
						else {
							objMAT.content[ii] = 0.0;
						}
					}
				}
				else if (objMAT.content.size() > 2 && ((signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_preNRandom)) {
					int isetzero = (i - (signed)objMAT.content.size()) / mapout_preNRandom;
					if (isetzero < 0) {
						isetzero = 0;
					}
					else if (isetzero >= (signed)objMAT.content.size()) {
						isetzero = (signed)objMAT.content.size() - 1;
					}
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						if (ii != isetzero) {
							objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						}
						else {
							objMAT.content[ii] = 0.0;
						}
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}
				else {
					long double norm = 0.0;
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
						norm += objMAT.content[ii];
					}
					for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
						objMAT.content[ii] /= norm;
					}
				}

				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					std::vector<long double> dv, mamssk;
					get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j, AveOnly);
					get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze, AveOnly);
					for (int k = 0; k < 6; k++) {
						xbox[j][k] = mamssk[k];
					}
					xbox_OM = 0.0;
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}

				for (int k = 0; k < (signed)bagtype->size(); k++) {
					std::vector<int> test_anades;
					std::vector<int> test_mamssk;
					std::vector<long double> test_anades_degree;
					bool test_OM = false;
					int test_CO = -100;
					long double test_OM_degree = 0.0;
					long double test_CO_degree = 0.0;
					bool post_if_exists = true;
					std::vector<long double> v, va;
					long double v_stdev;
					long double post_v = 0.0;
					stdev[k] += powl(calc_bag_chain(&v, &va, &v_stdev, regtype, &(*bagtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v) - ave[k], 2.0);
					std::vector<long double>().swap(v);
					std::vector<long double>().swap(va);
					std::vector<int>().swap(test_anades);
					std::vector<int>().swap(test_mamssk);
					std::vector<long double>().swap(test_anades_degree);
				}

			}

		}

		for (int ifunc = 0; ifunc < (signed)bagtype->size(); ifunc++) {
			stdev[ifunc] /= (long double)sample_count;
			stdev[ifunc] = sqrtl(stdev[ifunc]);
			zTargetVReg->push_back(((*TargetVReg)[ifunc] - ave[ifunc])/stdev[ifunc]);
			if (pr_id == 0) {
				writeoutput << (*bagtype)[ifunc].name << "\t[STDEV]\t" << stdev[ifunc] << "\t[TargetVReg]\t" << (*TargetVReg)[ifunc] << "\t[zTargetVReg]\t" << (*zTargetVReg)[ifunc] << "\n";
				writeoutput.flush();
			}
		}

		readmapout.close();
		readmapout.open(*MapoutFilename);
		for (int iline = 0; iline < 3; iline++) {
			std::getline(readmapout, sbuf);
		}

		for (int isystem = 0; isystem < system_count; isystem++) {

			bool specific_content = false;
			std::getline(readmapout, sbuf);
			MAT objMAT;
			std::stringstream lcstr;
			lcstr << sbuf;
			if (sbuf.find("*") == std::string::npos) {
				while (true) {
					lcstr >> sbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(0.0);
					if (lcstr.eof()) {
						break;
					}
				}
			}
			else {
				specific_content = true;
				lcstr >> sbuf;
				while (true) {
					long double ldbuf;
					lcstr >> sbuf >> ldbuf;
					objMAT.atomname.push_back(sbuf);
					objMAT.content.push_back(ldbuf);
					if (lcstr.eof()) {
						break;
					}
				}
			}

			designmate_iatom(&objMAT, atomdat);

			long double system_min = (*VInfinite);
			long double* v_at_system_min;
			long double* z_at_system_min;
			v_at_system_min = new long double[bagtype->size()];
			z_at_system_min = new long double[bagtype->size()];
			for (int j = 0; j < (signed)bagtype->size(); j++) {
				v_at_system_min[j] = 0.0;
				z_at_system_min[j] = 0.0;
			}
			long double* min_content;
			min_content = new long double [(signed)objMAT.content.size()];
			for (int j = 0; j < (signed)objMAT.content.size(); j++) {
				min_content[j] = 0.0;
			}

			int upto = 1;
			if (!specific_content) {
				upto = (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_NRandom + mapout_NRandom;
			}

			for (int i = 0; i < upto; i++) {

				if (!specific_content) {
					if (i < objMAT.content.size()) {
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							if (ii == i) {
								objMAT.content[ii] = 1.0;
							}
							else {
								objMAT.content[ii] = 0.0;
							}
						}
					}
					else if (objMAT.content.size() > 2 && (i < (signed)objMAT.content.size() + (signed)objMAT.content.size() * mapout_NRandom)) {
						int isetzero = (i - (signed)objMAT.content.size()) / mapout_NRandom;
						if (isetzero < 0) {
							isetzero = 0;
						}
						else if (isetzero >= (signed)objMAT.content.size()) {
							isetzero = (signed)objMAT.content.size() - 1;
						}
						long double norm = 0.0;
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							if (ii != isetzero) {
								objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
							}
							else {
								objMAT.content[ii] = 0.0;
							}
							norm += objMAT.content[ii];
						}
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							objMAT.content[ii] /= norm;
						}
					}
					else {
						long double norm = 0.0;
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							objMAT.content[ii] = folder((long double)rand() / RAND_MAX);
							norm += objMAT.content[ii];
						}
						for (int ii = 0; ii < (signed)objMAT.content.size(); ii++) {
							objMAT.content[ii] /= norm;
						}
					}
				}
				
				for (int j = 0; j < (signed)anadescriptor->size(); j++) {
					std::vector<long double> dv, mamssk;
					get_dv(&dv, &objMAT, atomdat, strdat, anadescriptor, &j, AveOnly);
					get_mamssk(&mamssk, &dv, VInfinite, &objMAT, DoNotAnalyze, AveOnly);
					for (int k = 0; k < 6; k++) {
						xbox[j][k] = mamssk[k];
					}
					xbox_OM = 0.0;
					std::vector<long double>().swap(dv);
					std::vector<long double>().swap(mamssk);
				}

				std::vector<long double> v, z;
				long double z_dist = 0.0;
				for (int k = 0; k < (signed)bagtype->size(); k++) {
					std::vector<int> test_anades;
					std::vector<int> test_mamssk;
					std::vector<long double> test_anades_degree;
					bool test_OM = false;
					int test_CO = -100;
					long double test_OM_degree = 0.0;
					long double test_CO_degree = 0.0;
					bool post_if_exists = true;
					std::vector<long double> v, va;
					long double v_stdev;
					long double post_v = 0.0;
					long double v_loc = calc_bag_chain(&v, &va, &v_stdev, regtype, &(*bagtype)[k], &objMAT, atomdat, strdat, anadescriptor, VInfinite, DoNotAnalyze, CationNames, &test_anades, &test_mamssk, &test_OM, &test_CO, &test_anades_degree, &test_OM_degree, &test_CO_degree, verbose, &post_if_exists, AveOnly, &post_v);
					long double z_loc = (v_loc - ave[k]) / stdev[k];
					v.push_back(v_loc);
					z.push_back(z_loc);
					z_dist += (*ZDistWeight)[k] * powl(z_loc - (*zTargetVReg)[k], 2.0);
					std::vector<long double>().swap(v);
					std::vector<long double>().swap(va);
					std::vector<int>().swap(test_anades);
					std::vector<int>().swap(test_mamssk);
					std::vector<long double>().swap(test_anades_degree);
				}

				if (mapout_Verbose == 1) {
					if (pr_id == 0) {
						writeoutput << isystem + 1 << "\t" << z_dist;
						writeoutput.flush();
					}
					int finite_comp = 0;
					for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
						if (objMAT.content[iatom] != 0.0) {
							finite_comp++;
							if (pr_id == 0) {
								writeoutput << "\t" << objMAT.atomname[iatom] << "\t" << objMAT.content[iatom];
								writeoutput.flush();
							}
						}
					}
					if (pr_id == 0) {
						for (int vatom = 0; vatom < max_comp - finite_comp; vatom++) {
							writeoutput << "\t\t";
							writeoutput.flush();
						}
						writeoutput << "\t...\t";
						writeoutput.flush();
						for (int k = 0; k < (signed)bagtype->size(); k++) {
							writeoutput << (*bagtype)[k].name << "\t[v]\t" << v[k] << "\t[z]\t" << z[k];
							writeoutput.flush();
							if (k != (signed)bagtype->size() - 1) {
								writeoutput << "\t";
								writeoutput.flush();
							}
							else {
								writeoutput << "\n";
								writeoutput.flush();
							}
						}
					}
				}

				if (z_dist < system_min) {
					system_min = z_dist;
					for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
						min_content[iatom] = objMAT.content[iatom];
					}
					for (int k = 0; k < (signed)bagtype->size(); k++) {
						v_at_system_min[k] = v[k];
						z_at_system_min[k] = z[k];
					}
				}

				std::vector<long double>().swap(v);
				std::vector<long double>().swap(z);

			}

			
			SMAT objSMATmin, objSMATmax;
			for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
				if (min_content[iatom] != 0.0) {
					objSMATmin.atomname.push_back(objMAT.atomname[iatom]);
					objSMATmin.content.push_back(min_content[iatom]);
				}
			}
			bool overlapmin = false;
			for (int ilist = (signed)objSMATLIST.smat_min.size() - 1; ilist >= 0; ilist--) {
				bool is_same = true;
				if (objSMATLIST.smat_min[ilist].atomname.size() == objSMATmin.atomname.size()) {
					for (int iele = 0; iele < (signed)objSMATmin.atomname.size(); iele++) {
						if (objSMATLIST.smat_min[ilist].atomname[iele] != objSMATmin.atomname[iele]) {
							is_same = false;
							break;
						}
					}
				}
				else {
					is_same = false;
				}
				if (is_same) {
					overlapmin = true;
					break;
				}
			}
			if (!overlapmin) {
				objSMATLIST.smat_min.push_back(objSMATmin);
				if (pr_id == 0) {
					writeoutput << isystem + 1 << "\t" << system_min;
					writeoutput.flush();
				}
				int finite_comp = 0;
				for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
					if (min_content[iatom] != 0.0) {
						finite_comp++;
						if (pr_id == 0) {
							writeoutput << "\t" << objMAT.atomname[iatom] << "\t" << min_content[iatom];
							writeoutput.flush();
						}
					}
				}
				if (pr_id == 0) {
					for (int vatom = 0; vatom < max_comp - finite_comp; vatom++) {
						writeoutput << "\t\t";
						writeoutput.flush();
					}
					writeoutput << "\t...\t";
					writeoutput.flush();
					for (int k = 0; k < (signed)bagtype->size(); k++) {
						writeoutput << (*bagtype)[k].name << "\t[v]\t" << v_at_system_min[k] << "\t[z]\t" << z_at_system_min[k];
						writeoutput.flush();
						if (k != (signed)bagtype->size() - 1) {
							writeoutput << "\t";
							writeoutput.flush();
						}
						else {
							writeoutput << "\n";
							writeoutput.flush();
						}
					}
				}
			}
			
			if (system_min < global_min) {
				global_min = system_min;
				std::stringstream globalstr;
				globalstr << "GlobalMIN:\t" << system_min << "\n(GlobalMIN_at):\t";
				for (int iatom = 0; iatom < (signed)objMAT.content.size(); iatom++) {
					if (min_content[iatom] != 0.0) {
						globalstr << "\t" << objMAT.atomname[iatom] << "\t" << min_content[iatom] << "\t";
					}
				}
				globalstr << "...\t";
				for (int k = 0; k < (signed)bagtype->size(); k++) {
					globalstr << (*bagtype)[k].name << "\t[v]\t" << v_at_system_min[k] << "\t[z]\t" << z_at_system_min[k];
					if (k != (signed)bagtype->size() - 1) {
						globalstr << "\t";
					}
				}
				global_min_str = globalstr.str();

			}
			
			std::vector<std::string>().swap(objMAT.atomname);
			std::vector<long double>().swap(objMAT.content);

			delete[] v_at_system_min;
			delete[] z_at_system_min;
			delete[] min_content;
			v_at_system_min = nullptr;
			z_at_system_min = nullptr;
			min_content = nullptr;

		}

		if (pr_id == 0) {
			writeoutput << "****************************************************************\n";
			writeoutput.flush();
			writeoutput << global_min_str << "\n";
			writeoutput.flush();
			writeoutput << "****************************************************************\n";
			writeoutput.flush();
		}

		std::vector<long double>().swap(ave);
		std::vector<long double>().swap(stdev);

	}

	readmapout.close();

}

int main(int argc, char* argv[]) {

	ierr = MPI_Init(&argc, &argv);
	if (ierr != 0)
	{
		std::cout << "  MPI_Init returned nonzero ierr.\n";
		exit(1);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &pr_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &pr_id);

	read_para();
	if (FaceFilename != "-") {
		read_face(&face, &FaceFeatureName, &FaceFilename);
	}
	readfold.open("fold.txt");
	FoldMode = false;
	bool TestModelStDevLimitFound = false;
	if (readfold.is_open()) {
		FoldMode = true;
		std::getline(readfold, sbuf);
		readfold >> NFolds >> TakeFold >> Train0Test1 >> Straight0Rotate1;
		if (!readfold.eof()) {
			readfold >> sbuf;
			if (sbuf.find("TestModelStDevLimit") != std::string::npos) {
				readfold >> TestModelStDevLimit;
				TestModelStDevLimitFound = true;
			}
		}
		if (pr_id == 0) {
			std::cout << "FoldMode = true\n";
			std::cout.flush();
			std::cout << "NFolds= " << NFolds << " , TakeFold= " << TakeFold << " , Train0Test1= " << Train0Test1 << " , Straight0Rotate1= " << Straight0Rotate1 << "\n";
			std::cout.flush();
			if (TestModelStDevLimitFound) {
				std::cout << "TestModelStDevLimit= " << TestModelStDevLimit << "\n";
				std::cout.flush();
			}
		}
	}
	readfold.close();

	if (pr_id == 0) {
		writeoutput.open(OutputFilename, std::ios::binary);
		writeoutput << std::fixed << std::setprecision(20);
	}
	if (pr_id == 0) {
		writeanalysis.open(AnalysisOutputFilename, std::ios::binary);
		writeanalysis << std::fixed << std::setprecision(20);
		writeanalysis_second.open("__R__." + AnalysisOutputFilename, std::ios::binary);
		writeanalysis_second << std::fixed << std::setprecision(20);
		writeanalysis << "No.\tName\tCompositionName\tStrcture\t";
		writeanalysis.flush();
		if (AnionName != "-") {
			writeanalysis << AnionName << "M\t";
			writeanalysis.flush();
			for (int k = 0; k < (signed)CationNames.size(); k++) {
				writeanalysis << CationNames[k] << AnionName << "\t";
				writeanalysis.flush();
			}
		}
		if (FaceFilename != "-") {
			for (int j = 0; j < (signed)FaceFeatureName.size(); j++) {
				writeanalysis << FaceFeatureName[j] << "\t";
				writeanalysis.flush();
			}
		}
		for (int j = 0; j < (signed)anadescriptor.size(); j++) {
			if (anadescriptor[j].ASB != "U") {
				writeanalysis << "Min" << anadescriptor[j].name << "\t";
				writeanalysis.flush();
				writeanalysis << "Ave" << anadescriptor[j].name << "\t";
				writeanalysis.flush();
				writeanalysis << "Max" << anadescriptor[j].name << "\t";
				writeanalysis.flush();
				writeanalysis << "StdDv" << anadescriptor[j].name << "\t";
				writeanalysis.flush();
				writeanalysis << "Skew" << anadescriptor[j].name << "\t";
				writeanalysis.flush();
				writeanalysis << "Kurto" << anadescriptor[j].name;
				writeanalysis.flush();
			}
			else {
				writeanalysis << anadescriptor[j].name;
				writeanalysis.flush();
			}
			if (j != (signed)anadescriptor.size() - 1) {
				writeanalysis << "\t";
				writeanalysis.flush();
			}
		}
	}

	read_atomdat(&atomdat, &descriptor, &anadescriptor, &AtomDatListFilename, &DoNotAnalyze);
	pick_neighbor(&atomdat, &NeighborSelectSpace, &NeighborZCrit, &NeighborCrit, &ExcludeNeighbors, &DoNotFindNeighbors, &DopeOnly, &descriptor);
	bool usestrname = true;
	if (UseStrName == 0) {
		usestrname = false;
	}
	read_strdat(&strdat, &anadescriptor, &sdescriptorname, &sitename, &StructureListFilename, &descriptor, &atomdat, &usestrname);
	bool AveOnly = false;
	if (SkipMinMaxKurto == 2) {
		AveOnly = true;
	}
	for (int i = 0; i < (signed)RegModelFilename.size(); i++) {
		regtype.push_back(read_reg_chain(&RegModelFilename[i], &anadescriptor, &AnionName, &CationNames, &verbose, &AveOnly));
		if (pr_id == 0) {
			writelatex.open(regtype[regtype.size() - 1].name + "." + std::to_string(i) + ".latex.txt");
			writelatex << std::fixed << std::setprecision(20);
			for (int j = 0; j < regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func.size(); j++) {
				writelatex << j + 1 << "\t" << to_times_e(regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func[j].coeff, 20) << "\t" << regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func[j].zcoeff << "\t" << regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func[j].latex << "\n";
				writelatex.flush();
			}
			writelatex << "INTERCEPT\t" << to_times_e(regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].intercept, 20);
			writelatex.flush();
			writelatex.close();
			writemathematica.open(regtype[regtype.size() - 1].name + "." + std::to_string(i) + ".mathematica.txt");
			writemathematica << std::fixed << std::setprecision(20);
		}
		std::string modname = "";
		for (char ch : regtype[regtype.size() - 1].name) {
			if (ch != '_') 
				modname += ch;
		} 
		std::string mathe = modname + std::to_string(i) + " = Simplify[ ";
		for (int j = 0; j < regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func.size(); j++) {
			std::string coeff_str = to_times10(regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func[j].coeff, 20);
			if (pr_id == 0) {
				writemathematica << j + 1 << "\t" << coeff_str << "\t" << regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func[j].zcoeff << "\t" << regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func[j].mathematica << "\n";
				writemathematica.flush();
			}
			if (j != 0) {
				mathe += " + ";
			}
			mathe += coeff_str + " " + regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].func[j].mathematica;
		}
		std::string intercept_str = to_times10(regtype[regtype.size() - 1].funcchain[regtype[regtype.size() - 1].funcchain.size() - 1].intercept, 20);
		if (pr_id == 0) {
			writemathematica << "INTERCEPT\t" << intercept_str;
			writemathematica.flush();
		}
		mathe += "+" + intercept_str + " ];";
		if (pr_id == 0) {
			writemathematica << "\n============================================\n" << mathe;
			writemathematica.flush();
			writemathematica.close();
		}
	}
	if (regtype.size() > 0) {
		for (int ireg = 0; ireg < (signed)regtype.size(); ireg++) {
			if (pr_id == 0) {
				writeanalysis_second << "reg" << ireg << "_" << regtype[ireg].name << "\t";
				writeanalysis_second.flush();
			}
			TERMANA objTERMANA;
			for (int iana = 0; iana < (signed)anadescriptor.size(); iana++) {
				objTERMANA.count = 0;
				objTERMANA.min_absz = VInfinite;
				objTERMANA.ave_absz = 0.0;
				objTERMANA.max_absz = -VInfinite;
				if (SkipMinMaxKurto != 2) {
					if (anadescriptor[iana].ASB != "U") {
						if (SkipMinMaxKurto == 0) {
							objTERMANA.name.clear();
							objTERMANA.name.push_back("Min" + anadescriptor[iana].name);
							regtype[ireg].termana.push_back(objTERMANA);
						}
						objTERMANA.name.clear();
						objTERMANA.name.push_back("Ave" + anadescriptor[iana].name);
						regtype[ireg].termana.push_back(objTERMANA);
						if (SkipMinMaxKurto == 0) {
							objTERMANA.name.clear();
							objTERMANA.name.push_back("Max" + anadescriptor[iana].name);
							regtype[ireg].termana.push_back(objTERMANA);
						}
						objTERMANA.name.clear();
						objTERMANA.name.push_back("StdDv" + anadescriptor[iana].name);
						regtype[ireg].termana.push_back(objTERMANA);
						if (SkipMinMaxKurto == 0) {
							objTERMANA.name.clear();
							objTERMANA.name.push_back("Kurto" + anadescriptor[iana].name);
							regtype[ireg].termana.push_back(objTERMANA);
						}
					}
					else {
						objTERMANA.name.clear();
						objTERMANA.name.push_back(anadescriptor[iana].name);
						regtype[ireg].termana.push_back(objTERMANA);
					}
				}
				else {
					objTERMANA.name.clear();
					objTERMANA.name.push_back(anadescriptor[iana].name);
					regtype[ireg].termana.push_back(objTERMANA);
				}
			}
			objTERMANA.name.clear();
			objTERMANA.name.push_back(AnionName + "M");
			regtype[ireg].termana.push_back(objTERMANA);
			for (int icat = 0; icat < (signed)CationNames.size(); icat++) {
				objTERMANA.name.clear();
				objTERMANA.name.push_back(CationNames[icat] + AnionName);
				regtype[ireg].termana.push_back(objTERMANA);
			}
		}
		if (NVariableLimit[0] > 1 && NVariableLimit[1] > 1) {
			Nfeatures = (signed)regtype[0].termana.size();
			estimate_jobsize(&xpr_ld_get, &jobsize, &jobsize_ld, &Nfeatures, NVariableLimit);
			set_pr_jobs(&xpr_ld_get, &jobsize, &jobsize_ld);
			initialize_occp(&arrange, &occup, &Nfeatures);
			run_through_models(&regtype, &xpr_ld_get, &Nfeatures, NVariableLimit);
		}
		bagtype.clear();
		for (int ireg = 0; ireg < (signed)regtype.size(); ireg++) {
			int saved_name = -1;
			for (int ibag = 0; ibag < (signed)bagtype.size(); ibag++) {
				if (bagtype[ibag].name == regtype[ireg].name) {
					saved_name = ibag;
					break;
				}
			}
			if (saved_name == -1) {
				BAGTYPE objBAGTYPE;
				objBAGTYPE.name = regtype[ireg].name;
				objBAGTYPE.which_regtype.push_back(ireg);
				objBAGTYPE.obs_v_max = -VInfinite;
				objBAGTYPE.obs_v_min = VInfinite;
				bagtype.push_back(objBAGTYPE);
			}
			else {
				bagtype[saved_name].which_regtype.push_back(ireg);
			}
		}
		for (int ibag = 0; ibag < (signed)bagtype.size(); ibag++) {
			if (pr_id == 0) {
				writeanalysis_second << "bag" << ibag << "_" << bagtype[ibag].name << "\tbag" << ibag << "_" << bagtype[ibag].name << "_stdev";
				writeanalysis_second.flush();
				if (ibag != (signed)bagtype.size() - 1) {
					writeanalysis_second << "\t";
					writeanalysis_second.flush();
				}
			}
			if (ibag < bUseAve.size()) {
				bagtype[ibag].user_defined_v_ave = bUseAve[ibag];
				if (bUseAve[ibag]) {
					bagtype[ibag].v_ave = UseAve[ibag];
				}
			}
			if (ibag < bUseStdDv.size()) {
				bagtype[ibag].user_defined_v_stddv = bUseStdDv[ibag];
				if (bUseStdDv[ibag]) {
					bagtype[ibag].v_stddv = UseStdDv[ibag];
				}
			}
			for (int it = 0; it < (signed)regtype[0].termana.size(); it++) {
				TERMANAFORBAG objTERMANAFORBAG;
				objTERMANAFORBAG.name = regtype[0].termana[it].name;
				objTERMANAFORBAG.min_absz = VInfinite;
				objTERMANAFORBAG.ave_absz = 0.0;
				objTERMANAFORBAG.max_absz = -VInfinite;
				objTERMANAFORBAG.metacount = 0;
				objTERMANAFORBAG.count = 0;
				bagtype[ibag].termanaforbag.push_back(objTERMANAFORBAG);
			}
			std::string secoeff_filename = bagtype[ibag].name + ".se.txt";
			std::ifstream readsecoeff;
			readsecoeff.open(secoeff_filename);
			if (readsecoeff.is_open()) {
				bagtype[ibag].stackingensemble.clear();
				if (pr_id == 0) {
					std::cout << "!! OVERWRITING Coefficients for stacking ensembles... " << bagtype[ibag].name << " from " << secoeff_filename << "\n";
					std::cout.flush();
				}
				std::getline(readsecoeff, sbuf);
				for (int kk = 0; kk < (signed)bagtype[ibag].which_regtype.size(); kk++) {
					STACKINGENSEMBLE objSTACKINGENSEMBLE;
					objSTACKINGENSEMBLE.coeff_sum = 0.0;
					for (int kkk = 0; kkk < kk + 2; kkk++) {
						if (kkk == 0 || kkk == 1) {
							readsecoeff >> sbuf;
						}
						else {
							readsecoeff >> sbuf >> sbuf;
						}
						long double ldbuf;
						readsecoeff >> ldbuf;
						objSTACKINGENSEMBLE.coeff.push_back(ldbuf);
						if (objSTACKINGENSEMBLE.coeff.size() > 1) {
							objSTACKINGENSEMBLE.coeff_sum += fabsl(ldbuf);
						}
					}
					std::getline(readsecoeff, sbuf);
					bagtype[ibag].stackingensemble.push_back(objSTACKINGENSEMBLE);
				}
			}
		}
		
	}
	if (PostTagNames.size() > 0) {
		for (int ipost = 0; ipost < (signed)PostTagNames.size(); ipost++) {
			if (pr_id == 0) {
				writeanalysis << "\t" << PostTagNames[ipost];
				writeanalysis.flush();
			}
		}
		if (regtype.size() > 0) {
			for (int ipost = 0; ipost < (signed)PostTagNames.size(); ipost++) {
				if (pr_id == 0) {
					writeanalysis_second << "\t" << PostTagNames[ipost];
					writeanalysis_second.flush();
				}
			}
		}
	}
	if (pr_id == 0) {
		writeanalysis << "\n";
		writeanalysis.flush();
		writeanalysis_second << "\n";
		writeanalysis_second.flush();
	}

	bool bUseOffset = false;
	if (UseOffset == 1) {
		bUseOffset = true;
	}

	bool headexists = true;
	if (HeadExists == 0) {
		headexists = false;
	}
	set_par(&MaterialsFilename, &headexists);
	bool collect_only = true;
	process_matdat(&zTargetVReg, &ZDistWeight, &TargetVReg, &TargetStdDvReg, &ModulateContentFactor, &DopeFactor, &MaterialsFilename, &regtype, &bagtype, &atomdat, &FaceFilename, &FaceFeatureName, &face, &strdat, &sitename, &anadescriptor, &AnionName, &Anions, &modulation_limit, &CationNames, &AnionValence, &VInfinite, &collect_only, &bUseOffset, &MaxIteration, &headexists, &usestrname, &PostTagNames, &DoNotAnalyze, &verbose, &AveOnly, &FoldMode, &NFolds, &TakeFold, &Train0Test1, &Straight0Rotate1, &TestModelStDevLimit);
	if (pr_id == 0) {
		writeanalysis.close();
	}
	if ((signed)regtype.size() > 0) {
		for (int ibag = 0; ibag < (signed)bagtype.size(); ibag++) {
			analyze_terms_zc(&bagtype[ibag], &ibag, &regtype, &RegModelFilename, &atomdat, &strdat, &anadescriptor, &VInfinite, &DoNotAnalyze, &AnionName, &CationNames, &verbose, &NCavailles, &AnalysisOutputFilename, &AnalysisMesh, &AveOnly);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (FurtherMaterialsDesign) {
			collect_only = false;
			process_matdat(&zTargetVReg, &ZDistWeight, &TargetVReg, &TargetStdDvReg, &ModulateContentFactor, &DopeFactor, &MaterialsFilename, &regtype, &bagtype, &atomdat, &FaceFilename, &FaceFeatureName, &face, &strdat, &sitename, &anadescriptor, &AnionName, &Anions, &modulation_limit, &CationNames, &AnionValence, &VInfinite, &collect_only, &bUseOffset, &MaxIteration, &headexists, &usestrname, &PostTagNames, &DoNotAnalyze, &verbose, &AveOnly, &FoldMode, &NFolds, &TakeFold, &Train0Test1, &Straight0Rotate1, &TestModelStDevLimit);
			MPI_Barrier(MPI_COMM_WORLD);
			if (pr_id == 0) {
				std::vector<std::string> sbufs;
				std::vector<std::string> sbufs_best;
				long double min_zdist = 100000000000000000000000000000000.0;
				long double current_zdist = 0.0;
				bool sbufs_update = false;
				for (int iprid = 0; iprid < pr_size; iprid++) {
					readoutputpr.open("MD." + std::to_string(iprid) + "." + OutputFilename);
					if (readoutputpr.is_open()) {
						sbufs.clear();
						while (true) {
							std::getline(readoutputpr, sbuf);
							sbufs.push_back(sbuf);
							if (sbuf.find("zDist=") != std::string::npos) {
								std::stringstream lcstr;
								lcstr.str("");
								lcstr.clear();
								lcstr << sbuf;
								while (true) {
									lcstr >> sbuf;
									if (sbuf == "zDist=") {
										lcstr >> current_zdist;
										if (current_zdist < min_zdist) {
											min_zdist = current_zdist;
											sbufs_update = true;
										}
										break;
									}
								}
							}
							if (sbuf == "=================================================================================") {
								for (int iline = 0; iline < (signed)sbufs.size(); iline++) {
									writeoutput << sbufs[iline] << "\n";
									writeoutput.flush();
								}
								
								if (sbufs_update) {
									sbufs_best.clear();
									for (int iline = 0; iline < (signed)sbufs.size(); iline++) {
										sbufs_best.push_back(sbufs[iline]);
									}
								}
								sbufs.clear();

								sbufs_update = false;
							}
							if (readoutputpr.eof()) {
								break;
							}
						}
					}
					for (int iline = 0; iline < (signed)sbufs_best.size(); iline++) {
						writeoutput << sbufs_best[iline] << "\n";
						writeoutput.flush();
					}
					readoutputpr.close();
					std::vector<std::string>().swap(sbufs);
					std::vector<std::string>().swap(sbufs_best);
				}
				char cmd[500];
				std::string masterrm = "del MD.*." + OutputFilename;
#ifdef _WIN32
				strcpy_s(cmd, sizeof(cmd), masterrm.c_str());
				system(cmd);
#elif _WIN64
				strcpy_s(cmd, sizeof(cmd), masterrm.c_str());
				system(cmd);
#else
				masterrm = "rm -f MD.*." + OutputFilename;
				strcpy(cmd, masterrm.c_str());
				system(cmd);
#endif
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	//if (FurtherMaterialsDesign && mapout) {
		//process_mapout(&zTargetVReg, &ZDistWeight, &TargetVReg, &MapOutSpaceFilename, &regtype, &bagtype, &atomdat, &strdat, &anadescriptor, &VInfinite, &DoNotAnalyze, &CationNames, &verbose, &AveOnly);
	//}
	
	if (pr_id == 0) {
		writeoutput.close();
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;

}
