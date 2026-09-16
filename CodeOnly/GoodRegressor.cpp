#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <complex>
#include <set>
#include <chrono>
#include <mpi.h>
#include "Eigen/Dense"

#define PI 3.14159265359

typedef struct {

	std::string setX;
	bool negative_exists;
	long double min, ave, max, stdev;

}GROUPX;

std::vector<GROUPX> GroupX;

typedef struct {

	std::string name;
	std::vector<long double> v;
	std::vector<int> yhit;

}DATASET;

std::vector<DATASET> dataset, testdataset;
std::vector<DATASET> specialdataset, specialdataset_quad;

typedef struct {

	Eigen::MatrixXd MatrixX;
	Eigen::MatrixXd MatrixY;
	Eigen::MatrixXd MatrixYexpected;
	Eigen::MatrixXd MatrixXTXi;
	Eigen::MatrixXd MatrixBeta;
	Eigen::MatrixXd MatrixBetaStandardized;
	Eigen::MatrixXd ShuffledMatrixBetaStandardized;
	std::vector<long double> fit_coeff, fit_coeff_standardized, fit_coeff_SE, fit_coeff_t, fit_coeff_p, fit_coeff_CI_start, fit_coeff_CI_end;
	long double SStot, SSres;
	long double MSE, RMSE, MAE;
	long double det_r2, det_radj2;
	int df;
	long double Fvalue, p_Fvalue;
	bool full_Frequentist;
	std::vector<int> fit_coeff_porder, fit_coeff_porder_on_arr;
	std::vector<long double> fit_coeff_pordered;

	Eigen::MatrixXd MatrixYB;
	Eigen::MatrixXd MatrixZB;
	Eigen::MatrixXd MatrixBetaB;
	Eigen::MatrixXd MatrixInitBetaB;
	Eigen::MatrixXd MatrixYexpectedB;
	Eigen::MatrixXd MatrixWB;
	Eigen::MatrixXd MatrixDB;
	Eigen::MatrixXd MatrixcB;
	Eigen::MatrixXd MatrixTB;
	Eigen::MatrixXd MatrixFisherKubetabetaB;
	Eigen::MatrixXd MatrixFisherKuphiphiB;
	long double fit_phiB, fit_phi_SEB, fit_phi_tB, fit_phi_pB, fit_phi_CI_startB, fit_phi_CI_endB;
	long double initphiB;
	long double logbetadenseB;
	std::vector<long double> fit_coeffB, fit_coeff_standardizedB, fit_coeff_SEB, fit_coeff_tB, fit_coeff_pB, fit_coeff_CI_startB, fit_coeff_CI_endB;
	long double SStotB, SSresB;
	long double det_r2_pseudo;
	std::vector<int> fit_coeff_porderB, fit_coeff_porder_on_arrB;
	std::vector<long double> fit_coeff_porderedB;

}REGELE;

REGELE RegEle, GoodRegEle, TestRegEle, TestGoodRegEle, specialRegEle, specialRegEle_quad, specialTestRegEle, specialTestRegEle_quad;

typedef struct {

	std::vector<std::string> userdefined_name;
	std::string FName;
	std::vector<long double> coeff;

}FUNCTYPE;

std::vector<FUNCTYPE> FuncType, PlainFuncType;

typedef struct {

	std::string setY;
	bool minlimitY, maxlimitY;
	long double minY, maxY;

}INSTALLY;

std::vector<INSTALLY> instally;

typedef struct {

	int NVariableLimit[2];
	bool strict_limit;
	int interact;

}INSTALLX;

std::vector<INSTALLX> installx;

typedef struct {

	std::string setX;
	bool minlimitX, maxlimitX;
	long double minX, maxX;

}SUBFILTERX;

std::vector<SUBFILTERX> subfilterx;

typedef struct {

	int master;
	std::vector<int> slave;

}HEGEL;

std::vector<HEGEL> Hegel;

int pr_size, ierr, pr_id;
signed long long int* pr_job_size;
signed long long int** pr_job;
long double* pr_job_size_ld;
long double** pr_job_ld;
signed long long int jobsize;
long double jobsize_ld;
bool xpr_ld_get;
bool finite_job;
auto start_t = std::chrono::system_clock::now();
int writeoutput_prid = 0;

std::vector<std::string> DataFilename;
std::string OutputFilename;
int OutputPrecision;
int tFTestMesh;
int NVariableLimit[2];
long double RegressionCoeffConfindenceIntervalLevel;
long double PickUpPLevel;
int Verbose;
long double TargetCPUTimeBeforeSwap;
long double TargetCPUTimeAfterSwap;
bool useTargetCPUTimeAfterSwap = false;
bool useExitMode = false;
bool useSomeShotMode = false;
bool useSuperFastSwap = false;
long double TrainRatio;
long double R2Threshold;
int MaxBFGSIterations;
int DuplicateX;
int NameDuplicateHandling;
long double Dmin, Dmax, Dincrement;
int NTestRobustness;
std::vector<std::string> betatype;

signed long long int* SigmaMinusXk;
long double* SigmaMinusXk_ld;
int Nfeatures;
std::vector<int> datadivisionflag;
std::vector<int> arrange;
std::vector<int> occup;
std::vector<signed long long int> jobsize_step;
std::vector<long double> jobsize_step_ld;
std::vector<std::string> column_str;
bool FullFrequentist_exists = false;
int good_occup;
std::vector<int> good_arrange;
long double max_R2 = -1000000.0;
long double max_testR2 = -1000000.0;
long double max_rR2 = -1000000.0;
bool strict_limit = false;
long double score_slope = 1.0;
long double score_intercept = 0.0;
int failure_count = 0;
int max_failure_count = 2;
bool use_dataname = false;
std::string userdefinedName;

bool ExistingCholskeySeries = false;
Eigen::MatrixXd ExistingShuffledStandardizedMatrixX;
Eigen::LLT<Eigen::MatrixXd> ExistingLLT;
std::vector<long double> Existing_x_average_series, Existing_x_stdev_series;
int CholskeyMonitor = 0;

std::string sbuf;
std::ifstream readpara, readdata, readinitial, readregele, readreport, readsub;
std::ofstream writeoutput, writeoutputpr, writeregele, writereport;

void write_matrixbin(Eigen::MatrixXd* matrix, FILE* writeregelebin) {

	int row = matrix->rows();
	fwrite(&row, sizeof(int), 1, writeregelebin);
	int col = matrix->cols();
	fwrite(&col, sizeof(int), 1, writeregelebin);
	for (int i = 0; i < matrix->rows(); i++) {
		for (int j = 0; j < matrix->cols(); j++) {
			long double ele = (*matrix)(i, j);
			fwrite(&ele, sizeof(long double), 1, writeregelebin);
		}
	}

}

void write_ldsbin(std::vector<long double>* ld, FILE* writeregelebin) {

	int s = (signed)ld->size();
	fwrite(&s, sizeof(int), 1, writeregelebin);
	for (int i = 0; i < ld->size(); i++) {
		fwrite(&(*ld)[i], sizeof(long double), 1, writeregelebin);
	}

}

void write_isbin(std::vector<int>* ib, FILE* writeregelebin) {

	int s = (signed)ib->size();
	fwrite(&s, sizeof(int), 1, writeregelebin);
	for (int i = 0; i < ib->size(); i++) {
		fwrite(&(*ib)[i], sizeof(int), 1, writeregelebin);
	}

}

void write_ldbin(long double* ld, FILE* writeregelebin) {

	fwrite(ld, sizeof(long double), 1, writeregelebin);

}

void write_ibin(int* i, FILE* writeregelebin) {

	fwrite(i, sizeof(int), 1, writeregelebin);

}

void write_bbin(const bool* b, FILE* f) {

	unsigned char tmp = *b ? 1 : 0;
	fwrite(&tmp, sizeof(tmp), 1, f);

}

void write_strbin(std::string* str, FILE* writeregelebin) {

	size_t size = str->size();
	fwrite(&size, sizeof(size), 1, writeregelebin);
	fwrite(str->data(), sizeof(char), size, writeregelebin);

}

void read_matrixbin(Eigen::MatrixXd* matrix, FILE* readregelebin) {

	int row;
	fread(&row, sizeof(int), 1, readregelebin);
	int col;
	fread(&col, sizeof(int), 1, readregelebin);

	*matrix = Eigen::MatrixXd::Random(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			long double elebuf;
			fread(&elebuf, sizeof(long double), 1, readregelebin);
			(*matrix)(i, j) = elebuf;
		}
	}

}

void read_ldsbin(std::vector<long double>* ld, FILE* readregelebin) {

	int s;
	fread(&s, sizeof(int), 1, readregelebin);
	std::vector<long double>().swap(*ld);
	ld->clear();
	for (int i = 0; i < s; i++) {
		long double ldbuf;
		fread(&ldbuf, sizeof(long double), 1, readregelebin);
		ld->push_back(ldbuf);
	}

}

void read_isbin(std::vector<int>* ib, FILE* readregelebin) {

	int s;
	fread(&s, sizeof(int), 1, readregelebin);
	std::vector<int>().swap(*ib);
	ib->clear();
	for (int i = 0; i < s; i++) {
		int ibuf;
		fread(&ibuf, sizeof(int), 1, readregelebin);
		ib->push_back(ibuf);
	}

}

void read_ldbin(long double* ld, FILE* readregelebin) {

	fread(ld, sizeof(long double), 1, readregelebin);

}

void read_ibin(int* i, FILE* readregelebin) {

	fread(i, sizeof(int), 1, readregelebin);

}

void read_bbin(bool* b, FILE* f) {

	unsigned char tmp;
	fread(&tmp, sizeof(tmp), 1, f);
	*b = (tmp != 0);

}

void read_strbin(std::string* str, FILE* readregelebin) {

	size_t size;
	fread(&size, sizeof(size), 1, readregelebin);
	str->resize(size);
	fread(&(*str)[0], 1, size, readregelebin);

}

void write_regelebin(REGELE* RegEle, std::string* filename) {

	FILE* writeregelebin = nullptr;
	const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
	errno_t error;
	error = fopen_s(&writeregelebin, ccfilename, "wb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&writeregelebin, ccfilename, "wb");
#else
	writeregelebin = fopen(ccfilename, "wb");
#endif

	write_matrixbin(&RegEle->MatrixX, writeregelebin);
	write_matrixbin(&RegEle->MatrixY, writeregelebin);
	write_matrixbin(&RegEle->MatrixYexpected, writeregelebin);
	write_matrixbin(&RegEle->MatrixXTXi, writeregelebin);
	write_matrixbin(&RegEle->MatrixBeta, writeregelebin);
	write_matrixbin(&RegEle->MatrixBetaStandardized, writeregelebin);
	write_matrixbin(&RegEle->ShuffledMatrixBetaStandardized, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_standardized, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_SE, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_t, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_p, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_CI_start, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_CI_end, writeregelebin);
	write_ldbin(&RegEle->SStot, writeregelebin);
	write_ldbin(&RegEle->SSres, writeregelebin);
	write_ldbin(&RegEle->MSE, writeregelebin);
	write_ldbin(&RegEle->RMSE, writeregelebin);
	write_ldbin(&RegEle->MAE, writeregelebin);
	write_ldbin(&RegEle->det_r2, writeregelebin);
	write_ldbin(&RegEle->det_radj2, writeregelebin);
	write_ibin(&RegEle->df, writeregelebin);
	write_ldbin(&RegEle->Fvalue, writeregelebin);
	write_ldbin(&RegEle->p_Fvalue, writeregelebin);
	write_bbin(&RegEle->full_Frequentist, writeregelebin);
	write_isbin(&RegEle->fit_coeff_porder, writeregelebin);
	write_isbin(&RegEle->fit_coeff_porder_on_arr, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_pordered, writeregelebin);

	write_matrixbin(&RegEle->MatrixYB, writeregelebin);
	write_matrixbin(&RegEle->MatrixZB, writeregelebin);
	write_matrixbin(&RegEle->MatrixBetaB, writeregelebin);
	write_matrixbin(&RegEle->MatrixInitBetaB, writeregelebin);
	write_matrixbin(&RegEle->MatrixYexpectedB, writeregelebin);
	write_matrixbin(&RegEle->MatrixWB, writeregelebin);
	write_matrixbin(&RegEle->MatrixDB, writeregelebin);
	write_matrixbin(&RegEle->MatrixcB, writeregelebin);
	write_matrixbin(&RegEle->MatrixTB, writeregelebin);
	write_matrixbin(&RegEle->MatrixFisherKubetabetaB, writeregelebin);
	write_matrixbin(&RegEle->MatrixFisherKuphiphiB, writeregelebin);
	write_ldbin(&RegEle->fit_phiB, writeregelebin);
	write_ldbin(&RegEle->fit_phi_SEB, writeregelebin);
	write_ldbin(&RegEle->fit_phi_tB, writeregelebin);
	write_ldbin(&RegEle->fit_phi_pB, writeregelebin);
	write_ldbin(&RegEle->fit_phi_CI_startB, writeregelebin);
	write_ldbin(&RegEle->fit_phi_CI_endB, writeregelebin);
	write_ldbin(&RegEle->initphiB, writeregelebin);
	write_ldbin(&RegEle->logbetadenseB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeffB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_standardizedB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_SEB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_tB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_pB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_CI_startB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_CI_endB, writeregelebin);
	write_ldbin(&RegEle->SStotB, writeregelebin);
	write_ldbin(&RegEle->SSresB, writeregelebin);
	write_ldbin(&RegEle->det_r2_pseudo, writeregelebin);
	write_isbin(&RegEle->fit_coeff_porderB, writeregelebin);
	write_isbin(&RegEle->fit_coeff_porder_on_arrB, writeregelebin);
	write_ldsbin(&RegEle->fit_coeff_porderedB, writeregelebin);

	fclose(writeregelebin);

}

bool read_regelebin(REGELE* RegEle, std::string* filename) {

	std::ifstream readtest;
	readtest.open(*filename);
	bool readable = false;
	if (readtest.is_open()) {
		readable = true;
	}
	readtest.close();

	if (readable) {


		FILE* readregelebin = nullptr;
		const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
		errno_t error;
		error = fopen_s(&readregelebin, ccfilename, "rb");
#elif _WIN64
		errno_t error;
		error = fopen_s(&readregelebin, ccfilename, "rb");
#else
		readregelebin = fopen(ccfilename, "rb");
#endif

		read_matrixbin(&RegEle->MatrixX, readregelebin);
		read_matrixbin(&RegEle->MatrixY, readregelebin);
		read_matrixbin(&RegEle->MatrixYexpected, readregelebin);
		read_matrixbin(&RegEle->MatrixXTXi, readregelebin);
		read_matrixbin(&RegEle->MatrixBeta, readregelebin);
		read_matrixbin(&RegEle->MatrixBetaStandardized, readregelebin);
		read_matrixbin(&RegEle->ShuffledMatrixBetaStandardized, readregelebin);
		read_ldsbin(&RegEle->fit_coeff, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_standardized, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_SE, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_t, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_p, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_CI_start, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_CI_end, readregelebin);
		read_ldbin(&RegEle->SStot, readregelebin);
		read_ldbin(&RegEle->SSres, readregelebin);
		read_ldbin(&RegEle->MSE, readregelebin);
		read_ldbin(&RegEle->RMSE, readregelebin);
		read_ldbin(&RegEle->MAE, readregelebin);
		read_ldbin(&RegEle->det_r2, readregelebin);
		read_ldbin(&RegEle->det_radj2, readregelebin);
		read_ibin(&RegEle->df, readregelebin);
		read_ldbin(&RegEle->Fvalue, readregelebin);
		read_ldbin(&RegEle->p_Fvalue, readregelebin);
		read_bbin(&RegEle->full_Frequentist, readregelebin);
		read_isbin(&RegEle->fit_coeff_porder, readregelebin);
		read_isbin(&RegEle->fit_coeff_porder_on_arr, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_pordered, readregelebin);

		read_matrixbin(&RegEle->MatrixYB, readregelebin);
		read_matrixbin(&RegEle->MatrixZB, readregelebin);
		read_matrixbin(&RegEle->MatrixBetaB, readregelebin);
		read_matrixbin(&RegEle->MatrixInitBetaB, readregelebin);
		read_matrixbin(&RegEle->MatrixYexpectedB, readregelebin);
		read_matrixbin(&RegEle->MatrixWB, readregelebin);
		read_matrixbin(&RegEle->MatrixDB, readregelebin);
		read_matrixbin(&RegEle->MatrixcB, readregelebin);
		read_matrixbin(&RegEle->MatrixTB, readregelebin);
		read_matrixbin(&RegEle->MatrixFisherKubetabetaB, readregelebin);
		read_matrixbin(&RegEle->MatrixFisherKuphiphiB, readregelebin);
		read_ldbin(&RegEle->fit_phiB, readregelebin);
		read_ldbin(&RegEle->fit_phi_SEB, readregelebin);
		read_ldbin(&RegEle->fit_phi_tB, readregelebin);
		read_ldbin(&RegEle->fit_phi_pB, readregelebin);
		read_ldbin(&RegEle->fit_phi_CI_startB, readregelebin);
		read_ldbin(&RegEle->fit_phi_CI_endB, readregelebin);
		read_ldbin(&RegEle->initphiB, readregelebin);
		read_ldbin(&RegEle->logbetadenseB, readregelebin);
		read_ldsbin(&RegEle->fit_coeffB, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_standardizedB, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_SEB, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_tB, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_pB, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_CI_startB, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_CI_endB, readregelebin);
		read_ldbin(&RegEle->SStotB, readregelebin);
		read_ldbin(&RegEle->SSresB, readregelebin);
		read_ldbin(&RegEle->det_r2_pseudo, readregelebin);
		read_isbin(&RegEle->fit_coeff_porderB, readregelebin);
		read_isbin(&RegEle->fit_coeff_porder_on_arrB, readregelebin);
		read_ldsbin(&RegEle->fit_coeff_porderedB, readregelebin);

		fclose(readregelebin);
	
	} 

	return readable;

}

void write_datasetbin(std::vector<DATASET>* dset, std::string* filename) {

	FILE* writedataset = nullptr;
	const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
	errno_t error;
	error = fopen_s(&writedataset, ccfilename, "wb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&writedataset, ccfilename, "wb");
#else
	writedataset = fopen(ccfilename, "wb");
#endif

	int s = dset->size();	
	write_ibin(&s, writedataset);
	for (int i = 0; i < s; i++) {
		write_strbin(&(*dset)[i].name, writedataset);
		write_ldsbin(&(*dset)[i].v, writedataset);
		write_isbin(&(*dset)[i].yhit, writedataset);
	}

	fclose(writedataset);

}

bool read_datasetbin(std::vector<DATASET>* dset, std::string* filename) {

	std::ifstream readtest;
	readtest.open(*filename);
	bool readable = false;
	if (readtest.is_open()) {
		readable = true;
	}
	readtest.close();

	if (readable) {
	
		FILE* readdataset = nullptr;
		const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
		errno_t error;
		error = fopen_s(&readdataset, ccfilename, "rb");
#elif _WIN64
		errno_t error;
		error = fopen_s(&readdataset, ccfilename, "rb");
#else
		readdataset = fopen(ccfilename, "rb");
#endif

		std::vector<DATASET>().swap(*dset);
		dset->clear();

		int s;
		read_ibin(&s, readdataset);
		for (int i = 0; i < s; i++) {
			DATASET objDATASET;
			read_strbin(&objDATASET.name, readdataset);
			read_ldsbin(&objDATASET.v, readdataset);
			read_isbin(&objDATASET.yhit, readdataset);
			dset->push_back(objDATASET);
		}

		fclose(readdataset);

	}

	return readable;

}

void write_arrangebin(std::vector<int>* arrange, std::string* filename) {

	FILE* writearrangebin = nullptr;
	const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
	errno_t error;
	error = fopen_s(&writearrangebin, ccfilename, "wb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&writearrangebin, ccfilename, "wb");
#else
	writearrangebin = fopen(ccfilename, "wb");
#endif

	write_isbin(arrange, writearrangebin);
	
	fclose(writearrangebin);

}

bool read_arrangebin(std::vector<int>* arrange, std::string* filename) {

	std::ifstream readtest;
	readtest.open(*filename);
	bool readable = false;
	if (readtest.is_open()) {
		readable = true;
	}
	readtest.close();

	if (readable) {

		FILE* readarrangebin = nullptr;
		const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
		errno_t error;
		error = fopen_s(&readarrangebin, ccfilename, "rb");
#elif _WIN64
		errno_t error;
		error = fopen_s(&readarrangebin, ccfilename, "rb");
#else
		readarrangebin = fopen(ccfilename, "rb");
#endif

		read_isbin(arrange, readarrangebin);

		fclose(readarrangebin);

	}

	return readable;

}

void write_Noccupbin(int* Noccup, std::string* filename) {

	FILE* writeNoccupbin = nullptr;
	const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
	errno_t error;
	error = fopen_s(&writeNoccupbin, ccfilename, "wb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&writeNoccupbin, ccfilename, "wb");
#else
	writeNoccupbin = fopen(ccfilename, "wb");
#endif

	write_ibin(Noccup, writeNoccupbin);

	fclose(writeNoccupbin);

}

bool read_Noccupbin(int* Noccup, std::string* filename) {

	std::ifstream readtest;
	readtest.open(*filename);
	bool readable = false;
	if (readtest.is_open()) {
		readable = true;
	}
	readtest.close();

	if (readable) {

		FILE* readNoccupbin = nullptr;
		const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
		errno_t error;
		error = fopen_s(&readNoccupbin, ccfilename, "rb");
#elif _WIN64
		errno_t error;
		error = fopen_s(&readNoccupbin, ccfilename, "rb");
#else
		readNoccupbin = fopen(ccfilename, "rb");
#endif

		read_ibin(Noccup, readNoccupbin);

		fclose(readNoccupbin);

	}

	return readable;

}

void write_maxR2bin(long double* max_R2, long double* max_testR2, long double* max_rR2, std::string* filename) {

	FILE* writemaxR2 = nullptr;
	const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
	errno_t error;
	error = fopen_s(&writemaxR2, ccfilename, "wb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&writemaxR2, ccfilename, "wb");
#else
	writemaxR2 = fopen(ccfilename, "wb");
#endif

	write_ldbin(max_R2, writemaxR2);
	write_ldbin(max_testR2, writemaxR2);
	write_ldbin(max_rR2, writemaxR2);

	fclose(writemaxR2);

}

bool read_maxR2bin(long double* max_R2, long double* max_testR2, long double* max_rR2, std::string* filename) {

	std::ifstream readtest;
	readtest.open(*filename);
	bool readable = false;
	if (readtest.is_open()) {
		readable = true;
	}
	readtest.close();

	if (readable) {

		FILE* readmaxR2 = nullptr;
		const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
		errno_t error;
		error = fopen_s(&readmaxR2, ccfilename, "rb");
#elif _WIN64
		errno_t error;
		error = fopen_s(&readmaxR2, ccfilename, "rb");
#else
		readmaxR2 = fopen(ccfilename, "rb");
#endif

		read_ldbin(max_R2, readmaxR2);
		read_ldbin(max_testR2, readmaxR2);
		read_ldbin(max_rR2, readmaxR2);

		fclose(readmaxR2);

	}

	return readable;

}

void write_NNendbin(int* NNend, std::string* filename) {

	FILE* writeNNendbin = nullptr;
	const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
	errno_t error;
	error = fopen_s(&writeNNendbin, ccfilename, "wb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&writeNNendbin, ccfilename, "wb");
#else
	writeNNendbin = fopen(ccfilename, "wb");
#endif

	write_ibin(NNend, writeNNendbin);

	fclose(writeNNendbin);

}

bool read_NNendbin(int* NNend, std::string* filename) {

	std::ifstream readtest;
	readtest.open(*filename);
	bool readable = false;
	if (readtest.is_open()) {
		readable = true;
	}
	readtest.close();

	if (readable) {

		FILE* readNNendbin = nullptr;
		const char* ccfilename = (*filename).c_str();

#ifdef _WIN32
		errno_t error;
		error = fopen_s(&readNNendbin, ccfilename, "rb");
#elif _WIN64
		errno_t error;
		error = fopen_s(&readNNendbin, ccfilename, "rb");
#else
		readNNendbin = fopen(ccfilename, "rb");
#endif

		read_ibin(NNend, readNNendbin);

		fclose(readNNendbin);

	}

	return readable;

}

bool is_number(const std::string* s) {
	std::istringstream iss(*s);
	double d;
	return iss >> std::noskipws >> d && iss.eof();
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

long double get_derv_special(std::string* ftype, long double* x) {

	long double dx = 1.0 / (long double)tFTestMesh;
	long double xpdx = *x + dx;
	return((get_v_special(ftype, &xpdx) - get_v_special(ftype, x)) / dx);

}

long double get_v(FUNCTYPE* functype, long double* x) {

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

long double get_trigamma(long double* v) {

	return(1.0 / *v + 1.0 / (2.0 * powl(*v, 2.0)) + 1.0 / (6.0 * powl(*v, 3.0)) - 1.0 / (30.0 * powl(*v, 5.0)));

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

	std::vector<int>().swap(up);
	std::vector<int>().swap(down);
	std::vector<bool>().swap(down_chk);

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

	std::vector<int>().swap(up);
	std::vector<int>().swap(down);
	std::vector<bool>().swap(down_chk);

	return bino;

}

void memcpy_for_vectorint(std::vector<int>* to_array, std::vector<int>* from_array) {

	std::vector<int>().swap((*to_array));
	for (int i = 0; i < (signed)from_array->size(); i++) {
		to_array->push_back((*from_array)[i]);
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

void kill_Jekyll(std::vector<int>* arrange, int* DuplicateX, std::vector<FUNCTYPE>* functype) {

	if (*DuplicateX > 1) {
		for (int i = 0; i < arrange->size(); i += *DuplicateX) {
			for (int j = *DuplicateX - 1; j >= 0; j--) {
				if ((*arrange)[i + j] != 0) {
					bool* funcchk;
					funcchk = new bool[functype->size()];
					for (int ifu = 0; ifu < functype->size(); ifu++) {
						funcchk[ifu] = false;
					}
					for (int k = 0; k < *DuplicateX; k++) {
						if (j != k && (*arrange)[i + k] != 0) {
							funcchk[(*arrange)[i + k] - 1] = true;
						}
					}
					if (funcchk[(*arrange)[i + j] - 1]) {
						for (int ifu = 0; ifu < functype->size(); ifu++) {
							if (funcchk[ifu] == false) {
								(*arrange)[i + j] = ifu + 1;
								break;
							}
						}
					}
					delete[] funcchk;
					funcchk = nullptr;
				}
			}
		}
	}

}

int find_column(std::string* target, std::vector<std::string>* column_str) {

	int itarget = -1;
	for (unsigned int i = 0; i < column_str->size(); i++) {
		if ((*column_str)[i] == *target) {
			itarget = (signed)i;
			//break;
		}
	}

	return itarget;

}

void read_para() {

	readpara.open("basic_properties_GoodRegressor.txt");

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	std::stringstream localstr;
	localstr << sbuf;
	while (true) {
		localstr >> sbuf;
		DataFilename.push_back(sbuf);
		if (localstr.eof()) {
			break;
		}
	}

	readpara >> sbuf >> OutputFilename;
	readpara >> sbuf >> OutputPrecision;
	readpara >> sbuf >> tFTestMesh;
	std::getline(readpara, sbuf);

	std::getline(readpara, sbuf);
	localstr.clear();
	localstr << sbuf;
	localstr >> sbuf;
	while (true) {
		INSTALLX objINSTALLX;
		objINSTALLX.interact = 0;
		localstr >> sbuf;
		std::size_t pos = sbuf.find('i');
		while (pos != std::string::npos) {
			if (installx.size() != 0) {
				objINSTALLX.interact = 1;
			}
			sbuf.erase(pos, 1);
			pos = sbuf.find('i');
		}
		pos = sbuf.find('I');
		while (pos != std::string::npos) {
			if (installx.size() != 0) {
				objINSTALLX.interact = 2;
			}
			sbuf.erase(pos, 1);
			pos = sbuf.find('I');
		}
		std::stringstream stoistream;
		stoistream << sbuf;
		stoistream >> objINSTALLX.NVariableLimit[1];
		objINSTALLX.strict_limit = false;
		if (objINSTALLX.NVariableLimit[1] < 0) {
			objINSTALLX.strict_limit = true;
			objINSTALLX.NVariableLimit[1] *= -1;
		}
		objINSTALLX.NVariableLimit[0] = objINSTALLX.NVariableLimit[1];
		installx.push_back(objINSTALLX);
		if (localstr.eof()) {
			break;
		}
	}
	/*readpara >> sbuf >> NVariableLimit[0] >> NVariableLimit[1];
	if (NVariableLimit[1] < 0) {
		strict_limit = true;
		NVariableLimit[1] *= -1;
	}*/

	readpara >> sbuf >> PickUpPLevel;
	RegressionCoeffConfindenceIntervalLevel = 1.0 - PickUpPLevel;
	readpara >> sbuf >> Verbose;
	std::getline(readpara, sbuf);

	std::getline(readpara, sbuf);
	localstr.clear();
	localstr << sbuf;
	localstr >> sbuf >> TargetCPUTimeBeforeSwap;
	if (TargetCPUTimeBeforeSwap < 0) {
		TargetCPUTimeBeforeSwap *= -1.0;
		useSuperFastSwap = true;
	}
	if (!localstr.eof()) {
		localstr >> TargetCPUTimeAfterSwap;
		if (TargetCPUTimeAfterSwap > 0) {
			useTargetCPUTimeAfterSwap = true;
		}
		if (!localstr.eof()) {
			localstr >> sbuf;
			if (sbuf == "!") {
				useExitMode = true;
				useSomeShotMode = false;
			}
			else if (sbuf == "!!") {
				useExitMode = false;
				useSomeShotMode = true;
			}
			else if (sbuf == "!!!") {
				useExitMode = true;
				useSomeShotMode = true;
			}
		}
	}

	readpara >> sbuf >> TrainRatio;
	readpara >> sbuf >> R2Threshold;
	readpara >> sbuf >> MaxBFGSIterations;
	readpara >> sbuf >> DuplicateX;
	readpara >> sbuf >> NameDuplicateHandling;
	if (NameDuplicateHandling == 2) {
		readpara >> Dmin >> Dmax >> Dincrement;
	}
	readpara >> sbuf >> NTestRobustness;
	std::getline(readpara, sbuf);

	std::getline(readpara, sbuf);
	localstr.clear();
	localstr << sbuf;
	localstr >> sbuf;
	while (true) {
		INSTALLY objINSTALLY;
		localstr >> objINSTALLY.setY;
		localstr >> sbuf;
		if (sbuf == "*") {
			objINSTALLY.minlimitY = false;
		}
		else {
			objINSTALLY.minlimitY = true;
			std::stringstream llocalstr;
			llocalstr << sbuf;
			llocalstr >> objINSTALLY.minY;
		}
		localstr >> sbuf;
		if (sbuf == "*") {
			objINSTALLY.maxlimitY = false;
		}
		else {
			objINSTALLY.maxlimitY = true;
			std::stringstream llocalstr;
			llocalstr << sbuf;
			llocalstr >> objINSTALLY.maxY;
		}
		instally.push_back(objINSTALLY);
		if (localstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	localstr.clear();
	localstr << sbuf;
	while (true) {
		GROUPX objGROUPX;
		localstr >> objGROUPX.setX;
		objGROUPX.negative_exists = false;
		objGROUPX.min = 1000000000000000000000000000000000000000.0;
		objGROUPX.ave = 0.0;
		objGROUPX.max = -objGROUPX.min;
		objGROUPX.stdev = 0.0;
		for (int ix = 0; ix < DuplicateX; ix++) {
			GroupX.push_back(objGROUPX);
		}
		if (localstr.eof()) {
			break;
		}
	}

	readpara >> sbuf;
	std::getline(readpara, sbuf);
	if (sbuf.find("-") == std::string::npos) {
		localstr.clear();
		localstr << sbuf;
		while (true) {
			SUBFILTERX objSUBFILTERX;
			localstr >> objSUBFILTERX.setX;
			localstr >> sbuf;
			if (sbuf == "*") {
				objSUBFILTERX.minlimitX = false;
			}
			else {
				objSUBFILTERX.minlimitX = true;
				std::stringstream llocalstr;
				llocalstr << sbuf;
				llocalstr >> objSUBFILTERX.minX;
			}
			localstr >> sbuf;
			if (sbuf == "*") {
				objSUBFILTERX.maxlimitX = false;
			}
			else {
				objSUBFILTERX.maxlimitX = true;
				std::stringstream llocalstr;
				llocalstr << sbuf;
				llocalstr >> objSUBFILTERX.maxX;
			}
			subfilterx.push_back(objSUBFILTERX);
			if (localstr.eof()) {
				break;
			}
		}
	}

	std::getline(readpara, sbuf);

	while (true) {
		std::getline(readpara, sbuf);
		if (sbuf.find("BetaRegression") == std::string::npos) {
			FUNCTYPE objFUNCTYPE;
			localstr.clear();
			localstr << sbuf;
			localstr >> sbuf >> sbuf;
			std::stringstream getname(sbuf);
			std::string seg;
			while (std::getline(getname, seg, '^')) {
				objFUNCTYPE.userdefined_name.push_back(seg);
			}
			localstr >> objFUNCTYPE.FName;
			int NCoeff = 0;
			if (objFUNCTYPE.FName == "pow" || objFUNCTYPE.FName == "abs") {
				NCoeff = 2;
			}
			else if (objFUNCTYPE.FName == "log" || objFUNCTYPE.FName == "erf" || objFUNCTYPE.FName == "erfc") {
				NCoeff = 3;
			}
			else if (objFUNCTYPE.FName == "sin" || objFUNCTYPE.FName == "cos" || objFUNCTYPE.FName == "tan" || objFUNCTYPE.FName == "sinh" || objFUNCTYPE.FName == "cosh" || objFUNCTYPE.FName == "tanh") {
				NCoeff = 4;
			}
			else if (objFUNCTYPE.FName == "exp") {
				NCoeff = 5;
			}
			for (int j = 0; j < NCoeff; j++) {
				long double ldbuf;
				localstr >> sbuf >> ldbuf;
				objFUNCTYPE.coeff.push_back(ldbuf);
			}
			FuncType.push_back(objFUNCTYPE);
		}
		else {
			break;
		}
	}

	FUNCTYPE objpFUNCTYPE;
	objpFUNCTYPE.userdefined_name.push_back("");
	objpFUNCTYPE.userdefined_name.push_back("1");
	objpFUNCTYPE.FName = "pow";
	objpFUNCTYPE.coeff.push_back(0);
	objpFUNCTYPE.coeff.push_back(1);
	PlainFuncType.push_back(objpFUNCTYPE);

	betatype.push_back("");
	localstr.clear();
	localstr << sbuf;
	localstr >> sbuf;
	while (true) {
		localstr >> sbuf;
		if (sbuf == "none") {
			break;
		}
		betatype.push_back(sbuf);
		if (localstr.eof()) {
			break;
		}
	}

	readpara.close();

}

void set_MatrixY(REGELE* RegEle, REGELE* TestRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, bool* for_test, bool* reset_slope_intercept_for_scores) {

	bool between_0_and_1 = true;
	long double min_y = 1000000000000000000000000000000000000000000000000000000000000.0;
	long double max_y = -1000000000000000000000000000000000000000000000000000000000000.0;
	RegEle->MatrixY = Eigen::MatrixXd::Random(dataset->size(), 1);
	RegEle->MatrixYB = Eigen::MatrixXd::Random(dataset->size(), 1);
	if (*for_test) {
		TestRegEle->MatrixY = Eigen::MatrixXd::Random(dataset->size(), 1);
		TestRegEle->MatrixYB = Eigen::MatrixXd::Random(dataset->size(), 1);
	}
	for (int i = 0; i < (signed)dataset->size(); i++) {
		(RegEle->MatrixY)(i, 0) = (*dataset)[i].v[0];
		if ((*dataset)[i].v[0] < 0.0 || (*dataset)[i].v[0] > 1.0) {
			between_0_and_1 = false;
		}
		if ((*dataset)[i].v[0] < min_y) {
			min_y = (*dataset)[i].v[0];
		}
		if ((*dataset)[i].v[0] > max_y) {
			max_y = (*dataset)[i].v[0];
		}
	}
	if (*for_test) {
		for (int i = 0; i < (signed)testdataset->size(); i++) {
			(TestRegEle->MatrixY)(i, 0) = (*testdataset)[i].v[0];
			if ((*testdataset)[i].v[0] < 0.0 || (*testdataset)[i].v[0] > 1.0) {
				between_0_and_1 = false;
			}
			if ((*testdataset)[i].v[0] < min_y) {
				min_y = (*testdataset)[i].v[0];
			}
			if ((*testdataset)[i].v[0] > max_y) {
				max_y = (*testdataset)[i].v[0];
			}
		}
	}

	if (!between_0_and_1 && *reset_slope_intercept_for_scores) {
		long double total_dataset_size = 0.0;
		if (*for_test) {
			total_dataset_size = (long double)dataset->size() + (long double)testdataset->size();
		}
		else {
			total_dataset_size = (long double)dataset->size();
		}
		score_slope = (1.0 / total_dataset_size) * (total_dataset_size - 1.0) / (max_y - min_y);
		score_intercept = (1.0 / total_dataset_size) * (-min_y * (total_dataset_size - 1.0) / (max_y - min_y) + 0.5);
	}

	for (int i = 0; i < (signed)dataset->size(); i++) {
		(RegEle->MatrixYB)(i, 0) = score_slope * (*dataset)[i].v[0] + score_intercept;
	}
	if (*for_test) {
		for (int i = 0; i < (signed)testdataset->size(); i++) {
			(TestRegEle->MatrixYB)(i, 0) = score_slope * (*testdataset)[i].v[0] + score_intercept;
		}
	}

}

void read_data(std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* datadivisionflag, bool* use_dataname, std::string* userdefinedName, long double* TrainRatio, std::vector<std::string>* DataFilename, bool* ReOutput, bool* ReInitial, std::string* InitialFilename, int* ibeta, std::string* setY, std::vector<GROUPX>* GroupX, bool* minlimitY, bool* maxlimitY, long double* minY, long double* maxY, std::vector<std::string>* column_str, int* DuplicateX, int* NameDuplicateHandling, long double* Dmin, long double* Dmax, long double* Dincrement, std::vector<SUBFILTERX>* subfilterx, int* Verbose) {

	readdata.open((*DataFilename)[0]);
	std::string Name = "Name";
	int histogramset = 0;
	long double histogram0 = 0.0;
	long double histogram1 = 0.0;
	if (*NameDuplicateHandling == 2) {
		histogramset = (int)((*Dmax - *Dmin) / *Dincrement) + 1;
		histogram0 = *Dmin - 0.5 * *Dincrement;
		histogram1 = *Dmax - 0.5 * *Dincrement;
	}

	if (*ReOutput) {
		GroupX->clear();
		if (*ReInitial) {
			readinitial.open(*InitialFilename);
			while (true) {
				std::getline(readinitial, sbuf);
				if (sbuf.find("DATASET") != std::string::npos || readinitial.eof()) {
					break;
				}
			}
		}
		int interaction_line_count = -1;
		while (true) {
			std::getline(readdata, sbuf);
			if (sbuf.find("INTERACTION") != std::string::npos) {
				interaction_line_count++;
				if (interaction_line_count == *ibeta) {
					break;
				}
			}
			if (readdata.eof()) {
				break;
			}
		}
	}

	std::stringstream localstr;
	unsigned int initial_column_str_count = 0;
	if (*ReOutput && *ReInitial) {
		std::getline(readinitial, sbuf);
		localstr.clear();
		localstr << sbuf;
		while (true) {
			localstr >> sbuf;
			column_str->push_back(sbuf);
			if (sbuf.find(Name) != std::string::npos) {
				*use_dataname = true;
				*userdefinedName = sbuf;
			}
			else if (sbuf != *setY) {
				GROUPX objGROUPX;
				objGROUPX.setX = sbuf;
				objGROUPX.negative_exists = false;
				objGROUPX.min = 1000000000000000000000000000000000000000.0;
				objGROUPX.ave = 0.0;
				objGROUPX.max = -objGROUPX.min;
				objGROUPX.stdev = 0.0;
				for (int ix = 0; ix < *DuplicateX; ix++) {
					GroupX->push_back(objGROUPX);
				}
			}
			if (localstr.eof()) {
				break;
			}
		}
		initial_column_str_count = column_str->size();
	}
	std::getline(readdata, sbuf);
	localstr.clear();
	localstr << sbuf;
	while (true) {
		localstr >> sbuf;
		column_str->push_back(sbuf);
		if (sbuf.find(Name) != std::string::npos) {
			*use_dataname = true;
			*userdefinedName = sbuf;
		}
		else if (*ReOutput && sbuf != *setY) {
			GROUPX objGROUPX;
			objGROUPX.setX = sbuf;
			objGROUPX.negative_exists = false;
			objGROUPX.min = 1000000000000000000000000000000000000000.0;
			objGROUPX.ave = 0.0;
			objGROUPX.max = -objGROUPX.min;
			objGROUPX.stdev = 0.0;
			for (int ix = 0; ix < *DuplicateX; ix++) {
				GroupX->push_back(objGROUPX);
			}
		}
		if (localstr.eof()) {
			break;
		}
	}

	int* WhichColumnSaved;
	int* WhichColumnSavedSubFilterX;
	WhichColumnSaved = new int[column_str->size()];
	WhichColumnSavedSubFilterX = new int[column_str->size()];
	for (unsigned int i = 0; i < column_str->size(); i++) {
		WhichColumnSaved[i] = -1;
		WhichColumnSavedSubFilterX[i] = -1;
	}
	WhichColumnSaved[find_column(setY, column_str)] = 0;
	unsigned int point_i = 0;
	for (unsigned int i = 0; i < GroupX->size(); i += *DuplicateX) {
		WhichColumnSaved[find_column(&(*GroupX)[i].setX, column_str)] = point_i + 1;
		point_i++;
	}
	if (*use_dataname) {
		WhichColumnSaved[find_column(userdefinedName, column_str)] = -2;
	}
	for (unsigned int i = 0; i < subfilterx->size(); i++) {
		int index_found = find_column(&(*subfilterx)[i].setX, column_str);
		if (index_found >= 0) {
			WhichColumnSavedSubFilterX[index_found] = i;
		}
	}

	readdata.close();
	if (*ReOutput && *ReInitial) {
		readinitial.close();
	}

	std::vector<DATASET> predataset;
	int wrongdatacount = 0;
	int duplicatecount = 0;
	int outofrangecount = 0;
	std::vector<bool> is_initial;
	for (int ifile = 0; ifile < (signed)DataFilename->size(); ifile++) {
		is_initial.push_back(false);
	}
	if (*ReOutput && *ReInitial) {
		DataFilename->push_back(*InitialFilename);
		is_initial.push_back(true);
	}

	for (int ifile = 0; ifile < (signed)DataFilename->size(); ifile++) {

		readdata.open((*DataFilename)[ifile]);

		int interaction_count = -1;
		if (*ReOutput) {
			if (!is_initial[ifile]) {
				while (true) {
					std::getline(readdata, sbuf);
					if (sbuf.find("INTERACTION") != std::string::npos) {
						interaction_count++;
						if (interaction_count == *ibeta) {
							break;
						}
					}
					if (readdata.eof()) {
						break;
					}
				}
			}
			else {
				while (true) {
					std::getline(readdata, sbuf);
					if (sbuf.find("DATASET") != std::string::npos || readdata.eof()) {
						break;
					}
				}
			}
		}

		std::getline(readdata, sbuf);

		if (!(*ReOutput && is_initial[ifile])) {

			int linecount = 1;
			while (true) {
				DATASET objDATASET;
				long double ldbuf;
				linecount++;
				bool full_number = true;
				bool in_minmaxYwindow = true;
				bool subfilter = true;
				int duplicate_with_existing_datum = -1;
				for (unsigned int i = 0; i < GroupX->size() + 1; i++) {
					objDATASET.v.push_back(0.0);
				}
				std::getline(readdata, sbuf);
				if (*ReOutput && sbuf == "******************************************") {
					break;
				}
				std::stringstream linestr;
				linestr << sbuf;
				for (unsigned int i = initial_column_str_count; i < column_str->size(); i++) {
					linestr >> sbuf;
					if (WhichColumnSaved[i] == -2) {
						objDATASET.name = sbuf;
						if (*NameDuplicateHandling == 0 || *NameDuplicateHandling == 1 || *NameDuplicateHandling == 2) {
							for (int ipredata = 0; ipredata < (signed)predataset.size(); ipredata++) {
								if (predataset[ipredata].name == sbuf) {
									duplicate_with_existing_datum = ipredata;
									break;
								}
							}
						}
					}
					else if (WhichColumnSaved[i] != -1) {
						if (is_number(&sbuf)) {
							std::stringstream lcstr;
							lcstr << sbuf;
							lcstr >> ldbuf;
							if (WhichColumnSaved[i] == 0) {
								objDATASET.v[0] = ldbuf;
							}
							else {
								int start = *DuplicateX * WhichColumnSaved[i] - (*DuplicateX - 1);
								int end = *DuplicateX * WhichColumnSaved[i];
								for (int iix = start; iix <= end; iix++) {
									objDATASET.v[iix] = ldbuf;
								}
							}
							if (WhichColumnSaved[i] == 0 && ((*minlimitY && ldbuf < *minY) || (*maxlimitY && ldbuf > *maxY))) {
								in_minmaxYwindow = false;
								outofrangecount++;
								if (pr_id == writeoutput_prid && *Verbose == 1) {
									writeoutput << "... Y is out of window!\n";
									writeoutput.flush();
									writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << ldbuf << "\n";
									writeoutput.flush();
								}
								break;
							}
						}
						else {
							full_number = false;
							wrongdatacount++;
							if (pr_id == writeoutput_prid && *Verbose == 1) {
								writeoutput << "... data are given wrong; this line will be ignored!\n";
								writeoutput.flush();
								writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << sbuf << "\n";
								writeoutput.flush();
							}
							break;
						}
					}
					if (WhichColumnSaved[i] != -2 && WhichColumnSavedSubFilterX[i] != -1) {
						if (is_number(&sbuf)) {
							std::stringstream lcstr;
							lcstr << sbuf;
							lcstr >> ldbuf;
							if (((*subfilterx)[WhichColumnSavedSubFilterX[i]].minlimitX && ldbuf < (*subfilterx)[WhichColumnSavedSubFilterX[i]].minX) || ((*subfilterx)[WhichColumnSavedSubFilterX[i]].maxlimitX && ldbuf > (*subfilterx)[WhichColumnSavedSubFilterX[i]].maxX)) {
								subfilter = false;
								outofrangecount++;
								if (pr_id == writeoutput_prid && *Verbose == 1) {
									writeoutput << "... X(" << (*subfilterx)[WhichColumnSavedSubFilterX[i]].setX << ") is out of window!\n";
									writeoutput.flush();
									writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << ldbuf << "\n";
									writeoutput.flush();
								}
								break;
							}
						}
						else {
							subfilter = false;
							wrongdatacount++;
							if (pr_id == writeoutput_prid && *Verbose == 1) {
								writeoutput << "... data are given wrong; this line will be ignored!\n";
								writeoutput.flush();
								writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << sbuf << "\n";
								writeoutput.flush();
							}
							break;
						}
					}
				}
				if (full_number && in_minmaxYwindow && subfilter) {
					if (duplicate_with_existing_datum == -1) {
						if (*NameDuplicateHandling == 2) {
							for (int ihit = 0; ihit < histogramset; ihit++) {
								objDATASET.yhit.push_back(0);
							}
							int whichhistogramset = (int)(((long double)histogramset - 1.0) * (objDATASET.v[0] - histogram0) / (histogram1 - histogram0));

							if (whichhistogramset >= 0 && whichhistogramset < (signed)objDATASET.yhit.size()) {
								objDATASET.yhit[whichhistogramset]++;
							}
						}
						predataset.push_back(objDATASET);

					}
					else {
						duplicatecount++;
						if ((*NameDuplicateHandling == 0 && objDATASET.v[0] < predataset[duplicate_with_existing_datum].v[0]) || (*NameDuplicateHandling == 1 && objDATASET.v[0] > predataset[duplicate_with_existing_datum].v[0])) {
							if (pr_id == writeoutput_prid && *Verbose == 1) {
								writeoutput << "... Datum was replaced due to duplicate issue!\n";
								writeoutput.flush();
								writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " ... " << objDATASET.name << ": " << objDATASET.v[0] << " <--- " << predataset[duplicate_with_existing_datum].v[0] << "\n";
								writeoutput.flush();
							}
							predataset[duplicate_with_existing_datum] = objDATASET;
						}
						else if ((*NameDuplicateHandling == 0 && objDATASET.v[0] > predataset[duplicate_with_existing_datum].v[0]) || (*NameDuplicateHandling == 1 && objDATASET.v[0] < predataset[duplicate_with_existing_datum].v[0])) {
							if (pr_id == writeoutput_prid && *Verbose == 1) {
								writeoutput << "... Datum was not replaced in spite of duplicate issue!\n";
								writeoutput.flush();
								writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " ... " << objDATASET.name << ": " << objDATASET.v[0] << " <-X- " << predataset[duplicate_with_existing_datum].v[0] << "\n";
								writeoutput.flush();
							}
						}
						else if (*NameDuplicateHandling == 2) {
							int whichhistogramset = (int)(((long double)histogramset - 1.0) * (objDATASET.v[0] - histogram0) / (histogram1 - histogram0));
							if (whichhistogramset >= 0 && whichhistogramset < (signed)objDATASET.yhit.size()) {
								predataset[duplicate_with_existing_datum].yhit[whichhistogramset]++;
								if (pr_id == writeoutput_prid && *Verbose == 1) {
									writeoutput << "... Datum duplicate!\n";
									writeoutput.flush();
									writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " ... " << objDATASET.name << "\n";
									writeoutput.flush();
								}
							}
						}
					}
				}
				if (readdata.eof()) {
					break;
				}
			}

		}
		else {

			int datasetcount = 0;
			while (true) {
				std::getline(readdata, sbuf);
				if (sbuf == "******************************************") {
					break;
				}
				std::stringstream linestr;
				linestr << sbuf;
				for (unsigned int i = 0; i < initial_column_str_count; i++) {
					linestr >> sbuf;
					if (WhichColumnSaved[i] >= 1) {
						if (is_number(&sbuf)) {
							std::stringstream lcstr;
							lcstr << sbuf;
							long double ldbuf;
							lcstr >> ldbuf;
							if (datasetcount < (signed)predataset.size()) {
								int start = *DuplicateX * WhichColumnSaved[i] - (*DuplicateX - 1);
								int end = *DuplicateX * WhichColumnSaved[i];
								for (int iix = start; iix <= end; iix++) {
									predataset[datasetcount].v[iix] = ldbuf;
								}
							}
							else {
								if (pr_id == writeoutput_prid && *Verbose == 1) {
									writeoutput << "... Something wrong in reading initial dataset!\n";
									writeoutput.flush();
									writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << datasetcount + 2 << "\n";
									writeoutput.flush();
								}
							}
						}
						else {
							int start = *DuplicateX * WhichColumnSaved[i] - (*DuplicateX - 1);
							int end = *DuplicateX * WhichColumnSaved[i];
							for (int iix = start; iix <= end; iix++) {
								predataset[datasetcount].v[iix] = -100000000000000;
							}
							if (pr_id == writeoutput_prid && *Verbose == 1) {
								writeoutput << "... Something wrong in reading initial dataset [non-value]!\n";
								writeoutput.flush();
								writeoutput << "check line: " << (*DataFilename)[ifile] << " -- " << datasetcount + 2 << " ... " << sbuf << "\n";
								writeoutput.flush();
							}
						}
					}
				}
				datasetcount++;
			}

		}

		readdata.close();

	}

	if (*NameDuplicateHandling == 2) {
		for (int idat = 0; idat < (signed)predataset.size(); idat++) {
			int maxhit = -100;
			int targethit = -1;
			for (int ihit = 0; ihit < (signed)predataset[idat].yhit.size(); ihit++) {
				if (predataset[idat].yhit[ihit] > maxhit) {
					maxhit = predataset[idat].yhit[ihit];
					targethit = ihit;
				}
			}
			if (targethit >= 0) {
				predataset[idat].v[0] = *Dmin + (long double)targethit * *Dincrement;
			}
		}
	}

	if (pr_id == writeoutput_prid) {
		writeoutput << "TOTAL DATA NO = " << predataset.size() << "\n";
		writeoutput.flush();
		writeoutput << "WRONG DATA NO = " << wrongdatacount << "\n";
		writeoutput.flush();
		writeoutput << "DUPLICATE DATA NO = " << duplicatecount << "\n";
		writeoutput.flush();
		writeoutput << "OUT-OF-RANGE DATA NO = " << outofrangecount << "\n";
		writeoutput.flush();
	}

	int Ntrain = (int)((long double)predataset.size() * *TrainRatio);
	int Ntest = (signed)predataset.size() - Ntrain;
	/*bool defacto_notest = false;
	if (Ntest == 0) {
		Ntest = 1;
		defacto_notest = true;
	}*/
	if (pr_id == writeoutput_prid) {
		writeoutput << "Ntrain = " << Ntrain << "\n";
		writeoutput.flush();
		writeoutput << "Ntest = " << Ntest << "\n";
		writeoutput.flush();
	}

	std::vector<int>().swap(*datadivisionflag);

	std::vector<int> iTest;
	for (int i = 0; i < Ntest; i++) {
		while (true) {
			int randomNumber = rand() % (signed)predataset.size();
			bool already = false;
			for (int j = 0; j < (signed)iTest.size(); j++) {
				if (iTest[j] == randomNumber) {
					already = true;
					break;
				}
			}
			if (!already) {
				iTest.push_back(randomNumber);
				break;
			}
		}
	}

	for (int i = 0; i < (signed)predataset.size(); i++) {
		bool for_test = false;
		for (int j = 0; j < (signed)iTest.size(); j++) {
			if (iTest[j] == i) {
				for_test = true;
				break;
			}
		}
		if (for_test) {
			testdataset->push_back(predataset[i]);
			datadivisionflag->push_back(-(signed)testdataset->size());
		}
		else {
			dataset->push_back(predataset[i]);
			datadivisionflag->push_back((signed)dataset->size());
		}
	}

	if (dataset->size() > 0) {

		for (int ix = 0; ix < (signed)GroupX->size(); ix++) {
			(*GroupX)[ix].min = 1000000000000000000000000000000000000000.0;
			(*GroupX)[ix].ave = 0.0;
			(*GroupX)[ix].max = -(*GroupX)[ix].min;
			(*GroupX)[ix].stdev = 0.0;
			(*GroupX)[ix].negative_exists = false;
		}
		for (int i = 0; i < (signed)dataset->size(); i++) {
			for (int j = 1; j < (signed)(*dataset)[0].v.size(); j++) {
				if ((*dataset)[i].v[j] < 0.0) {
					(*GroupX)[j - 1].negative_exists = true;
				}
				if ((*dataset)[i].v[j] < (*GroupX)[j - 1].min) {
					(*GroupX)[j - 1].min = (*dataset)[i].v[j];
				}
				if ((*dataset)[i].v[j] > (*GroupX)[j - 1].max) {
					(*GroupX)[j - 1].max = (*dataset)[i].v[j];
				}
				(*GroupX)[j - 1].ave += (*dataset)[i].v[j];
			}
		}
		for (int j = 1; j < (signed)(*dataset)[0].v.size(); j++) {
			(*GroupX)[j - 1].ave /= (long double)(signed)dataset->size();
		}
		for (int i = 0; i < (signed)dataset->size(); i++) {
			for (int j = 1; j < (signed)(*dataset)[0].v.size(); j++) {
				(*GroupX)[j - 1].stdev += powl((*dataset)[i].v[j] - (*GroupX)[j - 1].ave, 2.0);
			}
		}
		for (int j = 1; j < (signed)(*dataset)[0].v.size(); j++) {
			(*GroupX)[j - 1].stdev /= (long double)(signed)dataset->size();
			(*GroupX)[j - 1].stdev = sqrtl((*GroupX)[j - 1].stdev);
		}

	}

	if (pr_id == writeoutput_prid) {
		writeoutput << "******************************DATASET******************************\n";
		writeoutput.flush();
		if (*use_dataname) {
			writeoutput << *userdefinedName << "\t";
			writeoutput.flush();
		}
		writeoutput << *setY << "\t";
		writeoutput.flush();
		for (unsigned int i = 0; i < GroupX->size(); i++) {
			writeoutput << (*GroupX)[i].setX;
			writeoutput.flush();
			if (i != GroupX->size() - 1) {
				writeoutput << "\t";
				writeoutput.flush();
			}
			else {
				writeoutput << "\n";
				writeoutput.flush();
			}
		}
		for (unsigned int j = 0; j < predataset.size(); j++) {
			if (*use_dataname) {
				writeoutput << predataset[j].name << "\t";
				writeoutput.flush();
			}
			for (unsigned int i = 0; i < predataset[j].v.size(); i++) {
				writeoutput << predataset[j].v[i];
				writeoutput.flush();
				if (i != predataset[j].v.size() - 1) {
					writeoutput << "\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\n";
					writeoutput.flush();
				}
			}
		}
		writeoutput << "******************************************\n";
		writeoutput.flush();
		writeoutput << "******************************************\n";
		writeoutput.flush();
		writeoutput << "******************************************\n";
		writeoutput.flush();
	}

	if (*Verbose == 1) {
		if (*use_dataname) {
			writeoutputpr << *userdefinedName << "\t";
			writeoutputpr.flush();
		}
		writeoutputpr << *setY << "\t";
		writeoutputpr.flush();
		for (unsigned int i = 0; i < GroupX->size(); i++) {
			writeoutputpr << (*GroupX)[i].setX;
			writeoutputpr.flush();
			if (i != GroupX->size() - 1) {
				writeoutputpr << "\t";
				writeoutputpr.flush();
			}
			else {
				writeoutputpr << "\n";
				writeoutputpr.flush();
			}
		}
		for (unsigned int j = 0; j < predataset.size(); j++) {
			if (*use_dataname) {
				writeoutputpr << predataset[j].name << "\t";
				writeoutputpr.flush();
			}
			for (unsigned int i = 0; i < predataset[j].v.size(); i++) {
				writeoutputpr << predataset[j].v[i];
				writeoutputpr.flush();
				if (i != predataset[j].v.size() - 1) {
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
	
	if (pr_id == writeoutput_prid) {
		writeoutput << "Name\tNegativeExists\tMin\tAve\tMax\tStdev\n";
		writeoutput.flush();
		for (int ix = 0; ix < (signed)GroupX->size(); ix++) {
			writeoutput << (*GroupX)[ix].setX << "\t";
			writeoutput.flush();
			if ((*GroupX)[ix].negative_exists) {
				writeoutput << "true\t";
				writeoutput.flush();
			}
			else {
				writeoutput << "false\t";
				writeoutput.flush();
			}
			writeoutput << (*GroupX)[ix].min << "\t" << (*GroupX)[ix].ave << "\t" << (*GroupX)[ix].max << "\t" << (*GroupX)[ix].stdev << "\n";
			writeoutput.flush();
		}
	}
	if (*Verbose == 1) {
		writeoutputpr << "Name\tNegativeExists\tMin\tAve\tMax\tStdev\n";
		writeoutputpr.flush();
		for (int ix = 0; ix < (signed)GroupX->size(); ix++) {
			writeoutputpr << (*GroupX)[ix].setX << "\t";
			writeoutputpr.flush();
			if ((*GroupX)[ix].negative_exists) {
				writeoutputpr << "true\t";
				writeoutputpr.flush();
			}
			else {
				writeoutputpr << "false\t";
				writeoutputpr.flush();
			}
			writeoutputpr << (*GroupX)[ix].min << "\t" << (*GroupX)[ix].ave << "\t" << (*GroupX)[ix].max << "\t" << (*GroupX)[ix].stdev << "\n";
			writeoutputpr.flush();
		}
	}

	delete[] WhichColumnSaved;
	WhichColumnSaved = nullptr;
	delete[] WhichColumnSavedSubFilterX;
	WhichColumnSavedSubFilterX = nullptr;

	std::vector<DATASET>().swap(predataset);
	std::vector<int>().swap(iTest);
	std::vector<bool>().swap(is_initial);

}

void transpose_one_column_to_y(REGELE* new_RegEle, std::vector<DATASET>* new_dataset, REGELE* RegEle, int* data_size, int* variable_size, unsigned int* whichcolumn) {

	new_RegEle->MatrixY = Eigen::MatrixXd::Random(*data_size, 1);
	for (unsigned int iLine = 0; iLine < *data_size; iLine++) {
		DATASET objDATASET;
		objDATASET.v.push_back(0.0);
		for (unsigned int iColumn = 0; iColumn < *variable_size + 1; iColumn++) {
			if (iColumn != *whichcolumn - 1) {
				objDATASET.v.push_back(RegEle->MatrixX(iLine, iColumn + 1));
			}
			else {
				new_RegEle->MatrixY(iLine, 0) = RegEle->MatrixX(iLine, iColumn + 1);
				objDATASET.v[0] = RegEle->MatrixX(iLine, iColumn + 1);
			}
		}
		new_dataset->push_back(objDATASET);
	}

}

void estimate_jobsize(bool* xpr_ld_get, signed long long int* jobsize, long double* jobsize_ld, int* Nfeatures, int* NVariables, int* Verbose) {

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
		std::vector<int>().swap(recipe_jobcontrol);
	}

	if (!*xpr_ld_get && pr_id == writeoutput_prid) {
		for (int taken = NVariables[0]; taken <= NVariables[1]; taken++) {
			writeoutput << "jobsize up to Nvar = " << taken << " : " << jobsize_step[taken - NVariables[0]] << "\n";
			writeoutput.flush();
		}
	}
	else if (*xpr_ld_get && pr_id == writeoutput_prid) {
		writeoutput << "xpr_ld_get=\ttrue\n";
		writeoutput.flush();
		for (int taken = NVariables[0]; taken <= NVariables[1]; taken++) {
			writeoutput << "jobsize up to Nvar = " << taken << " : " << jobsize_step_ld[taken - NVariables[0]] << "\n";
			writeoutput.flush();
		}
	}

	if (*Verbose == 1) {
		if (!*xpr_ld_get) {
			for (int taken = NVariables[0]; taken <= NVariables[1]; taken++) {
				writeoutputpr << "jobsize up to Nvar = " << taken << " : " << jobsize_step[taken - NVariables[0]] << "\n";
				writeoutputpr.flush();
			}
		}
		else {
			writeoutputpr << "xpr_ld_get=\ttrue\n";
			writeoutputpr.flush();
			for (int taken = NVariables[0]; taken <= NVariables[1]; taken++) {
				writeoutputpr << "jobsize up to Nvar = " << taken << " : " << jobsize_step_ld[taken - NVariables[0]] << "\n";
				writeoutputpr.flush();
			}
		}
	}
	
}

void set_pr_jobs(bool* xpr_ld_get, signed long long int* jobsize, long double* jobsize_ld, int* Verbose) {

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

	if (pr_id == writeoutput_prid) {
		if (!*xpr_ld_get) {
			writeoutput << "jobsize=\t" << *jobsize << "\n";
			writeoutput.flush();
			for (int j = 0; j < pr_size; j++) {
				writeoutput << "pr_id\t" << j << "\t" << pr_job[j][0] + 1 << "\t" << pr_job[j][1] + 1 << "\t" << pr_job[j][1] - pr_job[j][0] + 1;
				writeoutput.flush();
				if (j != pr_size - 1) {
					writeoutput << "\n";
					writeoutput.flush();
				}
			}
		}
		else {
			writeoutput << "jobsize=\t" << *jobsize_ld << "\n";
			writeoutput.flush();
			for (int j = 0; j < pr_size; j++) {
				writeoutput << "pr_id\t" << j << "\t" << pr_job_ld[j][0] + 1 << "\t" << pr_job_ld[j][1] + 1 << "\t" << pr_job_ld[j][1] - pr_job_ld[j][0] + 1;
				writeoutput.flush();
				if (j != pr_size - 1) {
					writeoutput << "\n";
					writeoutput.flush();
				}
			}
		}
	}

	if (*Verbose == 1) {
		if (!*xpr_ld_get) {
			writeoutputpr << "jobsize=\t" << *jobsize << "\n";
			writeoutputpr.flush();
			for (int j = 0; j < pr_size; j++) {
				writeoutputpr << "pr_id\t" << j << "\t" << pr_job[j][0] + 1 << "\t" << pr_job[j][1] + 1 << "\t" << pr_job[j][1] - pr_job[j][0] + 1;
				writeoutputpr.flush();
				if (j != pr_size - 1) {
					writeoutputpr << "\n";
					writeoutputpr.flush();
				}
			}
		}
		else {
			writeoutputpr << "jobsize=\t" << *jobsize_ld << "\n";
			writeoutputpr.flush();
			for (int j = 0; j < pr_size; j++) {
				writeoutputpr << "pr_id\t" << j << "\t" << pr_job_ld[j][0] + 1 << "\t" << pr_job_ld[j][1] + 1 << "\t" << pr_job_ld[j][1] - pr_job_ld[j][0] + 1;
				writeoutputpr.flush();
				if (j != pr_size - 1) {
					writeoutputpr << "\n";
					writeoutputpr.flush();
				}
			}
		}
	}
	
	SigmaMinusXk = new signed long long int[2];
	SigmaMinusXk_ld = new long double[2];

	auto current_t = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_secondsfr = current_t - start_t;
	if (pr_id == writeoutput_prid) {
		writeoutput << "\nCPU_time_so_far=\t" << elapsed_secondsfr.count();
		writeoutput.flush();
	}
	if (*Verbose == 1) {
		writeoutputpr << "\nCPU_time_so_far=\t" << elapsed_secondsfr.count();
		writeoutputpr.flush();
	}
	else {
		writeoutputpr << "CPU_time_so_far=\t" << elapsed_secondsfr.count();
		writeoutputpr.flush();
	}

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

void write_arrange(std::ofstream* to_stream, std::vector<int>* arrange) {

	for (int iarr = 0; iarr < (signed)arrange->size(); iarr++) {
		*to_stream << (*arrange)[iarr];
		to_stream->flush();
		if (iarr != (signed)arrange->size() - 1) {
			*to_stream << "-";
			to_stream->flush();
		}
	}

}

void initialize_regression(REGELE* RegEle, std::vector<DATASET>* dataset, int* Noccup) {

	RegEle->MatrixX = Eigen::MatrixXd::Random(dataset->size(), *Noccup + 1);
	RegEle->MatrixYexpected = Eigen::MatrixXd::Random(dataset->size(), 1);
	RegEle->MatrixXTXi = Eigen::MatrixXd::Random(*Noccup + 1, *Noccup + 1);
	RegEle->MatrixBeta = Eigen::MatrixXd::Random(*Noccup + 1, 1);
	RegEle->MatrixBetaStandardized = Eigen::MatrixXd::Random(*Noccup + 1, 1);
	RegEle->ShuffledMatrixBetaStandardized = Eigen::MatrixXd::Random(*Noccup + 1, 1);

	std::vector<long double>().swap(RegEle->fit_coeff);
	std::vector<long double>().swap(RegEle->fit_coeff_standardized);
	std::vector<long double>().swap(RegEle->fit_coeff_SE);
	std::vector<long double>().swap(RegEle->fit_coeff_t);
	std::vector<long double>().swap(RegEle->fit_coeff_p);
	std::vector<long double>().swap(RegEle->fit_coeff_CI_start);
	std::vector<long double>().swap(RegEle->fit_coeff_CI_end);
	std::vector<int>().swap(RegEle->fit_coeff_porder);
	std::vector<int>().swap(RegEle->fit_coeff_porder_on_arr);
	std::vector<long double>().swap(RegEle->fit_coeff_pordered);

}

void set_data(std::vector<FUNCTYPE>* FuncType, REGELE* RegEle, std::vector<DATASET>* dataset, std::vector<int>* arrange, int* Noccup) {

	int* WhichColumnSaved;
	int* FuncList;
	WhichColumnSaved = new int[*Noccup];
	FuncList = new int[*Noccup];
	for (int i = 0; i < *Noccup; i++) {
		WhichColumnSaved[i] = -1;
		FuncList[i] = 0;
	}
	int ColumnCount = 0;
	for (int i = 0; i < arrange->size(); i++) {
		if ((*arrange)[i] >= 1) {
			WhichColumnSaved[ColumnCount] = i;
			FuncList[ColumnCount] = (*arrange)[i];
			ColumnCount++;
		}
	}

	for (int i = 0; i < (signed)dataset->size(); i++) {
		RegEle->MatrixX(i, 0) = 1.0;
		for (int j = 0; j < *Noccup; j++) {
			RegEle->MatrixX(i, j + 1) = get_v(&(*FuncType)[FuncList[j] - 1], &(*dataset)[i].v[WhichColumnSaved[j] + 1]);
		}
	}

	delete[] WhichColumnSaved;
	delete[] FuncList;
	WhichColumnSaved = nullptr;
	FuncList = nullptr;

}

void cholesky_update_one_column(Eigen::LLT<Eigen::MatrixXd>* llt, Eigen::MatrixXd* MatrixX, Eigen::VectorXd* NewVectorX, int* UpdatedColumnIndex) {

	Eigen::VectorXd OldVectorX = (*MatrixX).col(*UpdatedColumnIndex);
	Eigen::VectorXd d = *NewVectorX - OldVectorX;
	if (d.isZero(0)) {
		return;
	}
	Eigen::VectorXd u = (*MatrixX).transpose() * d;
	long double s = d.squaredNorm();
	int p = (*MatrixX).cols();
	Eigen::MatrixXd U(p, 2);
	U.setZero();
	U(*UpdatedColumnIndex, 0) = 1.0;
	U.col(1) = u;
	Eigen::Matrix2d C;
	C << s, 1.0,
		1.0, 0.0;
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(C);
	Eigen::Vector2d lambda = es.eigenvalues();
	Eigen::Matrix2d Q = es.eigenvectors();
	for (int k = 0; k < 2; ++k) {
		double alpha = lambda(k);
		if (std::abs(alpha) == 0.0) continue;
		Eigen::VectorXd v = U * Q.col(k);
		(*llt).rankUpdate(v, alpha);
	}
	(*MatrixX).col(*UpdatedColumnIndex) = *NewVectorX;
	
}

void cholesky_update_multiple_columns(Eigen::LLT<Eigen::MatrixXd>* llt, Eigen::MatrixXd* MatrixX, Eigen::MatrixXd* NewMatrixX) {
	
	std::vector<int> changed;
	changed.clear();
	for (int i = 0; i < (*MatrixX).cols(); i++) {
		if ((*MatrixX).col(i) != (*NewMatrixX).col(i)) {
			changed.push_back(i);
			///////FOR DEBUG///////
			/*
			if (pr_id == writeoutput_prid) {
				std::cout << i << " ";
				std::cout.flush();
			}
			*/
			//////////////////////
		}
	}
	///////FOR DEBUG///////
	/*
	if (pr_id == writeoutput_prid) {
		std::cout << " :\n";
		std::cout.flush();
	}
	*/
	//////////////////////

	if (!changed.empty()) {
		if ((!useSuperFastSwap && changed.size() > 2) || (useSuperFastSwap && changed.size() >= (*MatrixX).cols() - 1)) {
			*MatrixX = *NewMatrixX;
			Eigen::MatrixXd MatrixXTX = MatrixX->transpose() * (*MatrixX);
			llt->compute(MatrixXTX);
			return; 
		}
		///////FOR DEBUG///////
		/*
		if (pr_id == writeoutput_prid) {
			std::cout << "*******************:\n";
			std::cout.flush();
		}
		*/
		//////////////////////
		for (int idx = 0; idx < (int)changed.size(); idx++) {
			int j = changed[idx];
			Eigen::VectorXd NewVectorX = NewMatrixX->col(j);
			cholesky_update_one_column(llt, MatrixX, &NewVectorX, &j);
		}
	}

	std::vector<int>().swap(changed);

}

void fit_any(REGELE* RegEle, std::vector<DATASET>* dataset, int* Noccup, std::vector<long double>* x_stdev_series) {

	Eigen::MatrixXd MatrixStandardizedY = RegEle->MatrixY; 
	Eigen::MatrixXd MatrixStandardizedX = RegEle->MatrixX;
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
	for (int i = 0; i < (signed)dataset->size(); i++) {
		MatrixStandardizedY(i, 0) = (RegEle->MatrixY(i, 0) - y_average) / y_stdev;
	}
	std::vector<long double> x_average_series;
	x_stdev_series->clear();
	x_average_series.push_back(1.0);
	x_stdev_series->push_back(0.0);
	for (int i = 0; i < (signed)dataset->size(); i++) {
		MatrixStandardizedX(i, 0) = 1.0;
	}
	for (int j = 1; j < *Noccup + 1; j++) {
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
		x_average_series.push_back(x_average);
		x_stdev_series->push_back(x_stdev);
		for (int i = 0; i < (signed)dataset->size(); i++) {
			MatrixStandardizedX(i, j) = (RegEle->MatrixX(i, j) - x_average) / x_stdev;
		}
	}

	Eigen::MatrixXd MatrixStandardizedXTX, ShuffledMatrixStandardizedX;
	Eigen::VectorXd VectorRHS;
	Eigen::LLT<Eigen::MatrixXd> llt;
	std::vector<int> WhichCol_InOriginal, WhichCol_InNew;
	WhichCol_InOriginal.clear();
	WhichCol_InNew.clear();
	ShuffledMatrixStandardizedX = MatrixStandardizedX;
	if (ExistingCholskeySeries && ExistingShuffledStandardizedMatrixX.rows() == MatrixStandardizedX.rows() && ExistingShuffledStandardizedMatrixX.cols() == MatrixStandardizedX.cols()) {
		for (int i = 0; i < ExistingShuffledStandardizedMatrixX.cols(); i++) {
			WhichCol_InOriginal.push_back(-1);
			WhichCol_InNew.push_back(-1);
		}
		for (int i = 0; i < ExistingShuffledStandardizedMatrixX.cols(); i++) {
			for (int j = 0; j < MatrixStandardizedX.cols(); j++) {
				if (ExistingShuffledStandardizedMatrixX.col(i) == MatrixStandardizedX.col(j) && WhichCol_InOriginal[i] == -1 && WhichCol_InNew[j] == -1) {
					ShuffledMatrixStandardizedX.col(i) = MatrixStandardizedX.col(j);
					WhichCol_InOriginal[i] = j;
					WhichCol_InNew[j] = i;
					break;
				}
			}
		}
		for (int i = 0; i < ExistingShuffledStandardizedMatrixX.cols(); i++) {
			if (WhichCol_InOriginal[i] == -1) {
				for (int j = 0; j < MatrixStandardizedX.cols(); j++) {
					if (WhichCol_InNew[j] == -1) {
						ShuffledMatrixStandardizedX.col(i) = MatrixStandardizedX.col(j);
						WhichCol_InOriginal[i] = j;
						WhichCol_InNew[j] = i;
						break;
					}
				}
			}
		}
		bool found_all = true;
		for (int i = 0; i < ExistingShuffledStandardizedMatrixX.cols(); i++) {
			if (WhichCol_InOriginal[i] == -1 || WhichCol_InNew[i] == -1) {
				found_all = false;
				break;
			}
		}
		if (!found_all) {
			ShuffledMatrixStandardizedX = MatrixStandardizedX;
			for (int i = 0; i < ExistingShuffledStandardizedMatrixX.cols(); i++) {
				WhichCol_InOriginal[i] = i;
				WhichCol_InNew[i] = i;
			}
		}
		llt = ExistingLLT;
		cholesky_update_multiple_columns(&llt, &ExistingShuffledStandardizedMatrixX, &ShuffledMatrixStandardizedX);
	}
	else {
		for (int i = 0; i < *Noccup + 1; i++) {
			WhichCol_InOriginal.push_back(i);
			WhichCol_InNew.push_back(i);
		}
		MatrixStandardizedXTX = ShuffledMatrixStandardizedX.transpose() * ShuffledMatrixStandardizedX;
		llt.compute(MatrixStandardizedXTX);
	}
	VectorRHS = ShuffledMatrixStandardizedX.transpose() * MatrixStandardizedY;
	RegEle->ShuffledMatrixBetaStandardized = llt.solve(VectorRHS);
	RegEle->ShuffledMatrixBetaStandardized(0, 0) = 0.0;
	for (int i = 0; i < *Noccup + 1; i++) {
		RegEle->MatrixBetaStandardized(i, 0) = RegEle->ShuffledMatrixBetaStandardized(WhichCol_InNew[i], 0);
	}

	for (int j = 1; j < *Noccup + 1; j++) {
		RegEle->MatrixBeta(j, 0) = RegEle->MatrixBetaStandardized(j, 0) * y_stdev / (*x_stdev_series)[j];
	}
	RegEle->MatrixBeta(0, 0) = y_average;
	for (int j = 1; j < *Noccup + 1; j++) {
		RegEle->MatrixBeta(0, 0) -= RegEle->MatrixBeta(j, 0) * x_average_series[j];
	}
	Eigen::MatrixXd L = llt.matrixL();
	for (int j = 0; j < *Noccup + 1; j++) {
		Eigen::VectorXd e = Eigen::VectorXd::Zero(*Noccup + 1);
		e(j) = 1.0;
		Eigen::VectorXd z = L.triangularView<Eigen::Lower>().solve(e);
		RegEle->MatrixXTXi(j, j) = z.squaredNorm();
	}
	///////FOR DEBUG///////
	/*
	Eigen::MatrixXd MatrixBetaBuf = RegEle->MatrixBeta;
	MatrixStandardizedXTX = ShuffledMatrixStandardizedX.transpose() * ShuffledMatrixStandardizedX;
	Eigen::LLT<Eigen::MatrixXd> llt2;
	llt2.compute(MatrixStandardizedXTX);
	VectorRHS = ShuffledMatrixStandardizedX.transpose() * MatrixStandardizedY;
	Eigen::MatrixXd ShuffledMatrixBetaStandardized2 = llt2.solve(VectorRHS);
	Eigen::MatrixXd MatrixBetaStandardized2 = ShuffledMatrixBetaStandardized2;
	ShuffledMatrixBetaStandardized2(0, 0) = 0.0;
	Eigen::MatrixXd MatrixBeta2 = MatrixBetaBuf;
	for (int i = 0; i < ExistingShuffledStandardizedMatrixX.cols(); i++) {
		MatrixBetaStandardized2(i, 0) = ShuffledMatrixBetaStandardized2(WhichCol_InNew[i], 0);
	}
	for (int j = 1; j < *Noccup + 1; j++) {
		MatrixBeta2(j, 0) = MatrixBetaStandardized2(j, 0) * y_stdev / (*x_stdev_series)[j];
	}
	MatrixBeta2(0, 0) = y_average;
	for (int j = 1; j < *Noccup + 1; j++) {
		MatrixBeta2(0, 0) -= MatrixBeta2(j, 0) * x_average_series[j];
	}
	Eigen::VectorXd diff = MatrixBetaBuf - MatrixBeta2;
	if (pr_id == writeoutput_prid) {
		std::cout << "A " << MatrixBetaBuf.transpose() << "\n";
		std::cout.flush();
		std::cout << "B " << MatrixBeta2.transpose() << "\n";
		std::cout.flush();
		std::cout << "C " << diff.transpose() << "\n";
		std::cout.flush();
	}
	*/
	//////////////////////
	
	if (ExistingCholskeySeries && ExistingShuffledStandardizedMatrixX.rows() == MatrixStandardizedX.rows() && ExistingShuffledStandardizedMatrixX.cols() == MatrixStandardizedX.cols() && ((!useSuperFastSwap && CholskeyMonitor >= 1000) || (useSuperFastSwap && CholskeyMonitor >= 10))) {
		Eigen::MatrixXd MatrixBetaBuf = RegEle->MatrixBeta;
		MatrixStandardizedXTX = ShuffledMatrixStandardizedX.transpose() * ShuffledMatrixStandardizedX;
		llt.compute(MatrixStandardizedXTX);
		VectorRHS = ShuffledMatrixStandardizedX.transpose() * MatrixStandardizedY;
		RegEle->ShuffledMatrixBetaStandardized = llt.solve(VectorRHS);
		RegEle->ShuffledMatrixBetaStandardized(0, 0) = 0.0;
		for (int i = 0; i < ExistingShuffledStandardizedMatrixX.cols(); i++) {
			RegEle->MatrixBetaStandardized(i, 0) = RegEle->ShuffledMatrixBetaStandardized(WhichCol_InNew[i], 0);
		}
		for (int j = 1; j < *Noccup + 1; j++) {
			RegEle->MatrixBeta(j, 0) = RegEle->MatrixBetaStandardized(j, 0) * y_stdev / (*x_stdev_series)[j];
		}
		RegEle->MatrixBeta(0, 0) = y_average;
		for (int j = 1; j < *Noccup + 1; j++) {
			RegEle->MatrixBeta(0, 0) -= RegEle->MatrixBeta(j, 0) * x_average_series[j];
		}
		Eigen::VectorXd diff = MatrixBetaBuf - RegEle->MatrixBeta; 
		long double sum_abs_diff = diff.cwiseAbs().sum(); 
		if (sum_abs_diff > 0.000001) {
			ExistingCholskeySeries = false;
			if (!useSuperFastSwap) {
				std::cout << "[LOG] pr_id = " << pr_id << " : Reset Cholskey part\n";
				std::cout.flush();
			}
		}
		CholskeyMonitor = 0;
	}
		
	for (int j = 0; j < *Noccup + 1; j++) {
		RegEle->fit_coeff.push_back(RegEle->MatrixBeta(j, 0));
		RegEle->fit_coeff_standardized.push_back(RegEle->MatrixBetaStandardized(j, 0));
	}
	ExistingCholskeySeries = true;
	ExistingShuffledStandardizedMatrixX = ShuffledMatrixStandardizedX;
	ExistingLLT = llt;
	Existing_x_average_series = x_average_series;
	Existing_x_stdev_series = *x_stdev_series;
	std::vector<long double>().swap(x_average_series);
	std::vector<int>().swap(WhichCol_InOriginal);
	std::vector<int>().swap(WhichCol_InNew);
	CholskeyMonitor++;

}

void measure_SStot_SSres(REGELE* RegEle, std::vector<DATASET>* dataset, REGELE* trainedRegEle, bool* for_test) {

	if (!*for_test) {
		RegEle->MatrixYexpected = RegEle->MatrixX * RegEle->MatrixBeta;
	}
	else {
		RegEle->MatrixYexpected = RegEle->MatrixX * trainedRegEle->MatrixBeta;
	}

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

void test_robustness(long double* SStot, long double* SSres, long double* MSE, long double* RMSE, long double* MAE, long double* ybar, long double* R2, std::vector<long double>* y_regressed, std::vector<long double>* y_raw, std::vector<long double>* y_testregressed, std::vector<long double>* y_testraw, std::vector<int>* element_list) {

	*SStot = 0.0;
	*SSres = 0.0;
	*MSE = 0.0;
	*RMSE = 0.0;
	*MAE = 0.0;
	*ybar = 0.0;
	*R2 = 0.0;

	int icount = 0;
	for (int i = 0; i < (signed)element_list->size(); i++) {
		bool for_test = false;
		int target_index = (*element_list)[i];
		if (target_index >= (signed)y_regressed->size()) {
			for_test = true;
			target_index -= (signed)y_regressed->size();
		}
		if (for_test && target_index >= (signed)y_testregressed->size()) {
			break;
		}
		if (!for_test) {
			*ybar += (*y_raw)[target_index];
			icount++;
		}
		else {
			*ybar += (*y_testraw)[target_index];
			icount++;
		}
	}

	if (icount > 0) {

		*ybar /= (long double)icount;

		for (int i = 0; i < (signed)element_list->size(); i++) {
			bool for_test = false;
			int target_index = (*element_list)[i];
			if (target_index >= (signed)y_regressed->size()) {
				for_test = true;
				target_index -= (signed)y_regressed->size();
			}
			if (for_test && target_index >= (signed)y_testregressed->size()) {
				break;
			}
			if (!for_test) {
				*SStot += powl((*y_raw)[target_index] - *ybar, 2.0);
				*SSres += powl((*y_raw)[target_index] - (*y_regressed)[target_index], 2.0);
				*MAE += fabsl((*y_raw)[target_index] - (*y_regressed)[target_index]);
			}
			else {
				*SStot += powl((*y_testraw)[target_index] - *ybar, 2.0);
				*SSres += powl((*y_testraw)[target_index] - (*y_testregressed)[target_index], 2.0);
				*MAE += fabsl((*y_testraw)[target_index] - (*y_testregressed)[target_index]);
			}
		}

		*MSE = *SSres / (long double)icount;
		*RMSE = sqrtl(*MSE);
		*MAE /= (long double)icount;
		*R2 = 1.0 - *SSres / *SStot;

	}

}

void measure_t(REGELE* RegEle, int* Noccup, std::vector<long double>* x_stdev_series) {

	long double sigma2 = RegEle->SSres / (long double)RegEle->df;
	Eigen::MatrixXd XTX = RegEle->MatrixX.transpose() * RegEle->MatrixX;
	Eigen::LLT<Eigen::MatrixXd> llt_raw(XTX);
	Eigen::MatrixXd L_raw = llt_raw.matrixL();
	Eigen::VectorXd e0 = Eigen::VectorXd::Zero(*Noccup + 1);
	e0(0) = 1.0;
	Eigen::VectorXd z0 = L_raw.triangularView<Eigen::Lower>().solve(e0);
	long double inv00 = z0.squaredNorm();
	long double SE0 = sqrtl(sigma2 * inv00);
	RegEle->fit_coeff_SE.push_back(SE0);
	RegEle->fit_coeff_t.push_back(RegEle->MatrixBeta(0, 0) / SE0);

	Eigen::MatrixXd L = ExistingLLT.matrixL();
	for (int i = 1; i < *Noccup + 1; i++) {
		Eigen::VectorXd e = Eigen::VectorXd::Zero(*Noccup + 1);
		e(i) = 1.0;
		Eigen::VectorXd z = L.triangularView<Eigen::Lower>().solve(e);
		long double inv_std_jj = z.squaredNorm();
		long double sx = (*x_stdev_series)[i];
		long double SE = sqrtl(sigma2 * inv_std_jj / powl(sx, 2.0));
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

void order_p(REGELE* RegEle, std::vector<int>* arrange) {

	std::vector<int>().swap(RegEle->fit_coeff_porder);
	std::vector<int>().swap(RegEle->fit_coeff_porder_on_arr);
	std::vector<long double>().swap(RegEle->fit_coeff_pordered);

	bool* already_chk;
	already_chk = new bool[RegEle->fit_coeff_p.size() - 1];
	for (int i = 0; i < (signed)RegEle->fit_coeff_p.size() - 1; i++) {
		already_chk[i] = false;
	}

	for (int ii = 1; ii < RegEle->fit_coeff_p.size(); ii++) {
		long double find_max = -1000000.0;
		int target_index = -1;
		for (int i = 1; i < (signed)RegEle->fit_coeff_p.size(); i++) {
			if (already_chk[i - 1] == false && RegEle->fit_coeff_p[i] > find_max) {
				find_max = RegEle->fit_coeff_p[i];
				target_index = i - 1;
			}
		}
		if (target_index == -1) {
			for (int i = 1; i < (signed)RegEle->fit_coeff_p.size(); i++) {
				if (already_chk[i - 1] == false) {
					target_index = i - 1;
				}
			}
		}
		RegEle->fit_coeff_porder.push_back(target_index);
		RegEle->fit_coeff_pordered.push_back(RegEle->fit_coeff_p[target_index + 1]);
		already_chk[target_index] = true;
		int find_pos = -1;
		for (int find_in_arr = 0; find_in_arr < (signed)arrange->size(); find_in_arr++) {
			if ((*arrange)[find_in_arr] != 0) {
				find_pos++;
				if (find_pos == target_index) {
					RegEle->fit_coeff_porder_on_arr.push_back(find_in_arr);
					break;
				}
			}
		}

	}

	delete[] already_chk;
	already_chk = nullptr;

}

void initialize_B(REGELE* RegEle, std::vector<DATASET>* dataset, int* Noccup) {

	RegEle->MatrixX = Eigen::MatrixXd::Random(dataset->size(), *Noccup + 1);
	RegEle->MatrixXTXi = Eigen::MatrixXd::Random(*Noccup + 1, *Noccup + 1);
	RegEle->MatrixYexpectedB = Eigen::MatrixXd::Random(dataset->size(), 1);
	RegEle->MatrixZB = Eigen::MatrixXd::Random(dataset->size(), 1);
	RegEle->MatrixInitBetaB = Eigen::MatrixXd::Random(*Noccup + 1, 1);
	RegEle->MatrixBetaB = Eigen::MatrixXd::Random(*Noccup + 1, 1);
	RegEle->MatrixWB = Eigen::MatrixXd::Random(dataset->size(), dataset->size());
	RegEle->MatrixDB = Eigen::MatrixXd::Random(dataset->size(), dataset->size());
	RegEle->MatrixcB = Eigen::MatrixXd::Random(dataset->size(), 1);
	RegEle->MatrixTB = Eigen::MatrixXd::Random(dataset->size(), dataset->size());
	RegEle->MatrixFisherKubetabetaB = Eigen::MatrixXd::Random(*Noccup + 1, dataset->size());
	RegEle->MatrixFisherKuphiphiB = Eigen::MatrixXd::Random(1, 1);

	std::vector<long double>().swap(RegEle->fit_coeffB);
	std::vector<long double>().swap(RegEle->fit_coeff_standardizedB);
	std::vector<long double>().swap(RegEle->fit_coeff_SEB);
	std::vector<long double>().swap(RegEle->fit_coeff_tB);
	std::vector<long double>().swap(RegEle->fit_coeff_pB);
	std::vector<long double>().swap(RegEle->fit_coeff_CI_startB);
	std::vector<long double>().swap(RegEle->fit_coeff_CI_endB);

}

void initialize_parametersB(REGELE* RegEle, std::string* functype_gB, std::string* functype_giB, std::vector<DATASET>* dataset, int* Noccup) {

	for (int i = 0; i < (signed)dataset->size(); i++) {
		long double iv = RegEle->MatrixYB(i, 0);
		RegEle->MatrixZB(i, 0) = get_v_special(functype_gB, &iv);
	}
	RegEle->MatrixXTXi = (RegEle->MatrixX.transpose() * RegEle->MatrixX).inverse();
	RegEle->MatrixInitBetaB = RegEle->MatrixXTXi * RegEle->MatrixX.transpose() * RegEle->MatrixZB;
	RegEle->initphiB = 0.0;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		Eigen::MatrixXd MatrixxB = Eigen::MatrixXd::Random(1, *Noccup + 1);
		for (int j = 0; j < *Noccup + 1; j++) {
			MatrixxB(0, j) = RegEle->MatrixX(i, j);
		}
		long double v_here = (MatrixxB * RegEle->MatrixInitBetaB)(0, 0);
		long double myuB = get_v_special(functype_giB, &v_here);
		Eigen::MatrixXd MatrixeB = Eigen::MatrixXd::Random(dataset->size(), 1);
		MatrixeB = RegEle->MatrixZB - RegEle->MatrixX * RegEle->MatrixInitBetaB;
		long double sigma2B = (MatrixeB.transpose() * MatrixeB)(0, 0) / (((long double)dataset->size() - ((long double)*Noccup + 1.0)) * powl(get_derv_special(functype_gB, &myuB), 2.0));
		RegEle->initphiB += myuB * (1.0 - myuB) / sigma2B;
	}
	RegEle->initphiB /= (long double)dataset->size();
	RegEle->initphiB -= 1.0;

}

void expect_YB(REGELE* RegEle, std::vector<DATASET>* dataset, std::string* functype_giB, int* Noccup, REGELE* trainedRegEle, bool* for_test) {

	for (int i = 0; i < (signed)dataset->size(); i++) {
		long double v = 0;
		for (int j = 0; j < *Noccup + 1; j++) {
			if (!*for_test) {
				v += RegEle->MatrixX(i, j) * RegEle->MatrixBetaB(j, 0);
			}
			else {
				v += RegEle->MatrixX(i, j) * trainedRegEle->MatrixBetaB(j, 0);
			}
		}
		RegEle->MatrixYexpectedB(i, 0) = get_v_special(functype_giB, &v);
	}
}

void estimate_logbetadenseB(REGELE* RegEle, std::vector<DATASET>* dataset) {

	if (RegEle->fit_phiB > 0.0) {
		RegEle->logbetadenseB = 0.0;
		for (int i = 0; i < (signed)dataset->size(); i++) {
			long double x = RegEle->MatrixYB(i, 0);
			long double a = RegEle->MatrixYexpectedB(i, 0) * RegEle->fit_phiB;
			long double b = (1.0 - RegEle->MatrixYexpectedB(i, 0)) * RegEle->fit_phiB;
			RegEle->logbetadenseB += logl(get_beta_PDF(&x, &a, &b));
		}
	}
	else {
		RegEle->logbetadenseB = -100.0;
	}

}

void maximize_logbetadenseB(REGELE* RegEle, std::vector<DATASET>* dataset, std::string* functype_giB, int* Noccup, bool* for_test) {

	long double epsilon = 1.0e-10;
	long double epsilon_max = 1.0e-5;
	long double crit_max = 1.0e-2;

	long double test_mesh = epsilon;
	long double der_mesh = epsilon;
	int max_monitor = 10;
	int count = 0;

	while (true) {

		Eigen::MatrixXd MatrixBetaB_current = RegEle->MatrixBetaB;
		long double fit_phiB_current = RegEle->fit_phiB;
		expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
		estimate_logbetadenseB(RegEle, dataset);
		long double logbetadenseB_current = RegEle->logbetadenseB;
		long double logbetadenseB_start = RegEle->logbetadenseB;
		std::vector<long double> derB;
		long double derindB;
		long double metadelta_mesh = test_mesh;
		long double best_mesh = 0.0;
		long double previous_logbetadenseB = RegEle->logbetadenseB;
		int delta_monitoring = 0;

		for (int ii = 0; ii < 2; ii++) {

			der_mesh = epsilon;
			long double previous_derindB = 0.0;
			bool init_trial = true;
			while (true) {
				RegEle->fit_phiB = fit_phiB_current;
				RegEle->fit_phiB += der_mesh;
				expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
				estimate_logbetadenseB(RegEle, dataset);
				derindB = (RegEle->logbetadenseB - logbetadenseB_current) / der_mesh;
				if (init_trial) {
					previous_derindB = derindB;
				}
				if (!init_trial && (fabsl((previous_derindB - derindB) / derindB) < crit_max || der_mesh > epsilon_max)) {
					break;
				}
				else {
					init_trial = false;
					previous_derindB = derindB;
					der_mesh *= 10.0;
				}
			}
			test_mesh = der_mesh;

			metadelta_mesh = test_mesh;
			best_mesh = 0.0;
			previous_logbetadenseB = RegEle->logbetadenseB;
			delta_monitoring = 0;
			for (long double trial_mesh = test_mesh; true; trial_mesh += metadelta_mesh) {
				RegEle->MatrixBetaB = MatrixBetaB_current;
				RegEle->fit_phiB = fit_phiB_current;
				RegEle->fit_phiB += trial_mesh * derindB;
				expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
				estimate_logbetadenseB(RegEle, dataset);
				if (RegEle->logbetadenseB > previous_logbetadenseB) {
					previous_logbetadenseB = RegEle->logbetadenseB;
					best_mesh = trial_mesh;
					delta_monitoring++;
					if (delta_monitoring == 2) {
						metadelta_mesh *= 2.0;
						delta_monitoring = 0;
					}
				}
				else {
					break;
				}
			}
			RegEle->MatrixBetaB = MatrixBetaB_current;
			RegEle->fit_phiB = fit_phiB_current;
			RegEle->fit_phiB += best_mesh * derindB;
			expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
			estimate_logbetadenseB(RegEle, dataset);
			MatrixBetaB_current = RegEle->MatrixBetaB;
			fit_phiB_current = RegEle->fit_phiB;
			logbetadenseB_current = RegEle->logbetadenseB;
			previous_logbetadenseB = logbetadenseB_current;

		}

		for (int ii = 0; ii < 2; ii++) {

			long double previous_derindB = 0.0;
			bool init_trial = true;
			for (int i = 0; i < *Noccup + 1; i++) {
				der_mesh = epsilon;
				while (true) {
					RegEle->MatrixBetaB = MatrixBetaB_current;
					RegEle->MatrixBetaB(i, 0) += der_mesh;
					expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
					estimate_logbetadenseB(RegEle, dataset);
					long double vderB = (RegEle->logbetadenseB - logbetadenseB_current) / der_mesh;
					if (init_trial) {
						previous_derindB = vderB;
					}
					if (!init_trial && (fabsl((previous_derindB - vderB) / vderB) < crit_max || der_mesh > epsilon_max)) {
						derB.push_back(vderB);
						break;
					}
					else {
						init_trial = false;
						previous_derindB = vderB;
						der_mesh *= 10.0;
					}
				}

			}

			der_mesh = epsilon;
			previous_derindB = 0.0;
			init_trial = true;
			while (true) {
				RegEle->MatrixBetaB = MatrixBetaB_current;
				RegEle->fit_phiB = fit_phiB_current + der_mesh;
				expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
				estimate_logbetadenseB(RegEle, dataset);
				long double vderB = (RegEle->logbetadenseB - logbetadenseB_current) / der_mesh;
				if (init_trial) {
					previous_derindB = vderB;
				}
				if (!init_trial && (fabsl((previous_derindB - vderB) / vderB) < crit_max || der_mesh > epsilon_max)) {
					derB.push_back(vderB);
					break;
				}
				else {
					init_trial = false;
					previous_derindB = vderB;
					der_mesh *= 10.0;
				}
			}
			test_mesh = der_mesh;

			previous_logbetadenseB = RegEle->logbetadenseB;
			best_mesh = 0.0;
			metadelta_mesh = test_mesh;
			delta_monitoring = 0;
			for (long double trial_mesh = test_mesh; trial_mesh < true; trial_mesh += metadelta_mesh) {
				RegEle->MatrixBetaB = MatrixBetaB_current;
				RegEle->fit_phiB = fit_phiB_current;
				for (int i = 0; i < *Noccup + 1; i++) {
					RegEle->MatrixBetaB(i, 0) += trial_mesh * derB[i];
				}
				RegEle->fit_phiB += trial_mesh * derB[*Noccup + 1];
				expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
				estimate_logbetadenseB(RegEle, dataset);
				if (RegEle->logbetadenseB > previous_logbetadenseB) {
					previous_logbetadenseB = RegEle->logbetadenseB;
					best_mesh = trial_mesh;
					delta_monitoring++;
					if (delta_monitoring == 2) {
						metadelta_mesh *= 2.0;
						delta_monitoring = 0;
					}
				}
				else {
					break;
				}
			}

			RegEle->MatrixBetaB = MatrixBetaB_current;
			RegEle->fit_phiB = fit_phiB_current;
			for (int i = 0; i < *Noccup + 1; i++) {
				RegEle->MatrixBetaB(i, 0) += best_mesh * derB[i];
			}
			RegEle->fit_phiB += best_mesh * derB[*Noccup + 1];
			expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
			estimate_logbetadenseB(RegEle, dataset);

			MatrixBetaB_current = RegEle->MatrixBetaB;
			fit_phiB_current = RegEle->fit_phiB;
			logbetadenseB_current = RegEle->logbetadenseB;
			previous_logbetadenseB = logbetadenseB_current;

		}

		count++;
		if (RegEle->logbetadenseB < logbetadenseB_start || fabsl((RegEle->logbetadenseB - logbetadenseB_start) / logbetadenseB_start) < epsilon_max * fabsl(logbetadenseB_start) || best_mesh == 0.0 || count > max_monitor) {
			break;
		}

		std::vector<long double>().swap(derB);

	}

}

void maximize_logbetadense_BFGSB(REGELE* RegEle, std::vector<DATASET>* dataset, std::string* functype_giB, int* Noccup, int* max_monitor, bool* for_test) {

	long double epsilon, epsilon_max;

	/*inacuurate, too fast*/
	epsilon = 1.0e-8;
	epsilon_max = 1.0e-1;
	/*accurate, too slow*/
	epsilon = 1.0e-10;
	epsilon_max = 1.0e-2;

	long double crit_max = 1.0e-2;
	Eigen::MatrixXd MatrixHesseB = Eigen::MatrixXd::Identity(*Noccup + 2, *Noccup + 2);
	Eigen::MatrixXd MatrixPB = Eigen::MatrixXd::Random(*Noccup + 2, 1);
	Eigen::MatrixXd MatrixSB = Eigen::MatrixXd::Random(*Noccup + 2, 1);
	Eigen::MatrixXd MatrixYB = Eigen::MatrixXd::Random(*Noccup + 2, 1);
	Eigen::MatrixXd MatrixBetaB_current = RegEle->MatrixBetaB;
	expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
	estimate_logbetadenseB(RegEle, dataset);

	bool init_H = true;
	long double test_mesh = epsilon;
	long double der_mesh = epsilon;
	int count = 0;
	std::vector<long double> derB_before, derB_current;
	for (int i = 0; i < *Noccup + 2; i++) {
		derB_before.push_back(0.0);
		derB_current.push_back(0.0);
	}

	long double fit_phiB_current = RegEle->fit_phiB;
	der_mesh = epsilon;
	long double previous_derindB = 0.0;
	bool init_trial = true;
	long double logbetadenseB_current = RegEle->logbetadenseB;
	long double derindB;
	long double metadelta_mesh = test_mesh;
	long double best_mesh = 0.0;
	long double saved_best_mesh = 0.0;
	long double previous_logbetadenseB = RegEle->logbetadenseB;
	int delta_monitoring = 0;

	while (true) {
		RegEle->fit_phiB = fit_phiB_current;
		RegEle->fit_phiB += der_mesh;
		expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
		estimate_logbetadenseB(RegEle, dataset);
		derindB = (RegEle->logbetadenseB - logbetadenseB_current) / der_mesh;
		if (init_trial) {
			previous_derindB = derindB;
		}
		if (!init_trial && (fabsl((previous_derindB - derindB) / derindB) < crit_max || der_mesh > epsilon_max)) {
			break;
		}
		else {
			init_trial = false;
			previous_derindB = derindB;
			der_mesh *= 10.0;
		}
	}
	test_mesh = der_mesh;

	metadelta_mesh = test_mesh;
	best_mesh = 0.0;
	previous_logbetadenseB = RegEle->logbetadenseB;
	delta_monitoring = 0;
	for (long double trial_mesh = test_mesh; true; trial_mesh += metadelta_mesh) {
		RegEle->MatrixBetaB = MatrixBetaB_current;
		RegEle->fit_phiB = fit_phiB_current;
		RegEle->fit_phiB += trial_mesh * derindB;
		expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
		estimate_logbetadenseB(RegEle, dataset);
		if (RegEle->logbetadenseB > previous_logbetadenseB) {
			previous_logbetadenseB = RegEle->logbetadenseB;
			best_mesh = trial_mesh;
			delta_monitoring++;
			if (delta_monitoring == 2) {
				metadelta_mesh *= 2.0;
				delta_monitoring = 0;
			}
		}
		else {
			break;
		}
	}
	RegEle->MatrixBetaB = MatrixBetaB_current;
	RegEle->fit_phiB = fit_phiB_current;
	RegEle->fit_phiB += best_mesh * derindB;
	expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
	estimate_logbetadenseB(RegEle, dataset);
	MatrixBetaB_current = RegEle->MatrixBetaB;
	fit_phiB_current = RegEle->fit_phiB;
	logbetadenseB_current = RegEle->logbetadenseB;
	previous_logbetadenseB = logbetadenseB_current;

	while (true) {

		MatrixBetaB_current = RegEle->MatrixBetaB;
		fit_phiB_current = RegEle->fit_phiB;
		expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
		estimate_logbetadenseB(RegEle, dataset);
		logbetadenseB_current = RegEle->logbetadenseB;
		long double logbetadenseB_start = RegEle->logbetadenseB;

		metadelta_mesh = test_mesh;
		best_mesh = 0.0;
		saved_best_mesh = 0.0;
		previous_logbetadenseB = RegEle->logbetadenseB;
		delta_monitoring = 0;

		long double previous_derindB = 0.0;
		bool init_trial = true;
		for (int i = 0; i < *Noccup + 1; i++) {
			der_mesh = epsilon;
			while (true) {
				RegEle->MatrixBetaB = MatrixBetaB_current;
				RegEle->MatrixBetaB(i, 0) += der_mesh;
				expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
				estimate_logbetadenseB(RegEle, dataset);
				long double vderB = (RegEle->logbetadenseB - logbetadenseB_current) / der_mesh;
				if (init_trial) {
					previous_derindB = vderB;
				}
				if (!init_trial && (fabsl((previous_derindB - vderB) / vderB) < crit_max || der_mesh > epsilon_max)) {
					derB_current[i] = vderB;
					break;
				}
				else {
					init_trial = false;
					previous_derindB = vderB;
					der_mesh *= 10.0;
				}
			}
		}

		der_mesh = epsilon;
		previous_derindB = 0.0;
		init_trial = true;
		while (true) {
			RegEle->MatrixBetaB = MatrixBetaB_current;
			RegEle->fit_phiB = fit_phiB_current + der_mesh;
			expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
			estimate_logbetadenseB(RegEle, dataset);
			long double vderB = (RegEle->logbetadenseB - logbetadenseB_current) / der_mesh;
			if (init_trial) {
				previous_derindB = vderB;
			}
			if (!init_trial && (fabsl((previous_derindB - vderB) / vderB) < crit_max || der_mesh > epsilon_max)) {
				derB_current[*Noccup + 1] = vderB;
				break;
			}
			else {
				init_trial = false;
				previous_derindB = vderB;
				der_mesh *= 10.0;
			}
		}

		if (!init_H) {
			for (int i = 0; i < *Noccup + 2; i++) {
				MatrixYB(i, 0) = derB_current[i] - derB_before[i];
			}
			long double STYB = (MatrixSB.transpose() * MatrixYB)(0, 0);
			long double YTHY = (MatrixYB.transpose() * MatrixHesseB * MatrixYB)(0, 0);
			Eigen::MatrixXd MatrixSSTB = MatrixSB * MatrixSB.transpose();
			Eigen::MatrixXd MatrixHYSTB = MatrixHesseB * MatrixYB * MatrixSB.transpose();
			Eigen::MatrixXd MatrixSYTHB = MatrixSB * MatrixYB.transpose() * MatrixHesseB;
			MatrixHesseB = MatrixHesseB + ((STYB + YTHY) / powl(STYB, 2.0)) * MatrixSSTB - (1.0 / STYB) * (MatrixHYSTB + MatrixSYTHB);
		}

		for (int i = 0; i < *Noccup + 2; i++) {
			MatrixPB(i, 0) = 0.0;
			for (int j = 0; j < *Noccup + 2; j++) {
				MatrixPB(i, 0) += MatrixHesseB(i, j) * derB_current[j];
			}
		}
		init_H = false;

		test_mesh = der_mesh;
		previous_logbetadenseB = RegEle->logbetadenseB;
		best_mesh = 0.0;
		metadelta_mesh = test_mesh;
		for (int ii = 0; ii < 2; ii++) {
			delta_monitoring = 0;
			for (long double trial_mesh = saved_best_mesh + test_mesh; true; trial_mesh += metadelta_mesh) {
				RegEle->MatrixBetaB = MatrixBetaB_current;
				RegEle->fit_phiB = fit_phiB_current;
				for (int i = 0; i < *Noccup + 1; i++) {
					RegEle->MatrixBetaB(i, 0) += trial_mesh * derB_current[i];
				}
				RegEle->fit_phiB += trial_mesh * derB_current[*Noccup + 1];
				expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
				estimate_logbetadenseB(RegEle, dataset);
				if (RegEle->logbetadenseB > previous_logbetadenseB) {
					previous_logbetadenseB = RegEle->logbetadenseB;
					best_mesh = trial_mesh;
					delta_monitoring++;
					if (delta_monitoring == 2) {
						metadelta_mesh *= 2.0;
						delta_monitoring = 0;
					}
				}
				else {
					break;
				}
			}
			saved_best_mesh = best_mesh;
		}

		MatrixSB = best_mesh * MatrixPB;

		RegEle->MatrixBetaB = MatrixBetaB_current;
		RegEle->fit_phiB = fit_phiB_current;
		for (int i = 0; i < *Noccup + 1; i++) {
			RegEle->MatrixBetaB(i, 0) += best_mesh * derB_current[i];
		}
		RegEle->fit_phiB += best_mesh * derB_current[*Noccup + 1];
		expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
		estimate_logbetadenseB(RegEle, dataset);

		MatrixBetaB_current = RegEle->MatrixBetaB;
		fit_phiB_current = RegEle->fit_phiB;
		logbetadenseB_current = RegEle->logbetadenseB;
		previous_logbetadenseB = logbetadenseB_current;

		for (int i = 0; i < *Noccup + 2; i++) {
			derB_before[i] = derB_current[i];
		}

		count++;
		if (RegEle->logbetadenseB < logbetadenseB_start || fabsl((RegEle->logbetadenseB - logbetadenseB_start) / logbetadenseB_start) < epsilon * fabsl(logbetadenseB_start) || best_mesh == 0.0 || count > *max_monitor) {
			break;
		}

	}

	std::vector<long double>().swap(derB_before);
	std::vector<long double>().swap(derB_current);

}

void measure_SStot_SSresB(REGELE* RegEle, std::string* functype_gB, std::vector<DATASET>* dataset) {

	long double ybar = 0.0;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		long double v = RegEle->MatrixYB(i, 0);
		ybar += get_v_special(functype_gB, &v);
	}
	ybar /= (long double)(signed)dataset->size();

	RegEle->SStotB = 0.0;
	RegEle->SSresB = 0.0;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		long double v1 = RegEle->MatrixYB(i, 0);
		long double v2 = RegEle->MatrixYexpectedB(i, 0);
		RegEle->SStotB += powl(get_v_special(functype_gB, &v1) - ybar, 2.0);
		RegEle->SSresB += powl(get_v_special(functype_gB, &v1) - get_v_special(functype_gB, &v2), 2.0);
	}

	ybar = 0.0;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		ybar += RegEle->MatrixYB(i, 0);
	}
	ybar /= (long double)(signed)dataset->size();

	RegEle->SStot = 0.0;
	RegEle->SSres = 0.0;
	RegEle->MSE = 0.0;
	RegEle->RMSE = 0.0;
	RegEle->MAE = 0.0;
	if (dataset->size() > 0) {
		for (int i = 0; i < (signed)dataset->size(); i++) {
			RegEle->SStot += powl(RegEle->MatrixYB(i, 0) - ybar, 2.0);
			RegEle->SSres += powl(RegEle->MatrixYB(i, 0) - RegEle->MatrixYexpectedB(i, 0), 2.0);
			RegEle->MAE += fabsl(RegEle->MatrixYB(i, 0) - RegEle->MatrixYexpectedB(i, 0));
		}
		RegEle->MSE = RegEle->SSres / (long double)((signed)dataset->size());
		RegEle->RMSE = sqrtl(RegEle->MSE);
		RegEle->MAE /= (long double)((signed)dataset->size());
	}

}

void measuere_r2B(REGELE* RegEle, std::vector<DATASET>* dataset) {

	RegEle->det_r2_pseudo = 1.0 - RegEle->SSresB / RegEle->SStotB;
	RegEle->det_r2 = 1.0 - RegEle->SSres / RegEle->SStot;
	if (!(RegEle->det_r2_pseudo >= 0.0 && RegEle->det_r2_pseudo <= 1.0)) {
		RegEle->det_r2_pseudo = -1;
	}
	if (!(RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0)) {
		RegEle->det_r2 = -1;
	}

}

void construct_FisherInformationB(REGELE* RegEle, std::string* functype_gB, std::vector<DATASET>* dataset, int* Noccup, long double* CI_level, long double* critic_p) {

	std::vector<long double> trigamma_myu_phi, trigamma_oneminusmyu_phi, gdermyu;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		long double myu_phi = RegEle->MatrixYexpectedB(i, 0) * RegEle->fit_phiB;
		long double oneminusmyu_phi = (1.0 - RegEle->MatrixYexpectedB(i, 0)) * RegEle->fit_phiB;
		trigamma_myu_phi.push_back(get_trigamma(&myu_phi));
		trigamma_oneminusmyu_phi.push_back(get_trigamma(&oneminusmyu_phi));
		long double myu = RegEle->MatrixYexpectedB(i, 0);
		gdermyu.push_back(get_derv_special(functype_gB, &myu));
	}

	for (int i = 0; i < (signed)dataset->size(); i++) {
		for (int j = 0; j < (signed)dataset->size(); j++) {
			if (i == j) {
				RegEle->MatrixWB(i, j) = RegEle->fit_phiB * (trigamma_myu_phi[i] + trigamma_oneminusmyu_phi[i]) / powl(gdermyu[i], 2.0);
				RegEle->MatrixTB(i, j) = 1.0 / gdermyu[i];
				RegEle->MatrixDB(i, j) = trigamma_myu_phi[i] * powl(RegEle->MatrixYexpectedB(i, 0), 2.0) + trigamma_oneminusmyu_phi[i] * powl(1.0 - RegEle->MatrixYexpectedB(i, 0), 2.0) - get_trigamma(&RegEle->fit_phiB);
			}
			else {
				RegEle->MatrixWB(i, j) = 0.0;
				RegEle->MatrixTB(i, j) = 0.0;
				RegEle->MatrixDB(i, j) = 0.0;
			}
		}
	}

	for (int i = 0; i < (signed)dataset->size(); i++) {
		RegEle->MatrixcB(i, 0) = RegEle->fit_phiB * (trigamma_myu_phi[i] * RegEle->MatrixYexpectedB(i, 0) + trigamma_oneminusmyu_phi[i] * (1.0 - RegEle->MatrixYexpectedB(i, 0)));
	}

	Eigen::MatrixXd MatrixXTWXiB = (RegEle->MatrixX.transpose() * RegEle->MatrixWB * RegEle->MatrixX).inverse();
	Eigen::MatrixXd MatrixXdXTWXiB = RegEle->MatrixX * MatrixXTWXiB;
	long double gammaB = RegEle->MatrixDB.trace() - (1.0 / RegEle->fit_phiB) * (RegEle->MatrixcB.transpose() * RegEle->MatrixTB.transpose() * MatrixXdXTWXiB * RegEle->MatrixX.transpose() * RegEle->MatrixTB * RegEle->MatrixcB)(0, 0);

	RegEle->MatrixFisherKubetabetaB = (1.0 / RegEle->fit_phiB) * MatrixXTWXiB * (Eigen::MatrixXd::Identity(*Noccup + 1, *Noccup + 1) + (1.0 / (gammaB * RegEle->fit_phiB)) * RegEle->MatrixX.transpose() * RegEle->MatrixTB * RegEle->MatrixcB * RegEle->MatrixcB.transpose() * RegEle->MatrixTB.transpose() * MatrixXdXTWXiB);
	RegEle->MatrixFisherKuphiphiB(0, 0) = 1.0 / RegEle->MatrixDB.trace();

	std::vector<long double>().swap(trigamma_myu_phi);
	std::vector<long double>().swap(trigamma_oneminusmyu_phi);
	std::vector<long double>().swap(gdermyu);

	std::string str_probit = "probit";
	long double CI_level_reset = 1.0 - (1.0 - *CI_level) / 2.0;
	long double CI_factor = get_v_special(&str_probit, &CI_level_reset);
	for (int i = 0; i < (signed)*Noccup + 1; i++) {
		long double SEB = sqrtl(RegEle->MatrixFisherKubetabetaB(i, i));
		RegEle->fit_coeffB.push_back(RegEle->MatrixBetaB(i, 0));
		RegEle->fit_coeff_SEB.push_back(SEB);
		RegEle->fit_coeff_tB.push_back(RegEle->MatrixBetaB(i, 0) / SEB);
		RegEle->fit_coeff_CI_startB.push_back(RegEle->MatrixBetaB(i, 0) - CI_factor * SEB);
		RegEle->fit_coeff_CI_endB.push_back(RegEle->MatrixBetaB(i, 0) + CI_factor * SEB);
	}
	if (RegEle->MatrixFisherKuphiphiB(0, 0) > 0.0) {
		RegEle->fit_phi_SEB = sqrtl(RegEle->MatrixFisherKuphiphiB(0, 0));
	}
	else {
		RegEle->fit_phi_SEB = RegEle->fit_phiB * 100000000000000000000000000000000000000000.0;
	}
	RegEle->fit_phi_tB = RegEle->fit_phiB / RegEle->fit_phi_SEB;
	RegEle->fit_phi_CI_startB = RegEle->fit_phiB - CI_factor * RegEle->fit_phi_SEB;
	RegEle->fit_phi_CI_endB = RegEle->fit_phiB + CI_factor * RegEle->fit_phi_SEB;

	std::string str_cdf = "cdf";
	for (int i = 0; i < (signed)RegEle->fit_coeff_tB.size(); i++) {
		long double tB = -fabsl(RegEle->fit_coeff_tB[i]);
		RegEle->fit_coeff_pB.push_back(2.0 * get_v_special(&str_cdf, &tB));
	}
	long double phi_tB = -fabsl(RegEle->fit_phi_tB);
	RegEle->fit_phi_pB = 2.0 * get_v_special(&str_cdf, &phi_tB);

	for (int j = 0; j < *Noccup + 1; j++) {
		RegEle->fit_coeff_standardizedB.push_back(0.0);
	}

	RegEle->fit_coeff_standardizedB[0] = RegEle->fit_coeffB[0];
	for (int j = 1; j < *Noccup + 1; j++) {
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
		RegEle->fit_coeff_standardizedB[j] = x_stdev * RegEle->fit_coeffB[j];
		RegEle->fit_coeff_standardizedB[0] += x_average * RegEle->fit_coeffB[j];
	}

}

void check_FrequentistB(REGELE* RegEle, long double* critic_p) {

	RegEle->full_Frequentist = true;

	if (RegEle->fit_phi_pB > *critic_p) {
		RegEle->full_Frequentist = false;
	}

	if (RegEle->full_Frequentist) {
		for (int i = 0; i < (signed)RegEle->fit_coeff_pB.size(); i++) {
			if (RegEle->fit_coeff_pB[i] > *critic_p) {
				RegEle->full_Frequentist = false;
				break;
			}
		}
	}

}

void order_pB(REGELE* RegEle, std::vector<int>* arrange) {

	std::vector<int>().swap(RegEle->fit_coeff_porderB);
	std::vector<int>().swap(RegEle->fit_coeff_porder_on_arrB);
	std::vector<long double>().swap(RegEle->fit_coeff_porderedB);

	bool* already_chk;
	already_chk = new bool[RegEle->fit_coeff_pB.size() - 1];
	for (int i = 0; i < (signed)RegEle->fit_coeff_pB.size() - 1; i++) {
		already_chk[i] = false;
	}

	for (int ii = 1; ii < RegEle->fit_coeff_pB.size(); ii++) {
		long double find_max = -1000000.0;
		int target_index = -1;
		for (int i = 1; i < (signed)RegEle->fit_coeff_pB.size(); i++) {
			if (already_chk[i - 1] == false && RegEle->fit_coeff_pB[i] > find_max) {
				find_max = RegEle->fit_coeff_pB[i];
				target_index = i - 1;
			}
		}
		if (target_index == -1) {
			for (int i = 1; i < (signed)RegEle->fit_coeff_pB.size(); i++) {
				if (already_chk[i - 1] == false) {
					target_index = i - 1;
				}
			}
		}
		RegEle->fit_coeff_porderB.push_back(target_index);
		RegEle->fit_coeff_pordered.push_back(RegEle->fit_coeff_pB[target_index + 1]);
		already_chk[target_index] = true;
		int find_pos = -1;
		for (int find_in_arr = 0; find_in_arr < (signed)arrange->size(); find_in_arr++) {
			if ((*arrange)[find_in_arr] != 0) {
				find_pos++;
				if (find_pos == target_index) {
					RegEle->fit_coeff_porder_on_arrB.push_back(find_in_arr);
					break;
				}
			}
		}

	}

	delete[] already_chk;
	already_chk = nullptr;

}

void do_regression(std::vector<FUNCTYPE>* FuncType, REGELE* RegEle, std::vector<DATASET>* dataset, std::vector<int>* arrange, int* Noccup, long double* CI_level, long double* critic_p, REGELE* trainedRegEle, bool* for_test) {

	initialize_regression(RegEle, dataset, Noccup);
	set_data(FuncType, RegEle, dataset, arrange, Noccup);
	RegEle->df = (signed)dataset->size() - (*Noccup + 1);
	std::vector<long double> x_stdev_series;
	if (!*for_test) {
		fit_any(RegEle, dataset, Noccup, &x_stdev_series);
	}
	measure_SStot_SSres(RegEle, dataset, trainedRegEle, for_test);
	measuere_r2(RegEle, dataset);
	if (!*for_test) {
		measure_t(RegEle, Noccup, &x_stdev_series);
		measure_p(RegEle);
		measure_confidence_interval(RegEle, Noccup, CI_level);
		measure_F(RegEle, Noccup);
		check_Frequentist(RegEle, critic_p);
	}

	std::vector<long double>().swap(x_stdev_series);

}

void do_regressionB(std::vector<FUNCTYPE>* FuncType, std::string* functype_gB, std::string* functype_giB, REGELE* RegEle, std::vector<DATASET>* dataset, std::vector<int>* arrange, int* Noccup, long double* CI_level, long double* critic_p, bool* BFGS, REGELE* trainedRegEle, int* max_monitor, bool* for_test) {

	initialize_B(RegEle, dataset, Noccup);
	set_data(FuncType, RegEle, dataset, arrange, Noccup);
	initialize_parametersB(RegEle, functype_gB, functype_giB, dataset, Noccup);
	RegEle->df = (signed)dataset->size() - (*Noccup + 1);
	RegEle->MatrixBetaB = RegEle->MatrixInitBetaB;
	RegEle->fit_phiB = RegEle->initphiB;
	if (!*for_test) {
		if (*BFGS) {
			maximize_logbetadenseB(RegEle, dataset, functype_giB, Noccup, for_test);
			maximize_logbetadense_BFGSB(RegEle, dataset, functype_giB, Noccup, max_monitor, for_test);
		}
		expect_YB(RegEle, dataset, functype_giB, Noccup, RegEle, for_test);
		estimate_logbetadenseB(RegEle, dataset);
	}
	else {
		expect_YB(RegEle, dataset, functype_giB, Noccup, trainedRegEle, for_test);
	}
	measure_SStot_SSresB(RegEle, functype_gB, dataset);
	measuere_r2B(RegEle, dataset);
	if (!*for_test) {
		construct_FisherInformationB(RegEle, functype_gB, dataset, Noccup, CI_level, critic_p);
		check_FrequentistB(RegEle, critic_p);
	}

}

void run_through_models(bool* xpr_ld_get, int* Nfeatures, int* NVariables, std::string* betachannel, std::vector<FUNCTYPE>* functype, int* DuplicateX) {

	auto run_start_t = std::chrono::system_clock::now();
	signed long long int RunningIndexDelta = 1;
	long double RunningIndexDelta_ld = 1.0;
	max_R2 = -1000000.0;
	max_testR2 = -1000000.0;
	max_rR2 = -1000000.0;
	FullFrequentist_exists = false;
	failure_count = 0;
	bool initial = true;

	int Monitoring = 0;
	if (pr_id == writeoutput_prid) {
		if (*betachannel == "") {
			writeoutput << "TrainedModel\tNdescriptors\tR2\tR2test\tR2_adj\tR2test_adj\tdf\tp(F)\n";
			writeoutput.flush();
		}
		else {
			writeoutput << "TrainedModel\tNdescriptors\tBeta\tR2\tR2test\tR2_pseudo\tR2test_pseudo\tdf\n";
			writeoutput.flush();
		}
	}

	if (*betachannel == "") {
		writeoutputpr << "TrainedModel\tNdescriptors\tR2\tR2test\tR2_adj\tR2test_adj\tdf\tp(F)\n";
		writeoutputpr.flush();
	}
	else {
		writeoutputpr << "TrainedModel\tNdescriptors\tBeta\tR2\tR2test\tR2_pseudo\tR2test_pseudo\tdf\n";
		writeoutputpr.flush();
	}

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
			kill_Jekyll(&arrange, DuplicateX, functype);
			if (pr_id == writeoutput_prid && Verbose == 1) {
				write_arrange(&writeoutput, &arrange);
			}
			if (Verbose == 1) {
				write_arrange(&writeoutputpr, &arrange);
			}
			if (*betachannel == "") {
				bool for_test = false;
				do_regression(&FuncType, &RegEle, &dataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &RegEle, &for_test);
				bool testonly = true;
				do_regression(&FuncType, &TestRegEle, &testdataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &RegEle, &testonly);
			}
			else {
				std::string ibetachannel = "i" + *betachannel;
				bool BFGS = false;
				bool for_test = false;
				do_regressionB(&FuncType, betachannel, &ibetachannel, &RegEle, &dataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &RegEle, &MaxBFGSIterations, &for_test);
				bool testonly = true;
				do_regressionB(&FuncType, betachannel, &ibetachannel, &TestRegEle, &testdataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &RegEle, &MaxBFGSIterations, &testonly);
			}
			if (pr_id == writeoutput_prid && Verbose == 1) {
				if (*betachannel == "") {
					writeoutput << "\t" << occup[1] << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_radj2 << "\t" << TestRegEle.det_radj2 << "\t" << RegEle.df << "\t" << RegEle.p_Fvalue;
					writeoutput.flush();
					if (!RegEle.full_Frequentist) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\t___RT___\n";
					writeoutput.flush();
				}
				else {
					writeoutput << "\t" << occup[1] << "\t" << *betachannel << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_r2_pseudo << "\t" << TestRegEle.det_r2_pseudo << "\t" << RegEle.df;
					writeoutput.flush();
					if (!RegEle.full_Frequentist) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\t___RT___\n";
					writeoutput.flush();
				}
			}
			if (Verbose == 1) {
				if (*betachannel == "") {
					writeoutputpr << "\t" << occup[1] << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_radj2 << "\t" << TestRegEle.det_radj2 << "\t" << RegEle.df << "\t" << RegEle.p_Fvalue;
					writeoutputpr.flush();
					if (!RegEle.full_Frequentist) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\t___RT___\n";
					writeoutputpr.flush();
				}
				else {
					writeoutputpr << "\t" << occup[1] << "\t" << *betachannel << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_r2_pseudo << "\t" << TestRegEle.det_r2_pseudo << "\t" << RegEle.df;
					writeoutputpr.flush();
					if (!RegEle.full_Frequentist) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\t___RT___\n";
					writeoutputpr.flush();
				}
			}
			if (!FullFrequentist_exists && !RegEle.full_Frequentist && RegEle.det_r2 > max_R2 && TestRegEle.det_r2 > max_testR2 && TestRegEle.det_r2 / RegEle.det_r2 > max_rR2 && RegEle.det_r2 >= 0.0 && RegEle.det_r2 <= 1.0) {
				max_R2 = RegEle.det_r2;
				max_testR2 = TestRegEle.det_r2;
				max_rR2 = TestRegEle.det_r2 / RegEle.det_r2;
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
			}
			else if (!FullFrequentist_exists && RegEle.full_Frequentist && RegEle.det_r2 >= 0.0 && RegEle.det_r2 <= 1.0) {
				FullFrequentist_exists = true;
				max_R2 = RegEle.det_r2;
				max_testR2 = TestRegEle.det_r2;
				max_rR2 = TestRegEle.det_r2 / RegEle.det_r2;
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
			}
			else if (FullFrequentist_exists && RegEle.full_Frequentist && RegEle.det_r2 > max_R2 && RegEle.det_r2 >= 0.0 && TestRegEle.det_r2 > max_testR2 && TestRegEle.det_r2 / RegEle.det_r2 > max_rR2 && RegEle.det_r2 <= 1.0) {
				max_R2 = RegEle.det_r2;
				max_testR2 = TestRegEle.det_r2;
				max_rR2 = TestRegEle.det_r2 / RegEle.det_r2;
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
			}
			else if (initial) {
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
				initial = false;
			}
			if (Monitoring <= 20) {
				Monitoring++;
			}
			if (Monitoring == 20) {
				auto run_current_t = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
				long double onesteptime = (long double)elapsed_secondsfr.count() / 20.0;
				RunningIndexDelta = ((long double)pr_job_size[pr_id]) / (TargetCPUTimeBeforeSwap / onesteptime);
				if (RunningIndexDelta == 0) {
					RunningIndexDelta++;
				}
				if (RunningIndexDelta > 1 && pr_id == writeoutput_prid && Verbose == 1) {
					writeoutput << "\n!!!!!!RunningIndexDelta=\t" << RunningIndexDelta << "\n";
					writeoutput.flush();
				}
				if (RunningIndexDelta > 1 && Verbose == 1) {
					writeoutputpr << "\n!!!!!!RunningIndexDelta=\t" << RunningIndexDelta << "\n";
					writeoutputpr.flush();
				}
			}

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
				std::cout << pr_id << " !call_XPR_ld\n";
				std::cout.flush();
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
				}
			}
			kill_Jekyll(&arrange, DuplicateX, functype);
			if (pr_id == writeoutput_prid && Verbose == 1) {
				write_arrange(&writeoutput, &arrange);
			}
			if (Verbose == 1) {
				write_arrange(&writeoutputpr, &arrange);
			}
			if (*betachannel == "") {
				bool for_test = false;
				do_regression(&FuncType, &RegEle, &dataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &RegEle, &for_test);
				bool testonly = true;
				do_regression(&FuncType, &TestRegEle, &testdataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &RegEle, &testonly);
			}
			else {
				std::string ibetachannel = "i" + *betachannel;
				bool BFGS = false;
				bool for_test = false;
				do_regressionB(&FuncType, betachannel, &ibetachannel, &RegEle, &dataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &RegEle, &MaxBFGSIterations, &for_test);
				bool testonly = true;
				do_regressionB(&FuncType, betachannel, &ibetachannel, &TestRegEle, &testdataset, &arrange, &occup[1], &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &RegEle, &MaxBFGSIterations, &testonly);
			}
			if (pr_id == writeoutput_prid && Verbose == 1) {
				if (*betachannel == "") {
					writeoutput << "\t" << occup[1] << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_radj2 << "\t" << TestRegEle.det_radj2 << "\t" << RegEle.df << "\t" << RegEle.p_Fvalue;
					writeoutput.flush();
					if (!RegEle.full_Frequentist) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\t___RT___\n";
					writeoutput.flush();
				}
				else {
					writeoutput << "\t" << occup[1] << "\t" << *betachannel << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_r2_pseudo << "\t" << TestRegEle.det_r2_pseudo << "\t" << RegEle.df;
					writeoutput.flush();
					if (!RegEle.full_Frequentist) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\t___RT___\n";
					writeoutput.flush();
				}
			}
			if (Verbose == 1) {
				if (*betachannel == "") {
					writeoutputpr << "\t" << occup[1] << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_radj2 << "\t" << TestRegEle.det_radj2 << "\t" << RegEle.df << "\t" << RegEle.p_Fvalue;
					writeoutputpr.flush();
					if (!RegEle.full_Frequentist) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\t___RT___\n";
					writeoutputpr.flush();
				}
				else {
					writeoutputpr << "\t" << occup[1] << "\t" << *betachannel << "\t" << RegEle.det_r2 << "\t" << TestRegEle.det_r2 << "\t" << RegEle.det_r2_pseudo << "\t" << TestRegEle.det_r2_pseudo << "\t" << RegEle.df;
					writeoutputpr.flush();
					if (!RegEle.full_Frequentist) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\t___RT___\n";
					writeoutputpr.flush();
				}
			}
			if (!FullFrequentist_exists && !RegEle.full_Frequentist && RegEle.det_r2 > max_R2 && TestRegEle.det_r2 > max_testR2 && TestRegEle.det_r2 / RegEle.det_r2 > max_rR2 && RegEle.det_r2 >= 0.0 && RegEle.det_r2 <= 1.0) {
				max_R2 = RegEle.det_r2;
				max_testR2 = TestRegEle.det_r2;
				max_rR2 = TestRegEle.det_r2 / RegEle.det_r2;
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
			}
			else if (!FullFrequentist_exists && RegEle.full_Frequentist && RegEle.det_r2 >= 0.0 && RegEle.det_r2 <= 1.0) {
				FullFrequentist_exists = true;
				max_R2 = RegEle.det_r2;
				max_testR2 = TestRegEle.det_r2;
				max_rR2 = TestRegEle.det_r2 / RegEle.det_r2;
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
			}
			else if (FullFrequentist_exists && RegEle.full_Frequentist && RegEle.det_r2 > max_R2 && TestRegEle.det_r2 > max_testR2 && TestRegEle.det_r2 / RegEle.det_r2 > max_rR2 && RegEle.det_r2 >= 0.0 && RegEle.det_r2 <= 1.0) {
				max_R2 = RegEle.det_r2;
				max_testR2 = TestRegEle.det_r2;
				max_rR2 = TestRegEle.det_r2 / RegEle.det_r2;
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
			}
			else if (initial) {
				good_occup = occup[1];
				memcpy_for_vectorint(&good_arrange, &arrange);
				GoodRegEle = RegEle;
				TestGoodRegEle = TestRegEle;
				initial = false;
			}
			if (Monitoring <= 20) {
				Monitoring++;
			}
			if (Monitoring == 20) {
				auto run_current_t = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
				long double onesteptime = (long double)elapsed_secondsfr.count() / 20.0;
				RunningIndexDelta_ld = (pr_job_size_ld[pr_id]) / (TargetCPUTimeBeforeSwap / onesteptime);
				if (RunningIndexDelta_ld < 1.0) {
					RunningIndexDelta_ld = 1.0;
				}
				if (RunningIndexDelta_ld > 1.0 && pr_id == writeoutput_prid && Verbose == 1) {
					writeoutput << "\n!!!!!!RunningIndexDelta_ld=\t" << RunningIndexDelta_ld << "\n";
					writeoutput.flush();
				}
				if (RunningIndexDelta_ld > 1.0 && Verbose == 1) {
					writeoutputpr << "\n!!!!!!RunningIndexDelta_ld=\t" << RunningIndexDelta_ld << "\n";
					writeoutputpr.flush();
				}
			}

		}
	}

	if (*betachannel == "") {
		ExistingCholskeySeries = false;
		bool for_test = false;
		do_regression(&FuncType, &GoodRegEle, &dataset, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &GoodRegEle, &for_test);
		bool testonly = true;
		do_regression(&FuncType, &TestGoodRegEle, &testdataset, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &GoodRegEle, &testonly);
	}
	else {
		std::string ibetachannel = "i" + *betachannel;
		bool BFGS = false;
		bool for_test = false;
		do_regressionB(&FuncType, betachannel, &ibetachannel, &GoodRegEle, &dataset, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &GoodRegEle, &MaxBFGSIterations, &for_test);
		bool testonly = true;
		do_regressionB(&FuncType, betachannel, &ibetachannel, &TestGoodRegEle, &testdataset, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &GoodRegEle, &MaxBFGSIterations, &testonly);
	}

	if (pr_id == writeoutput_prid) {
		write_arrange(&writeoutput, &good_arrange);
		if (*betachannel == "") {
			writeoutput << "\t" << good_occup << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_radj2 << "\t" << TestGoodRegEle.det_radj2 << "\t" << GoodRegEle.df << "\t" << GoodRegEle.p_Fvalue;
			writeoutput.flush();
		}
		else {
			writeoutput << "\t" << good_occup << "\t" << *betachannel << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_r2_pseudo << "\t" << TestGoodRegEle.det_r2_pseudo << "\t" << GoodRegEle.df;
			writeoutput.flush();
		}
		if (!FullFrequentist_exists) {
			writeoutput << "\tNF";
			writeoutput.flush();
		}
		writeoutput << "\t***RT***\n";
		writeoutput.flush();
	}
	write_arrange(&writeoutputpr, &good_arrange);
	if (*betachannel == "") {
		writeoutputpr << "\t" << good_occup << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_radj2 << "\t" << TestGoodRegEle.det_radj2 << "\t" << GoodRegEle.df << "\t" << GoodRegEle.p_Fvalue;
		writeoutputpr.flush();
	}
	else {
		writeoutputpr << "\t" << good_occup << "\t" << *betachannel << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_r2_pseudo << "\t" << TestGoodRegEle.det_r2_pseudo << "\t" << GoodRegEle.df;
		writeoutputpr.flush();
	}
	if (!FullFrequentist_exists) {
		writeoutputpr << "\tNF";
		writeoutputpr.flush();
		failure_count++;
	}
	writeoutputpr << "\t***RT***\n";
	writeoutputpr.flush();

}

void pass_data_from_master(std::string* prefix, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* max_R2, long double* max_testR2, long double* max_rR2, int* master_prid) {

	std::string regelebin_filename, datasetbin_filename, arrangebin_filename, Noccup_filename, maxR2_filename;
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".regele.smt";
	write_regelebin(RegEle, &regelebin_filename);
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".testregele.smt";
	write_regelebin(TestRegEle, &regelebin_filename);
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".goodregele.smt";
	write_regelebin(GoodRegEle, &regelebin_filename);
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".testgoodregele.smt";
	write_regelebin(TestGoodRegEle, &regelebin_filename);
	datasetbin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".dataset.smt";
	write_datasetbin(dataset, &datasetbin_filename);
	datasetbin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".testdataset.smt";
	write_datasetbin(testdataset, &datasetbin_filename);
	arrangebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".arrange.smt";
	write_arrangebin(arrange, &arrangebin_filename);
	arrangebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".goodarrange.smt";
	write_arrangebin(good_arrange, &arrangebin_filename);
	Noccup_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".Noccup.smt";
	write_Noccupbin(Noccup, &Noccup_filename);
	maxR2_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".maxR2.smt";
	write_maxR2bin(max_R2, max_testR2, max_rR2, &maxR2_filename);

}

void pass_data_to_slaves(std::string* prefix, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* max_R2, long double* max_testR2, long double* max_rR2, int* master_prid) {

	std::string regelebin_filename, datasetbin_filename, arrangebin_filename, Noccup_filename, maxR2_filename;
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".regele.smt";
	read_regelebin(RegEle, &regelebin_filename);
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".testregele.smt";
	read_regelebin(TestRegEle, &regelebin_filename);
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".goodregele.smt";
	read_regelebin(GoodRegEle, &regelebin_filename);
	regelebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".testgoodregele.smt";
	read_regelebin(TestGoodRegEle, &regelebin_filename);
	datasetbin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".dataset.smt";
	read_datasetbin(dataset, &datasetbin_filename);
	datasetbin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".testdataset.smt";
	read_datasetbin(testdataset, &datasetbin_filename);
	arrangebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".arrange.smt";
	read_arrangebin(arrange, &arrangebin_filename);
	arrangebin_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".goodarrange.smt";
	read_arrangebin(good_arrange, &arrangebin_filename);
	Noccup_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".Noccup.smt";
	read_Noccupbin(Noccup, &Noccup_filename);
	maxR2_filename = *prefix + "_M_" + std::to_string(*master_prid) + ".maxR2.smt";
	read_maxR2bin(max_R2, max_testR2, max_rR2, &maxR2_filename);

}

void divide_slaves(std::vector<int>* slave_prid_workcount, std::vector<int>* slave_prid_workstart, std::vector<int>* slave_prid_workend, int* current_pridposition, std::vector<int>* target_jobs, std::vector<int>* slave_prids) {

	slave_prid_workcount->clear();
	slave_prid_workstart->clear();
	slave_prid_workend->clear();
	for (int i = 0; i < (signed)slave_prids->size(); i++) {
		slave_prid_workcount->push_back(0);
		slave_prid_workstart->push_back(-1);
		slave_prid_workend->push_back(-1);
	}
	int pridpointer = 0;
	for (int j = 0; j < (signed)target_jobs->size(); j++) {
		(*slave_prid_workcount)[pridpointer]++;
		pridpointer++;
		if (pridpointer >= (signed)slave_prids->size()) {
			pridpointer = 0;
		}
	}
	(*slave_prid_workstart)[0] = 0;
	for (int i = 0; i < (signed)slave_prid_workcount->size(); i++) {
		if (i != 0) {
			(*slave_prid_workstart)[i] = (*slave_prid_workend)[i - 1] + 1;
		}
		if ((*slave_prid_workstart)[i] >= (signed)target_jobs->size()) {
			(*slave_prid_workstart)[i] = (signed)target_jobs->size() - 1;
		}
		if ((*slave_prid_workcount)[i] > 0) {
			(*slave_prid_workend)[i] = (*slave_prid_workstart)[i] + (*slave_prid_workcount)[i] - 1;
		}
		else {
			(*slave_prid_workend)[i] = (*slave_prid_workstart)[i];
		}
	}
	*current_pridposition = -1;
	for (int i = 0; i < (signed)slave_prids->size(); i++) {
		if ((*slave_prids)[i] == pr_id) {
			*current_pridposition = i;
			break;
		}
	}

}

void rm_smt_all(std::string* prefix, int* master_prid) {

	char cmd[500];
	std::string masterrm = "del " + *prefix + "_M_" + std::to_string(*master_prid) + ".*.smt";
#ifdef _WIN32
	strcpy_s(cmd, sizeof(cmd), masterrm.c_str());
	system(cmd);
#elif _WIN64
	strcpy_s(cmd, sizeof(cmd), masterrm.c_str());
	system(cmd);
#else
	masterrm = "rm -f " + *prefix + "_M_" + std::to_string(*master_prid) + ".*.smt";
	strcpy(cmd, masterrm.c_str());
	system(cmd);
#endif

}

void rm_smt_regele(std::string* prefix, int* master_prid) {

	char cmd[500];
	std::string mastertrainrm = "del " + *prefix + "_S*_M_" + std::to_string(*master_prid) + ".train.smt 2>nul";
	std::string mastertestrm = "del " + *prefix + "_S*_M_" + std::to_string(*master_prid) + ".test.smt 2>nul";
#ifdef _WIN32
	strcpy_s(cmd, sizeof(cmd), mastertrainrm.c_str());
	system(cmd);
	strcpy_s(cmd, sizeof(cmd), mastertestrm.c_str());
	system(cmd);
#elif _WIN64
	strcpy_s(cmd, sizeof(cmd), mastertrainrm.c_str());
	system(cmd);
	strcpy_s(cmd, sizeof(cmd), mastertestrm.c_str());
	system(cmd);
#else
	mastertrainrm = "rm -f " + *prefix + "_S*_M_" + std::to_string(*master_prid) + ".train.smt";
	mastertestrm = "rm -f " + *prefix + "_S*_M_" + std::to_string(*master_prid) + ".test.smt";
	strcpy(cmd, mastertrainrm.c_str());
	system(cmd);
	strcpy(cmd, mastertestrm.c_str());
	system(cmd);
#endif

}

bool add_descriptors(std::vector<FUNCTYPE>* FuncType, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* CI_level, long double* critic_p, long double* R2Threshold, std::string* betachannel, int* max_montior, bool* for_test, int* DuplicateX, std::vector<int>* slave_prids, int* master_prid) {

	ExistingCholskeySeries = false;
	bool gnew_GoodRegEle_found = false;
	auto run_start_t = std::chrono::system_clock::now();

	bool communicate = true;
	if (slave_prids->size() == 1 && (*slave_prids)[0] == *master_prid) {
		communicate = false;
	}
	
	std::string ms_prefix;
#ifdef _WIN32
	ms_prefix = "smt\\add_descriptors";
#elif _WIN64
	ms_prefix = "smt\\add_descriptors";
#else
	ms_prefix = "smt//add_descriptors";
#endif

	std::string regelebin_filename;
	if (communicate && pr_id == *master_prid) {
		pass_data_from_master(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (communicate) {
		pass_data_to_slaves(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<int> slave_prid_workcount, slave_prid_workstart, slave_prid_workend;
	int current_pridposition;
	divide_slaves(&slave_prid_workcount, &slave_prid_workstart, &slave_prid_workend, &current_pridposition, good_arrange, slave_prids);
	MPI_Barrier(MPI_COMM_WORLD);
	if (communicate && pr_id == *master_prid) {
		rm_smt_all(&ms_prefix, master_prid);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	REGELE new_GoodRegEle, new_TestGoodRegEle;
	std::vector<int> new_good_arrange;
	bool new_GoodRegEle_found = false;
	if (slave_prid_workcount[current_pridposition] > 0) {
		for (int j = slave_prid_workstart[current_pridposition]; j <= slave_prid_workend[current_pridposition]; j++) {

			int to = j;
			memcpy_for_vectorint(arrange, good_arrange);
			if ((*good_arrange)[j] == 0) {
				(*arrange)[to] = 1;
				kill_Jekyll(arrange, DuplicateX, FuncType);
				if (*betachannel == "") {
					do_regression(FuncType, RegEle, dataset, arrange, Noccup, CI_level, critic_p, RegEle, for_test);
					bool testonly = true;
					do_regression(FuncType, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, RegEle, &testonly);
				}
				else {
					std::string ibetachannel = "i" + *betachannel;
					bool BFGS = false;
					do_regressionB(FuncType, betachannel, &ibetachannel, RegEle, dataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, for_test);
					bool testonly = true;
					do_regressionB(FuncType, betachannel, &ibetachannel, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, &testonly);
				}
				if (pr_id == writeoutput_prid && Verbose == 1) {
					write_arrange(&writeoutput, arrange);
					if (*betachannel == "") {
						writeoutput << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
						writeoutput.flush();
					}
					else {
						writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
						writeoutput.flush();
					}
					if (!RegEle->full_Frequentist) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\t___AD___";
					writeoutput.flush();
					if (*master_prid != pr_id) {
						writeoutput << "\tM" << *master_prid;
						writeoutput.flush();
					}
					writeoutput << "\n";
					writeoutput.flush();
				}
				if (Verbose == 1) {
					write_arrange(&writeoutputpr, arrange);
					if (*betachannel == "") {
						writeoutputpr << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
						writeoutputpr.flush();
					}
					else {
						writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
						writeoutputpr.flush();
					}
					if (!RegEle->full_Frequentist) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\t___AD___";
					writeoutputpr.flush();
					if (*master_prid != pr_id) {
						writeoutputpr << "\tM" << *master_prid;
						writeoutputpr.flush();
					}
					writeoutputpr << "\n";
					writeoutputpr.flush();
				}

				regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".train.smt";
				write_regelebin(RegEle, &regelebin_filename);
				regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".test.smt";
				write_regelebin(TestRegEle, &regelebin_filename);

			}

		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (pr_id == *master_prid) {
		for (int j = 0; j < (signed)good_arrange->size(); j++) {
			int to = j;
			memcpy_for_vectorint(arrange, good_arrange);
			if ((*good_arrange)[j] == 0) {
				(*arrange)[to] = 1;
				kill_Jekyll(arrange, DuplicateX, FuncType);
				regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".train.smt";
				if (read_regelebin(RegEle, &regelebin_filename)) {
					regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".test.smt";
					if (read_regelebin(TestRegEle, &regelebin_filename)) {
						if (!FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && RegEle->det_r2 - max_R2 > *R2Threshold && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
							max_R2 = RegEle->det_r2;
							max_testR2 = TestRegEle->det_r2;
							max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							good_occup = *Noccup;
							memcpy_for_vectorint(&new_good_arrange, arrange);
							new_GoodRegEle = *RegEle;
							new_TestGoodRegEle = *TestRegEle;
							new_GoodRegEle_found = true;
						}
						else if (!FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
							FullFrequentist_exists = true;
							max_R2 = RegEle->det_r2;
							max_testR2 = TestRegEle->det_r2;
							max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							good_occup = *Noccup;
							memcpy_for_vectorint(&new_good_arrange, arrange);
							new_GoodRegEle = *RegEle;
							new_TestGoodRegEle = *TestRegEle;
							new_GoodRegEle_found = true;
						}
						else if (FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && RegEle->det_r2 - max_R2 > *R2Threshold && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
							max_R2 = RegEle->det_r2;
							max_testR2 = TestRegEle->det_r2;
							max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							good_occup = *Noccup;
							memcpy_for_vectorint(&new_good_arrange, arrange);
							new_GoodRegEle = *RegEle;
							new_TestGoodRegEle = *TestRegEle;
							new_GoodRegEle_found = true;
						}
					}
				}
			}
		}
		rm_smt_regele(&ms_prefix, master_prid);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (pr_id == *master_prid) {
		if (new_GoodRegEle_found) {
			gnew_GoodRegEle_found = true;
			good_occup = *Noccup;
			memcpy_for_vectorint(good_arrange, &new_good_arrange);
			*GoodRegEle = new_GoodRegEle;
			*TestGoodRegEle = new_TestGoodRegEle;
			if (pr_id == writeoutput_prid && Verbose == 1) {
				write_arrange(&writeoutput, good_arrange);
				if (*betachannel == "") {
					writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
					writeoutput.flush();
				}
				else {
					writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
					writeoutput.flush();
				}
				if (!FullFrequentist_exists) {
					writeoutput << "\tNF";
					writeoutput.flush();
				}
				writeoutput << "\n";
				writeoutput.flush();
			}
			if (Verbose == 1) {
				write_arrange(&writeoutputpr, good_arrange);
				if (*betachannel == "") {
					writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
					writeoutputpr.flush();
				}
				else {
					writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
					writeoutputpr.flush();
				}
				if (!FullFrequentist_exists) {
					writeoutputpr << "\tNF";
					writeoutputpr.flush();
				}
				writeoutputpr << "\n";
				writeoutputpr.flush();
			}
			std::vector<int>().swap(new_good_arrange);
		}

		if (*betachannel == "") {
			ExistingCholskeySeries = false;
			do_regression(FuncType, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, for_test);
			bool testonly = true;
			do_regression(FuncType, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, &testonly);
		}
		else {
			std::string ibetachannel = "i" + *betachannel;
			bool BFGS = false;
			do_regressionB(FuncType, betachannel, &ibetachannel, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, for_test);
			bool testonly = true;
			do_regressionB(FuncType, betachannel, &ibetachannel, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, &testonly);
		}

		if (gnew_GoodRegEle_found && pr_id == writeoutput_prid) {
			write_arrange(&writeoutput, good_arrange);
			if (*betachannel == "") {
				writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutput.flush();
			}
			else {
				writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutput.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutput << "\tNF";
				writeoutput.flush();
			}
			writeoutput << "\t***AD***\n";
			writeoutput.flush();
		}
		if (gnew_GoodRegEle_found) {
			write_arrange(&writeoutputpr, good_arrange);
			if (*betachannel == "") {
				writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutputpr.flush();
			}
			else {
				writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutputpr.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutputpr << "\tNF";
				writeoutputpr.flush();
				failure_count++;
			}
			writeoutputpr << "\t***AD***\n";
			writeoutputpr.flush();
		}
	}

	std::vector<int>().swap(slave_prid_workcount);
	std::vector<int>().swap(slave_prid_workstart);
	std::vector<int>().swap(slave_prid_workend);

	MPI_Barrier(MPI_COMM_WORLD);
	return gnew_GoodRegEle_found;

}

bool add_descriptors_Luhmann(std::vector<FUNCTYPE>* FuncType, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* CI_level, long double* critic_p, long double* R2Threshold, std::string* betachannel, int* max_montior, bool* for_test, int* DuplicateX) {

	ExistingCholskeySeries = false;
	bool gnew_GoodRegEle_found = false;
	auto run_start_t = std::chrono::system_clock::now();

	REGELE new_GoodRegEle, new_TestGoodRegEle;
	std::vector<int> new_good_arrange;
	bool new_GoodRegEle_found = false;
	for (int j = 0; j < (signed)good_arrange->size(); j++) {

		int to = j;
		memcpy_for_vectorint(arrange, good_arrange);
		if ((*good_arrange)[j] == 0) {
			(*arrange)[to] = 1;
			kill_Jekyll(arrange, DuplicateX, FuncType);
			if (*betachannel == "") {
				do_regression(FuncType, RegEle, dataset, arrange, Noccup, CI_level, critic_p, RegEle, for_test);
				bool testonly = true;
				do_regression(FuncType, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, RegEle, &testonly);
			}
			else {
				std::string ibetachannel = "i" + *betachannel;
				bool BFGS = false;
				do_regressionB(FuncType, betachannel, &ibetachannel, RegEle, dataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, for_test);
				bool testonly = true;
				do_regressionB(FuncType, betachannel, &ibetachannel, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, &testonly);
			}
			if (pr_id == writeoutput_prid && Verbose == 1) {
				write_arrange(&writeoutput, arrange);
				if (*betachannel == "") {
					writeoutput << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
					writeoutput.flush();
				}
				else {
					writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
					writeoutput.flush();
				}
				if (!RegEle->full_Frequentist) {
					writeoutput << "\tNF";
					writeoutput.flush();
				}
				writeoutput << "\t___AD___\n";
				writeoutput.flush();
			}
			if (Verbose == 1) {
				write_arrange(&writeoutputpr, arrange);
				if (*betachannel == "") {
					writeoutputpr << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
					writeoutputpr.flush();
				}
				else {
					writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
					writeoutputpr.flush();
				}
				if (!RegEle->full_Frequentist) {
					writeoutputpr << "\tNF";
					writeoutputpr.flush();
				}
				writeoutputpr << "\t___AD___\n";
				writeoutputpr.flush();
			}
			if (!FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && RegEle->det_r2 - max_R2 > *R2Threshold && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
				max_R2 = RegEle->det_r2;
				max_testR2 = TestRegEle->det_r2;
				max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
				good_occup = *Noccup;
				memcpy_for_vectorint(&new_good_arrange, arrange);
				new_GoodRegEle = *RegEle;
				new_TestGoodRegEle = *TestRegEle;
				new_GoodRegEle_found = true;
			}
			else if (!FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
				FullFrequentist_exists = true;
				max_R2 = RegEle->det_r2;
				max_testR2 = TestRegEle->det_r2;
				max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
				good_occup = *Noccup;
				memcpy_for_vectorint(&new_good_arrange, arrange);
				new_GoodRegEle = *RegEle;
				new_TestGoodRegEle = *TestRegEle;
				new_GoodRegEle_found = true;
			}
			else if (FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && RegEle->det_r2 - max_R2 > *R2Threshold && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
				max_R2 = RegEle->det_r2;
				max_testR2 = TestRegEle->det_r2;
				max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
				good_occup = *Noccup;
				memcpy_for_vectorint(&new_good_arrange, arrange);
				new_GoodRegEle = *RegEle;
				new_TestGoodRegEle = *TestRegEle;
				new_GoodRegEle_found = true;
			}
		}

		auto run_current_t = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
		long double duration = (long double)elapsed_secondsfr.count();
		if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
			std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during add_descriptors.\n";
			std::cout.flush();
			break;
		}

	}

	if (new_GoodRegEle_found) {
		gnew_GoodRegEle_found = true;
		good_occup = *Noccup;
		memcpy_for_vectorint(good_arrange, &new_good_arrange);
		*GoodRegEle = new_GoodRegEle;
		*TestGoodRegEle = new_TestGoodRegEle;
		if (pr_id == writeoutput_prid && Verbose == 1) {
			write_arrange(&writeoutput, good_arrange);
			if (*betachannel == "") {
				writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutput.flush();
			}
			else {
				writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutput.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutput << "\tNF";
				writeoutput.flush();
			}
			writeoutput << "\n";
			writeoutput.flush();
		}
		if (Verbose == 1) {
			write_arrange(&writeoutputpr, good_arrange);
			if (*betachannel == "") {
				writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutputpr.flush();
			}
			else {
				writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutputpr.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutputpr << "\tNF";
				writeoutputpr.flush();
			}
			writeoutputpr << "\n";
			writeoutputpr.flush();
		}
		if (new_good_arrange.size() > 0) {
			std::vector<int>().swap(new_good_arrange);
		}
	}

	if (*betachannel == "") {
		ExistingCholskeySeries = false;
		do_regression(FuncType, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, for_test);
		bool testonly = true;
		do_regression(FuncType, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, &testonly);
	}
	else {
		std::string ibetachannel = "i" + *betachannel;
		bool BFGS = false;
		do_regressionB(FuncType, betachannel, &ibetachannel, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, for_test);
		bool testonly = true;
		do_regressionB(FuncType, betachannel, &ibetachannel, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, &testonly);
	}

	if (gnew_GoodRegEle_found && pr_id == writeoutput_prid) {
		write_arrange(&writeoutput, good_arrange);
		if (*betachannel == "") {
			writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
			writeoutput.flush();
		}
		else {
			writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
			writeoutput.flush();
		}
		if (!FullFrequentist_exists) {
			writeoutput << "\tNF";
			writeoutput.flush();
		}
		writeoutput << "\t***AD***\n";
		writeoutput.flush();
	}
	if (gnew_GoodRegEle_found) {
		write_arrange(&writeoutputpr, good_arrange);
		if (*betachannel == "") {
			writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
			writeoutputpr.flush();
		}
		else {
			writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
			writeoutputpr.flush();
		}
		if (!FullFrequentist_exists) {
			writeoutputpr << "\tNF";
			writeoutputpr.flush();
			failure_count++;
		}
		writeoutputpr << "\t***AD***\n";
		writeoutputpr.flush();
	}

	return gnew_GoodRegEle_found;

}

bool exit_mode(std::string* macro_filename_prefix) {

	int the_others_yet = 0;
	for (int iprid = 0; iprid < pr_size; iprid++) {
		if (iprid != pr_id) {
			std::string testfilename = *macro_filename_prefix + "." + std::to_string(iprid) + ".reg";
			std::ifstream readtestfile;
			readtestfile.open(testfilename);
			if (!readtestfile.is_open()) {
				the_others_yet++;
			}
			readtestfile.close();
		}
	}
	if ((long double)the_others_yet / (long double)pr_size < 0.01) {
		return true;
	}
	else {
		return false;
	}

}

bool swap_descriptors(std::vector<FUNCTYPE>* FuncType, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* CI_level, long double* critic_p, std::string* betachannel, int* max_montior, bool* for_test, int* DuplicateX, std::string* macro_filename_prefix, std::vector<int>* slave_prids, int* master_prid) {

	ExistingCholskeySeries = false;
	order_p(GoodRegEle, good_arrange);
	std::vector<int> from_list;
	memcpy_for_vectorint(&from_list, &GoodRegEle->fit_coeff_porder_on_arr);

	std::string ms_prefix;
#ifdef _WIN32
	ms_prefix = "smt\\swap_descriptors";
#elif _WIN64
	ms_prefix = "smt\\swap_descriptors";
#else
	ms_prefix = "smt//swap_descriptors";
#endif

	std::string regelebin_filename;
	std::vector<int> slave_prid_workcount, slave_prid_workstart, slave_prid_workend;
	int current_pridposition;
	bool communicate = true;
	if (slave_prids->size() == 1 && (*slave_prids)[0] == *master_prid) {
		communicate = false;
	}
	if (communicate && pr_id == *master_prid) {
		pass_data_from_master(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	int NNend = (signed)from_list.size() - 1;
	if (communicate) {
		pass_data_to_slaves(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
		order_p(GoodRegEle, good_arrange);
		memcpy_for_vectorint(&from_list, &GoodRegEle->fit_coeff_porder_on_arr);
	}
	if (useSomeShotMode) {
		NNend = 0;
		for (int j = 0; j < (signed)(*GoodRegEle).fit_coeff_pordered.size(); j++) {
			if ((*GoodRegEle).fit_coeff_pordered[j] >= *critic_p) {
				NNend = j;
			}
			else {
				break;
			}
		}
		if (pr_id == writeoutput_prid) {
			writeoutput << "swap_descriptors: useSomeShotMode: 0 --> NNend= " << NNend << "\n";
			writeoutput.flush();
		}
		writeoutputpr << "swap_descriptors: useSomeShotMode: 0 --> NNend= " << NNend << "\n";
		writeoutputpr.flush();
	}
	int Nend = (signed)from_list.size() - 1;
	MPI_Barrier(MPI_COMM_WORLD);

	bool gnew_GoodRegEle_found = false;
	auto run_start_t = std::chrono::system_clock::now();

	bool chance_added = false;
	for (int i = 0; i <= Nend; i++) {
		
		if (communicate && pr_id == *master_prid && i <= NNend) {
			pass_data_from_master(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate && i <= NNend) {
			pass_data_to_slaves(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (i <= NNend) {
			divide_slaves(&slave_prid_workcount, &slave_prid_workstart, &slave_prid_workend, &current_pridposition, good_arrange, slave_prids);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate && pr_id == *master_prid && i <= NNend) {
			rm_smt_all(&ms_prefix, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		int from = from_list[i];
		REGELE new_GoodRegEle, new_TestGoodRegEle;
		std::vector<int> new_good_arrange;
		bool new_GoodRegEle_found = false;
		
		if (slave_prid_workcount[current_pridposition] > 0 && i <= NNend) {

			bool local_FullFrequentist_exists = FullFrequentist_exists;
			long double local_max_R2 = max_R2;
			long double local_max_testR2 = max_testR2;
			long double local_max_rR2 = max_rR2;
			int target_j = -1;
			REGELE target_regele, target_testregele;

			for (int j = slave_prid_workstart[current_pridposition]; j <= slave_prid_workend[current_pridposition]; j++) {

				int to = j;
				memcpy_for_vectorint(arrange, good_arrange);
				if ((*good_arrange)[j] == 0) {
					(*arrange)[from] = 0;
					(*arrange)[to] = 1;
					kill_Jekyll(arrange, DuplicateX, FuncType);
					if (*betachannel == "") {
						do_regression(FuncType, RegEle, dataset, arrange, Noccup, CI_level, critic_p, RegEle, for_test);
						bool testonly = true;
						do_regression(FuncType, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, RegEle, &testonly);
					}
					else {
						std::string ibetachannel = "i" + *betachannel;
						bool BFGS = false;
						do_regressionB(FuncType, betachannel, &ibetachannel, RegEle, dataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, for_test);
						bool testonly = true;
						do_regressionB(FuncType, betachannel, &ibetachannel, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, &testonly);
					}
					if (pr_id == writeoutput_prid && Verbose == 1) {
						write_arrange(&writeoutput, arrange);
						if (*betachannel == "") {
							writeoutput << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
							writeoutput.flush();
						}
						else {
							writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
							writeoutput.flush();
						}
						if (!RegEle->full_Frequentist) {
							writeoutput << "\tNF";
							writeoutput.flush();
						}
						writeoutput << "\t___SD___";
						writeoutput.flush();
						if (*master_prid != pr_id) {
							writeoutput << "\tM" << *master_prid;
							writeoutput.flush();
						}
						writeoutput << "\n";
						writeoutput.flush();
					}
					if (Verbose == 1) {
						write_arrange(&writeoutputpr, arrange);
						if (*betachannel == "") {
							writeoutputpr << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
							writeoutputpr.flush();
						}
						else {
							writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
							writeoutputpr.flush();
						}
						if (!RegEle->full_Frequentist) {
							writeoutputpr << "\tNF";
							writeoutputpr.flush();
						}
						writeoutputpr << "\t___SD___";
						writeoutputpr.flush();
						if (*master_prid != pr_id) {
							writeoutputpr << "\tM" << *master_prid;
							writeoutputpr.flush();
						}
						writeoutputpr << "\n";
						writeoutputpr.flush();
					}

					if (!local_FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > local_max_R2 && TestRegEle->det_r2 > local_max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > local_max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
						local_max_R2 = RegEle->det_r2;
						local_max_testR2 = TestRegEle->det_r2;
						local_max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
						target_j = j;
						target_regele = *RegEle;
						target_testregele = *TestRegEle;
					}
					else if (!local_FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
						local_FullFrequentist_exists = true;
						local_max_R2 = RegEle->det_r2;
						local_max_testR2 = TestRegEle->det_r2;
						local_max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
						target_j = j;
						target_regele = *RegEle;
						target_testregele = *TestRegEle;
					}
					else if (local_FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > local_max_R2 && TestRegEle->det_r2 > local_max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > local_max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
						local_max_R2 = RegEle->det_r2;
						local_max_testR2 = TestRegEle->det_r2;
						local_max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
						target_j = j;
						target_regele = *RegEle;
						target_testregele = *TestRegEle;
					}

				}

				if (useExitMode && exit_mode(macro_filename_prefix)) {
					break;
				}
				auto run_current_t = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
				long double duration = (long double)elapsed_secondsfr.count();
				if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
					std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swap_descriptors.\n";
					std::cout.flush();
					break;
				}

			}

			if (target_j != -1) {
				regelebin_filename = ms_prefix + "_S" + std::to_string(target_j) + "_M_" + std::to_string(*master_prid) + ".train.smt";
				write_regelebin(&target_regele, &regelebin_filename);
				regelebin_filename = ms_prefix + "_S" + std::to_string(target_j) + "_M_" + std::to_string(*master_prid) + ".test.smt";
				write_regelebin(&target_testregele, &regelebin_filename);
			}

		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (pr_id == *master_prid && i <= NNend) {
			for (int j = 0; j < (signed)good_arrange->size(); j++) {
				int to = j;
				memcpy_for_vectorint(arrange, good_arrange);
				if ((*good_arrange)[j] == 0) {
					(*arrange)[from] = 0;
					(*arrange)[to] = 1;
					kill_Jekyll(arrange, DuplicateX, FuncType);
					regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".train.smt";
					if (read_regelebin(RegEle, &regelebin_filename)) {
						regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".test.smt";
						if (read_regelebin(TestRegEle, &regelebin_filename)) {
							if (!FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
								max_R2 = RegEle->det_r2;
								max_testR2 = TestRegEle->det_r2;
								max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
								good_occup = *Noccup;
								memcpy_for_vectorint(&new_good_arrange, arrange);
								new_GoodRegEle = *RegEle;
								new_TestGoodRegEle = *TestRegEle;
								new_GoodRegEle_found = true;
							}
							else if (!FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
								FullFrequentist_exists = true;
								max_R2 = RegEle->det_r2;
								max_testR2 = TestRegEle->det_r2;
								max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
								good_occup = *Noccup;
								memcpy_for_vectorint(&new_good_arrange, arrange);
								new_GoodRegEle = *RegEle;
								new_TestGoodRegEle = *TestRegEle;
								new_GoodRegEle_found = true;
							}
							else if (FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
								max_R2 = RegEle->det_r2;
								max_testR2 = TestRegEle->det_r2;
								max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
								good_occup = *Noccup;
								memcpy_for_vectorint(&new_good_arrange, arrange);
								new_GoodRegEle = *RegEle;
								new_TestGoodRegEle = *TestRegEle;
								new_GoodRegEle_found = true;
							}
						}
					}
				}
			}
			rm_smt_regele(&ms_prefix, master_prid);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (pr_id == *master_prid && i <= NNend) {
			if (new_GoodRegEle_found) {
				if (i == NNend && i < Nend) {
					NNend++;
				}
				gnew_GoodRegEle_found = true;
				good_occup = *Noccup;
				memcpy_for_vectorint(good_arrange, &new_good_arrange);
				*GoodRegEle = new_GoodRegEle;
				*TestGoodRegEle = new_TestGoodRegEle;
				if (pr_id == writeoutput_prid && Verbose == 1) {
					write_arrange(&writeoutput, good_arrange);
					if (*betachannel == "") {
						writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutput.flush();
					}
					else {
						writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutput.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\n";
					writeoutput.flush();
				}
				if (Verbose == 1) {
					write_arrange(&writeoutputpr, good_arrange);
					if (*betachannel == "") {
						writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutputpr.flush();
					}
					else {
						writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutputpr.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\n";
					writeoutputpr.flush();
				}
				std::vector<int>().swap(new_good_arrange);
			}
			else {
				if (!chance_added && i < Nend) {
					NNend = i + 1;
					chance_added = true;
				}
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);

		std::string NNendbin_filename = ms_prefix + "_M_" + std::to_string(*master_prid) + ".NNend.smt";
		if (communicate && pr_id == *master_prid) {
			write_NNendbin(&NNend, &NNendbin_filename);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate) {
			read_NNendbin(&NNend, &NNendbin_filename);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate && pr_id == *master_prid) {
			rm_smt_all(&ms_prefix, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (i <= NNend) {
			if (pr_id == writeoutput_prid) {
				writeoutput << "swap_descriptors: useSomeShotMode: 0 --> NNendUpdated = " << NNend << "\n";
				writeoutput.flush();
			}
			writeoutputpr << "swap_descriptors: useSomeShotMode: 0 --> NNendUpdated = " << NNend << "\n";
			writeoutputpr.flush();
		}
		
		if (useExitMode && exit_mode(macro_filename_prefix)) {
			break;
		}
		auto run_current_t = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
		long double duration = (long double)elapsed_secondsfr.count();
		if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
			std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swap_descriptors.\n";
			std::cout.flush();
			break;
		}
		
	}

	if (pr_id == *master_prid) {
		if (*betachannel == "") {
			ExistingCholskeySeries = false;
			do_regression(FuncType, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, for_test);
			bool testonly = true;
			do_regression(FuncType, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, &testonly);
		}
		else {
			std::string ibetachannel = "i" + *betachannel;
			bool BFGS = false;
			do_regressionB(FuncType, betachannel, &ibetachannel, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, for_test);
			bool testonly = true;
			do_regressionB(FuncType, betachannel, &ibetachannel, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, &testonly);
		}

		if (gnew_GoodRegEle_found && pr_id == writeoutput_prid) {
			write_arrange(&writeoutput, good_arrange);
			if (*betachannel == "") {
				writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutput.flush();
			}
			else {
				writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutput.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutput << "\tNF";
				writeoutput.flush();
			}
			writeoutput << "\t***SD***\n";
			writeoutput.flush();
		}

		if (gnew_GoodRegEle_found) {
			write_arrange(&writeoutputpr, good_arrange);
			if (*betachannel == "") {
				writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutputpr.flush();
			}
			else {
				writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutputpr.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutputpr << "\tNF";
				writeoutputpr.flush();
				failure_count++;
			}
			writeoutputpr << "\t***SD***\n";
			writeoutputpr.flush();
		}
	}

	std::vector<int>().swap(slave_prid_workcount);
	std::vector<int>().swap(slave_prid_workstart);
	std::vector<int>().swap(slave_prid_workend);
	std::vector<int>().swap(from_list);

	return gnew_GoodRegEle_found;

}

bool swap_descriptors_Luhmann(std::vector<FUNCTYPE>* FuncType, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* CI_level, long double* critic_p, std::string* betachannel, int* max_montior, bool* for_test, int* DuplicateX, std::string* macro_filename_prefix) {

	ExistingCholskeySeries = false;
	order_p(GoodRegEle, good_arrange);
	std::vector<int> from_list;
	memcpy_for_vectorint(&from_list, &GoodRegEle->fit_coeff_porder_on_arr);

	bool gnew_GoodRegEle_found = false;
	auto run_start_t = std::chrono::system_clock::now();

	int Nend = (signed)from_list.size() - 1;
	int NNend = (signed)from_list.size() - 1;
	if (useSomeShotMode) {
		NNend = 0;
		for (int j = 0; j < (signed)(*GoodRegEle).fit_coeff_pordered.size(); j++) {
			if ((*GoodRegEle).fit_coeff_pordered[j] >= *critic_p) {
				NNend = j;
			}
			else {
				break;
			}
		}
		if (pr_id == writeoutput_prid) {
			writeoutput << "swap_descriptors: useSomeShotMode: 0 --> NNend= " << NNend << "\n";
			writeoutput.flush();
		}
		writeoutputpr << "swap_descriptors: useSomeShotMode: 0 --> NNend= " << NNend << "\n";
		writeoutputpr.flush();
	}

	bool chance_added = false;
	for (int i = 0; i <= Nend; i++) {

		int from = from_list[i];
		REGELE new_GoodRegEle, new_TestGoodRegEle;
		std::vector<int> new_good_arrange;
		bool new_GoodRegEle_found = false;

		if (i <= NNend) {

			for (int j = 0; j < (signed)good_arrange->size(); j++) {

				int to = j;
				memcpy_for_vectorint(arrange, good_arrange);
				if ((*good_arrange)[j] == 0) {
					(*arrange)[from] = 0;
					(*arrange)[to] = 1;
					kill_Jekyll(arrange, DuplicateX, FuncType);
					if (*betachannel == "") {
						do_regression(FuncType, RegEle, dataset, arrange, Noccup, CI_level, critic_p, RegEle, for_test);
						bool testonly = true;
						do_regression(FuncType, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, RegEle, &testonly);
					}
					else {
						std::string ibetachannel = "i" + *betachannel;
						bool BFGS = false;
						do_regressionB(FuncType, betachannel, &ibetachannel, RegEle, dataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, for_test);
						bool testonly = true;
						do_regressionB(FuncType, betachannel, &ibetachannel, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, &testonly);
					}
					if (pr_id == writeoutput_prid && Verbose == 1) {
						write_arrange(&writeoutput, arrange);
						if (*betachannel == "") {
							writeoutput << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
							writeoutput.flush();
						}
						else {
							writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
							writeoutput.flush();
						}
						if (!RegEle->full_Frequentist) {
							writeoutput << "\tNF";
							writeoutput.flush();
						}
						writeoutput << "\t___SD___\n";
						writeoutput.flush();
					}
					if (Verbose == 1) {
						write_arrange(&writeoutputpr, arrange);
						if (*betachannel == "") {
							writeoutputpr << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
							writeoutputpr.flush();
						}
						else {
							writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
							writeoutputpr.flush();
						}
						if (!RegEle->full_Frequentist) {
							writeoutputpr << "\tNF";
							writeoutputpr.flush();
						}
						writeoutputpr << "\t___SD___\n";
						writeoutputpr.flush();
					}
					if (!FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
						max_R2 = RegEle->det_r2;
						max_testR2 = TestRegEle->det_r2;
						max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
						good_occup = *Noccup;
						memcpy_for_vectorint(&new_good_arrange, arrange);
						new_GoodRegEle = *RegEle;
						new_TestGoodRegEle = *TestRegEle;
						new_GoodRegEle_found = true;
					}
					else if (!FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
						FullFrequentist_exists = true;
						max_R2 = RegEle->det_r2;
						max_testR2 = TestRegEle->det_r2;
						max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
						good_occup = *Noccup;
						memcpy_for_vectorint(&new_good_arrange, arrange);
						new_GoodRegEle = *RegEle;
						new_TestGoodRegEle = *TestRegEle;
						new_GoodRegEle_found = true;
					}
					else if (FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
						max_R2 = RegEle->det_r2;
						max_testR2 = TestRegEle->det_r2;
						max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
						good_occup = *Noccup;
						memcpy_for_vectorint(&new_good_arrange, arrange);
						new_GoodRegEle = *RegEle;
						new_TestGoodRegEle = *TestRegEle;
						new_GoodRegEle_found = true;
					}
				}

				if (useExitMode && exit_mode(macro_filename_prefix)) {
					break;
				}
				auto run_current_t = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
				long double duration = (long double)elapsed_secondsfr.count();
				if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
					std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swap_descriptors.\n";
					std::cout.flush();
					break;
				}

			}

			if (new_GoodRegEle_found) {
				if (i == NNend && i < Nend) {
					NNend++;
					if (pr_id == writeoutput_prid) {
						writeoutput << "swap_descriptors: useSomeShotMode: 0 --> NNendUpdated = " << NNend << "\n";
						writeoutput.flush();
					}
					writeoutputpr << "swap_descriptors: useSomeShotMode: 0 --> NNendUpdated = " << NNend << "\n";
					writeoutputpr.flush();
				}
				gnew_GoodRegEle_found = true;
				good_occup = *Noccup;
				memcpy_for_vectorint(good_arrange, &new_good_arrange);
				*GoodRegEle = new_GoodRegEle;
				*TestGoodRegEle = new_TestGoodRegEle;
				if (pr_id == writeoutput_prid && Verbose == 1) {
					write_arrange(&writeoutput, good_arrange);
					if (*betachannel == "") {
						writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutput.flush();
					}
					else {
						writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutput.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\n";
					writeoutput.flush();
				}
				if (Verbose == 1) {
					write_arrange(&writeoutputpr, good_arrange);
					if (*betachannel == "") {
						writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutputpr.flush();
					}
					else {
						writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutputpr.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\n";
					writeoutputpr.flush();
				}
				if (new_good_arrange.size() > 0) {
					std::vector<int>().swap(new_good_arrange);
				}
			}
			else {
				if (!chance_added && i < Nend) {
					NNend = i + 1;
					chance_added = true;
					if (pr_id == writeoutput_prid) {
						writeoutput << "swap_descriptors: useSomeShotMode: 0 --> NNendUpdated = " << NNend << "\n";
						writeoutput.flush();
					}
					writeoutputpr << "swap_descriptors: useSomeShotMode: 0 --> NNendUpdated = " << NNend << "\n";
					writeoutputpr.flush();
				}
			}

		}
		
		if (useExitMode && exit_mode(macro_filename_prefix)) {
			break;
		}
		auto run_current_t = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
		long double duration = (long double)elapsed_secondsfr.count();
		if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
			std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swap_descriptors.\n";
			std::cout.flush();
			break;
		}

	}

	if (*betachannel == "") {
		ExistingCholskeySeries = false;
		do_regression(FuncType, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, for_test);
		bool testonly = true;
		do_regression(FuncType, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, &testonly);
	}
	else {
		std::string ibetachannel = "i" + *betachannel;
		bool BFGS = false;
		do_regressionB(FuncType, betachannel, &ibetachannel, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, for_test);
		bool testonly = true;
		do_regressionB(FuncType, betachannel, &ibetachannel, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, &testonly);
	}

	if (gnew_GoodRegEle_found && pr_id == writeoutput_prid) {
		write_arrange(&writeoutput, good_arrange);
		if (*betachannel == "") {
			writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
			writeoutput.flush();
		}
		else {
			writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
			writeoutput.flush();
		}
		if (!FullFrequentist_exists) {
			writeoutput << "\tNF";
			writeoutput.flush();
		}
		writeoutput << "\t***SD***\n";
		writeoutput.flush();
	}

	if (gnew_GoodRegEle_found) {
		write_arrange(&writeoutputpr, good_arrange);
		if (*betachannel == "") {
			writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
			writeoutputpr.flush();
		}
		else {
			writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
			writeoutputpr.flush();
		}
		if (!FullFrequentist_exists) {
			writeoutputpr << "\tNF";
			writeoutputpr.flush();
			failure_count++;
		}
		writeoutputpr << "\t***SD***\n";
		writeoutputpr.flush();
	}

	return gnew_GoodRegEle_found;

}

bool switch_functype(std::vector<FUNCTYPE>* FuncType, std::vector<GROUPX>* GroupX, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* CI_level, long double* critic_p, std::string* betachannel, int* max_montior, bool* for_test, int* DuplicateX, std::string* macro_filename_prefix, std::vector<int>* slave_prids, int* master_prid) {

	ExistingCholskeySeries = false;
	order_p(GoodRegEle, good_arrange);
	std::vector<int> from_list;
	memcpy_for_vectorint(&from_list, &GoodRegEle->fit_coeff_porder_on_arr);

	std::string ms_prefix;
#ifdef _WIN32
	ms_prefix = "smt\\switch_functype";
#elif _WIN64
	ms_prefix = "smt\\switch_functype";
#else
	ms_prefix = "smt//switch_functype";
#endif

	std::string regelebin_filename;
	std::vector<int> slave_prid_workcount, slave_prid_workstart, slave_prid_workend;
	int current_pridposition;
	bool communicate = true;
	if (slave_prids->size() == 1 && (*slave_prids)[0] == *master_prid) {
		communicate = false;
	}
	if (communicate && pr_id == *master_prid) {
		pass_data_from_master(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (communicate) {
		pass_data_to_slaves(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
		order_p(GoodRegEle, good_arrange);
		memcpy_for_vectorint(&from_list, &GoodRegEle->fit_coeff_porder_on_arr);
	}

	int NNend = 0;
	if (useSomeShotMode) {
		NNend = (signed)(*GoodRegEle).fit_coeff_pordered.size() - 1;
		for (int j = (signed)(*GoodRegEle).fit_coeff_pordered.size() - 1; j >= 0; j--) {
			if ((*GoodRegEle).fit_coeff_pordered[j] < *critic_p) {
				NNend = j;
			}
			else {
				break;
			}
		}
		if (pr_id == writeoutput_prid) {
			writeoutput << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNend = " << NNend << "\n";
			writeoutput.flush();
		}
		writeoutputpr << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNend = " << NNend << "\n";
		writeoutputpr.flush();
	}
	int Nend = 0;
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<int> iFuncType;
	for (int i = 0; i < (signed)FuncType->size(); i++) {
		iFuncType.push_back(i);
	}

	bool gnew_GoodRegEle_found = false;
	auto run_start_t = std::chrono::system_clock::now();

	bool chance_added = false;
	for (int i = (signed)from_list.size() - 1; i >= Nend; i--) {

		if (communicate && pr_id == *master_prid && i >= NNend) {
			pass_data_from_master(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate && i >= NNend) {
			pass_data_to_slaves(&ms_prefix, RegEle, TestRegEle, GoodRegEle, TestGoodRegEle, dataset, testdataset, arrange, good_arrange, Noccup, &max_R2, &max_testR2, &max_rR2, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (i >= NNend) {
			divide_slaves(&slave_prid_workcount, &slave_prid_workstart, &slave_prid_workend, &current_pridposition, &iFuncType, slave_prids);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate && pr_id == *master_prid && i >= NNend) {
			rm_smt_all(&ms_prefix, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		int from = from_list[i];
		REGELE new_GoodRegEle, new_TestGoodRegEle;
		std::vector<int> new_good_arrange;
		bool new_GoodRegEle_found = false;

		if (slave_prid_workcount[current_pridposition] > 0 && i >= NNend) {

			bool local_FullFrequentist_exists = FullFrequentist_exists;
			long double local_max_R2 = max_R2;
			long double local_max_testR2 = max_testR2;
			long double local_max_rR2 = max_rR2;
			int target_j = -1;
			REGELE target_regele, target_testregele;

			for (int j = slave_prid_workstart[current_pridposition]; j <= slave_prid_workend[current_pridposition]; j++) {

				int to = j;
				memcpy_for_vectorint(arrange, good_arrange);
				if ((*good_arrange)[from] != j + 1) {
					(*arrange)[from] = j + 1;
					bool go_on = true;
					if ((*GroupX)[from].negative_exists && (*FuncType)[j].FName == "log") {
						go_on = false;
					}
					if (go_on) {
						kill_Jekyll(arrange, DuplicateX, FuncType);
						if (*betachannel == "") {
							do_regression(FuncType, RegEle, dataset, arrange, Noccup, CI_level, critic_p, RegEle, for_test);
							bool testonly = true;
							do_regression(FuncType, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, RegEle, &testonly);
						}
						else {
							std::string ibetachannel = "i" + *betachannel;
							bool BFGS = false;
							do_regressionB(FuncType, betachannel, &ibetachannel, RegEle, dataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, for_test);
							bool testonly = true;
							do_regressionB(FuncType, betachannel, &ibetachannel, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, &testonly);
						}
						if (pr_id == writeoutput_prid && Verbose == 1) {
							write_arrange(&writeoutput, arrange);
							if (*betachannel == "") {
								writeoutput << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
								writeoutput.flush();
							}
							else {
								writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
								writeoutput.flush();
							}
							if (!RegEle->full_Frequentist) {
								writeoutput << "\tNF";
								writeoutput.flush();
							}
							writeoutput << "\t___SF___";
							writeoutput.flush();
							if (*master_prid != pr_id) {
								writeoutput << "\tM" << *master_prid;
								writeoutput.flush();
							}
							writeoutput << "\n";
							writeoutput.flush();
						}
						if (Verbose == 1) {
							write_arrange(&writeoutputpr, arrange);
							if (*betachannel == "") {
								writeoutputpr << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
								writeoutputpr.flush();
							}
							else {
								writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
								writeoutputpr.flush();
							}
							if (!RegEle->full_Frequentist) {
								writeoutputpr << "\tNF";
								writeoutputpr.flush();
							}
							writeoutputpr << "\t___SF___";
							writeoutputpr.flush();
							if (*master_prid != pr_id) {
								writeoutputpr << "\tM" << *master_prid;
								writeoutputpr.flush();
							}
							writeoutputpr << "\n";
							writeoutputpr.flush();
						}

						if (!local_FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > local_max_R2 && TestRegEle->det_r2 > local_max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > local_max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
							local_max_R2 = RegEle->det_r2;
							local_max_testR2 = TestRegEle->det_r2;
							local_max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							target_j = j;
							target_regele = *RegEle;
							target_testregele = *TestRegEle;
						}
						else if (!local_FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
							local_FullFrequentist_exists = true;
							local_max_R2 = RegEle->det_r2;
							local_max_testR2 = TestRegEle->det_r2;
							local_max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							target_j = j;
							target_regele = *RegEle;
							target_testregele = *TestRegEle;
						}
						else if (local_FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > local_max_R2 && TestRegEle->det_r2 > local_max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > local_max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
							local_max_R2 = RegEle->det_r2;
							local_max_testR2 = TestRegEle->det_r2;
							local_max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							target_j = j;
							target_regele = *RegEle;
							target_testregele = *TestRegEle;
						}

					}
				}

				regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".train.smt";
				write_regelebin(RegEle, &regelebin_filename);
				regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".test.smt";
				write_regelebin(TestRegEle, &regelebin_filename);

				if (useExitMode && exit_mode(macro_filename_prefix)) {
					break;
				}
				auto run_current_t = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
				long double duration = (long double)elapsed_secondsfr.count();
				if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
					std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swtich_functype.\n";
					std::cout.flush();
					break;
				}

			}

			if (target_j != -1) {
				regelebin_filename = ms_prefix + "_S" + std::to_string(target_j) + "_M_" + std::to_string(*master_prid) + ".train.smt";
				write_regelebin(&target_regele, &regelebin_filename);
				regelebin_filename = ms_prefix + "_S" + std::to_string(target_j) + "_M_" + std::to_string(*master_prid) + ".test.smt";
				write_regelebin(&target_testregele, &regelebin_filename);
			}

		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (pr_id == *master_prid && i >= NNend) {
			for (int j = 0; j < (signed)FuncType->size(); j++) {
				int to = j;
				memcpy_for_vectorint(arrange, good_arrange);
				if ((*good_arrange)[from] != j + 1) {
					(*arrange)[from] = j + 1;
					bool go_on = true;
					if ((*GroupX)[from].negative_exists && (*FuncType)[j].FName == "log") {
						go_on = false;
					}
					if (go_on) {
						kill_Jekyll(arrange, DuplicateX, FuncType);
					}
					regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".train.smt";
					if (read_regelebin(RegEle, &regelebin_filename)) {
						regelebin_filename = ms_prefix + "_S" + std::to_string(j) + "_M_" + std::to_string(*master_prid) + ".test.smt";
						if (read_regelebin(TestRegEle, &regelebin_filename)) {
							if (!FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
								max_R2 = RegEle->det_r2;
								max_testR2 = TestRegEle->det_r2;
								max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
								good_occup = *Noccup;
								memcpy_for_vectorint(&new_good_arrange, arrange);
								new_GoodRegEle = *RegEle;
								new_TestGoodRegEle = *TestRegEle;
								new_GoodRegEle_found = true;
							}
							else if (!FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
								FullFrequentist_exists = true;
								max_R2 = RegEle->det_r2;
								max_testR2 = TestRegEle->det_r2;
								max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
								good_occup = *Noccup;
								memcpy_for_vectorint(&new_good_arrange, arrange);
								new_GoodRegEle = *RegEle;
								new_TestGoodRegEle = *TestRegEle;
								new_GoodRegEle_found = true;
							}
							else if (FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0 && TestRegEle->det_r2 >= 0.0 && TestRegEle->det_r2 <= 1.0) {
								max_R2 = RegEle->det_r2;
								max_testR2 = TestRegEle->det_r2;
								max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
								good_occup = *Noccup;
								memcpy_for_vectorint(&new_good_arrange, arrange);
								new_GoodRegEle = *RegEle;
								new_TestGoodRegEle = *TestRegEle;
								new_GoodRegEle_found = true;
							}
						}
					}
				}
			}
			rm_smt_regele(&ms_prefix, master_prid);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (pr_id == *master_prid && i >= NNend) {
			if (new_GoodRegEle_found) {
				if (i == NNend && i > Nend) {
					NNend--;
				}
				gnew_GoodRegEle_found = true;
				good_occup = *Noccup;
				memcpy_for_vectorint(good_arrange, &new_good_arrange);
				*GoodRegEle = new_GoodRegEle;
				*TestGoodRegEle = new_TestGoodRegEle;
				if (pr_id == writeoutput_prid && Verbose == 1) {
					write_arrange(&writeoutput, good_arrange);
					if (*betachannel == "") {
						writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutput.flush();
					}
					else {
						writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutput.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\n";
					writeoutput.flush();
				}
				if (Verbose == 1) {
					write_arrange(&writeoutputpr, good_arrange);
					if (*betachannel == "") {
						writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutputpr.flush();
					}
					else {
						writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutputpr.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\n";
					writeoutputpr.flush();
				}
				std::vector<int>().swap(new_good_arrange);
			}
			else {
				if (!chance_added && i > Nend) {
					NNend = i - 1;
					chance_added = true;
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		std::string NNendbin_filename = ms_prefix + "_M_" + std::to_string(*master_prid) + ".NNend.smt";
		if (communicate && pr_id == *master_prid) {
			write_NNendbin(&NNend, &NNendbin_filename);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate) {
			read_NNendbin(&NNend, &NNendbin_filename);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (communicate && pr_id == *master_prid) {
			rm_smt_all(&ms_prefix, master_prid);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (i >= NNend) {
			if (pr_id == writeoutput_prid) {
				writeoutput << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNendUpdated = " << NNend << "\n";
				writeoutput.flush();
			}
			writeoutputpr << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNendUpdated = " << NNend << "\n";
			writeoutputpr.flush();
		}
		
		if (useExitMode && exit_mode(macro_filename_prefix)) {
			break;
		}
		auto run_current_t = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
		long double duration = (long double)elapsed_secondsfr.count();
		if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
			std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swtich_functype.\n";
			std::cout.flush();
			break;
		}

	}

	if (pr_id == *master_prid) {
		if (*betachannel == "") {
			ExistingCholskeySeries = false;
			do_regression(FuncType, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, for_test);
			bool testonly = true;
			do_regression(FuncType, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, &testonly);
		}
		else {
			std::string ibetachannel = "i" + *betachannel;
			bool BFGS = false;
			do_regressionB(FuncType, betachannel, &ibetachannel, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, for_test);
			bool testonly = true;
			do_regressionB(FuncType, betachannel, &ibetachannel, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, &testonly);
		}

		if (gnew_GoodRegEle_found && pr_id == writeoutput_prid) {
			write_arrange(&writeoutput, good_arrange);
			if (*betachannel == "") {
				writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutput.flush();
			}
			else {
				writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutput.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutput << "\tNF";
				writeoutput.flush();
			}
			writeoutput << "\t***SF***\n";
			writeoutput.flush();
		}

		if (gnew_GoodRegEle_found) {
			write_arrange(&writeoutputpr, good_arrange);
			if (*betachannel == "") {
				writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
				writeoutputpr.flush();
			}
			else {
				writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
				writeoutputpr.flush();
			}
			if (!FullFrequentist_exists) {
				writeoutputpr << "\tNF";
				writeoutputpr.flush();
				failure_count++;
			}
			writeoutputpr << "\t***SF***\n";
			writeoutputpr.flush();
		}
	}

	std::vector<int>().swap(slave_prid_workcount);
	std::vector<int>().swap(slave_prid_workstart);
	std::vector<int>().swap(slave_prid_workend);
	std::vector<int>().swap(iFuncType);
	std::vector<int>().swap(from_list);

	return gnew_GoodRegEle_found;

}

bool switch_functype_Luhmann(std::vector<FUNCTYPE>* FuncType, std::vector<GROUPX>* GroupX, REGELE* RegEle, REGELE* TestRegEle, REGELE* GoodRegEle, REGELE* TestGoodRegEle, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, std::vector<int>* arrange, std::vector<int>* good_arrange, int* Noccup, long double* CI_level, long double* critic_p, std::string* betachannel, int* max_montior, bool* for_test, int* DuplicateX, std::string* macro_filename_prefix) {

	ExistingCholskeySeries = false;
	order_p(GoodRegEle, good_arrange);
	std::vector<int> from_list;
	memcpy_for_vectorint(&from_list, &GoodRegEle->fit_coeff_porder_on_arr);

	bool gnew_GoodRegEle_found = false;
	auto run_start_t = std::chrono::system_clock::now();

	int Nend = 0;
	int NNend = 0;
	if (useSomeShotMode) {
		NNend = (signed)(*GoodRegEle).fit_coeff_pordered.size() - 1;
		for (int j = (signed)(*GoodRegEle).fit_coeff_pordered.size() - 1; j >= 0; j--) {
			if ((*GoodRegEle).fit_coeff_pordered[j] < *critic_p) {
				NNend = j;
			}
			else {
				break;
			}
		}
		if (pr_id == writeoutput_prid) {
			writeoutput << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNend = " << NNend << "\n";
			writeoutput.flush();
		}
		writeoutputpr << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNend = " << NNend << "\n";
		writeoutputpr.flush();
	}

	bool chance_added = false;
	for (int i = (signed)from_list.size() - 1; i >= Nend; i--) {

		int from = from_list[i];
		REGELE new_GoodRegEle, new_TestGoodRegEle;
		std::vector<int> new_good_arrange;
		bool new_GoodRegEle_found = false;

		if (i >= NNend) {

			for (int j = 0; j < (signed)FuncType->size(); j++) {

				int to = j;
				memcpy_for_vectorint(arrange, good_arrange);
				if ((*good_arrange)[from] != j + 1) {
					(*arrange)[from] = j + 1;
					bool go_on = true;
					if ((*GroupX)[from].negative_exists && (*FuncType)[j].FName == "log") {
						go_on = false;
					}
					if (go_on) {
						kill_Jekyll(arrange, DuplicateX, FuncType);
						if (*betachannel == "") {
							do_regression(FuncType, RegEle, dataset, arrange, Noccup, CI_level, critic_p, RegEle, for_test);
							bool testonly = true;
							do_regression(FuncType, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, RegEle, &testonly);
						}
						else {
							std::string ibetachannel = "i" + *betachannel;
							bool BFGS = false;
							do_regressionB(FuncType, betachannel, &ibetachannel, RegEle, dataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, for_test);
							bool testonly = true;
							do_regressionB(FuncType, betachannel, &ibetachannel, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, &testonly);
						}
						if (pr_id == writeoutput_prid && Verbose == 1) {
							write_arrange(&writeoutput, arrange);
							if (*betachannel == "") {
								writeoutput << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
								writeoutput.flush();
							}
							else {
								writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
								writeoutput.flush();
							}
							if (!RegEle->full_Frequentist) {
								writeoutput << "\tNF";
								writeoutput.flush();
							}
							writeoutput << "\t___SF___\n";
							writeoutput.flush();
						}
						if (Verbose == 1) {
							write_arrange(&writeoutputpr, arrange);
							if (*betachannel == "") {
								writeoutputpr << "\t" << *Noccup << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_radj2 << "\t" << TestRegEle->det_radj2 << "\t" << RegEle->df << "\t" << RegEle->p_Fvalue;
								writeoutputpr.flush();
							}
							else {
								writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << RegEle->det_r2 << "\t" << TestRegEle->det_r2 << "\t" << RegEle->det_r2_pseudo << "\t" << TestRegEle->det_r2_pseudo << "\t" << RegEle->df;
								writeoutputpr.flush();
							}
							if (!RegEle->full_Frequentist) {
								writeoutputpr << "\tNF";
								writeoutputpr.flush();
							}
							writeoutputpr << "\t___SF___\n";
							writeoutputpr.flush();
						}
						if (!FullFrequentist_exists && !RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
							max_R2 = RegEle->det_r2;
							max_testR2 = TestRegEle->det_r2;
							max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							good_occup = *Noccup;
							memcpy_for_vectorint(&new_good_arrange, arrange);
							new_GoodRegEle = *RegEle;
							new_TestGoodRegEle = *TestRegEle;
							new_GoodRegEle_found = true;
						}
						else if (!FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
							FullFrequentist_exists = true;
							max_R2 = RegEle->det_r2;
							max_testR2 = TestRegEle->det_r2;
							max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							good_occup = *Noccup;
							memcpy_for_vectorint(&new_good_arrange, arrange);
							new_GoodRegEle = *RegEle;
							new_TestGoodRegEle = *TestRegEle;
							new_GoodRegEle_found = true;
						}
						else if (FullFrequentist_exists && RegEle->full_Frequentist && RegEle->det_r2 > max_R2 && TestRegEle->det_r2 > max_testR2 && TestRegEle->det_r2 / RegEle->det_r2 > max_rR2 && RegEle->det_r2 >= 0.0 && RegEle->det_r2 <= 1.0) {
							max_R2 = RegEle->det_r2;
							max_testR2 = TestRegEle->det_r2;
							max_rR2 = TestRegEle->det_r2 / RegEle->det_r2;
							good_occup = *Noccup;
							memcpy_for_vectorint(&new_good_arrange, arrange);
							new_GoodRegEle = *RegEle;
							new_TestGoodRegEle = *TestRegEle;
							new_GoodRegEle_found = true;
						}
					}
				}

				if (useExitMode && exit_mode(macro_filename_prefix)) {
					break;
				}
				auto run_current_t = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
				long double duration = (long double)elapsed_secondsfr.count();
				if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
					std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swtich_functype.\n";
					std::cout.flush();
					break;
				}

			}

			if (new_GoodRegEle_found) {
				if (i == NNend && i > Nend) {
					NNend--;
					if (pr_id == writeoutput_prid) {
						writeoutput << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNendUpdated = " << NNend << "\n";
						writeoutput.flush();
					}
					writeoutputpr << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNendUpdated = " << NNend << "\n";
					writeoutputpr.flush();
				}
				gnew_GoodRegEle_found = true;
				good_occup = *Noccup;
				memcpy_for_vectorint(good_arrange, &new_good_arrange);
				*GoodRegEle = new_GoodRegEle;
				*TestGoodRegEle = new_TestGoodRegEle;
				if (pr_id == writeoutput_prid && Verbose == 1) {
					write_arrange(&writeoutput, good_arrange);
					if (*betachannel == "") {
						writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutput.flush();
					}
					else {
						writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutput.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutput << "\tNF";
						writeoutput.flush();
					}
					writeoutput << "\n";
					writeoutput.flush();
				}
				if (Verbose == 1) {
					write_arrange(&writeoutputpr, good_arrange);
					if (*betachannel == "") {
						writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
						writeoutputpr.flush();
					}
					else {
						writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
						writeoutputpr.flush();
					}
					if (!FullFrequentist_exists) {
						writeoutputpr << "\tNF";
						writeoutputpr.flush();
					}
					writeoutputpr << "\n";
					writeoutputpr.flush();
				}
				if (new_good_arrange.size() > 0) {
					std::vector<int>().swap(new_good_arrange);
				}
			}
			else {
				if (!chance_added && i > Nend) {
					NNend = i - 1;
					chance_added = false;
					if (pr_id == writeoutput_prid) {
						writeoutput << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNendUpdated = " << NNend << "\n";
						writeoutput.flush();
					}
					writeoutputpr << "switch_functype: useSomeShotMode: " << (signed)from_list.size() - 1 << " --> NNendUpdated = " << NNend << "\n";
					writeoutputpr.flush();
				}
			}

		}
		
		if (useExitMode && exit_mode(macro_filename_prefix)) {
			break;
		}
		auto run_current_t = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
		long double duration = (long double)elapsed_secondsfr.count();
		if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
			std::cout << "[WARNING] pr_id = " << pr_id << " : over TargetCPUTimeAfterSwap during swtich_functype.\n";
			std::cout.flush();
			break;
		}

	}

	if (*betachannel == "") {
		ExistingCholskeySeries = false;
		do_regression(FuncType, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, for_test);
		bool testonly = true;
		do_regression(FuncType, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, GoodRegEle, &testonly);
	}
	else {
		std::string ibetachannel = "i" + *betachannel;
		bool BFGS = false;
		do_regressionB(FuncType, betachannel, &ibetachannel, GoodRegEle, dataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, for_test);
		bool testonly = true;
		do_regressionB(FuncType, betachannel, &ibetachannel, TestGoodRegEle, testdataset, good_arrange, &good_occup, CI_level, critic_p, &BFGS, GoodRegEle, max_montior, &testonly);
	}

	if (gnew_GoodRegEle_found && pr_id == writeoutput_prid) {
		write_arrange(&writeoutput, good_arrange);
		if (*betachannel == "") {
			writeoutput << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
			writeoutput.flush();
		}
		else {
			writeoutput << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
			writeoutput.flush();
		}
		if (!FullFrequentist_exists) {
			writeoutput << "\tNF";
			writeoutput.flush();
		}
		writeoutput << "\t***SF***\n";
		writeoutput.flush();
	}

	if (gnew_GoodRegEle_found) {
		write_arrange(&writeoutputpr, good_arrange);
		if (*betachannel == "") {
			writeoutputpr << "\t" << *Noccup << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_radj2 << "\t" << TestGoodRegEle->det_radj2 << "\t" << GoodRegEle->df << "\t" << GoodRegEle->p_Fvalue;
			writeoutputpr.flush();
		}
		else {
			writeoutputpr << "\t" << *Noccup << "\t" << *betachannel << "\t" << GoodRegEle->det_r2 << "\t" << TestGoodRegEle->det_r2 << "\t" << GoodRegEle->det_r2_pseudo << "\t" << TestGoodRegEle->det_r2_pseudo << "\t" << GoodRegEle->df;
			writeoutputpr.flush();
		}
		if (!FullFrequentist_exists) {
			writeoutputpr << "\tNF";
			writeoutputpr.flush();
			failure_count++;
		}
		writeoutputpr << "\t***SF***\n";
		writeoutputpr.flush();
	}

	return gnew_GoodRegEle_found;

}

void do_special_regression(REGELE* sRegEle, std::vector<DATASET>* sdataset, std::vector<FUNCTYPE>* pfunctype, REGELE* RegEle, std::vector<DATASET>* dataset, int* power, long double* CI_level, long double* critic_p, REGELE* trainedRegEle, bool* for_test) {

	ExistingCholskeySeries = false;
	std::vector<DATASET>().swap(*sdataset);
	for (int i = 0; i < (signed)dataset->size(); i++) {
		DATASET objDATASET;
		objDATASET.v.push_back((RegEle->MatrixYB(i, 0) - score_intercept) / score_slope);
		for (int j = 1; j <= *power; j++) {
			objDATASET.v.push_back(powl((RegEle->MatrixYexpectedB(i, 0) - score_intercept) / score_slope, (long double)j));
		}
		sdataset->push_back(objDATASET);
	}

	bool for_test_for_set_MatrixY = false;
	bool reset_slope_intercept_for_scores = false;
	set_MatrixY(sRegEle, sRegEle, sdataset, sdataset, &for_test_for_set_MatrixY, &reset_slope_intercept_for_scores);

	std::vector<int> sarrange;
	for (int si = 0; si < *power; si++) {
		sarrange.push_back(1);
	}
	if (!*for_test) {
		do_regression(pfunctype, sRegEle, sdataset, &sarrange, power, CI_level, critic_p, sRegEle, for_test);
	}
	else {
		do_regression(pfunctype, sRegEle, sdataset, &sarrange, power, CI_level, critic_p, trainedRegEle, for_test);
	}

	std::vector<int>().swap(sarrange);

}

bool write_regele(std::string* filename, REGELE* RegEle, REGELE* TestRegEle, std::vector<int>* datadivisionflag, std::vector<FUNCTYPE>* functype, std::vector<FUNCTYPE>* pfunctype, std::vector<int>* arrange, int* Noccup, std::vector<GROUPX>* GroupX, std::string* setY, std::vector<DATASET>* dataset, std::vector<DATASET>* testdataset, long double* CI_level, long double* critic_p, std::string* betachannel, bool* use_dataname, std::string* userdefinedName, int* max_montior, int* NTestRobustness, int* Verbose) {

	ExistingCholskeySeries = false;
	writeregele.open(*filename, std::ios::binary);
	writeregele << std::fixed << std::setprecision(OutputPrecision);

	std::vector<int> pos;
	std::vector<int> funct;
	pos.push_back(-1);
	funct.push_back(0);

	writeregele << "ModelDF\t" << *Noccup;
	writeregele.flush();
	if (*betachannel != "") {
		writeregele << "\tReg" << *setY << "\t=\t" << (1.0 / score_slope) << "\t*\t" << *betachannel << "_inverse[f(x)]\t+\t" << (-score_intercept / score_slope);
		writeregele.flush();
	}
	writeregele << "\n";
	writeregele.flush();
	writeregele << "Model\t";
	writeregele.flush();
	for (int iarr = 0; iarr < (signed)arrange->size(); iarr++) {
		writeregele << (*arrange)[iarr];
		writeregele.flush();
		if (iarr != (signed)arrange->size() - 1) {
			writeregele << "-";
			writeregele.flush();
		}
	}
	writeregele << "\t";
	writeregele.flush();
	for (int iarr = 0; iarr < (signed)arrange->size(); iarr++) {
		writeregele << (*arrange)[iarr];
		writeregele.flush();
		if ((*arrange)[iarr] != 0) {
			pos.push_back(iarr);
			funct.push_back((*arrange)[iarr]);
		}
		if (iarr != (signed)arrange->size() - 1) {
			writeregele << "\t";
			writeregele.flush();
		}
		else {
			writeregele << "\n";
			writeregele.flush();
		}
	}

	writeregele << "SStot\t" << RegEle->SStot;
	writeregele.flush();
	if (*betachannel != "") {
		writeregele << "\tSStotB\t" << RegEle->SStotB;
		writeregele.flush();
	}
	writeregele << "\n";
	writeregele.flush();
	writeregele << "SSres\t" << RegEle->SSres;
	writeregele.flush();
	if (*betachannel != "") {
		writeregele << "\tSSresB\t" << RegEle->SSresB;
		writeregele.flush();
	}
	writeregele << "\n";
	writeregele.flush();
	writeregele << "det_r2\t" << RegEle->det_r2 << "\n";
	writeregele.flush();
	if (*betachannel == "") {
		writeregele << "det_radj2\t" << RegEle->det_radj2 << "\n";
		writeregele.flush();
	}
	else {
		int power = 1;
		bool for_test = false;
		do_special_regression(&specialRegEle, &specialdataset, pfunctype, RegEle, dataset, &power, CI_level, critic_p, &specialRegEle, &for_test);
		power = 2;
		do_special_regression(&specialRegEle_quad, &specialdataset_quad, pfunctype, RegEle, dataset, &power, CI_level, critic_p, &specialRegEle_quad, &for_test);
		writeregele << "det_r2_pseudo\t" << RegEle->det_r2_pseudo << "\n";
		writeregele.flush();
		writeregele << "det_r2_revised\t" << specialRegEle.det_r2 << "\t:\tRevisedReg" << *setY << "\t=\t" << specialRegEle.fit_coeff[0] << "\t+\t" << specialRegEle.fit_coeff[1] << "\tReg" << *setY << "\n";
		writeregele.flush();
		writeregele << "det_r2_revised_quad\t" << specialRegEle_quad.det_r2 << "\t:\tQuadRevisedReg" << *setY << "\t=\t" << specialRegEle_quad.fit_coeff[0] << "\t+\t" << specialRegEle_quad.fit_coeff[1] << "\tReg" << *setY << "\t+\t" << specialRegEle_quad.fit_coeff[2] << "\tReg" << *setY << "^2\n";
		writeregele.flush();
	}
	writeregele << "MSE\t" << RegEle->MSE << "\n";
	writeregele.flush();
	writeregele << "RMSE\t" << RegEle->RMSE << "\n";
	writeregele.flush();
	writeregele << "MAE\t" << RegEle->MAE << "\n";
	writeregele.flush();
	if (*betachannel != "") {
		writeregele << "MSE_revised\t" << specialRegEle.MSE << "\n";
		writeregele.flush();
		writeregele << "RMSE_revised\t" << specialRegEle.RMSE << "\n";
		writeregele.flush();
		writeregele << "MAE_revised\t" << specialRegEle.MAE << "\n";
		writeregele.flush();
		writeregele << "MSE_revised_quad\t" << specialRegEle_quad.MSE << "\n";
		writeregele.flush();
		writeregele << "RMSE_revised_quad\t" << specialRegEle_quad.RMSE << "\n";
		writeregele.flush();
		writeregele << "MAE_revised_quad\t" << specialRegEle_quad.MAE << "\n";
		writeregele.flush();
	}
	writeregele << "ErrorDF\t" << RegEle->df << "\n";
	writeregele.flush();
	if (*betachannel == "") {
		writeregele << "Fvalue\t" << RegEle->Fvalue << "\n";
		writeregele.flush();
		writeregele << "p_Fvalue\t" << RegEle->p_Fvalue << "\n";
		writeregele.flush();
	}
	writeregele << "full_Frequentist\t";
	writeregele.flush();
	if (RegEle->full_Frequentist) {
		writeregele << "1\n";
		writeregele.flush();
	}
	else {
		writeregele << "0\n";
		writeregele.flush();
	}
	writeregele << "fit_variable\tfit_coeff\tfit_coeff_standardized\tfit_coeff_SE\tfit_coeff_t\tfit_coeff_p\tfit_coeff_CI_start\tfit_coeff_CI_end\tVIF\n";
	writeregele.flush();
	std::vector<std::string> variable_name_list;

	int coeff_size;
	if (*betachannel == "") {
		coeff_size = RegEle->fit_coeff.size();
	}
	else {
		coeff_size = RegEle->fit_coeffB.size();
	}
	long double* VIF;
	VIF = new long double[coeff_size];
	std::string* variable_name;
	variable_name = new std::string[coeff_size];

	for (unsigned int i = 0; i < coeff_size; i++) {
		if (i == 0) {
			variable_name[i] = "intercept";
			VIF[i] = 0.0;
		}
		else {
			if ((*functype)[funct[i] - 1].userdefined_name.size() == 1) {
				variable_name[i] = (*functype)[funct[i] - 1].userdefined_name[0] + (*GroupX)[pos[i]].setX;
			}
			else if ((*functype)[funct[i] - 1].userdefined_name.size() > 1) {
				variable_name[i] = (*functype)[funct[i] - 1].userdefined_name[0] + (*GroupX)[pos[i]].setX + "^" + (*functype)[funct[i] - 1].userdefined_name[1];
			}
			else {
				variable_name[i] = "F?";
			}
			variable_name_list.push_back(variable_name[i]);
			REGELE newRegEle;
			std::vector<DATASET> newdataset;
			std::vector<int> newarrange;
			int newNoccup = *Noccup - 1;
			for (int newi = 0; newi < newNoccup; newi++) {
				newarrange.push_back(1);
			}
			int data_size = (signed)dataset->size();
			transpose_one_column_to_y(&newRegEle, &newdataset, RegEle, &data_size, &newNoccup, &i);
			bool for_test = false;
			do_regression(pfunctype, &newRegEle, &newdataset, &newarrange, &newNoccup, CI_level, critic_p, &newRegEle, &for_test);
			VIF[i] = 1.0 / (1.0 - newRegEle.det_r2);
		}

	}

	int* order_z;
	int* order_z_except_for_intercept;
	bool* already_check;
	order_z = new int[coeff_size];
	order_z_except_for_intercept = new int[coeff_size - 1];
	already_check = new bool[coeff_size];
	for (unsigned int i = 0; i < coeff_size; i++) {
		order_z[i] = -1;
		already_check[i] = false;
	}
	int intercept_target = -1;
	for (unsigned int i = 0; i < coeff_size; i++) {
		long double find_max_abs = -10000000000000000000000.0;
		int target_index = -1;
		for (unsigned int j = 0; j < coeff_size; j++) {
			if (*betachannel == "") {
				if (!already_check[j] && fabsl(RegEle->fit_coeff_standardized[j]) > find_max_abs) {
					find_max_abs = fabsl(RegEle->fit_coeff_standardized[j]);
					target_index = j;
				}
			}
			else {
				if (!already_check[j] && fabsl(RegEle->fit_coeff_standardizedB[j]) > find_max_abs) {
					find_max_abs = fabsl(RegEle->fit_coeff_standardizedB[j]);
					target_index = j;
				}
			}
		}
		if (target_index == -1) {
			for (unsigned int j = 0; j < coeff_size; j++) {
				if (!already_check[j]) {
					target_index = j;
				}
			}
		}
		order_z[i] = target_index;
		if (target_index == 0) {
			intercept_target = i;
		}
		already_check[target_index] = true;
	}
	if (order_z[0] == -1) {
		writeregele.close();
		delete[] VIF;
		VIF = nullptr;
		delete[] variable_name;
		std::vector<int>().swap(pos);
		std::vector<int>().swap(funct);
		std::vector<std::string>().swap(variable_name_list);
		delete[] order_z;
		order_z = nullptr;
		return false;
	}

	for (unsigned int i = 0; i < coeff_size; i++) {
		if (i < intercept_target) {
			order_z_except_for_intercept[i] = order_z[i] - 1;
		}
		else if (i > intercept_target) {
			order_z_except_for_intercept[i - 1] = order_z[i] - 1;
		}
	}
	for (unsigned int i = 0; i < coeff_size; i++) {
		if (*betachannel == "") {
			writeregele << variable_name[order_z[i]] << "\t" << RegEle->fit_coeff[order_z[i]] << "\t" << RegEle->fit_coeff_standardized[order_z[i]] << "\t" << RegEle->fit_coeff_SE[order_z[i]] << "\t" << RegEle->fit_coeff_t[order_z[i]] << "\t" << RegEle->fit_coeff_p[order_z[i]] << "\t" << RegEle->fit_coeff_CI_start[order_z[i]] << "\t" << RegEle->fit_coeff_CI_end[order_z[i]] << "\t" << VIF[order_z[i]];
			writeregele.flush();
		}
		else {
			writeregele << variable_name[order_z[i]] << "\t" << RegEle->fit_coeffB[order_z[i]] << "\t" << RegEle->fit_coeff_standardizedB[order_z[i]] << "\t" << RegEle->fit_coeff_SEB[order_z[i]] << "\t" << RegEle->fit_coeff_tB[order_z[i]] << "\t" << RegEle->fit_coeff_pB[order_z[i]] << "\t" << RegEle->fit_coeff_CI_startB[order_z[i]] << "\t" << RegEle->fit_coeff_CI_endB[order_z[i]] << "\t" << VIF[order_z[i]];
			writeregele.flush();
		}
		if (i != coeff_size - 1) {
			writeregele << "\n";
			writeregele.flush();
		}
	}
	delete[] VIF;
	VIF = nullptr;
	delete[] variable_name;
	delete[] already_check;
	already_check = nullptr;

	std::vector<long double> y_regressed, y_raw;
	std::vector<long double> y_testregressed, y_testraw;

	std::vector<long double> ObsEtas;
	long double ObsEta_absmax = -100000000.0;
	std::string betachanneli;
	if (*betachannel != "") {
		betachanneli = "i" + *betachannel;
		for (int i = 0; i < (signed)dataset->size(); i++) {
			long double v = RegEle->MatrixYexpectedB(i, 0);
			long double gv = get_v_special(betachannel, &v);
			ObsEtas.push_back(gv);
			if (fabsl(gv) > ObsEta_absmax) {
				ObsEta_absmax = fabsl(gv);
			}
		}
	}

	if (*betachannel == "") {
		writeregele << "\ndataset.size()\t" << (signed)dataset->size() << "\n";
		writeregele.flush();
		if (*use_dataname) {
			writeregele << *userdefinedName << "\t";
			writeregele.flush();
		}
		writeregele << "Reg" << *setY << "\t" << *setY << "\t";
		writeregele.flush();
	}
	else {
		writeregele << "\ndataset.size()\t" << (signed)dataset->size() << "\nScore=\t" << score_slope << "\t" << *setY << "\t+\t" << score_intercept << "\nEta\tMyu\t";
		writeregele.flush();
		if (*use_dataname) {
			writeregele << *userdefinedName << "\t";
			writeregele.flush();
		}
		writeregele << "ObsEta\tRegScore" << *setY << "\tScore" << *setY << "\tReg" << *setY << "\tRevisedReg" << *setY << "\tQuadRevisedReg" << *setY << "\t" << *setY << "\t";
		writeregele.flush();
	}
	for (unsigned int i = 0; i < variable_name_list.size(); i++) {
		writeregele << variable_name_list[order_z_except_for_intercept[i]];
		writeregele.flush();
		if (i != variable_name_list.size() - 1) {
			writeregele << "\t";
			writeregele.flush();
		}
		else {
			writeregele << "\n";
			writeregele.flush();
		}
	}
	ObsEta_absmax = 10.0;
	long double eta_mesh = 0.1;
	long double eta_pos = -ObsEta_absmax;
	for (int i = 0; i < (signed)dataset->size(); i++) {
		if (*betachannel == "") {
			if (*use_dataname) {
				writeregele << (*dataset)[i].name << "\t";
				writeregele.flush();
			}
			writeregele << RegEle->MatrixYexpected(i, 0) << "\t" << RegEle->MatrixY(i, 0) << "\t";
			writeregele.flush();
			y_regressed.push_back(RegEle->MatrixYexpected(i, 0));
			y_raw.push_back(RegEle->MatrixY(i, 0));
		}
		else {
			writeregele << eta_pos << "\t" << get_v_special(&betachanneli, &eta_pos) << "\t";
			writeregele.flush();
			if (*use_dataname) {
				writeregele << (*dataset)[i].name << "\t";
				writeregele.flush();
			}
			writeregele << ObsEtas[i] << "\t" << RegEle->MatrixYexpectedB(i, 0) << "\t" << RegEle->MatrixYB(i, 0) << "\t" << (RegEle->MatrixYexpectedB(i, 0) - score_intercept) / score_slope << "\t" << specialRegEle.MatrixYexpected(i, 0) << "\t" << specialRegEle_quad.MatrixYexpected(i, 0) << "\t" << (RegEle->MatrixYB(i, 0) - score_intercept) / score_slope << "\t";
			writeregele.flush();
			y_regressed.push_back(specialRegEle.MatrixYexpected(i, 0));
			y_raw.push_back((RegEle->MatrixYB(i, 0) - score_intercept) / score_slope);
			eta_pos += eta_mesh;
			if (eta_pos > ObsEta_absmax) {
				eta_pos = ObsEta_absmax;
			}
		}
		for (unsigned int j = 0; j < variable_name_list.size(); j++) {
			writeregele << RegEle->MatrixX(i, order_z_except_for_intercept[j] + 1);
			writeregele.flush();
			if (j != variable_name_list.size() - 1) {
				writeregele << "\t";
				writeregele.flush();
			}
			else if (i != (signed)dataset->size() - 1) {
				writeregele << "\n";
				writeregele.flush();
			}
		}
	}
	if (eta_pos < ObsEta_absmax && *betachannel != "") {
		writeregele << "\n";
		writeregele.flush();
		while (true) {
			writeregele << eta_pos << "\t" << get_v_special(&betachanneli, &eta_pos);
			writeregele.flush();
			eta_pos += eta_mesh;
			if (eta_pos > ObsEta_absmax) {
				break;
			}
			else {
				writeregele << "\n";
				writeregele.flush();
			}
		}
	}

	if (testdataset->size() > 0) {

		bool for_test = true;
		if (*betachannel == "") {
			do_regression(functype, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, RegEle, &for_test);
		}
		else {
			bool BFGS = false;
			std::string ibetachannel = "i" + *betachannel;
			do_regressionB(functype, betachannel, &ibetachannel, TestRegEle, testdataset, arrange, Noccup, CI_level, critic_p, &BFGS, RegEle, max_montior, &for_test);
			ObsEtas.clear();
			ObsEta_absmax = -100000000.0;
			betachanneli = "i" + *betachannel;
			for (int i = 0; i < (signed)testdataset->size(); i++) {
				long double v = TestRegEle->MatrixYexpectedB(i, 0);
				long double gv = get_v_special(betachannel, &v);
				ObsEtas.push_back(gv);
				if (fabsl(gv) > ObsEta_absmax) {
					ObsEta_absmax = fabsl(gv);
				}
			}
			ObsEta_absmax = 10.0;
			eta_mesh = 0.1;
			eta_pos = -ObsEta_absmax;
			int power = 1;
			do_special_regression(&specialTestRegEle, &specialdataset, pfunctype, TestRegEle, testdataset, &power, CI_level, critic_p, &specialRegEle, &for_test);
			power = 2;
			do_special_regression(&specialTestRegEle_quad, &specialdataset_quad, pfunctype, TestRegEle, testdataset, &power, CI_level, critic_p, &specialRegEle_quad, &for_test);
		}
		writeregele << "\n******************************TEST******************************\n";
		writeregele.flush();
		writeregele << "det_r2\t" << TestRegEle->det_r2 << "\n";
		writeregele.flush();
		writeregele << "det_radj2\t" << TestRegEle->det_radj2 << "\n";
		writeregele.flush();
		writeregele << "MSE\t" << TestRegEle->MSE << "\n";
		writeregele.flush();
		writeregele << "RMSE\t" << TestRegEle->RMSE << "\n";
		writeregele.flush();
		writeregele << "MAE\t" << TestRegEle->MAE << "\n";
		writeregele.flush();
		if (*betachannel != "") {
			writeregele << "det_r2_revised\t" << specialTestRegEle.det_r2 << "\n";
			writeregele.flush();
			writeregele << "det_radj2_revised\t" << specialTestRegEle.det_radj2 << "\n";
			writeregele.flush();
			writeregele << "MSE_revised\t" << specialTestRegEle.MSE << "\n";
			writeregele.flush();
			writeregele << "RMSE_revised\t" << specialTestRegEle.RMSE << "\n";
			writeregele.flush();
			writeregele << "MAE_revised\t" << specialTestRegEle.MAE << "\n";
			writeregele.flush();
			writeregele << "det_r2_revised_quad\t" << specialTestRegEle_quad.det_r2 << "\n";
			writeregele.flush();
			writeregele << "det_radj2_revised_quad\t" << specialTestRegEle_quad.det_radj2 << "\n";
			writeregele.flush();
			writeregele << "MSE_revised_quad\t" << specialTestRegEle_quad.MSE << "\n";
			writeregele.flush();
			writeregele << "RMSE_revised_quad\t" << specialTestRegEle_quad.RMSE << "\n";
			writeregele.flush();
			writeregele << "MAE_revised_quad\t" << specialTestRegEle_quad.MAE << "\n";
			writeregele.flush();
		}
		writeregele << "testdataset.size()\t" << (signed)testdataset->size() << "\n";
		writeregele.flush();
		if (*betachannel == "") {
			if (*use_dataname) {
				writeregele << *userdefinedName << "\t";
				writeregele.flush();
			}
			writeregele << "Reg" << *setY << "\t" << *setY << "\t";
			writeregele.flush();
		}
		else {
			writeregele << "Eta\tMyu\t";
			writeregele.flush();
			if (*use_dataname) {
				writeregele << *userdefinedName << "\t";
				writeregele.flush();
			}
			writeregele << "ObsEta\tRegScore" << *setY << "\tScore" << *setY << "\tReg" << *setY << "\tRevisedReg" << *setY << "\tQuadRevisedReg" << *setY << "\t" << *setY << "\t";
			writeregele.flush();
		}
		for (unsigned int i = 0; i < variable_name_list.size(); i++) {
			writeregele << variable_name_list[order_z_except_for_intercept[i]];
			writeregele.flush();
			if (i != variable_name_list.size() - 1) {
				writeregele << "\t";
				writeregele.flush();
			}
			else {
				writeregele << "\n";
				writeregele.flush();
			}
		}

		for (int i = 0; i < (signed)testdataset->size(); i++) {
			if (*betachannel == "") {
				if (*use_dataname) {
					writeregele << (*testdataset)[i].name << "\t";
					writeregele.flush();
				}
				writeregele << TestRegEle->MatrixYexpected(i, 0) << "\t" << TestRegEle->MatrixY(i, 0) << "\t";
				writeregele.flush();
				y_testregressed.push_back(TestRegEle->MatrixYexpected(i, 0));
				y_testraw.push_back(TestRegEle->MatrixY(i, 0));
			}
			else {
				writeregele << eta_pos << "\t" << get_v_special(&betachanneli, &eta_pos) << "\t";
				writeregele.flush();
				if (*use_dataname) {
					writeregele << (*testdataset)[i].name << "\t";
					writeregele.flush();
				}
				writeregele << ObsEtas[i] << "\t" << TestRegEle->MatrixYexpectedB(i, 0) << "\t" << TestRegEle->MatrixYB(i, 0) << "\t" << (TestRegEle->MatrixYexpectedB(i, 0) - score_intercept) / score_slope << "\t" << specialTestRegEle.MatrixYexpected(i, 0) << "\t" << specialTestRegEle_quad.MatrixYexpected(i, 0) << "\t" << (TestRegEle->MatrixYB(i, 0) - score_intercept) / score_slope << "\t";
				writeregele.flush();
				y_testregressed.push_back(specialTestRegEle.MatrixYexpected(i, 0));
				y_testraw.push_back((TestRegEle->MatrixYB(i, 0) - score_intercept) / score_slope);
				eta_pos += eta_mesh;
				if (eta_pos > ObsEta_absmax) {
					eta_pos = ObsEta_absmax;
				}
			}
			for (unsigned int j = 0; j < variable_name_list.size(); j++) {
				writeregele << TestRegEle->MatrixX(i, order_z_except_for_intercept[j] + 1);
				writeregele.flush();
				if (j != variable_name_list.size() - 1) {
					writeregele << "\t";
					writeregele.flush();
				}
				else if (i != (signed)testdataset->size() - 1) {
					writeregele << "\n";
					writeregele.flush();
				}
			}
		}


		writeregele << "\n******************************ROBUSTNESS******************************\n";
		writeregele.flush();
		long double SStot = 0.0;
		long double SSres = 0.0;
		long double MSE = 0.0;
		long double RMSE = 0.0;
		long double MAE = 0.0;
		long double ybar = 0.0;
		long double R2;
		int total_dataset = (signed)dataset->size() + (signed)testdataset->size();
		std::vector<long double> LMSE, LRMSE, LMAE, LR2;
		if (*Verbose == 1) {
			writeregele << "i\tMSE\tRMSE\tMAE\tR2\n";
			writeregele.flush();
		}
		for (int irob = 0; irob < *NTestRobustness; irob++) {
			std::vector<int> iTest;
			for (int iirob = 0; iirob < (signed)testdataset->size(); iirob++) {
				while (true) {
					int randomNumber = rand() % total_dataset;
					bool already = false;
					for (int j = 0; j < (signed)iTest.size(); j++) {
						if (iTest[j] == randomNumber) {
							already = true;
							break;
						}
					}
					if (!already) {
						iTest.push_back(randomNumber);
						break;
					}
				}
			}
			test_robustness(&SStot, &SSres, &MSE, &RMSE, &MAE, &ybar, &R2, &y_regressed, &y_raw, &y_testregressed, &y_testraw, &iTest);
			LMSE.push_back(MSE);
			LRMSE.push_back(RMSE);
			LMAE.push_back(MAE);
			LR2.push_back(R2);
			std::vector<int>().swap(iTest);
			if (*Verbose == 1) {
				writeregele << irob + 1 << "\t" << MSE << "\t" << RMSE << "\t" << MAE << "\t" << R2 << "\n";
				writeregele.flush();
			}
		}

		long double minMSE = 100000000000000.0;
		long double minRMSE = 100000000000000.0;
		long double minMAE = 100000000000000.0;
		long double minR2 = 100000000000000.0;
		long double aveMSE = 0.0;
		long double aveRMSE = 0.0;
		long double aveMAE = 0.0;
		long double aveR2 = 0.0;
		long double maxMSE = -100000000000000.0;
		long double maxRMSE = -100000000000000.0;
		long double maxMAE = -100000000000000.0;
		long double maxR2 = -100000000000000.0;
		long double stdevMSE = 0.0;
		long double stdevRMSE = 0.0;
		long double stdevMAE = 0.0;
		long double stdevR2 = 0.0;
		for (int irob = 0; irob < (signed)LMSE.size(); irob++) {
			if (LMSE[irob] < minMSE) {
				minMSE = LMSE[irob];
			}
			aveMSE += LMSE[irob];
			if (LMSE[irob] > maxMSE) {
				maxMSE = LMSE[irob];
			}
			if (LRMSE[irob] < minRMSE) {
				minRMSE = LRMSE[irob];
			}
			aveRMSE += LRMSE[irob];
			if (LRMSE[irob] > maxRMSE) {
				maxRMSE = LRMSE[irob];
			}
			if (LMAE[irob] < minMAE) {
				minMAE = LMAE[irob];
			}
			aveMAE += LMAE[irob];
			if (LMAE[irob] > maxMAE) {
				maxMAE = LMAE[irob];
			}
			if (LR2[irob] < minR2) {
				minR2 = LR2[irob];
			}
			aveR2 += LR2[irob];
			if (LR2[irob] > maxR2) {
				maxR2 = LR2[irob];
			}
		}

		aveMSE /= (long double)LMSE.size();
		aveRMSE /= (long double)LRMSE.size();
		aveMAE /= (long double)LMAE.size();
		aveR2 /= (long double)LR2.size();

		for (int irob = 0; irob < (signed)LMSE.size(); irob++) {
			stdevMSE += powl(LMSE[irob] - aveMSE, 2.0);
			stdevRMSE += powl(LRMSE[irob] - aveRMSE, 2.0);
			stdevMAE += powl(LMAE[irob] - aveMAE, 2.0);
			stdevR2 += powl(LR2[irob] - aveR2, 2.0);
		}

		stdevMSE /= (long double)LMSE.size();
		stdevRMSE /= (long double)LRMSE.size();
		stdevMAE /= (long double)LMAE.size();
		stdevR2 /= (long double)LR2.size();

		stdevMSE = sqrtl(stdevMSE);
		stdevRMSE = sqrtl(stdevRMSE);
		stdevMAE = sqrtl(stdevMAE);
		stdevR2 = sqrtl(stdevR2);

		writeregele << "Tested_dataset_size\t" << testdataset->size() << "\n";
		writeregele.flush();
		writeregele << "NRandom\t" << LMSE.size() << "\n";
		writeregele.flush();
		writeregele << "Metric\tMin\tAve\tMax\tStDev\n";
		writeregele.flush();
		writeregele << "MSE\t" << minMSE << "\t" << aveMSE << "\t" << maxMSE << "\t" << stdevMSE << "\n";
		writeregele.flush();
		writeregele << "RMSE\t" << minRMSE << "\t" << aveRMSE << "\t" << maxRMSE << "\t" << stdevRMSE << "\n";
		writeregele.flush();
		writeregele << "MAE\t" << minMAE << "\t" << aveMAE << "\t" << maxMAE << "\t" << stdevMAE << "\n";
		writeregele.flush();
		writeregele << "R2\t" << minR2 << "\t" << aveR2 << "\t" << maxR2 << "\t" << stdevR2;
		writeregele.flush();

		std::vector<long double>().swap(LMSE);
		std::vector<long double>().swap(LRMSE);
		std::vector<long double>().swap(LMAE);
		std::vector<long double>().swap(LR2);

	}

	if (variable_name_list.size() > 1) {

		std::vector<bool> divisible;
		for (unsigned int j2 = 0; j2 < variable_name_list.size(); j2++) {
			bool divisible_ele = true;
			for (int idat = 0; idat < (signed)datadivisionflag->size(); idat++) {
				if ((*datadivisionflag)[idat] > 0) {
					int dataindex = (*datadivisionflag)[idat] - 1;
					if (RegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1) == 0.0 || RegEle->MatrixX(0, order_z_except_for_intercept[j2] + 1) * RegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1) < 0) {
						divisible_ele = false;
						break;
					}
				}
				else {
					int dataindex = -(*datadivisionflag)[idat] - 1;
					if (TestRegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1) == 0.0 || RegEle->MatrixX(0, order_z_except_for_intercept[j2] + 1) * TestRegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1) < 0) {
						divisible_ele = false;
						break;
					}
				}
			}
			divisible.push_back(divisible_ele);
		}

		writeregele << "\n******************************INTERACTION******************************\n";
		writeregele.flush();
		if (*use_dataname) {
			writeregele << *userdefinedName << "\t";
			writeregele.flush();
		}
		writeregele << *setY << "\t";
		writeregele.flush();
		for (unsigned int i = 0; i < variable_name_list.size(); i++) {
			writeregele << "(" << variable_name_list[order_z_except_for_intercept[i]] << ")\t";
			writeregele.flush();
		}
		for (unsigned int i1 = 0; i1 < variable_name_list.size(); i1++) {
			for (unsigned int i2 = i1 + 1; i2 < variable_name_list.size(); i2++) {
				writeregele << "(" << variable_name_list[order_z_except_for_intercept[i1]] << "*" << variable_name_list[order_z_except_for_intercept[i2]] << ")";
				writeregele.flush();
				if (divisible[i1] && divisible[i2]) {
					writeregele << "\t(" << variable_name_list[order_z_except_for_intercept[i1]] << "/" << variable_name_list[order_z_except_for_intercept[i2]] << ")";
					writeregele.flush();
				}
				if (!(i1 == variable_name_list.size() - 2 && i2 == variable_name_list.size() - 1)) {
					writeregele << "\t";
					writeregele.flush();
				}
				else {
					writeregele << "\n";
					writeregele.flush();
				}
			}
		}

		for (int idat = 0; idat < (signed)datadivisionflag->size(); idat++) {
			if ((*datadivisionflag)[idat] > 0) {
				int dataindex = (*datadivisionflag)[idat] - 1;
				if (*use_dataname) {
					writeregele << (*dataset)[dataindex].name << "\t";
					writeregele.flush();
				}
				writeregele << RegEle->MatrixY(dataindex, 0) << "\t";
				writeregele.flush();
				for (unsigned int j = 0; j < variable_name_list.size(); j++) {
					writeregele << RegEle->MatrixX(dataindex, order_z_except_for_intercept[j] + 1) << "\t";
					writeregele.flush();
				}
				for (unsigned int j1 = 0; j1 < variable_name_list.size(); j1++) {
					for (unsigned int j2 = j1 + 1; j2 < variable_name_list.size(); j2++) {
						writeregele << RegEle->MatrixX(dataindex, order_z_except_for_intercept[j1] + 1) * RegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1);
						writeregele.flush();
						if (divisible[j1] && divisible[j2]) {
							writeregele << "\t" << RegEle->MatrixX(dataindex, order_z_except_for_intercept[j1] + 1) / RegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1);
							writeregele.flush();
						}
						if (!(j1 == variable_name_list.size() - 2 && j2 == variable_name_list.size() - 1)) {
							writeregele << "\t";
							writeregele.flush();
						}
						else {
							if (idat != (signed)datadivisionflag->size() - 1) {
								writeregele << "\n";
								writeregele.flush();
							}
						}
					}
				}
			}
			else {
				int dataindex = -(*datadivisionflag)[idat] - 1;
				if (*use_dataname) {
					writeregele << (*testdataset)[dataindex].name << "\t";
					writeregele.flush();
				}
				writeregele << TestRegEle->MatrixY(dataindex, 0) << "\t";
				writeregele.flush();
				for (unsigned int j = 0; j < variable_name_list.size(); j++) {
					writeregele << TestRegEle->MatrixX(dataindex, order_z_except_for_intercept[j] + 1) << "\t";
					writeregele.flush();
				}
				for (unsigned int j1 = 0; j1 < variable_name_list.size(); j1++) {
					for (unsigned int j2 = j1 + 1; j2 < variable_name_list.size(); j2++) {
						writeregele << TestRegEle->MatrixX(dataindex, order_z_except_for_intercept[j1] + 1) * TestRegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1);
						writeregele.flush();
						if (divisible[j1] && divisible[j2]) {
							writeregele << "\t" << TestRegEle->MatrixX(dataindex, order_z_except_for_intercept[j1] + 1) / TestRegEle->MatrixX(dataindex, order_z_except_for_intercept[j2] + 1);
							writeregele.flush();
						}
						if (!(j1 == variable_name_list.size() - 2 && j2 == variable_name_list.size() - 1)) {
							writeregele << "\t";
							writeregele.flush();
						}
						else {
							if (idat != (signed)datadivisionflag->size() - 1) {
								writeregele << "\n";
								writeregele.flush();
							}
						}
					}
				}
			}
		}

		std::vector<bool>().swap(divisible);

	}

	writeregele.close();

	std::vector<int>().swap(pos);
	std::vector<int>().swap(funct);
	std::vector<std::string>().swap(variable_name_list);
	std::vector<long double>().swap(ObsEtas);
	std::vector<long double>().swap(y_regressed);
	std::vector<long double>().swap(y_raw);
	std::vector<long double>().swap(y_testregressed);
	std::vector<long double>().swap(y_testraw);

	delete[] order_z;
	order_z = nullptr;
	delete[] order_z_except_for_intercept;
	order_z_except_for_intercept = nullptr;
	return true;

}

long double get_r2_by_reading(std::string* filename) {

	readregele.open(*filename);

	while (true) {
		readregele >> sbuf;
		if (sbuf == "det_r2" || readregele.eof()) {
			break;
		}
	}
	long double ldbuf;
	readregele >> ldbuf;
	readregele.close();

	return(ldbuf);

}

long double get_tr2_by_reading(std::string* filename) {

	readregele.open(*filename);

	int count = 0;
	while (true) {
		readregele >> sbuf;
		if (sbuf == "det_r2") {
			count++;
		}
		if (count == 2 || readregele.eof()) {
			break;
		}
	}
	long double ldbuf;
	readregele >> ldbuf;
	readregele.close();

	return(ldbuf);

}

long double get_aver2_by_reading(std::string* filename) {

	readregele.open(*filename);
	
	while (true) {
		readregele >> sbuf;
		if (sbuf == "R2" || readregele.eof()) {
			break;
		}
	}
	long double ldbuf;
	readregele >> ldbuf >> ldbuf;
	readregele.close();

	return(ldbuf);

}

bool get_full_Frequentist_by_reading(std::string* filename) {

	readregele.open(*filename);

	while (true) {
		readregele >> sbuf;
		if (sbuf == "full_Frequentist") {
			break;
		}
	}
	int ibuf;
	readregele >> ibuf;

	readregele.close();

	if (ibuf == 0) {
		return false;
	}
	else {
		return true;
	}

}

void rewrite_regele(std::ofstream* to_stream, std::string* from_filename, std::string* setY, std::string* betatype, std::string* macro_filename_prefix, bool* failed_case) {

	readregele.open(*from_filename);

	while (true) {
		std::getline(readregele, sbuf);
		*to_stream << sbuf << "\n";
		to_stream->flush();
		if (readregele.eof()) {
			break;
		}
	}

	readregele.close();

	readregele.open(*from_filename);

	int ModelDF = 0;
	long double R2 = 0.0;
	long double tR2 = 0.0;
	long double MSE = 0.0;
	long double tMSE = 0.0;
	long double RMSE = 0.0;
	long double tRMSE = 0.0;
	long double MAE = 0.0;
	long double tMAE = 0.0;
	long double aveR2 = 0.0;
	int cR2 = 0;
	int ctR2 = 0;
	int cMSE = 0;
	int ctMSE = 0;
	int cRMSE = 0;
	int ctRMSE = 0;
	int cMAE = 0;
	int ctMAE = 0;
	std::vector<std::string> descriptors;
	bool finalized = false;
	bool process_r2 = true;
	bool process_MSE = true;
	bool process_RMSE = true;
	bool process_MAE = true;

	while (!*failed_case) {

		readregele >> sbuf;

		if (sbuf == "ModelDF") {
			readregele >> ModelDF;
		}

		if (process_r2 && *betatype == "" && cR2 == 0 && sbuf == "det_r2") {
			readregele >> R2;
			cR2++;
		}
		else if (process_r2 && *betatype != "" && cR2 == 0 && sbuf == "det_r2_revised") {
			readregele >> R2;
			cR2++;
		}
		else if (process_r2 && *betatype == "" && cR2 == 1 && sbuf == "det_r2") {
			readregele >> tR2;
			process_r2 = false;
		}
		else if (process_r2 && *betatype != "" && cR2 == 1 && sbuf == "det_r2_revised") {
			readregele >> tR2;
			process_r2 = false;
		}

		if (process_MSE && *betatype == "" && cMSE == 0 && sbuf == "MSE") {
			readregele >> MSE;
			cMSE++;
		}
		else if (process_MSE && *betatype != "" && cMSE == 0 && sbuf == "MSE_revised") {
			readregele >> MSE;
			cMSE++;
		}
		else if (process_MSE && *betatype == "" && cMSE == 1 && sbuf == "MSE") {
			readregele >> tMSE;
			process_MSE = false;
		}
		else if (process_MSE && *betatype != "" && cMSE == 1 && sbuf == "MSE_revised") {
			readregele >> tMSE;
			process_MSE = false;
		}

		if (process_RMSE && *betatype == "" && cRMSE == 0 && sbuf == "RMSE") {
			readregele >> RMSE;
			cRMSE++;
		}
		else if (process_RMSE && *betatype != "" && cRMSE == 0 && sbuf == "RMSE_revised") {
			readregele >> RMSE;
			cRMSE++;
		}
		else if (process_RMSE && *betatype == "" && cRMSE == 1 && sbuf == "RMSE") {
			readregele >> tRMSE;
			process_RMSE = false;
		}
		else if (process_RMSE && *betatype != "" && cRMSE == 1 && sbuf == "RMSE_revised") {
			readregele >> tRMSE;
			process_RMSE = false;
		}

		if (process_MAE && *betatype == "" && cMAE == 0 && sbuf == "MAE") {
			readregele >> MAE;
			cMAE++;
		}
		else if (process_MAE && *betatype != "" && cMAE == 0 && sbuf == "MAE_revised") {
			readregele >> MAE;
			cMAE++;
		}
		else if (process_MAE && *betatype == "" && cMAE == 1 && sbuf == "MAE") {
			readregele >> tMAE;
			process_MAE = false;
			//finalized = true;
		}
		else if (process_MAE && *betatype != "" && cMAE == 1 && sbuf == "MAE_revised") {
			readregele >> tMAE;
			process_MAE = false;
			//finalized = true;
		}

		if (sbuf == "R2") {
			readregele >> aveR2 >> aveR2;
			finalized = true;
		}

		if (sbuf == "fit_variable") {
			std::getline(readregele, sbuf);
			while (true) {
				readregele >> sbuf;
				if (sbuf != "intercept" && sbuf != "dataset.size()") {
					descriptors.push_back(sbuf);
					std::getline(readregele, sbuf);
				}
				else if (sbuf == "intercept") {
					std::getline(readregele, sbuf);
				}
				else {
					break;
				}
			}
		}

		if (readregele.eof() || finalized) {
			break;
		}

	}

	writereport << *setY << "\t" << ModelDF << "\t";
	writereport.flush();
	if (*betatype == "") {
		writereport << "NULL\t";
		writereport.flush();
	}
	else {
		writereport << *betatype << "\t";
		writereport.flush();
	}
	writereport << *macro_filename_prefix << ".txt\t" << aveR2 << "\t" << R2 << "\t" << tR2 << "\t" << MSE << "\t" << tMSE << "\t" << RMSE << "\t" << tRMSE << "\t" << MAE << "\t" << tMAE << "\t";
	writereport.flush();

	for (int i = 0; i < (signed)descriptors.size(); i++) {
		writereport << descriptors[i];
		writereport.flush();
		if (i != (signed)descriptors.size() - 1) {
			writereport << "\t";
			writereport.flush();
		}
		else {
			writereport << "\n";
			writereport.flush();
		}
	}
	if (*failed_case) {
		writereport << "\n";
		writereport.flush();
	}

	std::vector<std::string>().swap(descriptors);
	readregele.close();

}

void search_globally(std::string* setY, std::string* betatype, std::string* macro_filename_prefix) {

	std::string target_filename;
	FullFrequentist_exists = false;
	max_R2 = -1000000000.0;
	max_testR2 = -1000000000.0;
	max_rR2 = -1000000000.0;
	long double max_aveR2 = -1000000000.0;
	bool found = false;
	for (int iprid = 0; iprid < pr_size; iprid++) {
		std::string filename = *macro_filename_prefix + "." + std::to_string(iprid) + ".reg";
		long double r2 = get_r2_by_reading(&filename);
		long double tr2 = get_tr2_by_reading(&filename);
		long double aver2 = get_aver2_by_reading(&filename);
		bool FF = get_full_Frequentist_by_reading(&filename);
		if (!FullFrequentist_exists && !FF && r2 > max_R2 && tr2 > max_testR2 && tr2 / r2 > max_rR2 && r2 >= 0.0 && r2 <= 1.0 && tr2 >= 0.0 && tr2 <= 1.0) {
			max_R2 = r2;
			max_testR2 = tr2;
			max_rR2 = tr2 / r2;
			max_aveR2 = aver2;
			target_filename = filename;
			found = true;
		}
		else if (!FullFrequentist_exists && FF && r2 >= 0.0 && r2 <= 1.0 && tr2 >= 0.0 && tr2 <= 1.0) {
			FullFrequentist_exists = true;
			max_R2 = r2;
			max_testR2 = tr2;
			max_rR2 = tr2 / r2;
			max_aveR2 = aver2;
			target_filename = filename;
			found = true;
		}
		else if (FullFrequentist_exists && FF && aver2 >= max_aveR2 && r2 >= 0.0 && r2 <= 1.0 && tr2 >= 0.0 && tr2 <= 1.0) {
			max_R2 = r2;
			max_testR2 = tr2;
			max_rR2 = tr2 / r2;
			max_aveR2 = aver2;
			target_filename = filename;
			found = true;
		}
	}
	if (found) {
		if (!FullFrequentist_exists) {
			writeoutput << "[WARNING] NF\n";
			writeoutput.flush();
		}
		bool failed_case = false;
		rewrite_regele(&writeoutput, &target_filename, setY, betatype, macro_filename_prefix, &failed_case);
	}
	else {
		std::string filename = *macro_filename_prefix + "." + std::to_string(writeoutput_prid) + ".reg";
		long double r2 = get_r2_by_reading(&filename);
		long double tr2 = get_tr2_by_reading(&filename);
		bool FF = get_full_Frequentist_by_reading(&filename);
		max_R2 = r2;
		max_testR2 = tr2;
		max_rR2 = tr2 / r2;
		target_filename = filename;
		bool failed_case = true;
		rewrite_regele(&writeoutput, &target_filename, setY, betatype, macro_filename_prefix, &failed_case);
	}

}

void remove_regele() {

	char cmd[100];
#ifdef _WIN32
	strcpy_s(cmd, "del *.reg");
	system(cmd);
#elif _WIN64
	strcpy_s(cmd, "del *.reg");
	system(cmd);
#else
	strcpy(cmd, "rm -f *.reg");
	system(cmd);
#endif

}

void remove_prtxt() {

	char cmd[100];
#ifdef _WIN32
	strcpy_s(cmd, "del *.pr.txt");
	system(cmd);
#elif _WIN64
	strcpy_s(cmd, "del *.pr.txt");
	system(cmd);
#else
	strcpy(cmd, "rm -f *.pr.txt");
	system(cmd);
#endif

}

void mkdir_smt() {

	char cmd[100];
#ifdef _WIN32
	strcpy_s(cmd, "mkdir smt");
	system(cmd);
#elif _WIN64
	strcpy_s(cmd, "mkdir smt");
	system(cmd);
#else
	strcpy(cmd, "mkdir smt");
	system(cmd);
#endif

}

void remove_smt() {

	char cmd[100];
#ifdef _WIN32
	strcpy_s(cmd, "rmdir /s /q smt");
	system(cmd);
#elif _WIN64
	strcpy_s(cmd, "rmdir /s /q smt");
	system(cmd);
#else
	strcpy(cmd, "rm -rf smt");
	system(cmd);
#endif

}

void write_finalreport_per_core(int* ibeta, std::string* setY, std::string* macro_filename_prefix) {

	if (betatype[*ibeta] != "") {
		std::string ibetachannel = "i" + betatype[*ibeta];
		bool BFGS = true;
		bool for_test_last = false;
		do_regressionB(&FuncType, &betatype[*ibeta], &ibetachannel, &GoodRegEle, &dataset, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &GoodRegEle, &MaxBFGSIterations, &for_test_last);
		bool testonly = true;
		do_regressionB(&FuncType, &betatype[*ibeta], &ibetachannel, &TestGoodRegEle, &testdataset, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &BFGS, &GoodRegEle, &MaxBFGSIterations, &testonly);
	}

	if (pr_id == writeoutput_prid) {
		write_arrange(&writeoutput, &good_arrange);
		if (betatype[*ibeta] == "") {
			writeoutput << "\t" << good_occup << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_radj2 << "\t" << TestGoodRegEle.det_radj2 << "\t" << GoodRegEle.df << "\t" << GoodRegEle.p_Fvalue;
			writeoutput.flush();
		}
		else {
			writeoutput << "\t" << good_occup << "\t" << betatype[*ibeta] << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_r2_pseudo << "\t" << TestGoodRegEle.det_r2_pseudo << "\t" << GoodRegEle.df;
			writeoutput.flush();
		}
		if (!GoodRegEle.full_Frequentist) {
			writeoutput << "\tNF";
			writeoutput.flush();
		}
		writeoutput << "\t********\n";
		writeoutput.flush();
	}
	write_arrange(&writeoutputpr, &good_arrange);
	if (betatype[*ibeta] == "") {
		writeoutputpr << "\t" << good_occup << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_radj2 << "\t" << TestGoodRegEle.det_radj2 << "\t" << GoodRegEle.df << "\t" << GoodRegEle.p_Fvalue;
		writeoutput.flush();
	}
	else {
		writeoutputpr << "\t" << good_occup << "\t" << betatype[*ibeta] << "\t" << GoodRegEle.det_r2 << "\t" << TestGoodRegEle.det_r2 << "\t" << GoodRegEle.det_r2_pseudo << "\t" << TestGoodRegEle.det_r2_pseudo << "\t" << GoodRegEle.df;
		writeoutput.flush();
	}
	if (!GoodRegEle.full_Frequentist) {
		writeoutputpr << "\tNF";
		writeoutputpr.flush();
	}
	writeoutputpr << "\t********\n";
	writeoutputpr.flush();

	std::string filename = *macro_filename_prefix + "." + std::to_string(pr_id) + ".reg";
	write_regele(&filename, &GoodRegEle, &TestRegEle, &datadivisionflag, &FuncType, &PlainFuncType, &good_arrange, &good_occup, &GroupX, setY, &dataset, &testdataset, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &betatype[*ibeta], &use_dataname, &userdefinedName, &MaxBFGSIterations, &NTestRobustness, &Verbose);

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
	int successful_for_self;
	int* gssuccessful_for_self;
	gssuccessful_for_self = new int[pr_size];

	if (pr_id == writeoutput_prid) {
		remove_regele();
		remove_prtxt();
		remove_smt();
		mkdir_smt();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	read_para();
	if (pr_id == writeoutput_prid) {
		writereport.open(OutputFilename + ".txt", std::ios::binary);
		writereport << std::fixed << std::setprecision(OutputPrecision);
		writereport << "*** In cases of linear regression, Beta will be printed as NULL.\n";
		writereport.flush();
		writereport << "*** In cases of beta regression, det_r2_revised will be printed as R2.\n";
		writereport.flush();
		writereport << "Y\tNdescriptors\tBeta\tLogFilename\t<R2>\tR2\tR2test\tMSE\tMSEtest\tRMSE\tRMSEtest\tMAE\tMAEtest\tDescriptors\n";
		writereport.flush();
	}

	if (FuncType.size() < DuplicateX) {
		DuplicateX = FuncType.size();
	}

	std::string prefilename = "";

	for (int iy = 0; iy < (signed)instally.size(); iy++) {

		std::string InitialFilename = "";
		for (int ix = 0; ix < (signed)installx.size(); ix++) {

			std::string setY = instally[iy].setY;
			long double minY = instally[iy].minY;
			long double maxY = instally[iy].maxY;
			bool limitminY = instally[iy].minlimitY;
			bool limitmaxY = instally[iy].maxlimitY;
			long double max_R2 = -1000000.0;
			long double max_testR2 = -1000000.0;
			long double max_rR2 = -1000000.0;
			NVariableLimit[0] = installx[ix].NVariableLimit[0];
			NVariableLimit[1] = installx[ix].NVariableLimit[1];
			if (NVariableLimit[1] > GroupX.size()) {
				break;
			}

			strict_limit = installx[ix].strict_limit;
			std::string postfix_x = "";
			if (strict_limit) {
				postfix_x = "L";
			}
			postfix_x += std::to_string(NVariableLimit[1]);
			std::string macro_filename_prefix = OutputFilename + "_" + setY + "_" + std::to_string(ix + 1) + postfix_x;
			if (installx[ix].interact == 1) {
				macro_filename_prefix += "i";
			}
			else if (installx[ix].interact == 2) {
				macro_filename_prefix += "I";
			}

			if (ix == 0) {
				InitialFilename = macro_filename_prefix + ".txt";
			}
			if (pr_id == writeoutput_prid) {
				writeoutput.open(macro_filename_prefix + ".txt", std::ios::binary);
				writeoutput << std::fixed << std::setprecision(OutputPrecision);
			}
			writeoutputpr.open(macro_filename_prefix + "." + std::to_string(pr_id) + ".pr.txt", std::ios::binary);
			writeoutputpr << std::fixed << std::setprecision(OutputPrecision);

			std::vector<std::string> InputDataFilename;
			InputDataFilename.clear();
			bool ReOutput = false;
			bool ReInitial = false;
			if (installx[ix].interact == 0) {
				InputDataFilename = DataFilename;
			}
			else {
				InputDataFilename.push_back(prefilename);
				ReOutput = true;
				if (installx[ix].interact == 2) {
					ReInitial = true;
				}
			}

			int cbeta = 0;
			read_data(&dataset, &testdataset, &datadivisionflag, &use_dataname, &userdefinedName, &TrainRatio, &InputDataFilename, &ReOutput, &ReInitial, &InitialFilename, &cbeta, &setY, &GroupX, &limitminY, &limitmaxY, &minY, &maxY, &column_str, &DuplicateX, &NameDuplicateHandling, &Dmin, &Dmax, &Dincrement, &subfilterx, &Verbose);
			bool for_test = true;
			bool reset_slope_intercept_for_scores = true;
			set_MatrixY(&RegEle, &TestRegEle, &dataset, &testdataset, &for_test, &reset_slope_intercept_for_scores);

			Nfeatures = GroupX.size();
			estimate_jobsize(&xpr_ld_get, &jobsize, &jobsize_ld, &Nfeatures, NVariableLimit, &Verbose);
			set_pr_jobs(&xpr_ld_get, &jobsize, &jobsize_ld, &Verbose);

			writeoutput << "\n";
			writeoutput.flush();
			writeoutputpr << "\n";
			writeoutputpr.flush();

			initialize_occp(&arrange, &occup, &Nfeatures);
			MPI_Barrier(MPI_COMM_WORLD);

			for (int ibeta = 0; ibeta < (signed)betatype.size(); ibeta++) {

				if (ibeta != 0 && installx[ix].interact != 0) {
					read_data(&dataset, &testdataset, &datadivisionflag, &use_dataname, &userdefinedName, &TrainRatio, &InputDataFilename, &ReOutput, &ReInitial, &InitialFilename, &ibeta, &setY, &GroupX, &limitminY, &limitmaxY, &minY, &maxY, &column_str, &DuplicateX, &NameDuplicateHandling, &Dmin, &Dmax, &Dincrement, &subfilterx, &Verbose);
					bool for_test = true;
					bool reset_slope_intercept_for_scores = true;
					set_MatrixY(&RegEle, &TestRegEle, &dataset, &testdataset, &for_test, &reset_slope_intercept_for_scores);
					Nfeatures = GroupX.size();
					estimate_jobsize(&xpr_ld_get, &jobsize, &jobsize_ld, &Nfeatures, NVariableLimit, &Verbose);
					set_pr_jobs(&xpr_ld_get, &jobsize, &jobsize_ld, &Verbose);
					writeoutput << "\n";
					writeoutput.flush();
					writeoutputpr << "\n";
					writeoutputpr.flush();
					initialize_occp(&arrange, &occup, &Nfeatures);
					MPI_Barrier(MPI_COMM_WORLD);
				}

				run_through_models(&xpr_ld_get, &Nfeatures, NVariableLimit, &betatype[ibeta], &FuncType, &DuplicateX);
				bool init_search = true;
				if (failure_count >= 0) {

					successful_for_self = 0;
					for (int iprid = 0; iprid < pr_size; iprid++) {
						gssuccessful_for_self[iprid] = 0;
					}
					int current_Hegel_pos = -1;
					int turn = 0;

					while (true) {

						MPI_Barrier(MPI_COMM_WORLD);
						MPI_Allgather(&successful_for_self, 1, MPI_INT, gssuccessful_for_self, 1, MPI_INT, MPI_COMM_WORLD);
						MPI_Barrier(MPI_COMM_WORLD);
						bool all_successful_for_self = true;
						bool all_yet_for_self = true;
						for (int iprid = 0; iprid < pr_size; iprid++) {
							if (gssuccessful_for_self[iprid] == 0) {
								all_successful_for_self = false;
							}
							else {
								all_yet_for_self = false;
							}
						}
						if (all_successful_for_self) {
							break;
						}
						else {
							turn++;
							std::vector<HEGEL>().swap(Hegel);
							Hegel.clear();
							for (int iprid = 0; iprid < pr_size; iprid++) {
								if (gssuccessful_for_self[iprid] == 0) {
									HEGEL objHEGEL;
									objHEGEL.master = iprid;
									objHEGEL.slave.push_back(iprid);
									Hegel.push_back(objHEGEL);
								}
							}
							int HegelPointer = 0;
							for (int iprid = 0; iprid < pr_size; iprid++) {
								if (gssuccessful_for_self[iprid] == 1) {
									Hegel[HegelPointer].slave.push_back(iprid);
									HegelPointer++;
									if (HegelPointer >= (signed)Hegel.size()) {
										HegelPointer = 0;
									}
								}
							}
							for (int iH = 0; iH < (signed)Hegel.size(); iH++) {
								for (int jH = 0; jH < (signed)Hegel[iH].slave.size(); jH++) {
									if (Hegel[iH].slave[jH] == pr_id) {
										current_Hegel_pos = iH;
									}
								}
							}
							if (pr_id == writeoutput_prid) {
								if (!all_yet_for_self) {
									writeoutput << "HEGEL @" << macro_filename_prefix << " turn= " << turn << "\n";
									writeoutput.flush();
									for (int iH = 0; iH < (signed)Hegel.size(); iH++) {
										writeoutput << "Master " << Hegel[iH].master << " ... Slave(s)";
										writeoutput.flush();
										for (int jH = 0; jH < (signed)Hegel[iH].slave.size(); jH++) {
											writeoutput << " " << Hegel[iH].slave[jH];
											writeoutput.flush();
											if (jH == (signed)Hegel[iH].slave.size() - 1) {
												writeoutput << "\n";
												writeoutput.flush();
											}
										}
									}
									writeoutputpr << "HEGEL @" << macro_filename_prefix << " turn= " << turn << "\n";
									writeoutputpr.flush();
									for (int iH = 0; iH < (signed)Hegel.size(); iH++) {
										writeoutputpr << "Master " << Hegel[iH].master << " ... Slave(s)";
										writeoutputpr.flush();
										for (int jH = 0; jH < (signed)Hegel[iH].slave.size(); jH++) {
											writeoutputpr << " " << Hegel[iH].slave[jH];
											writeoutputpr.flush();
											if (jH == (signed)Hegel[iH].slave.size() - 1) {
												writeoutputpr << "\n";
												writeoutputpr.flush();
											}
										}
									}
								}
								else {
									writeoutput << "LUHMANN-AUTOPOIESES @" << macro_filename_prefix << " turn= " << turn << "\n";
									writeoutput.flush();
									writeoutputpr << "LUHMANN-AUTOPOIESES @" << macro_filename_prefix << " turn= " << turn << "\n";
									writeoutputpr.flush();
								}
							}
						}
						MPI_Barrier(MPI_COMM_WORLD);

						bool add_updated = true;
						bool swap_updated = true;
						bool switch_updated = true;
						bool for_test = false;
						auto run_start_t = std::chrono::system_clock::now();
						if (!init_search) {
							int add_occup = good_occup + 1;
							if (!(strict_limit && add_occup > NVariableLimit[1])) {
								if (all_yet_for_self) {
									add_updated = add_descriptors_Luhmann(&FuncType, &RegEle, &TestRegEle, &GoodRegEle, &TestGoodRegEle, &dataset, &testdataset, &arrange, &good_arrange, &add_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &R2Threshold, &betatype[ibeta], &MaxBFGSIterations, &for_test, &DuplicateX);
								}
								else {
									add_updated = add_descriptors(&FuncType, &RegEle, &TestRegEle, &GoodRegEle, &TestGoodRegEle, &dataset, &testdataset, &arrange, &good_arrange, &add_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &R2Threshold, &betatype[ibeta], &MaxBFGSIterations, &for_test, &DuplicateX, &Hegel[current_Hegel_pos].slave, &Hegel[current_Hegel_pos].master);
								}
							}
							else {
								add_updated = false;
							}
						}
						if (useExitMode && exit_mode(&macro_filename_prefix)) {
							std::cout << "[WARNING] pr_id = " << pr_id << " @ " << instally[iy].setY << " @ " << std::to_string(ix + 1) << postfix_x << " : use exit_mode before swap_updated.\n";
							std::cout.flush();
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						auto run_current_t = std::chrono::system_clock::now();
						std::chrono::duration<double> elapsed_secondsfr = run_current_t - run_start_t;
						long double duration = (long double)elapsed_secondsfr.count();
						if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
							std::cout << "[WARNING] pr_id = " << pr_id << " @ " << instally[iy].setY << " @ " << std::to_string(ix + 1) << postfix_x << " : over TargetCPUTimeAfterSwap before swap_updated.\n";
							std::cout.flush();
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						run_start_t = std::chrono::system_clock::now();
						if (all_yet_for_self) {
							swap_updated = swap_descriptors_Luhmann(&FuncType, &RegEle, &TestRegEle, &GoodRegEle, &TestGoodRegEle, &dataset, &testdataset, &arrange, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &betatype[ibeta], &MaxBFGSIterations, &for_test, &DuplicateX, &macro_filename_prefix);
						}
						else {
							swap_updated = swap_descriptors(&FuncType, &RegEle, &TestRegEle, &GoodRegEle, &TestGoodRegEle, &dataset, &testdataset, &arrange, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &betatype[ibeta], &MaxBFGSIterations, &for_test, &DuplicateX, &macro_filename_prefix, &Hegel[current_Hegel_pos].slave, &Hegel[current_Hegel_pos].master);
						}
						if (useExitMode && exit_mode(&macro_filename_prefix)) {
							std::cout << "[WARNING] pr_id = " << pr_id << " @ " << instally[iy].setY << " @ " << std::to_string(ix + 1) << postfix_x << " : use exit_mode before switch_updated.\n";
							std::cout.flush();
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix); 
								successful_for_self = 1;
							}
							//break;
						}
						run_current_t = std::chrono::system_clock::now();
						elapsed_secondsfr = run_current_t - run_start_t;
						duration = (long double)elapsed_secondsfr.count();
						if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
							std::cout << "[WARNING] pr_id = " << pr_id << " @ " << instally[iy].setY << " @ " << std::to_string(ix + 1) << postfix_x << " : over TargetCPUTimeAfterSwap before switch_updated.\n";
							std::cout.flush();
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						run_start_t = std::chrono::system_clock::now();
						if (all_yet_for_self) {
							switch_updated = switch_functype_Luhmann(&FuncType, &GroupX, &RegEle, &TestRegEle, &GoodRegEle, &TestGoodRegEle, &dataset, &testdataset, &arrange, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &betatype[ibeta], &MaxBFGSIterations, &for_test, &DuplicateX, &macro_filename_prefix);

						}
						else {
							switch_updated = switch_functype(&FuncType, &GroupX, &RegEle, &TestRegEle, &GoodRegEle, &TestGoodRegEle, &dataset, &testdataset, &arrange, &good_arrange, &good_occup, &RegressionCoeffConfindenceIntervalLevel, &PickUpPLevel, &betatype[ibeta], &MaxBFGSIterations, &for_test, &DuplicateX, &macro_filename_prefix, &Hegel[current_Hegel_pos].slave, &Hegel[current_Hegel_pos].master);
						}
						init_search = false;
						if (!init_search && !add_updated && !swap_updated && !switch_updated) {
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						if (good_occup + 1 > (signed)GroupX.size()) {
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						if (failure_count > max_failure_count) {
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						if (useExitMode && exit_mode(&macro_filename_prefix)) {
							std::cout << "[WARNING] pr_id = " << pr_id << " @ " << instally[iy].setY << " @ " << std::to_string(ix + 1) << postfix_x << " : use exit_mode after swap_updated and switch_updated.\n";
							std::cout.flush();
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						run_current_t = std::chrono::system_clock::now();
						elapsed_secondsfr = run_current_t - run_start_t;
						duration = (long double)elapsed_secondsfr.count();
						if (useTargetCPUTimeAfterSwap && duration > TargetCPUTimeAfterSwap) {
							std::cout << "[WARNING] pr_id = " << pr_id << " @ " << instally[iy].setY << " @ " << std::to_string(ix + 1) << postfix_x << " : over TargetCPUTimeAfterSwap after swap_updated and switch_updated.\n";
							std::cout.flush();
							if (successful_for_self == 0) {
								write_finalreport_per_core(&ibeta, &setY, &macro_filename_prefix);
								successful_for_self = 1;
							}
							//break;
						}
						
						MPI_Barrier(MPI_COMM_WORLD);

					}

				}

				MPI_Barrier(MPI_COMM_WORLD);
				if (pr_id == writeoutput_prid) {
					auto current_t = std::chrono::system_clock::now();
					std::chrono::duration<double> elapsed_secondsfr = current_t - start_t;
					writeoutput << "CPU_time_so_far=\t" << elapsed_secondsfr.count() << "\n";
					writeoutput.flush();
					writeoutput << "******************************************\n";
					writeoutput.flush();
					writeoutput << "***************FINAL_ANSWER***************\n";
					writeoutput.flush();
					writeoutput << "******************************************\n";
					writeoutput.flush();
					search_globally(&setY, &betatype[ibeta], &macro_filename_prefix);
					writeoutput << "******************************************\n";
					writeoutput.flush();
					writeoutput << "******************************************\n";
					writeoutput.flush();
					writeoutput << "******************************************\n";
					writeoutput.flush();
					remove_regele();
				}
				auto current_t_pr = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr_pr = current_t_pr - start_t;
				writeoutputpr << "CPU_time_so_far=\t" << elapsed_secondsfr_pr.count() << "\n";
				writeoutputpr.flush();
				writeoutputpr << "******************************************\n";
				writeoutputpr.flush();
				writeoutputpr << "***************FINAL_ANSWER***************\n";
				writeoutputpr.flush();
				writeoutputpr << "******************************************\n";
				writeoutputpr.flush();
				writeoutputpr << "******************************************\n";
				writeoutputpr.flush();
				writeoutputpr << "******************************************\n";
				writeoutputpr.flush();
				writeoutputpr << "******************************************\n";
				writeoutputpr.flush();
				MPI_Barrier(MPI_COMM_WORLD);

			}

			writeoutputpr.close();
			if (pr_id == writeoutput_prid) {
				writeoutput.close();
				remove_prtxt();
			}

			dataset.clear();
			testdataset.clear();
			column_str.clear();
			std::vector<std::string>().swap(InputDataFilename);
			prefilename = macro_filename_prefix + ".txt";
			MPI_Barrier(MPI_COMM_WORLD);

		}
	}

	writereport.close();

	if (pr_id == writeoutput_prid) {

		readreport.open(OutputFilename + ".txt");
		std::getline(readreport, sbuf);
		std::getline(readreport, sbuf);
		std::getline(readreport, sbuf);
		std::vector<std::string> total_recommendedlog;

		for (int iy = 0; iy < (signed)instally.size(); iy++) {
			
			long double aveR2; 
			long double max_aveR2 = -100000000000.0;
			std::string rlog;

			for (int ix = 0; ix < (signed)installx.size() * (signed)betatype.size(); ix++) {
				int ibuf, ldbuf;
				std::string logbuf;
				std::getline(readreport, logbuf);
				std::stringstream bootstrap;
				bootstrap.clear();
				bootstrap << logbuf;
				bootstrap >> sbuf >> ibuf >> sbuf >> sbuf >> aveR2;
				if (aveR2 >= max_aveR2) {
					max_aveR2 = aveR2;
					rlog = logbuf;
				}
			}

			total_recommendedlog.push_back(rlog);

		}

		readreport.close();

		writereport.open(OutputFilename + ".txt", std::ios::binary | std::ios::app);
		writereport << std::fixed << std::setprecision(OutputPrecision);
		writereport << "==========RECOMMEND==========\n";
		writereport.flush();
		writereport << "Y\tNdescriptors\tBeta\tLogFilename\t<R2>\tR2\tR2test\tMSE\tMSEtest\tRMSE\tRMSEtest\tMAE\tMAEtest\tDescriptors\n";
		writereport.flush();
		std::vector<int> great_betaindex;
		std::vector<std::string> great_filename;
		for (int itlog = 0; itlog < (signed)total_recommendedlog.size(); itlog++) {
			writereport << total_recommendedlog[itlog] << "\n";
			writereport.flush();
			std::stringstream findfilename;
			findfilename << total_recommendedlog[itlog];
			findfilename >> sbuf >> sbuf >> sbuf;
			if (sbuf == "NULL") {
				sbuf = "";
			}
			int ibeta = -1;
			for (int ib = 0; ib < (signed)betatype.size(); ib++) {
				if (sbuf == betatype[ib]) {
					ibeta = ib;
					break;
				}
			}
			great_betaindex.push_back(ibeta);
			findfilename >> sbuf;
			great_filename.push_back(sbuf);
		}
		writereport << "==========RECOMMEND===========\n";
		writereport.flush();

		for (int isub = 0; isub < (signed)great_filename.size(); isub++) {
			readsub.open(great_filename[isub]);
			int ibetacount = -1;
			while (true) {
				std::getline(readsub, sbuf);
				if (sbuf == "***************FINAL_ANSWER***************") {
					ibetacount++;
					if (ibetacount == great_betaindex[isub]) {
						break;
					}
				}
			}
			std::getline(readsub, sbuf);
			writereport << "******************************************" << "\n";
			writereport.flush();
			writereport << "***************FINAL_ANSWER***************" << "\n";
			writereport.flush();
			writereport << "******************************************" << "\n";
			writereport.flush();
			while (true) {
				std::getline(readsub, sbuf);
				writereport << sbuf << "\n";
				writereport.flush();
				if (sbuf == "******************************************") {
					break;
				}
			}
			writereport << "******************************************" << "\n";
			writereport.flush();
			writereport << "******************************************" << "\n";
			writereport.flush();
			readsub.close();
		}

		writereport.close();

		std::vector<std::string>().swap(total_recommendedlog);
		std::vector<int>().swap(great_betaindex);
		std::vector<std::string>().swap(great_filename);

	}

	if (pr_id == writeoutput_prid) {
		remove_smt();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;

}
