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

typedef struct {

	std::string name;
	std::vector<long double> v;
	std::vector<int> yhit;

}DATASET;

std::vector<DATASET> predataset;

typedef struct {

	std::string setX;
	bool negative_exists;
	long double min, ave, max, stdev;
	bool interact;

}GROUPX;

std::vector<GROUPX> GroupX;

typedef struct {

	std::string setX;
	bool minlimitX, maxlimitX;
	long double minX, maxX;

}SUBFILTERX;

std::vector<SUBFILTERX> subfilterx;

typedef struct {

	std::string setY;
	bool minlimitY, maxlimitY;
	long double minY, maxY;

}INSTALLY;

std::vector<INSTALLY> instally;

typedef struct {

	std::string postfix;
	int SamplingAction;
	long double SamplingFactor;
	int pick_count;

}SAMPLINGFACTOR;

std::vector<SAMPLINGFACTOR> SamplingFactor;

typedef struct {

	std::string name;
	int multi0div1;
	int whichx1, whichx2;

}INTERFEATURES;

std::vector<INTERFEATURES> Interfeatures;

std::string DataFilenamePrefix;
std::string LogOutputFilename;
int OutputPrecision;
int NameDuplicateHandling;
long double Dmin, Dmax, Dincrement;
std::vector<std::string> XSpace;
std::vector<std::string> InteractXSpace;
std::string ParseData;
int ParseGroup;

bool use_dataname = false;
std::string userdefinedName;
std::vector<std::string> column_str;

std::string sbuf;
std::ifstream readpara, readdata, readparse;
std::ofstream writelog, writeparse;

int writeoutput_prid = 0;
int pr_id = 0;

bool is_number(const std::string* s) {
	std::istringstream iss(*s);
	double d;
	return iss >> std::noskipws >> d && iss.eof();
}

void read_para() {

	readpara.open("basic_properties_GoodCurator.txt");

	readpara >> sbuf >> DataFilenamePrefix;
	readpara >> sbuf >> LogOutputFilename;
	readpara >> sbuf >> OutputPrecision;
	readpara >> sbuf >> NameDuplicateHandling;
	if (NameDuplicateHandling == 2) {
		readpara >> Dmin >> Dmax >> Dincrement;
	}
	std::getline(readpara, sbuf);
	
	std::stringstream localstr;
	std::getline(readpara, sbuf);
	localstr << sbuf;
	localstr >> sbuf;
	while (true) {
		localstr >> sbuf;
		XSpace.push_back(sbuf);
		if (localstr.eof()) {
			break;
		}
	}

	std::getline(readpara, sbuf);
	localstr.clear();
	localstr << sbuf;
	localstr >> sbuf;
	while (true) {
		localstr >> sbuf;
		if (sbuf == "-") {
			break;
		}
		InteractXSpace.push_back(sbuf);
		if (localstr.eof()) {
			break;
		}
	}

	readpara >> sbuf >> ParseData;
	readpara >> sbuf >> ParseGroup;
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
		if (sbuf.find("STOP") != std::string::npos || readpara.eof()) {
			break;
		}
		else {
			localstr.clear();
			localstr << sbuf;
			SAMPLINGFACTOR objSAMPLINGFACTOR;
			localstr >> objSAMPLINGFACTOR.postfix >> objSAMPLINGFACTOR.SamplingAction >> objSAMPLINGFACTOR.SamplingFactor;
			objSAMPLINGFACTOR.pick_count = 0;
			SamplingFactor.push_back(objSAMPLINGFACTOR);
		}
	}

	readpara.close();

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

bool pick_data(int* testindex, SAMPLINGFACTOR* SamplingFactor, std::vector<DATASET>* predataset, std::string* ParseData, int* ParseGroup) {

	std::vector<std::string> atomname1, atomname2;
	std::vector<long double> natom1, natom2;

	readparse.open(*ParseData + ".temporary.txt");

	atomname1.clear();
	natom1.clear();

	for (int iline = 0; iline < *testindex; iline++) {
		std::getline(readparse, sbuf);
	}
	int barcount = 0;
	while (true) {
		readparse >> sbuf;
		if (sbuf == "||") {
			barcount++;
			if (barcount == 2) {
				break;
			}
		}
	}
	barcount = 0;
	while (true) {
		readparse >> sbuf;
		if (sbuf == "||") {
			barcount++;
			if (barcount == *ParseGroup) {
				break;
			}
		}
		else {
			atomname1.push_back(sbuf);
			long double ldbuf;
			readparse >> ldbuf;
			natom1.push_back(ldbuf);
		}
	}
	
	readparse.close();

	readparse.open(*ParseData + ".temporary.txt");

	bool pick = true;
	int linecount = 0;
	while (true) {
		atomname2.clear();
		natom2.clear();
		barcount = 0;
		while (true) {
			readparse >> sbuf;
			if (readparse.eof()) {
				break;
			}
			if (sbuf == "||") {
				barcount++;
				if (barcount == 2) {
					break;
				}
			}
		}
		if (readparse.eof()) {
			break;
		}
		barcount = 0;
		while (true) {
			readparse >> sbuf;
			if (sbuf == "||") {
				barcount++;
				if (barcount == *ParseGroup) {
					break;
				}
			}
			else {
				atomname2.push_back(sbuf);
				long double ldbuf;
				readparse >> ldbuf;
				natom2.push_back(ldbuf);
			}
		}
		if (atomname1.size() == atomname2.size()) {
			bool atomname_match = true;
			for (int iatom = 0; iatom < (signed)atomname1.size(); iatom++) {
				if (atomname1[iatom] != atomname2[iatom]) {
					atomname_match = false;
					break;
				}
			}
			if (atomname_match) {
				long double natomsum1 = 0.0;
				long double natomsum2 = 0.0;
				for (int iatom = 0; iatom < (signed)natom1.size(); iatom++) {
					natomsum1 += natom1[iatom];
					natomsum2 += natom2[iatom];
				}
				long double atom_dist = 0.0;
				for (int iatom = 0; iatom < (signed)natom1.size(); iatom++) {
					atom_dist += powl(natom1[iatom] / natomsum1 - natom2[iatom] / natomsum2, 2.0);
				}
				atom_dist = sqrtl(atom_dist);
				if (atom_dist <= SamplingFactor->SamplingFactor) {
					if (SamplingFactor->SamplingAction == 0) {
						if ((*predataset)[*testindex].v[0] > (*predataset)[linecount].v[0]) {
							pick = false;
						}
					}
					else {
						if ((*predataset)[*testindex].v[0] < (*predataset)[linecount].v[0]) {
							pick = false;
						}
					}
				}
			}
		}
		if (!pick) {
			break;
		}
		else {
			std::getline(readparse, sbuf);
			linecount++;
			if (readparse.eof() || linecount >= predataset->size()) {
				break;
			}
		}
	}

	readparse.close();

	std::vector<std::string>().swap(atomname1);
	std::vector<std::string>().swap(atomname2);
	std::vector<long double>().swap(natom1);
	std::vector<long double>().swap(natom2);

	return pick;

}

void read_data(std::vector<DATASET>* predataset, bool* use_dataname, std::string* userdefinedName, std::vector<INTERFEATURES>* Interfeatures, std::vector<std::string>* DataFilename, std::string* setY, std::vector<std::string>* XSpace, std::vector<std::string>* InteractXSpace, std::vector<GROUPX>* GroupX, bool* minlimitY, bool* maxlimitY, long double* minY, long double* maxY, std::vector<std::string>* column_str, int* DuplicateX, int* NameDuplicateHandling, long double* Dmin, long double* Dmax, long double* Dincrement, std::vector<SUBFILTERX>* subfilterx, std::string* ParseData, int* ParseGroup, std::vector<SAMPLINGFACTOR>* SamplingFactor) {

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

	std::stringstream localstr;
	unsigned int initial_column_str_count = 0;
	Interfeatures->clear();

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
		else if (sbuf != *setY) {
			bool add = false;
			for (int ixspace = 0; ixspace < (signed)XSpace->size(); ixspace++) {
				if (sbuf.find((*XSpace)[ixspace]) != std::string::npos) {
					add = true;
					break;
				}
			}
			if (add) {
				GROUPX objGROUPX;
				objGROUPX.setX = sbuf;
				objGROUPX.negative_exists = false;
				objGROUPX.min = 1000000000000000000000000000000000000000.0;
				objGROUPX.ave = 0.0;
				objGROUPX.max = -objGROUPX.min;
				objGROUPX.stdev = 0.0;
				objGROUPX.interact = false;
				for (int iint = 0; iint < (signed)InteractXSpace->size(); iint++) {
					if (objGROUPX.setX.find((*InteractXSpace)[iint]) != std::string::npos) {
						objGROUPX.interact = true;
						break;
					}
				}
				for (int ix = 0; ix < *DuplicateX; ix++) {
					GroupX->push_back(objGROUPX);
				}
			}
			bool conisder = false;

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
		WhichColumnSavedSubFilterX[find_column(&(*subfilterx)[i].setX, column_str)] = i;
	}

	readdata.close();

	int wrongdatacount = 0;
	int duplicatecount = 0;
	int outofrangecount = 0;

	readparse.open(*ParseData);
	writeparse.open(*ParseData + ".temporary.txt", std::ios::binary);
	std::string linebuf;

	for (int ifile = 0; ifile < (signed)DataFilename->size(); ifile++) {

		readdata.open((*DataFilename)[ifile]);
		std::getline(readdata, sbuf);

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
			std::getline(readparse, linebuf);
			std::stringstream linestr;
			linestr << sbuf;
			for (unsigned int i = initial_column_str_count; i < column_str->size(); i++) {
				linestr >> sbuf;
				if (WhichColumnSaved[i] == -2) {
					objDATASET.name = sbuf;
					if (*NameDuplicateHandling == 0 || *NameDuplicateHandling == 1 || *NameDuplicateHandling == 2) {
						for (int ipredata = 0; ipredata < (signed)predataset->size(); ipredata++) {
							if ((*predataset)[ipredata].name == sbuf) {
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
							if (pr_id == writeoutput_prid) {
								writelog << "... Y is out of window!\n";
								writelog.flush();
								writelog << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << ldbuf << "\n";
								writelog.flush();
							}
							break;
						}
					}
					else {
						full_number = false;
						wrongdatacount++;
						if (pr_id == writeoutput_prid) {
							writelog << "... data are given wrong; this line will be ignored!\n";
							writelog.flush();
							writelog << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << sbuf << "\n";
							writelog.flush();
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
							if (pr_id == writeoutput_prid) {
								writelog << "... X(" << (*subfilterx)[WhichColumnSavedSubFilterX[i]].setX << ") is out of window!\n";
								writelog.flush();
								writelog << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << ldbuf << "\n";
								writelog.flush();
							}
							break;
						}
					}
					else {
						subfilter = false;
						wrongdatacount++;
						if (pr_id == writeoutput_prid) {
							writelog << "... data are given wrong; this line will be ignored!\n";
							writelog.flush();
							writelog << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " column: " << i << " " << (*column_str)[i] << " ... " << sbuf << "\n";
							writelog.flush();
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
					predataset->push_back(objDATASET);
					writeparse << linebuf << "\n";
					writeparse.flush();
				}
				else {
					duplicatecount++;
					if ((*NameDuplicateHandling == 0 && objDATASET.v[0] < (*predataset)[duplicate_with_existing_datum].v[0]) || (*NameDuplicateHandling == 1 && objDATASET.v[0] > (*predataset)[duplicate_with_existing_datum].v[0])) {
						if (pr_id == writeoutput_prid) {
							writelog << "... Datum was replaced due to duplicate issue!\n";
							writelog.flush();
							writelog << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " ... " << objDATASET.name << ": " << objDATASET.v[0] << " <--- " << (*predataset)[duplicate_with_existing_datum].v[0] << "\n";
							writelog.flush();
						}
						(*predataset)[duplicate_with_existing_datum] = objDATASET;
					}
					else if ((*NameDuplicateHandling == 0 && objDATASET.v[0] > (*predataset)[duplicate_with_existing_datum].v[0]) || (*NameDuplicateHandling == 1 && objDATASET.v[0] < (*predataset)[duplicate_with_existing_datum].v[0])) {
						if (pr_id == writeoutput_prid) {
							writelog << "... Datum was not replaced in spite of duplicate issue!\n";
							writelog.flush();
							writelog << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " ... " << objDATASET.name << ": " << objDATASET.v[0] << " <-X- " << (*predataset)[duplicate_with_existing_datum].v[0] << "\n";
							writelog.flush();
						}
					}
					else if (*NameDuplicateHandling == 2) {
						int whichhistogramset = (int)(((long double)histogramset - 1.0) * (objDATASET.v[0] - histogram0) / (histogram1 - histogram0));
						if (whichhistogramset >= 0 && whichhistogramset < (signed)objDATASET.yhit.size()) {
							(*predataset)[duplicate_with_existing_datum].yhit[whichhistogramset]++;
							if (pr_id == writeoutput_prid) {
								writelog << "... Datum duplicate!\n";
								writelog.flush();
								writelog << "check line: " << (*DataFilename)[ifile] << " -- " << linecount << " ... " << objDATASET.name << "\n";
								writelog.flush();
							}
						}
					}
				}
			}
			if (readdata.eof()) {
				break;
			}
		}

		readdata.close();

	}

	readparse.close();
	writeparse.close();

	if (*NameDuplicateHandling == 2) {
		for (int idat = 0; idat < (signed)predataset->size(); idat++) {
			int maxhit = -100;
			int targethit = -1;
			for (int ihit = 0; ihit < (signed)(*predataset)[idat].yhit.size(); ihit++) {
				if ((*predataset)[idat].yhit[ihit] > maxhit) {
					maxhit = (*predataset)[idat].yhit[ihit];
					targethit = ihit;
				}
			}
			if (targethit >= 0) {
				(*predataset)[idat].v[0] = *Dmin + (long double)targethit * *Dincrement;
			}
		}
	}

	if (pr_id == writeoutput_prid) {
		writelog << "TOTAL DATA NO = " << predataset->size() << "\n";
		writelog.flush();
		writelog << "WRONG DATA NO = " << wrongdatacount << "\n";
		writelog.flush();
		writelog << "DUPLICATE DATA NO = " << duplicatecount << "\n";
		writelog.flush();
		writelog << "OUT-OF-RANGE DATA NO = " << outofrangecount << "\n";
		writelog.flush();
	}

	if (predataset->size() > 0) {

		for (int ix = 0; ix < (signed)GroupX->size(); ix++) {
			(*GroupX)[ix].min = 1000000000000000000000000000000000000000.0;
			(*GroupX)[ix].ave = 0.0;
			(*GroupX)[ix].max = -(*GroupX)[ix].min;
			(*GroupX)[ix].stdev = 0.0;
			(*GroupX)[ix].negative_exists = false;
		}
		for (int i = 0; i < (signed)predataset->size(); i++) {
			for (int j = 1; j < (signed)(*predataset)[0].v.size(); j++) {
				if ((*predataset)[i].v[j] < 0.0) {
					(*GroupX)[j - 1].negative_exists = true;
				}
				if ((*predataset)[i].v[j] < (*GroupX)[j - 1].min) {
					(*GroupX)[j - 1].min = (*predataset)[i].v[j];
				}
				if ((*predataset)[i].v[j] > (*GroupX)[j - 1].max) {
					(*GroupX)[j - 1].max = (*predataset)[i].v[j];
				}
				(*GroupX)[j - 1].ave += (*predataset)[i].v[j];
			}
		}
		for (int j = 1; j < (signed)(*predataset)[0].v.size(); j++) {
			(*GroupX)[j - 1].ave /= (long double)(signed)predataset->size();
		}
		for (int i = 0; i < (signed)predataset->size(); i++) {
			for (int j = 1; j < (signed)(*predataset)[0].v.size(); j++) {
				(*GroupX)[j - 1].stdev += powl((*predataset)[i].v[j] - (*GroupX)[j - 1].ave, 2.0);
			}
		}
		for (int j = 1; j < (signed)(*predataset)[0].v.size(); j++) {
			(*GroupX)[j - 1].stdev /= (long double)(signed)predataset->size();
			(*GroupX)[j - 1].stdev = sqrtl((*GroupX)[j - 1].stdev);
		}

	}

	if (pr_id == writeoutput_prid) {
		writelog << "******************************DATASET******************************\n";
		writelog.flush();
		if (*use_dataname) {
			writelog << *userdefinedName << "\t";
			writelog.flush();
		}
		writelog << *setY << "\t";
		writelog.flush();
		for (unsigned int i = 0; i < GroupX->size(); i++) {
			writelog << (*GroupX)[i].setX;
			writelog.flush();
			if (i != GroupX->size() - 1) {
				writelog << "\t";
				writelog.flush();
			}
		}
		if (InteractXSpace->size() > 0) {
			for (unsigned int i = 0; i < GroupX->size(); i++) {
				for (unsigned int j = i + 1; j < GroupX->size(); j++) {
					if ((*GroupX)[i].interact && (*GroupX)[j].interact) {
						INTERFEATURES objINTERFEATURES;
						objINTERFEATURES.name = "(" + (*GroupX)[i].setX + "*" + (*GroupX)[j].setX + ")";
						objINTERFEATURES.multi0div1 = 0;
						objINTERFEATURES.whichx1 = i;
						objINTERFEATURES.whichx2 = j;
						Interfeatures->push_back(objINTERFEATURES);
						if ((((*GroupX)[i].min < 0 && (*GroupX)[i].max < 0) || ((*GroupX)[i].min > 0 && (*GroupX)[i].max > 0)) && (((*GroupX)[j].min < 0 && (*GroupX)[j].max < 0) || ((*GroupX)[j].min > 0 && (*GroupX)[j].max > 0))) {
							objINTERFEATURES.name = "(" + (*GroupX)[i].setX + "/" + (*GroupX)[j].setX + ")";
							objINTERFEATURES.multi0div1 = 1;
							objINTERFEATURES.whichx1 = i;
							objINTERFEATURES.whichx2 = j;
							Interfeatures->push_back(objINTERFEATURES);
						}
					}
				}
			}
			if (Interfeatures->size() > 0) {
				writelog << "\t";
				writelog.flush();
				for (unsigned int i = 0; i < Interfeatures->size(); i++) {
					writelog << (*Interfeatures)[i].name;
					writelog.flush();
					if (i != Interfeatures->size() - 1) {
						writelog << "\t";
						writelog.flush();
					}
				}
			}
		}
		if (*ParseData != "-" && SamplingFactor->size() > 0) {
			writelog << "\t";
			writelog.flush();
			for (int isampling = 0; isampling < (signed)SamplingFactor->size(); isampling++) {
				writelog << "Pick" << (*SamplingFactor)[isampling].postfix;
				writelog.flush();
				if (isampling != (signed)SamplingFactor->size() - 1) {
					writelog << "\t";
					writelog.flush();
				}
			}
		}
		writelog << "\n";
		writelog.flush();
		for (unsigned int j = 0; j < predataset->size(); j++) {
			if (*use_dataname) {
				writelog << (*predataset)[j].name << "\t";
				writelog.flush();
			}
			for (unsigned int i = 0; i < (*predataset)[j].v.size(); i++) {
				writelog << (*predataset)[j].v[i];
				writelog.flush();
				if (i != (*predataset)[j].v.size() - 1) {
					writelog << "\t";
					writelog.flush();
				}
			}
			if (Interfeatures->size() > 0) {
				writelog << "\t";
				writelog.flush();
				for (unsigned int i = 0; i < Interfeatures->size(); i++) {
					int target_x1 = (*Interfeatures)[i].whichx1;
					int target_x2 = (*Interfeatures)[i].whichx2;
					int mul0div1 = (*Interfeatures)[i].multi0div1;
					long double val = (*predataset)[j].v[target_x1 + 1];
					if (mul0div1 == 0) {
						val *= (*predataset)[j].v[target_x2 + 1];
					}
					else {
						val /= (*predataset)[j].v[target_x2 + 1];
					}
					writelog << val;
					writelog.flush();
					if (i != Interfeatures->size() - 1) {
						writelog << "\t";
						writelog.flush();
					}
				}
			}
			if (*ParseData != "-" && SamplingFactor->size() > 0) {
				writelog << "\t";
				writelog.flush();
				for (int isampling = 0; isampling < (signed)SamplingFactor->size(); isampling++) {
					int signedj = (signed)j;
					if (pick_data(&signedj, &(*SamplingFactor)[isampling], predataset, ParseData, ParseGroup)) {
						(*SamplingFactor)[isampling].pick_count++;
						writelog << "1";
						writelog.flush();
					}
					else {
						writelog << "0";
						writelog.flush();
					}
					if (isampling != (signed)SamplingFactor->size() - 1) {
						writelog << "\t";
						writelog.flush();
					}
				}
			}
			writelog << "\n";
			writelog.flush();
		}
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();
	}

	if (pr_id == writeoutput_prid) {

		for (int ix = 0; ix < (signed)GroupX->size(); ix++) {
			writelog << (*GroupX)[ix].setX;
			writelog.flush();
			if (ix != (signed)GroupX->size() - 1) {
				writelog << "\t";
				writelog.flush();
			}
			else {
				writelog << "\n";
				writelog.flush();
			}
		}
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "Name\tNegativeExists\tMin\tAve\tMax\tStdev\n";
		writelog.flush();
		for (int ix = 0; ix < (signed)GroupX->size(); ix++) {
			writelog << (*GroupX)[ix].setX << "\t";
			writelog.flush();
			if ((*GroupX)[ix].negative_exists) {
				writelog << "true\t";
				writelog.flush();
			}
			else {
				writelog << "false\t";
				writelog.flush();
			}
			writelog << (*GroupX)[ix].min << "\t" << (*GroupX)[ix].ave << "\t" << (*GroupX)[ix].max << "\t" << (*GroupX)[ix].stdev << "\n";
			writelog.flush();
		}

		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();
		for (int isampling = 0; isampling < (signed)SamplingFactor->size(); isampling++) {
			writelog << (*SamplingFactor)[isampling].postfix << "\t---->\t" << (*SamplingFactor)[isampling].pick_count << "\t" << 100.0 * (long double)(*SamplingFactor)[isampling].pick_count / (long double)predataset->size() << "%\n";
			writelog.flush();
		}
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();
		writelog << "******************************************\n";
		writelog.flush();

	}

	delete[] WhichColumnSaved;
	WhichColumnSaved = nullptr;
	delete[] WhichColumnSavedSubFilterX;
	WhichColumnSavedSubFilterX = nullptr;

	char cmd[100];
#ifdef _WIN32
	strcpy_s(cmd, "del *.temporary.txt");
	system(cmd);
#elif _WIN64
	strcpy_s(cmd, "del *.temporary.txt");
	system(cmd);
#else
	strcpy(cmd, "rm *.temporary.txt");
	system(cmd);
#endif

}



int main() {
	
	read_para();
	
	writelog.open(LogOutputFilename, std::ios::binary);
	writelog << std::fixed << std::setprecision(OutputPrecision);

	for (int iy = 0; iy < (signed)instally.size(); iy++) {
		predataset.clear();
		std::string setY = instally[iy].setY;
		long double minY = instally[iy].minY;
		long double maxY = instally[iy].maxY;
		bool limitminY = instally[iy].minlimitY;
		bool limitmaxY = instally[iy].maxlimitY;
		std::vector<std::string> InputDataFilename;
		InputDataFilename.clear();
		InputDataFilename.push_back(DataFilenamePrefix + ".txt");
		int DuplicateX = 1;
		read_data(&predataset, &use_dataname, &userdefinedName, &Interfeatures, &InputDataFilename, &setY, &XSpace, &InteractXSpace, &GroupX, &limitminY, &limitmaxY, &minY, &maxY, &column_str, &DuplicateX, &NameDuplicateHandling, &Dmin, &Dmax, &Dincrement, &subfilterx, &ParseData, &ParseGroup, &SamplingFactor);
	}

	writelog.close();

	return 0;

}
