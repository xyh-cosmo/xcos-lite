#include <imcmc/parser++.hpp>

#include "Cosmology.hpp"
#include "HST.hpp"

using namespace std;
using namespace imcmc::parser;

namespace Data{

Data_HST::Data_HST(){

	my_rank			= MPI::COMM_WORLD.Get_rank();

    HST_H0          = 73.8;
    HST_sigma_H0    = 2.4;
    Params.AddParam("H0", true);
	Params.UpdateParams("H0", 70);
}

void Data_HST::Read_H0(std::string paramfile){
	if( Read::Has_Key_in_File(paramfile, "HST_H0") ){
		if( Read::Has_Value(paramfile, "HST_H0", "double") )
			HST_H0		= Read::Read_Double_from_File(paramfile, "HST_H0");
	}
	else
		HST_H0			= 73.8;
}

void Data_HST::Read_Sigma_H0(std::string paramfile){
	if( Read::Has_Key_in_File(paramfile, "HST_sigma_H0") ){
		if( Read::Has_Value(paramfile, "HST_sigma_H0", "double") )
			HST_sigma_H0= Read::Read_Double_from_File(paramfile, "HST_sigma_H0");
	}
	else
		HST_sigma_H0	= 2.4;

	MPI_cout("HST_sigma_H0	= " + Read::DoubleToString(HST_sigma_H0), my_rank);
}

void Data_HST::Init(std::string paramfile){
	Read_H0(paramfile);
	Read_Sigma_H0(paramfile);
	Params.UpdateParams("H0", 70);
}


double Data_HST::H0(){
	return Params["H0"];
}

}
