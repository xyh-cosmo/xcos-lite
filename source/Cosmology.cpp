#include "Cosmology.hpp"
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "SNE.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
// using namespace Cosmology;

namespace Cosmology{

    CosmologyWorkspace::CosmologyWorkspace(){
    	my_rank    = MPI::COMM_WORLD.Get_rank();
        background = NULL;
    }

    CosmologyWorkspace::~CosmologyWorkspace(){
    	if( background != NULL ){
			background->~Background();		//	this will cause double free() error
	        background = NULL;
    	}
    }

	void CosmologyWorkspace::AddParams( ParamVector& par ){
		par.CopyInto(this->Params);
	}

	bool CosmologyWorkspace::Init( std::string paramfile ){

		background      = new Background;

	//  now intialize from paramfile
		if( background->Init(paramfile) == false ){
            cout << "@Error: failed to initialize the background.\n";
            return false;
        }

    //  copy background parameters into CosmologyWorkspace::Params;
		background->Params.CopyInto(this->Params);

        MPI_cout("CosmologyWorkspace initialization completed!", my_rank);

        return true;
	}


	bool CosmologyWorkspace::Update( imcmc_double param ){

	    background->Params.UpdateParams(param);

        if( background->Update() == false ){
            cout << "## ==> failed to update background\n";
            return false;
        }

        return true;
	}

}
