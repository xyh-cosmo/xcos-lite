#include "CosmoParams.hpp"

//	define some macros to contral the outputs

#ifndef _max_width_
    #define _max_parname_width_ 20
	#define _max_value_width_	15
	#define _mpw_				_max_parname_width_
	#define _mvw_				_max_value_width_
#endif

#ifndef _set_parname_width_
	#define _set_parname_width_	std::setw(_mpw_)
	#define _spw_ 				_set_parname_width_
#endif

#ifndef _set_value_width_
	#define _set_value_width_ 	std::setw(_mvw_)
	#define _svw_ 				_set_value_width_
#endif

#ifndef _set_precision_
	#define _set_precision_ 	std::setprecision(10)
	#define _sp_ 				_set_precision_
#endif

void Print_imcmc_double( icos_double par, std::string pname ){
    std::cout << "\t" + pname << " = " << par[pname] << "\n";
}

//  ============================================================================
//  ============================================================================
void Print_ParamValue( std::string mesg ){
	if( mesg.size() > 0 )
		std::cout << "# " + mesg << std::endl;
}

void Print_ParamValue( std::string pname, std::string pvalue, std::string mesg ){
    
	std::cout << std::left << _spw_ << pname << " = " << _svw_ << _sp_ << pvalue << "\t"; 
	
	if( mesg.size() >0 )
		std::cout << "# " + mesg << std::endl;
}

void Print_ParamValue( std::string pname, bool pvalue, std::string mesg ){
    
    if( pvalue == true ){
    	std::cout << std::left << _spw_ << pname << " = " << _svw_ << "true\t";
    } else {
    	std::cout << std::left << _spw_ << pname << " = " << _svw_ << "false\t";
    }
	
	if( mesg.size() > 0 )
		std::cout << "# " + mesg << std::endl;
	else
		std::cout << std::endl;
}

void Print_ParamValue( std::string pname, int pvalue, std::string mesg ){
    
	std::cout << std::left << _spw_ << pname << " = " << _svw_ << _sp_ << pvalue << "\t";
	
	if( mesg.size() > 0 )
		std::cout << "# " + mesg << std::endl;
	else
		std::cout << std::endl;
}

void Print_ParamValue( std::string pname, double pvalue, std::string mesg ){
    
	std::cout << std::left << _spw_ << pname << " = " << _svw_ << _sp_ << pvalue << "\t";
	
	if( mesg.size() > 0 )
		std::cout << "# " + mesg << std::endl;
	else
		std::cout << std::endl;
}


void Print_ParamValue( std::ofstream& os, std::string mesg ){
	if( mesg.size() > 0 )
		os << "# " + mesg << std::endl;
}

void Print_ParamValue( std::ofstream& os, std::string pname, std::string pvalue, std::string mesg ){
    
	os << std::left << _spw_ << pname << " = " << _svw_ << _sp_ << pvalue << "\t";
	
	if( mesg.size() > 0 )
		os << "# " + mesg << std::endl;
	else
		os << std::endl;
}

void Print_ParamValue( std::ofstream& os, std::string pname, bool pvalue, std::string mesg ){
    
    if( pvalue == true ){
    	os << std::left << _spw_ << pname << " = " << _svw_ << "true\t";
    } else {
    	os << std::left << _spw_ << pname << " = " << _svw_ << "false\t";
    }
	
	if( mesg.size() > 0 )
		os << "# " + mesg << std::endl;
	else
		os << std::endl;
}

void Print_ParamValue( std::ofstream& os, std::string pname, int pvalue, std::string mesg ){
    
	os << std::left << _spw_ << pname << " = " << _svw_ << _sp_ << pvalue << "\t";
	
	if( mesg.size() > 0 ) {
		os << "# " + mesg << std::endl;
	} else {
		os << std::endl;
	}
}

void Print_ParamValue( std::ofstream& os, std::string pname, double pvalue, std::string mesg ){
    
	os << std::left << _spw_ << pname  << " = " << _svw_ << _sp_ << pvalue << "\t";
	
	if( mesg.size() > 0 ) {
		os << "# " + mesg << std::endl;
	} else {
		os << std::endl;
	}
}

//  ================================================================================================


bool ParamVector::Do_MCMC( std::string pname ){ //  check whether pname is a mcmc sampling parameter
    
    bool do_mcmc = false;

    icos_string_iterator itx = MCMC_Params.begin();

    while( itx != MCMC_Params.end() ){

        if( *itx == pname ){
            do_mcmc = true;
            break;
        } else {
            ++itx;
        }
    }

    return do_mcmc;
}

//  check whether the name passed into functions like Update(...) has already existed
bool ParamVector::HasName( std::string name ){
    
    bool has = false;
    if( Value.count(name) == 1 ) { //  Ok, found one
        has = true;
    } else {
        if( Value.count(name) == 0 ) {//  this is also fone ...
            has = false;
        } else if( Value.count(name) > 1 ) { //  This should NEVER happen
            std::string err = "ParamVector::HasName ==> found duplicate of: " + name;
            throw std::runtime_error(err);
        }
    }
    return has;
}

void ParamVector::AddParam( std::string PName, bool do_mcmc ){
    
    if( !HasName(PName) ){
        Value[PName] = -9999;

        if( do_mcmc ) 
            MCMC_Params.push_back(PName);
    }
    else{
        std::string err = "ParamVector::HasName ==> already has: " + PName + ", do not add again.";
        throw std::runtime_error(err);
    }
}

void ParamVector::AddParam_with_Value( std::string PName, double PValue, bool do_mcmc ){
    
    if( !HasName(PName) ){
        Value[PName] = PValue;

        if( do_mcmc )
            MCMC_Params.push_back(PName);
    }
    else{
        std::string err = "ParamVector::HasName ==> already has: " + PName + ", do not add again.";
        throw std::runtime_error(err);
    }
}

void ParamVector::UpdateParams( std::string PName, double PValue ){
    if( HasName(PName) ){
        Value[PName] = PValue;
    }
}

void ParamVector::UpdateParams( icos_double value ){
    
    it = value.begin();
    while( it != value.end() ){
        UpdateParams( it->first, it->second );
        ++it;
    }
}


double ParamVector::operator()( std::string ParamName ){
    
#ifndef _DEBUG_
    return Value[ParamName];
#elif defined(_DEBUG_)  //  debug mode is more safer, will check the exsitence
                        //  of ParamName, and this will be helpful when debugging
                        //  likelihoods etc.
    double val;
    if( Value.count(ParamName) == 1 )
        val = Value[ParamName];
    else{
        std::string err = "struct ParamVector::operator() ==> no \"" + ParamName + "\"";
        throw std::runtime_error(err);
    }
    return val;
#endif
}

double ParamVector::operator[]( std::string ParamName ){
    
#ifndef _DEBUG_
    return Value[ParamName];
#elif defined(_DEBUG_)  //  debug mode is more safer, will check the exsitence
                        //  of ParamName, and this will be helpful when debugging
                        //  likelihoods etc.
    double val;
    if( Value.count(ParamName) == 1 )
        val = Value[ParamName];
    else{
        std::string err = "struct ParamVector::operator[] ==> no \"" + ParamName + "\"";
        throw std::runtime_error(err);
    }
    return val;
#endif
}

void ParamVector::PrintValues(){
    
    std::cout << "\n****************  ParamVector Values  ****************\n";
    std::cout   << " = Note =\n"
                << " Parameters stored in ParamVector are mainly to be\n"
                << " sampled by a MCMC sampler, but it's still OK to put\n"
                << " non-MCMC variables into ParamVector.\n";

    std::cout << " *  "
              << std::setw(12)
              << "Param Name"
              << " ="
              << std::setw(15)
              << std::setprecision(10)
              << "Param Value"
              << std::setw(20)
              << "\tis mcmc variable? [yes/no] "
              << std::endl;

    it = Value.begin();
    while( it != Value.end() ){
    
        std::cout << " *  "
                  << std::setw(12)
                  << it->first
                  << " ="
                  << std::setw(15)
                  << std::setprecision(10)
                  << it->second ;

        if( Do_MCMC(it->first) )
            std::cout << "\t\tyes \n";
        else
            std::cout << "\t\tno \n";

        ++it;
    }

    std::cout << "\n";
}

void ParamVector::CopyInto(struct ParamVector& pv ){
    
    it = Value.begin();
    while( it != Value.end() ){
        if( !pv.HasName(it->first) ){
            pv.AddParam_with_Value(it->first, it->second);
//std::cout << " adding " << it->first << "\t = " << it->second << "\n";
        }
        ++it;
    }

    icos_string_iterator itx = MCMC_Params.begin();
    while( itx != MCMC_Params.end() ){
        if( !pv.Do_MCMC(*itx) ){    // if Do_MCMC() is true, then *itx has already in MCMC_Params
            pv.MCMC_Params.push_back(*itx);
//std::cout << " adding MCMC parameter: " << *itx << "\n";
        }
        ++itx;
    }
}
