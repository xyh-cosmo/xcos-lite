//  this header defines the data types to be used throught the project icos

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdlib>

#include <imcmc/imcmc.hpp>

#ifndef __COSMOPARAMS__
#define __COSMOPARAMS__

typedef std::map<std::string, double>           icos_double;
typedef std::map<std::string, double>::iterator icos_double_iterator;
typedef std::vector<std::string>                icos_string;
typedef std::vector<std::string>::iterator      icos_string_iterator;

void Print_imcmc_double( icos_double par, std::string pname );

void Print_ParamValue( std::string mesg );
void Print_ParamValue( std::string pname, std::string pvalue, std::string mesg="" );
void Print_ParamValue( std::string pname, bool pvalue, std::string mesg="" );
void Print_ParamValue( std::string pname, int pvalue, std::string mesg="" );
void Print_ParamValue( std::string pname, double pvalue, std::string mesg="" );

void Print_ParamValue( std::ofstream& os, std::string mesg );
void Print_ParamValue( std::ofstream& os, std::string pname, std::string pvalue, std::string mesg="" );
void Print_ParamValue( std::ofstream& os, std::string pname, bool pvalue, std::string mesg="" );
void Print_ParamValue( std::ofstream& os, std::string pname, int pvalue, std::string mesg="" );
void Print_ParamValue( std::ofstream& os, std::string pname, double pvalue, std::string mesg="" );


#ifndef __do_mcmc__
    #define __do_mcmc__ true
#endif

#ifndef __undo_mcmc__
    #define __undo_mcmc__ false
#endif

struct ParamVector{ //  double is used by default

    icos_string             MCMC_Params;    //  passed to imcmc::emcee_workspace::add_likelihood(..)
    icos_double             Value;          //  actually has both name and value
    icos_double_iterator    it;

    bool HasName( std::string name );
    bool Do_MCMC( std::string pname );
    void AddParam( std::string PName, bool do_mcmc=false );
    void AddParam_with_Value( std::string PName, double PValue, bool do_mcmc=false );
    void UpdateParams( std::string PName, double PValue );
    void UpdateParams( icos_double value );
    double operator()( std::string ParamName );
    double operator[]( std::string ParamName );

    void PrintValues();
    void CopyInto(struct ParamVector& pv );
};

typedef ParamVector PVector;

#endif  //_COSMOPARAMS_
