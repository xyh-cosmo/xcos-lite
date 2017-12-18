#include "CosmoParams.hpp"
#include "Background.hpp"
#include "Verbose.hpp"
#include "GSL_Math.hpp"
#include "Constants.hpp"

#ifndef __COSMOLOGY__
#define __COSMOLOGY__


namespace Cosmology{

    struct CosmologyWorkspace{

        int my_rank;

        ParamVector   Params;           //  used to store all relevant parameters in different modules.
        Background    *background;

        bool    Init( std::string paramfile );
        bool    Update( imcmc::imcmc_double param );
        void    AddParams( ParamVector& par );    //  used to add parameters coming from module other than cosmology, for example,
                                                        //  the SNe Lightcureve parameters.
        CosmologyWorkspace();
        ~CosmologyWorkspace();
    };
}

#endif
