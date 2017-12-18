#ifndef _ICOSMO_CONSTANTS_
#define _ICOSMO_CONSTANTS_


namespace Constants{

#ifndef __DEFINED_CONSTS__
#define __DEFINED_CONSTS__

//  Mathematical constants taken from GSL
    #define _PI_                3.14159265358979323846264338328     // pi 
    #define _SQRT_PI_           1.77245385090551602729816748334     // sqrt(pi) 
    #define _E_                 2.71828182845904523536028747135     // e 
    #define LOG2_E_             1.44269504088896340735992468100     // log_2 (e) 
    #define LOG10_E_            0.43429448190325182765112891892     // log_10 (e) 
    #define LN_2_               0.69314718055994530941723212146     // ln(2) 
    #define LN_10_              2.30258509299404568401799145468     // ln(10) 
    #define LN_PI_              1.14472988584940017414342735135     // ln(pi) 
    #define SQRT_2_             1.41421356237309504880168872421     // sqrt(2) 
    #define SQRT_1_2_           0.70710678118654752440084436210     // sqrt(1/2) 
    #define SQRT_3_             1.73205080756887729352744634151     // sqrt(3) 


//  These are taken from class
    #define _Mpc_over_m_        3.085677581282e22                   /**< conversion factor from meters to megaparsecs */
                                                                    /* remark: CAMB uses 3.085678e22: good to know if you want to compare  with high accuracy */

    #define _Gyr_over_Mpc_      3.06601394e2                        /**< conversion factor from megaparsecs to gigayears
				                                                        (c=1 units, Julian years of 365.25 days) */
    #define _c_                 2.99792458e8                        /**< c in m/s */
	#define _c_in_km_s_			2.99792458e5						/**< c in km/s */
    #define _G_                 6.67428e-11                         /**< Newton constant in m^3/Kg/s^2 */
    #define _eV_                1.602176487e-19                     /**< 1 eV expressed in J */

    /* parameters entering in Stefan-Boltzmann constant sigma_B */
    #define _k_B_ 1.3806504e-23
    #define _h_P_ 6.62606896e-34
    /* remark: sigma_B = 2 pi^5 k_B^4 / (15h^3c^2) = 5.670400e-8
                       = Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

    #define _sigma_B_ (2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2))

#endif  //  __DEFINED_CONSTS__
}


#endif
