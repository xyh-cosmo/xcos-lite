/*
    == parameters ==
    stop_on_A : if ture, the calclulation will stop after saved 'mu_zs2zp.txt'
    LSST_SN_density_file : photoz number density, used to calibrate predictions of mu(zp)
    LSST_SN_A_nrow  : # of rows of matrix A
    LSST_SN_A_ncol  : # of cols of matrix A
    LSST_SN_dataset : # LSST SN data filename
*/


#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "Cosmology.hpp"
#include "SNE.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
using namespace arma;

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

namespace Data{

    SNeIa_LSST::SNeIa_LSST(){
        sn_num  = 0;
        sn_z    = NULL;
        sn_mu   = NULL;
        sn_err  = NULL;
        sn_zerr = NULL;

        acc_ns      = NULL;
        spline_ns   = NULL;
        zs = NULL;
        zp = NULL;

        acc_muzperr = NULL;
        spline_muzperr  = NULL;
        muzperr_z = NULL;
        muzperr = NULL;
    }

    SNeIa_LSST::~SNeIa_LSST(){

        if( sn_z != NULL ){
            delete[] sn_z;
            sn_z = NULL;
        }

        if( sn_mu != NULL ){
            delete[] sn_mu;
            sn_mu = NULL;
        }

        if( sn_err != NULL ){
            delete[] sn_err;
            sn_err = NULL;
        }

        if( sn_zerr != NULL ){
            delete[] sn_zerr;
            sn_zerr = NULL;
        }

        if( acc_ns != NULL )
            gsl_interp_accel_free(acc_ns);

        if( spline_ns != NULL )
            gsl_spline_free(spline_ns);

        if( zp != NULL )
            delete[] zp;

        if( zs != NULL )
            delete[] zs;

        if( acc_muzperr != NULL ) gsl_interp_accel_free(acc_muzperr);
        if( spline_muzperr != NULL ) gsl_spline_free(spline_muzperr);
        if( muzperr_z != NULL ) delete[] muzperr_z;
        if( muzperr != NULL ) delete[] muzperr;
    }

    void SNeIa_LSST::init_matrix_A( string& paramfile ){

        bool stop_on_A = false;
        stop_on_A = Read::Read_Bool_from_File(paramfile,"stop_on_A");

        string SN_density_file;
        if( Read::Has_Key_in_File(paramfile,"LSST_SN_density_file") ){
            SN_density_file = Read::Read_String_from_File(paramfile,"LSST_SN_density_file");
        }
        else{
            throw runtime_error("# cannot find key LSST_SN_density_file!\n");
        }

        arma::mat sn_density;
        sn_density.load(SN_density_file, arma::auto_detect);

        if( sn_density.n_rows <= 0 || sn_density.n_cols != 2 ){
            throw runtime_error("# error in reading SN_density_file!\n");
        }

        int ns_size = sn_density.n_rows;

        if( Read::Has_Key_in_File(paramfile,"LSST_SN_A_nrow") ){
            A_nrow = Read::Read_Int_from_File(paramfile,"LSST_SN_A_nrow");
            if( A_nrow < 100 ){
                throw runtime_error("## LSST_SN_A_nrow is too small (<100), choose a larger number!\n");
            }
        }else{
            A_nrow = 100;
            cout << "==> LSST_SN_A_nrow set to default value 100\n";
        }

        if( Read::Has_Key_in_File(paramfile,"LSST_SN_A_ncol") ){
            A_ncol = Read::Read_Int_from_File(paramfile,"LSST_SN_A_ncol");
            if( A_ncol < 200 ){
                throw runtime_error("## LSST_SN_A_ncol is too small (<200), choose a lager number!\n");
            }
        }else{
            A_ncol = 200;
            cout << "==> LSST_SN_A_ncol  set to default value 200\n";
        }

        A = arma::zeros(A_nrow,A_ncol);

        double *z = new double[ns_size];
        double *ns= new double[ns_size];

        double ns_norm = 0.0;
        zs_min = zp_min = 10;
        zs_max = zp_max = 0.0;

        for( int i=0; i<ns_size; ++i ){
            z[i] = sn_density(i,0);
            ns[i]= sn_density(i,1);
            ns_norm += ns[i];

            if( zs_min > z[i] ) zs_min = z[i];
            if( zs_max < z[i] ) zs_max = z[i];

            if( zp_min > z[i] ) zp_min = z[i];
            if( zp_max < z[i] ) zp_max = z[i];
        }

        for( int i=0; i<ns_size; ++i )
            ns[i] /= ns_norm;

        acc_ns = gsl_interp_accel_alloc();
        spline_ns = gsl_spline_alloc(gsl_interp_cspline, ns_size);

    //  initialize spline_ns and pass it into photoz_mu_integrand(**)
        gsl_spline_init(spline_ns, z, ns, ns_size);
        // exit(0);

        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        gsl_function F;
        F.function = photoz_mu_integrand;
        F.params = this;

    //  #############################################
    //  NOET the small difference between dzp and dzs
    //  #############################################
        double dzp = (zp_max-zp_min)/(A_nrow-1.);
        double dzs = (zs_max-zs_min)/A_ncol;
        double result, error;

        zp = new double[A_nrow];
        zs = new double[A_ncol];

        for( int i=0; i<A_nrow; ++i ){

            zp_current = zp_min + i*dzp;
            zp[i] = zp_current;

            double norm = 0.0;
            for( int j=0; j<A_ncol; ++j ){
                double zs1 = zs_min + j*dzs;
                double zs2 = zs_min + (j+1)*dzs;
                // cout << "zs1 = " << zs1 << "\tzs2 = " << zs2 << endl;
                gsl_integration_qags( &F, zs1, zs2, 1E-8, 1E-8, 1000, w, &result, &error );
                A(i,j) = result;
                norm += result;
            }

        //  normalization
            for( int j=0; j<A_ncol; ++j )
                A(i,j) /= norm;
        }

        // initialize zs[], they are central values of each small redshift bin.
        for( int i=0; i<A_ncol; ++i ){
            zs[i] = zs_min+(i+0.5)*dzs;
        }

        delete[] z;
        delete[] ns;

        gsl_integration_workspace_free(w);

        if( MPI::COMM_WORLD.Get_rank() == 0 ){
            ofstream of("mu_zs2zp.txt");    // improvement is needed here !!!
            A.save(of, arma::raw_ascii);
            of.close();
            if(stop_on_A){
                cout << "\n\n==> You requested to stop after having saved mu_zs2zp.txt\n";
                exit(0);
            }
        }
    }


    bool SNeIa_LSST::ReadData( string& paramfile ){

        init_matrix_A(paramfile);

        bool init_status = false;
        arma::mat sn_data_temp;

        if( Read::Has_Key_in_File(paramfile,"LSST_SN_dataset") ){
            string sn_datafile = Read::Read_String_from_File(paramfile,"LSST_SN_dataset");
            sn_data_temp.load(sn_datafile, arma::auto_detect);
        }

        if( sn_data_temp.n_cols < 3 ){
            throw runtime_error(" # of cols of your SNeIa data file is less than 3, check your data file or modify code here, catch me!");
        }

        if( sn_data_temp.n_rows > 0 ){

            sn_num  = sn_data_temp.n_rows;
            sn_z    = new double[sn_num];
            sn_mu   = new double[sn_num];
            sn_err  = new double[sn_num];

            for( int i=0; i<sn_num; ++i ){
                sn_z[i]     = sn_data_temp(i,0);
                sn_mu[i]    = sn_data_temp(i,1);
                sn_err[i]   = sn_data_temp(i,2);
            }

            if( sn_data_temp.n_cols == 4){
                sn_zerr = new double[sn_num];

                for( int i=0; i<sn_num; ++i )
                    sn_zerr[i] = sn_data_temp(i,3);
            }

            init_status = true;
        }

        if( Read::Has_Key_in_File(paramfile,"LSST_SN_dispersion") ){
            string sn_disp = Read::Read_String_from_File(paramfile,"LSST_SN_dispersion");
            arma::mat sd;
            sd.load(sn_disp,arma::auto_detect);

            if( sd.n_rows <= 0 || sd.n_cols < 2 ){
                throw runtime_error("# error when reading LSST_SN_dispersion file.\n");
            }
            else{
                muzperr_z = new double[sd.n_rows];
                muzperr   = new double[sd.n_rows];
                muzperr_zmin = 10;
                muzperr_zmax = -10;

                for( size_t i=0; i<sd.n_rows; ++i ){
                    muzperr_z[i] = sd(i,0);
                    muzperr[i] = sd(i,1);

                    if( muzperr_zmin > muzperr_z[i] ) muzperr_zmin = muzperr_z[i];
                    if( muzperr_zmax < muzperr_z[i] ) muzperr_zmax = muzperr_z[i];
                }

                acc_muzperr = gsl_interp_accel_alloc();
                spline_muzperr = gsl_spline_alloc(gsl_interp_cspline,sd.n_rows);
                gsl_spline_init(spline_muzperr,muzperr_z,muzperr,sd.n_rows);
            }
        }

        return init_status;
    }

    double photoz_mu_integrand( double z, void *param ){
        SNeIa_LSST *sn = static_cast<SNeIa_LSST*>(param);
        double zp = sn->zp_current;
        double dz = z-zp;
        double ns = gsl_spline_eval(sn->spline_ns, z, sn->acc_ns);
        // double ns = 1.0;
        // cout << "ns = " << ns << endl;
        double sgz = 0.02*(1.0+z);
        // double sgz = 0.01*(1.0+z);
        return ns*exp(-0.5*dz*dz/sgz/sgz);
    }

}
