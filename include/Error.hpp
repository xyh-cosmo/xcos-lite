#ifndef __ICOSMO_ERROR__
#define __ICOSMO_ERROR__

	#ifndef icosmo_error
		#define print_icosmo_error  {std::cout << "\n****    ICOSMO ERROR    ****\n"\
		<< "\n--- File Name: " << __FILE__ << "\n"\
		<< "--- Line #: " << __LINE__ << "\n"\
		<< "--- Function Name: " << __FUNCTION__ << "\n";}
	#endif	//	icosmo_error

#endif	//	__ICOSMO_ERROR__


#ifndef __ICOSMO_WARNING__
#define __ICOSMO_WARNING__

	#ifndef icosmo_warning
		#define print_icosmo_warning {std::cout << "\n****    ICOSMO WARNING    ****\n"\
		<< "--- File Name: " << __FILE__ << "\n"\
		<< "--- Line #: " << __LINE__ << "\n"\
		<< "--- Function Name: " << __FUNCTION__ << "\n";}
	#endif	//	icosmo_warning


#endif	//	__ICOSMO_WARNING__
