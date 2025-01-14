# file commons.pxd
#
# this is a copy of functions/definitions that we want to see (and hence use) in Python
#
# 2019-01-24: cleaned up version
# 2021-01-19: removed "bins" functions
# 2021-12-19: forked out of entropy.pxd

# import definitions from C : 
#cdef extern from "library_commons.h":
#	void ANN_set_verbosity    		(int level)	
	

cdef extern from "verbosity.h":
	int lib_verbosity          
	int lib_warning_level     


cdef extern from "ANN_threads.h":
	int 	USE_PTHREAD
	void 	set_multithreading_state(int do_mp)
	int  	get_cores_number  		(int get_what)
	void 	set_cores_number  		(int n)
	int 	adapt_cores_number		(int npts_eff)


