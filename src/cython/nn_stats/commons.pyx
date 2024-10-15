#cython: language_level=3, boundscheck=False, cython.wraparound=False
# Cython wrappers for C functions
#
# 2018-04-10: enhanced memory mapping from C to Python
# 2021-12-02: multi-threaded entropy algorithm

# 2022-05-17: OK
def set_verbosity(int level=1):
    """
    sets the verbosity level of the library
    
    :param verbosity: an integer that indicates the verbosity level of the library (default=1)
    :returns: no output
    
    verbosity can be:
      - <0 : no messages, even if an error is encountered (not recommended!)
      -  0 : messages only if an error is encountered
      -  1 : messages for errors and important warnings only (default)
      -  2 : more warnings and/or more details on the warnings
      -  ...
        
    """
    commons.lib_verbosity=level
    
def get_verbosity():
    """
    gets the current verbosity level of the library
    
    :param none:
    :returns: no output values, but a message indicating the verbosity level is printed in the console.
    
    verbosity level explanation:
      - <0 : no messages, even if an error is encountered (not recommended!)
      -  0 : messages only if an error is encountered
      -  1 : messages for errors and important warnings only (default)
      -  2 : more warnings and/or more details on the warnings
      -  ...
        
    """
    print("verbosity level", commons.lib_verbosity)




# 2021-12-02
def multithreading(do_what="info", int nb_cores=0):
     """ multithreading(do_what=["info", "auto", "single"], n_cores=0)

     Selects the multi-threading schemes and eventually sets the number of threads to use
     
     :param do_what: either "info" or "auto" or "single", see below.
     :param nb_cores: an integer; if specified and positive, then the provided number nb_cores of threads will be used. (default=0 for "auto") 
     :returns: no output.
     
     The parameter do_what can be chosen as follows:
        - "info": nothing is done, but informations on the current multithreading state are displayed.
        - "auto": then the optimal (self-adapted to largest) number of threads will be used) (default).
        - "single": then the algorithms will run single-threaded (no multithreading) 
     
     if nb_cores (a positive number) is specified, then the provided number nb_cores of threads will be used.
     """
     
     if (do_what=="info"): 
        print("currently using %d out of %d cores available"
                    %(commons.get_cores_number(0x0001), commons.get_cores_number(0x0010)), end="")
        if (commons.get_cores_number(0x0001)==-1): 
                print(" (-1 means largest number available, so %d here)" %commons.get_cores_number(0x0010))
        else:   print()
        return
     elif (do_what=="auto"):                    # autodetect/automatic (default)
        commons.set_multithreading_state(2)
        commons.set_cores_number(-1)
     elif (do_what=="single"):                  # single threading
        commons.set_multithreading_state(0)     # no multithreading
        commons.set_cores_number(1)
     elif ( isinstance(do_what, (int, PNP.int32)) and (do_what>0) ):      # imposed nb of threads
# see https://stackoverflow.com/questions/21851985/difference-between-np-int-np-int-int-and-np-int-t-in-cython
        commons.set_multithreading_state(1) 
        commons.set_cores_number(do_what)
     else: raise ValueError("invalid parameter value")
     return


# 2021-12-20
def get_threads_number():
     """
     returns the current number of threads.
     
     :param none:
     :returns: an integer, the current number of threads used.
     """
     
     if (commons.USE_PTHREAD>0): return(commons.get_cores_number(0x0001))
     return(0)

