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
    
def get_verbosity(verb=1):
    """
    gets the current verbosity level of the library
    
    :param verb: local verbosity of this very function:
      - 0 for silent running; 
      - 1 for printing the verbosity level in the console
                 
    :returns: the verbosity level.
    
    verbosity level explanation:
      - <0 : no messages, even if an error is encountered (not recommended!)
      -  0 : messages only if an error is encountered
      -  1 : messages for errors and important warnings only (default)
      -  2 : more warnings and/or more details on the warnings
      -  ...
        
    """
    if (verb): print("verbosity level", commons.lib_verbosity)
    return(commons.lib_verbosity)



def set_kernel(int kernel_type=0, double scale=1.0):
   """
   sets the kernel ot use for local averagings, and its scale
    
    :param kernel_type: an integer that indicates the kernel type (default: no kernel)
    :param scale: a positive real number that defines the (observation) scale for the kernel
    :returns: no output
    
    kernel_type can be:
      -  0 : no kernel (should be the faster method)
      -  1 : brickwall, all points have the same weight (same as no kernel)
      -  2 : Gaussian
      -  3 : triangle
      -  4 : quartique
      -  5 : Epanechnikov
      -  6 : triweight
      -  7 : tricube
      -  8 : cosine
      -  9 : exponential
      - 10 : logistic
      - 11 : sigmoid
      - 12 : Silverman
    """
   commons.select_kernel(kernel_type, scale)

def get_kernel(verb=1):
    """
    gets the current kernel type used for averaging, and its associated (observation) scale
    
    :param verb: local verbosity of this very function:
      - 0 for silent running (no output in the console); 
      - 1 for printing the kernel type in the console
                 
    :returns: the kernel type.
    
    kernel_type explanation: see the help of the function "set_kernel()"
    """
    if (verb): 
      print("kernel type", commons.current_kernel_type, end=" ")
      if   (commons.current_kernel_type==0):  print("(no kernel)", end=" ")
      elif (commons.current_kernel_type==1):  print("(brickwall)", end=" ")
      elif (commons.current_kernel_type==2):  print("(Gaussian)", end=" ")
      elif (commons.current_kernel_type==3):  print("(triangle)", end=" ")
      elif (commons.current_kernel_type==4):  print("(quartique)", end=" ")
      elif (commons.current_kernel_type==5):  print("(Epanechnikov)", end=" ")
      elif (commons.current_kernel_type==6):  print("(triweight)", end=" ")
      elif (commons.current_kernel_type==7):  print("(tricube)", end=" ")
      elif (commons.current_kernel_type==8):  print("(cosine)", end=" ")
      elif (commons.current_kernel_type==9):  print("(exponential)", end=" ")
      elif (commons.current_kernel_type==10): print("(logistic)", end=" ")
      elif (commons.current_kernel_type==11): print("(sigmoid)", end=" ")
      elif (commons.current_kernel_type==12): print("(Silverman)", end=" ")

      if   (commons.current_kernel_type<2):  print("(scale unused)")
      else: print("with scale %2.2f" %commons.current_obs_scale)
    return(commons.current_kernel_type, commons.current_obs_scale)



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

