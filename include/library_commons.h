/*
 *  library_commons.h
 *  
 *
 *  Created by Nicolas Garnier on 2020/03/05.
 *  Copyright 2012-2022 ENS-Lyon - CNRS. All rights reserved.
 *
 *  
 */
#ifndef _LIBRARY_COMMONS_H
#define _LIBRARY_COMMONS_H
/**********************************************************************************************************/
#ifdef NAN
    #define my_NAN NAN
#else
    #define my_NAN 0.0
#endif

/**********************************************************************************************************/
extern const int k_default;
extern int tree_k_max;      // for safety checks
#endif

