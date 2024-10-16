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
struct dimension_parameters {
    int nx;
    int ny;
    int npts;
    int mx;
    int my;
    int mz;
};
typedef struct dimension_parameters dim_param;

struct embedding_parameters {
    int mx;
    int my;
    int mz;
    int stride_x;
    int stride_y;
    int stride_z;
    int lag;
};
typedef struct embedding_parameters embed_param; 

/**********************************************************************************************************/
extern const int k_default;
extern int tree_k_max;      // for safety checks
#endif

