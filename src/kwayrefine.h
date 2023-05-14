/**
 * @file kwayrefine.h
 * @brief KWay refinement function prototypes.
 * @author Dominique LaSalle <mtmetis@domnet.org>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-19
 */




#ifndef MTMETIS_KWAYREFINE_H
#define MTMETIS_KWAYREFINE_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"
#include "kwinfo.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define par_kwayrefine MTMETIS_par_kwayrefine
#define par_kwayrefine_chunk MTMETIS_par_kwayrefine_chunk
/**
* @brief Parallel kway-refinement.
*
* @param ctrl The control strucutre.
* @param graph The graph who's partition to refine.
* @param kwinfo The uncoarsening information structure. 
*
* @return Total of moved vertices.
*/
vtx_type par_kwayrefine(
    ctrl_type * const ctrl, 
    graph_type * const graph,
    kwinfo_type * kwinfo);

vtx_type par_kwayrefine_chunk(
    ctrl_type * const ctrl, 
    graph_type * const graph,
    kwinfo_type * kwinfo);




#endif
