/**
 * @file GraphMatrixReader.cpp
 * @brief Implementation of GraphMatrixReader.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "GraphMatrixReader.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


GraphMatrixReader::GraphMatrixReader(
    std::unique_ptr<IGraphReader>& reader) :
  m_reader(std::move(reader))
{
  // do nothing
}


GraphMatrixReader::~GraphMatrixReader()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphMatrixReader::getInfo(
    dim_t & nrows,
    dim_t & ncols,
    ind_t & nnz)
{
  int nvwgt;
  bool ewgts;
  m_reader->getInfo(nrows, nnz, nvwgt, ewgts);

  ncols = nrows;
}


void GraphMatrixReader::read(
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval,
    double * progress)
{
  m_reader->read(rowptr, rowind, nullptr, rowval, progress);
}




}
