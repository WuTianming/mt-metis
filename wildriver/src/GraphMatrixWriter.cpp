/**
 * @file GraphMatrixWriter.cpp
 * @brief Implementation of GraphMatrixWriter.
 * @author Dominique LaSalle <dominique@solidlake.com>
 * Copyright 2017, Solid Lake LLC
 * @version 1
 * @date 2017-07-25
 */




#include "GraphMatrixWriter.hpp"
#include "Exception.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


GraphMatrixWriter::GraphMatrixWriter(
    std::unique_ptr<IGraphWriter>& writer) :
  m_writer(std::move(writer))
{
  // do nothing
}


GraphMatrixWriter::~GraphMatrixWriter()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphMatrixWriter::setInfo(
    dim_t const numRows,
    dim_t const numCols,
    ind_t const nnz)
{
  if (numRows != numCols) {
    throw BadParameterException("Only square matrices can be written as " \
        "graphs");
  }

  m_writer->setInfo(numRows, nnz, 0, true);
}


void GraphMatrixWriter::write(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  m_writer->write(rowptr, rowind, nullptr, rowval);
}




}
