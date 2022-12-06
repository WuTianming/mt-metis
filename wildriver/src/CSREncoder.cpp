/**
 * @file CSREncoder.cpp
 * @brief Implementation of CSREncoder.
 * @author Dominique LaSalle <dominique@solidlake.com>
 * Copyright 2017, Solid Lake LLC
 * @version 1
 * @date 2017-07-26
 */




#include <cassert>


#include "CSREncoder.hpp"




namespace WildRiver
{



/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


CSREncoder::CSREncoder(
    IRowMatrixWriter * const writer) :
  m_infoSet(false),
  m_numRows(NULL_DIM),
  m_numCols(NULL_DIM),
  m_nnz(NULL_IND),
  m_writer(writer)
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void CSREncoder::setInfo(
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz)
{
	m_numRows = nrows;
  m_numCols = ncols;
  m_nnz = nnz;

  assert(m_infoSet == false);

  assert(m_numRows != NULL_DIM); 
  assert(m_numCols != NULL_DIM); 
  assert(m_nnz != NULL_IND);

  m_infoSet = true;

  m_writer->writeHeader(m_numRows, m_numCols, m_nnz);
}


void CSREncoder::write(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  for (dim_t i = 0; i < m_numRows; ++i) {
    ind_t const size = rowptr[i+1] - rowptr[i];
    m_writer->setNextRow(size, rowind+rowptr[i], rowval+rowptr[i]);
  }
}


}
