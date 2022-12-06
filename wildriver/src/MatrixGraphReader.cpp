/**
 * @file MatrixGraphReader.cpp
 * @brief An adapter class for treating matrices as graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "MatrixGraphReader.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixGraphReader::MatrixGraphReader(
    std::unique_ptr<IMatrixReader>& reader) :
  m_reader(std::move(reader))
{
  // do nothing
}


MatrixGraphReader::~MatrixGraphReader()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/

void MatrixGraphReader::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  dim_t nrows, ncols;
  ind_t nnz;

  m_reader->getInfo(nrows,ncols,nnz);

  nvtxs = nrows;
  nedges = nnz;
  
  // matrices never have vertex weights and always have edge weights
  nvwgt = 0;
  ewgts = true;
}


void MatrixGraphReader::read(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt,
    double * progress)
{
  dim_t nrows, ncols;
  ind_t nnz;

  m_reader->getInfo(nrows,ncols,nnz);

  // fill vertex weights with ones if requested
  if (vwgt) {
    for (dim_t i=0; i<nrows; ++i) {
      vwgt[i] = 1;
    }
  }

  m_reader->read(xadj,adjncy,adjwgt,progress);
}




}
