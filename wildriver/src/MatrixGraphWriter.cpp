/**
 * @file MatrixGraphWriter.cpp
 * @brief An adapter class for treating matrices as graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "MatrixGraphWriter.hpp"
#include "Exception.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixGraphWriter::MatrixGraphWriter(
    std::unique_ptr<IMatrixWriter>& writer) :
  m_writer(std::move(writer))
{
  // do nothing
}


MatrixGraphWriter::~MatrixGraphWriter()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixGraphWriter::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int const nvwgt,
    bool const ewgts)
{
  if (nvwgt != 0) {
    throw BadParameterException("MatrixGraphWriter cannot handle vertex " \
        "weights");
  }

  // TODO: add an internal state for edge weights, such that they can be
  // surpessed via this function.

  m_writer->setInfo(nvtxs,nvtxs,nedges);
}


void MatrixGraphWriter::write(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  m_writer->write(xadj,adjncy,adjwgt);
}




}
