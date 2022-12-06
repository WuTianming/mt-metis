/**
 * @file VectorOutHandle.cpp
 * @brief Implementation of VectorOutHandle class.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2015-10-25
 */




#include "VectorOutHandle.hpp"
#include "VectorWriterFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


VectorOutHandle::VectorOutHandle(
    std::string const & name) :
  m_writer(VectorWriterFactory::make(name))
{
  // do nothing
}


VectorOutHandle::~VectorOutHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void VectorOutHandle::write(
    val_t const * const vals,
    ind_t const size,
    double * const progress)
{
  m_writer->setSize(size);
  m_writer->write(vals,progress);
}




}
