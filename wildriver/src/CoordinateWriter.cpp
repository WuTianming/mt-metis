/**
* @file CoordinateWriter.cpp
* @brief Implementation of the CoordinateWriter class.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2018
* @version 1
* @date 2018-04-07
*/




#include "CoordinateWriter.hpp"
#include "TextFile.hpp"



namespace WildRiver
{



/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/

CoordinateWriter::CoordinateWriter(
    TextFile * const file) :
  m_numWrittenRows(0),
  m_file(file)
{
  // do nothing
}



/******************************************************************************
* PUBLIC METHODS **************************************************************
******************************************************************************/

void CoordinateWriter::writeHeader(
      dim_t,
      dim_t,
      ind_t)
{
  // no op
}



void CoordinateWriter::setNextRow(
    dim_t const numNonZeros,
    dim_t const * const columns,
    val_t const * const values)
{
  std::string const rowStr(std::to_string(m_numWrittenRows) + \
      std::string(" "));
  for (ind_t j = 0; j < numNonZeros; ++j) {
    if (values != nullptr) {
      m_file->setNextLine(rowStr + " " + std::to_string(columns[j]) + " " + \
          std::to_string(values[j]));
    } else {
      m_file->setNextLine(rowStr + " " + std::to_string(columns[j]));
    }
  }

  ++m_numWrittenRows;
}




}
