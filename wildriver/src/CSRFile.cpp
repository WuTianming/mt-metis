/**
 * @file CSRFile.cpp
 * @brief Class functions for reading and writing metis graph files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include <cstddef>
#include <cassert>


#include <sstream>


#include "CSRFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


namespace
{


static size_t const BUFFER_SIZE = 4096;


}


std::string const CSRFile::NAME = "CSR";




/******************************************************************************
* HELPER FUNCTIONS ************************************************************
******************************************************************************/


namespace
{


/**
* @brief Determine the given line is a comment.
*
* @param line The line.
*
* @return True if the line is a comment.
*/
bool isComment(
    std::string const & line) noexcept
{
  return line.size() > 0 && line[0] == '#';
}




}




/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/


bool CSRFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".csr");

  return TextFile::matchExtension(f,extensions);
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


bool CSRFile::nextNoncommentLine(
    std::string & line)
{
  do {
    if (!m_file.nextLine(line)) {
      return false;
    }
  } while (isComment(line));

  return true;
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


CSRFile::CSRFile(
    std::string const & fname) :
  m_oneBased(false),
  m_line(BUFFER_SIZE,'\0'),
  m_file(fname),
  m_decoder(nullptr),
  m_encoder(nullptr)
{
  // do nothing
}


CSRFile::~CSRFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void CSRFile::getInfo(
    dim_t & nrows,
    dim_t & ncols,
    ind_t & nnz)
{
  if (m_decoder.get() == nullptr) {
    // make our decoder if we haven't yet
    m_decoder.reset(new CSRDecoder(this));
  }

  m_decoder->getInfo(nrows, ncols, nnz);
}


void CSRFile::read(
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval,
    double * progress)
{
  if (m_decoder.get() == nullptr) {
    throw UnsetInfoException("Cannot call read() before calling getInfo()");
  }

  m_decoder->read(rowptr, rowind, rowval, progress);
}


void CSRFile::setInfo(
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz)
{
  if (m_encoder.get() == nullptr) {
    // make our encoder if we haven't yet
    m_encoder.reset(new CSREncoder(this));
  }

  m_encoder->setInfo(nrows, ncols, nnz);
}


void CSRFile::write(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  if (m_encoder.get() == nullptr) {
    throw UnsetInfoException("Cannot call write() before calling setInfo()");
  }

  m_encoder->write(rowptr, rowind, rowval);
}


void CSRFile::readHeader(
    dim_t & numRows,
    dim_t & numCols,
    ind_t & nnz)
{
  // open our file for reading
  m_file.openRead();

  // go to the start of the file
  m_file.resetStream();

  dim_t minColumn = NULL_DIM;
  dim_t maxColumn = 0;

  numRows = 0;
  nnz = 0;

  while (nextNoncommentLine(m_line)) {
    char const * const lineEnd = m_line.data() + m_line.length();

    char * sptr;
    char * eptr = (char*)m_line.data();

    // parse the line
    dim_t degree = 0;
    while (true) {
      sptr = eptr;
      dim_t col = static_cast<dim_t>(std::strtoull(sptr,&eptr,10));
      if (eptr == sptr) {
        // nothing left to read
        break;
      }

      if (col > maxColumn) {
        maxColumn = col;
      }
      if (minColumn == NULL_DIM || col < minColumn) {
        minColumn = col;
      }

      // skip value without converting
      do {
        ++eptr;
      } while (eptr < lineEnd && *eptr != ' ');
      ++degree;
    }

    nnz += degree;
    ++numRows;
  }

  // decide where we started counting
  if (minColumn > 0) {
    // 1-based
    numCols = maxColumn;
    m_oneBased = true;
  } else {
    // 0-based
    m_oneBased = false;
    numCols = maxColumn+1;
  }

  // go back to the start
  m_file.resetStream();
}


void CSRFile::getNextRow(
    dim_t * const numNonZeros,
    dim_t * const columns,
    val_t * const values)
{
  if (!nextNoncommentLine(m_line)) {
    throw BadFileException(std::string("Unexcepted end of file at line") + \
        std::to_string(m_file.getCurrentLine()));
  }

  // TOOD: don't cast constness away
  char * sptr;
  char * eptr = (char*)m_line.data();

  const dim_t offset = m_oneBased ? 1 : 0;

  dim_t degree = 0;
  dim_t col;
  val_t val;

  // Loop through row until we streamed to the end
  while (true) {
    sptr = eptr;
    col = static_cast<dim_t>(std::strtoull(sptr,&eptr,10));
    if (eptr == sptr) {
      // nothing left to read
      break;
    }

    sptr = eptr;

    val = std::strtod(sptr,&eptr);
    if (eptr == sptr) {
      throw BadFileException(std::string("Failed to read column on "
            "line ") + std::to_string(m_file.getCurrentLine()));
    }

    columns[degree] = col - offset;
    if (values) {
      values[degree] = val;
    }
    ++degree;
  }

  *numNonZeros = degree;
}


void CSRFile::writeHeader(
    dim_t,
    dim_t,
    ind_t)
{
  // open for writing
  m_file.openWrite();
}


void CSRFile::setNextRow(
    dim_t numNonZeros,
    dim_t const * const columns,
    val_t const * const values)
{
  std::stringstream streamBuffer;

  const dim_t offset = m_oneBased ? 1 : 0;

  if (numNonZeros > 0) {
    const size_t last = numNonZeros - 1;
    for (size_t i = 0; i < last; ++i) {
      streamBuffer << (columns[i]+offset) << " " << values[i] << " ";
    }
    streamBuffer << (columns[last]+offset) << " " << values[last];
  }

  m_file.setNextLine(streamBuffer.str());
}




}
