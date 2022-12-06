/**
 * @file PlainVectorFile.cpp
 * @brief Implemenation of PlainVectorFile class.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-08
 */




#include <typeinfo>


#include "PlainVectorFile.hpp"




namespace WildRiver
{



/******************************************************************************
* HELPER FUNCTIONS ************************************************************
******************************************************************************/


namespace
{


bool isComment(
    std::string const & line) noexcept
{
  if (line.size() > 0) {
    switch (line[0]) {
      case '#':
      case '%':
      case '/':
        return true;
      default:
        return false;
    }
  } else {
    return false;
  }
}




}




/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/


bool PlainVectorFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".txt");
  extensions.push_back(".vec");
  extensions.push_back(".perm");
  extensions.push_back(".part");
  extensions.push_back(".cluster");

  return TextFile::matchExtension(f,extensions);
}



/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


bool PlainVectorFile::nextNoncommentLine(
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
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


PlainVectorFile::PlainVectorFile(
    std::string const & name) :
  m_file(name),
  m_size(NULL_DIM),
  m_buffer()
{
  // do nothing
}


PlainVectorFile::~PlainVectorFile()
{
  // do nothing
}


/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void PlainVectorFile::setSize(
    const ind_t size)
{
  m_size = size;
}



ind_t PlainVectorFile::getSize()
{
  if (m_size == NULL_DIM) {
    size_t nlines = 0;
    std::string m_buffer;

    if (!m_file.isOpenRead()) {
      m_file.openRead();
    }

    // cout non-comment lines
    while (nextNoncommentLine(m_buffer)) {
      ++nlines;
    }

    m_size = nlines;

    // go back to the beginning
    m_file.resetStream();
  }

  return m_size;
}


void PlainVectorFile::read(
    val_t * const vals,
    double * progress)
{
  if (!m_file.isOpenRead()) {
    m_file.openRead();
  }

  const size_t n = getSize();

  m_file.resetStream();

  ind_t const interval = n > 100 ? n / 100 : 1;
  double const increment = 1.0 / 100.0;
  
  for (size_t i = 0; i < n; ++i) {
    if (!nextNoncommentLine(m_buffer)) {
      throw EOFException("Hit end of file before getting next value");
    }

    if (typeid(val_t) == typeid(double) || typeid(val_t) == typeid(float)) {
      vals[i] = static_cast<val_t>(std::stod(m_buffer));
    } else {
      vals[i] = static_cast<val_t>(std::stoll(m_buffer));
    }

    if (progress != nullptr && i % interval == 0) {
      *progress += increment;
    }
  }
}


void PlainVectorFile::write(
    val_t const * const vals,
    double  * progress)
{
  if (m_size == NULL_DIM) {
    throw UnsetInfoException("Size of vector is not set before call to " \
        "write()");
  }

  if (!m_file.isOpenWrite()) {
    m_file.openWrite();
  }

  const size_t n = getSize();
  for (size_t i = 0; i < n; ++i) {
    m_file.setNextLine(std::to_string(vals[i]));
  }
}




}
