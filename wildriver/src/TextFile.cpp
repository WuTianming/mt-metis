/**
 * @file TextFile.cpp
 * @brief Abstract class for text files (both graphs and matrices).
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "TextFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


enum {
  FILE_STATE_UNOPENED,
  FILE_STATE_READ,
  FILE_STATE_WRITE
};




/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/


bool TextFile::matchExtension(
    std::string const & f,
    std::vector<std::string> const & extensions)
{
  for (std::string const & ext : extensions) {
    size_t len = ext.size();

    if (f.compare(f.size()-len,len,ext) == 0) {
      return true;
    }
  }

  return false;
}



/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


TextFile::TextFile(
    std::string const & name) :
  m_state(FILE_STATE_UNOPENED),
  m_currentLine(0),
  m_name(name),
  m_stream()
{
  // throw exceptions when things go wrong
  m_stream.exceptions(std::fstream::failbit | std::fstream::badbit);
}


TextFile::~TextFile()
{
  m_stream.close();
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void TextFile::openWrite()
{
  if (m_state != FILE_STATE_UNOPENED) {
    throw BadFileStateException(std::string("Attempting to re-open file '") + \
        m_name + std::string("' for writing."));
  }

  m_stream.open(getFilename(),std::fstream::out | std::fstream::trunc);

  if (!m_stream.good()) {
    throw BadFileException(std::string("Failed to open file '") + \
        m_name + std::string("'"));
  }

  m_state = FILE_STATE_WRITE;
}


void TextFile::openRead()
{
  if (m_state != FILE_STATE_UNOPENED) {
    throw BadFileStateException(std::string("Attempting to re-open file '") + \
        m_name + std::string("' for reading."));
  }

  m_stream.open(getFilename(),std::fstream::in);

  if (!m_stream.good()) {
    throw BadFileException(std::string("Failed to open file '") + \
        m_name + std::string("'"));
  }

  m_state = FILE_STATE_READ;
}


bool TextFile::isOpenWrite() const
{
  return m_state == FILE_STATE_WRITE;
}


bool TextFile::isOpenRead() const
{
  return m_state == FILE_STATE_READ;
}


void TextFile::resetStream()
{
  m_stream.clear();
  m_stream.seekg(0,std::ifstream::beg);
  m_currentLine = 0;
}


bool TextFile::nextLine(
    std::string & line)
{
  try {
    std::getline(m_stream,line);
  } catch (std::fstream::failure const & e) {
    // at eof
    return false;
  }

  if (line.size() > 0 || !m_stream.eof()) {
    ++m_currentLine;
    return true;
  } else {
    return false;
  }
}


void TextFile::setNextLine(
    std::string const & line)
{
  m_stream << line << "\n";
}




}
