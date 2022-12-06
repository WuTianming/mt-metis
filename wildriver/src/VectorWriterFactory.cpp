/**
 * @file VectorReaderFactory.cpp
 * @brief Implemnetation of class for 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-07
 */




#include "Exception.hpp"
#include "VectorWriterFactory.hpp"
#include "PlainVectorFile.hpp"




namespace WildRiver
{



/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::unique_ptr<IVectorWriter> VectorWriterFactory::make(
    std::string const & name)
{
  std::unique_ptr<IVectorWriter> file;

  // determine what type of reader to instantiate based on extension
  if (PlainVectorFile::hasExtension(name)) {
    file.reset(new PlainVectorFile(name));
  } else {
    throw UnknownExtensionException(std::string("Unknown filetype: ") + name);
  }
 
  return file;
}




}
