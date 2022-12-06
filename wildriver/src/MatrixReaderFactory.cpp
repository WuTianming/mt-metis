/**
 * @file MatrixReaderFactory.cpp
 * @brief Implementation of the MatrixReaderFactory class.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-05
 */




#include "MatrixReaderFactory.hpp"
#include "CSRFile.hpp"
#include "MatrixMarketFile.hpp"
#include "MetisFile.hpp"
#include "SNAPFile.hpp"
#include "GraphReaderFactory.hpp"
#include "GraphMatrixReader.hpp"




namespace WildRiver
{


/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::unique_ptr<IMatrixReader> MatrixReaderFactory::make(
    std::string const & name)
{
  std::unique_ptr<IMatrixReader> file;

  // determine what type of reader to instantiate based on extension
  if (CSRFile::hasExtension(name)) {
    file.reset(new CSRFile(name));
  } else if (MatrixMarketFile::hasExtension(name)) {
    file.reset(new MatrixMarketFile(name));
  } else if (MetisFile::hasExtension(name)) {
    std::unique_ptr<IGraphReader> graphPtr(GraphReaderFactory::make(name));
    file.reset(new GraphMatrixReader(graphPtr));
  } else if (SNAPFile::hasExtension(name)) {
    std::unique_ptr<IGraphReader> graphPtr(GraphReaderFactory::make(name));
    file.reset(new GraphMatrixReader(graphPtr));
  } else {
    throw UnknownExtensionException(std::string("Unknown filetype: ") + name);
  }
 
  return file;
}




}
