/**
 * @file MatrixWriterFactory.cpp
 * @brief Implementation of the MatrixWriterFactory class.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-05
 */




#include "MatrixWriterFactory.hpp"
#include "CSRFile.hpp"
#include "MatrixMarketFile.hpp"
#include "MetisFile.hpp"
#include "GraphMatrixWriter.hpp"
#include "GraphWriterFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::unique_ptr<IMatrixWriter> MatrixWriterFactory::make(
    std::string const & name,
    bool useAdapter)
{
  std::unique_ptr<IMatrixWriter> file;

  // determine what type of reader to instantiate based on extension
  if (CSRFile::hasExtension(name)) {
    file.reset(new CSRFile(name));
  } else if (MatrixMarketFile::hasExtension(name)) {
    file.reset(new MatrixMarketFile(name));
  } else if (useAdapter) {
    std::unique_ptr<IGraphWriter> graphPtr(GraphWriterFactory::make(name, \
        false));
    file.reset(new GraphMatrixWriter(graphPtr));
  } else {
    throw UnknownExtensionException(std::string("Unknown filetype: ") + name);
  }
 
  return file;
}




}
