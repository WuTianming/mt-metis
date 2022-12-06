/**
 * @file GraphFactory.cpp
 * @brief Implementation of GraphFactory for instantiating graph files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-05
 */



#include "GraphReaderFactory.hpp"
#include "MetisFile.hpp"
#include "SNAPFile.hpp"
#include "MatrixReaderFactory.hpp"
#include "MatrixGraphReader.hpp"




namespace WildRiver
{


/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::unique_ptr<IGraphReader> GraphReaderFactory::make(
    std::string const & name)
{
  std::unique_ptr<IGraphReader> file;
  // determine what type of reader to instantiate based on extension
  if (MetisFile::hasExtension(name)) {
    file.reset(new MetisFile(name));
  } else if (SNAPFile::hasExtension(name)) {
    file.reset(new SNAPFile(name));
  } else {
    // need to wrap it with an adapter
    std::unique_ptr<IMatrixReader> matPtr(MatrixReaderFactory::make(name));
    file.reset(new MatrixGraphReader(matPtr));
  }
 
  return file;
}




}
