/**
 * @file VectorWriterFactory.hpp
 * @brief Class for instantiating vector writers. 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-07
 */




#ifndef WILDRIVER_VECTORWRITERFACTORY_HPP
#define WILDRIVER_VECTORWRITERFACTORY_HPP




#include <memory>
#include <string>

#include "IVectorFile.hpp"




namespace WildRiver
{


class VectorWriterFactory
{
  public:
    /**
     * @brief Allocate a new vector writer subclass based on teh file
     * extension.
     *
     * @param name The filename/path to open.
     *
     * @return The newly opened vector writer.
     */
    static std::unique_ptr<IVectorWriter> make(
        std::string const & name);


};




}




#endif
