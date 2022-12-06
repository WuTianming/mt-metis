/**
 * @file VectorReaderFactory.hpp
 * @brief Class for instantiating vector readers. 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-07
 */




#ifndef WILDRIVER_VECTORREADERFACTORY_HPP
#define WILDRIVER_VECTORREADERFACTORY_HPP




#include <memory>
#include <string>

#include "IVectorReader.hpp"




namespace WildRiver
{


class VectorReaderFactory
{
  public:
    /**
     * @brief Allocate a new vector readder subclass based on teh file
     * extension.
     *
     * @param name The filename/path to open.
     *
     * @return The newly opened vector reader.
     */
    static std::unique_ptr<IVectorReader> make(
        std::string const & name);


};




}




#endif
