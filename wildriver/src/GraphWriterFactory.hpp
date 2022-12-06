/**
 * @file GraphWriterFactory.hpp
 * @brief Class for instantiating graph writers.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-05
 */




#ifndef WILDRIVER_GRAPHWRITERFACTORY_HPP
#define WILDRIVER_GRAPHWRITERFACTORY_HPP



#include <memory>

#include "IGraphWriter.hpp"



namespace WildRiver
{


class GraphWriterFactory
{
  public:
    /**
     * @brief Allocate a new graph file subclass, based on the file extension.
     * The returned pointer must be delete'd by the caller.
     *
     * @param fname The filename/path.
     * @param useAdapter Whether or not to use an adapter for a matrix file. 
     *
     * @return A pointer to the new graph file instantion.
     *
     * @throw UnknownExtensionException If no class reading the specified file
     * type can be found.
     */
    static std::unique_ptr<IGraphWriter> make(
        std::string const & name,
        bool useAdapter = true);




};




}




#endif
