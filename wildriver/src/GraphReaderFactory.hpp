/**
 * @file GraphFactory.hpp
 * @brief Class for instantiating graph files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-05
 */




#ifndef WILDRIVER_GRAPHREADERFACTORY_HPP
#define WILDRIVER_GRAPHREADERFACTORY_HPP



#include <memory>

#include "IGraphReader.hpp"



namespace WildRiver
{


class GraphReaderFactory
{
  public:
    /**
     * @brief Allocate a new graph file subclass, based on the file extension.
     * The returned pointer must be delete'd by the caller.
     *
     * @param fname The filename/path.
     *
     * @return A pointer to the new graph file instantion.
     *
     * @throw UnknownExtensionException If no class reading the specified file
     * type can be found.
     */
    static std::unique_ptr<IGraphReader> make(
        std::string const & name);




};




}




#endif
