/**
 * @file MatrixWriterFactory.hpp
 * @brief Class for instantiating matrix writers.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-05
 */




#ifndef WILDRIVER_MATRIXWRITERFACTORY_HPP
#define WILDRIVER_MATRIXWRITERFACTORY_HPP




#include <memory>
#include <string>
#include "IMatrixWriter.hpp"



namespace WildRiver
{


class MatrixWriterFactory
{
  public:
    /**
     * @brief Instantiate a matrix file.
     *
     * @param fname The name of the file to open.
     * @param useAdapter Whether or not to use an adapter for a graph file.
     *
     * @return The instantiated object.
     *
     * @throw UnknownExtensionException If no class reading the specified file
     * type can be found.
     */
    static std::unique_ptr<IMatrixWriter> make(
        std::string const & name,
        bool useAdapter = true);
};




}




#endif
