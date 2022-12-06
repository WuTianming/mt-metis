/**
 * @file MatrixReaderFactory.hpp
 * @brief Class for instantiating matrix readers.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-05
 */




#ifndef WILDRIVER_MATRIXREADERFACTORY_HPP
#define WILDRIVER_MATRIXREADERFACTORY_HPP




#include <memory>
#include <string>
#include "IMatrixReader.hpp"



namespace WildRiver
{


class MatrixReaderFactory
{
  public:
    /**
     * @brief Instantiate a matrix file.
     *
     * @param fname The name of the file to open.
     *
     * @return The instantiated object.
     *
     * @throw UnknownExtensionException If no class reading the specified file
     * type can be found.
     */
    static std::unique_ptr<IMatrixReader> make(
        std::string const & name);
};




}




#endif
