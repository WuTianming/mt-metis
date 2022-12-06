/**
* @file GraphMatrixWriter.hpp
* @brief Adapter for writing graphs as matrices.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2017, Solid Lake LLC
* @version 1
* @date 2017-07-24
*/




#ifndef WILDRIVER_GRAPHMATRIXWRITER_HPP
#define WILDRIVER_GRAPHMATRIXWRITER_HPP




#include <vector>
#include <string>
#include <memory>

#include "IGraphWriter.hpp"
#include "IMatrixWriter.hpp"




namespace WildRiver
{


class GraphMatrixWriter :
    public IMatrixWriter
{
  public:
    /**
    * @brief Create an adapter.
    *
    * @param writer The writer to adapt.
    */
    GraphMatrixWriter(
        std::unique_ptr<IGraphWriter>& writer);


    /**
    * @brief Destructor.
    */
    virtual ~GraphMatrixWriter();


    /**
     * @brief Set the matrix information for this file.
     *
     * @param nrows The number of rows in the matrix.
     * @param ncols The number of columns in the matrix.
     * @param nnz The number of non-zeroes in the matrix.
     */
    virtual void setInfo(
        dim_t nrows,
        dim_t ncols,
        ind_t nnz) override;


    /**
     * @brief Write the given CSR structure to the file. The information for
     * the matrix must already be set.
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     */
    virtual void write(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval) override;




  private:
    std::unique_ptr<IGraphWriter> m_writer;




};




}




#endif
