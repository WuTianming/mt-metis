/**
* @file MatrixMarketFile.hpp
* @brief Class for reading/writing matrix market files. 
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2017
* @version 1
* @date 2017-07-28
*/




#ifndef WILDRIVER_MATRIXMARKETFILE_HPP
#define WILDRIVER_MATRIXMARKETFILE_HPP




#include <memory>


#include "IMatrixReader.hpp"
#include "IMatrixWriter.hpp"
#include "TextFile.hpp"




namespace WildRiver
{


/**
* @brief A class for reading and writing matrices stored in MatrixMarket
* format (https://math.nist.gov/MatrixMarket/formats.html).
*/
class MatrixMarketFile :
    public IMatrixReader,
    public IMatrixWriter
{
  public:
    /**
     * @brief Check if the given filename matches an extension for this
     * filetype.
     *
     * @param f The filename.
     *
     * @return True if the extension matches this filetype.
     */
    static bool hasExtension(
        std::string const & f);


    /**
     * @brief Create a new MatrixMarketFile for reading and writing.
     *
     * @param fname The filename/path.
     */
    MatrixMarketFile(
        std::string const & fname);


    /**
     * @brief Close file and free any memory.
     */
    virtual ~MatrixMarketFile();


    /**
     * @brief Get the number of rows, columns, and non-zeros in the matrix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz THe number of non-zeros.
     */
    virtual void getInfo(
        dim_t & nrows,
        dim_t & ncols,
        ind_t & nnz) override;


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
     * @brief Get the sparse matrix in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix
     *
     * |rowptr| = nrows + 1
     * |rowind| = nnz
     * |rowval| = nnz
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     * @param progress The variable to update as the matrix is loaded (may be
     * null).
     */
    virtual void read(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress) override;


    /**
     * @brief Write the given CSR structure to teh file. The information for
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


    /**
    * @brief Read the header data from the file.
    */
    virtual void readHeader();


    /**
     * @brief Write the header data to the file.
     */
    virtual void writeHeader();


    /**
     * @brief Read in the matrix in general coordinate format.
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     * @param progress The variable to update as the matrix is loaded (may be
     * null).
    */
    virtual void readCoordinates(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress);


    /**
     * @brief Read in the matrix in symmetric coordinate format.
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     * @param progress The variable to update as the matrix is loaded (may be
     * null).
    */
    virtual void readSymmetricCoordinates(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress);


    /**
    * @brief Read in the matrix in array format.
    */
    virtual void readArray();


    /**
     * @brief Write the matrix out in coordinate format. 
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     */
    virtual void writeCoordinates(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval);


  private:
    /**
     * @brief Whether the info has been written/read.
     */
    bool m_infoSet;

    /**
    * @brief The number of rows in the underyling matrix.
    */
    dim_t m_nrows;

    /**
    * @brief The number of cols in the underlying matrix.
    */
    dim_t m_ncols;

    /**
    * @brief The number of non-zeros in the underlying matrix.
    */
    ind_t m_nnz;

    /**
     * @brief The last read line from the text file.
     */
    std::string m_line;

    /**
     * @brief The underlying text file.
     */
    TextFile m_file;

    /**
     * @brief What type of structure is stored (matrix, vector, etc.).
     */
    int m_entity;

    /**
     * @brief Whether the martix is coordinate or array stored.
     */
    int m_format;

    /**
     * @brief What type of values are in the matrix (real, complex, integer,
     * pattern).
     */
    int m_type;

    /**
     * @brief Whether or not the matrix is symmetric.
     */
    bool m_symmetric;


    /**
    * @brief Get the next non-comment line from the file.
    *
    * @param line The line to fill.
    *
    * @return True if the line was filled.
    */
    bool nextNoncommentLine(
        std::string & line);




};




}




#endif
