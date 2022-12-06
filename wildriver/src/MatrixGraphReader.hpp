/**
 * @file MatrixGraphReader.hpp
 * @brief An adapter class for treating matrices as graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_MATRIXGRAPHREADER_HPP
#define WILDRIVER_MATRIXGRAPHREADER_HPP




#include <vector>
#include <string>
#include <memory>

#include "IGraphReader.hpp"
#include "IMatrixReader.hpp"




namespace WildRiver
{


class MatrixGraphReader :
  public IGraphReader
{
  public:
    /**
     * @brief Create an adapter so that a matrix can be treated as a graph. The
     * supplied pointer will be free'd when this object is destructed.
     *
     * @param mfile A pointer the MatrixFile to adapt.
     */
    MatrixGraphReader(
        std::unique_ptr<IMatrixReader>& file);


    /**
     * @brief Frees the underlying matrix file.
     */
    ~MatrixGraphReader();


    /**
     * @brief Get information about the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (directed).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are specified.
     */
    void getInfo(
        dim_t & nvtxs,
        ind_t & nedges,
        int & nvwgt,
        bool & ewgts) override;


    /**
     * @brief Read the CSR structure of the graph.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights (length nvtxs*nvwgt). This may be NULL in
     * order to ignore vertex weights. If it is specified and the file does not
     * contain vertex weights, it will be filled with ones.
     * @param adjwgt The edge weights (length nedges). This may be NULL in
     * order to ignore edge weights. If it is specified and the file does not
     * contain edge weights, it will be filled with ones.
     * @param progress The variable to update as the graph is loaded (may be
     * null).
     */
    void read(
        ind_t * xadj,
        dim_t * adjncy,
        val_t * vwgt,
        val_t * adjwgt,
        double * progress) override;




  private:
    std::unique_ptr<IMatrixReader> m_reader;


    // disable copying
    MatrixGraphReader(
        MatrixGraphReader const & rhs);
    MatrixGraphReader& operator=(
        MatrixGraphReader const & rhs);



};


}




#endif
