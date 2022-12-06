/**
* @file MatrixGraphWriter.hpp
* @brief Adapter class for writing matrix data as a graph.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2017
* @version 1
* @date 2017-07-24
*/




#ifndef WILDRIVER_MATRIXGRAPHWRITER_HPP
#define WILDRIVER_MATRIXGRAPHWRITER_HPP




#include <vector>
#include <string>
#include <memory>

#include "IGraphWriter.hpp"
#include "IMatrixWriter.hpp"




namespace WildRiver
{


class MatrixGraphWriter :
    public IGraphWriter
{
  public:
    /**
     * @brief Create an adapter so that a matrix can be treated as a graph. The
     * supplied pointer will be free'd when this object is destructed.
     *
     * @param mfile A pointer the MatrixFile to adapt.
     */
    MatrixGraphWriter(
        std::unique_ptr<IMatrixWriter>& file);


    /**
     * @brief Frees the underlying matrix file.
     */
    ~MatrixGraphWriter();


    /**
     * @brief Set the information for the graph. This must be called before
     * writing the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (an undirected edge counts as two).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are present.
     */
    void setInfo(
        dim_t nvtxs,
        ind_t nedges,
        int nvwgt,
        bool ewgts) override;


    /**
     * @brief Write a graph file from the given CSR structure.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights. 
     * @param adjwgt The edge weights.
     */
    void write(
        ind_t const * xadj,
        dim_t const * adjncy,
        val_t const * vwgt,
        val_t const * adjwgt) override;




  private:
    std::unique_ptr<IMatrixWriter> m_writer;


    // disable copying
    MatrixGraphWriter(
        MatrixGraphWriter const & rhs);
    MatrixGraphWriter& operator=(
        MatrixGraphWriter const & rhs);



};


}




#endif
