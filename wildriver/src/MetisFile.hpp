/**
 * @file MetisFile.hpp
 * @brief Class for reading/writing metis files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_METISFILE_HPP
#define WILDRIVER_METISFILE_HPP




#include "IGraphReader.hpp"
#include "IGraphWriter.hpp"
#include "TextFile.hpp"
#include "MatrixEntry.hpp"




namespace WildRiver {


class MetisFile : 
  public IGraphReader,
  public IGraphWriter
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
     * @brief Close file and free any memory.
     */
    virtual ~MetisFile();


    /**
     * @brief Create a new MetisFile for reading and writing.
     *
     * @param fname The filename/path.
     */
    MetisFile(
        std::string const & fname);


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
    virtual void read(
        ind_t * xadj,
        dim_t * adjncy,
        val_t * vwgt,
        val_t * adjwgt,
        double * progress) override;

  
    /**
     * @brief Get information about the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (directed).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are specified.
     */
    virtual void getInfo(
        dim_t & nvtxs,
        ind_t & nedges,
        int & nvwgt,
        bool & ewgts) override;


    /**
     * @brief Write a graph file from the given CSR structure.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights. 
     * @param adjwgt The edge weights.
     */
    virtual void write(
        ind_t const * xadj,
        dim_t const * adjncy,
        val_t const * vwgt,
        val_t const * adjwgt) override;


    /**
     * @brief Set the information for the graph. This must be called before
     * writing the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (an undirected edge counts as two).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are present.
     */
    virtual void setInfo(
        dim_t nvtxs,
        ind_t nedges,
        int nvwgt,
        bool ewgts) override;




  private:
    /**
     * @brief Whether or not the graph information has been set.
     */
    bool m_infoSet;


    /**
     * @brief The number of vertices in the graph. 
     */
    dim_t m_numVertices;


    /**
     * @brief The number of non-zeros in the matrix.
     */
    ind_t m_numEdges;


    /**
    * @brief The current row being processed.
    */
    dim_t m_currentVertex;


    /**
    * @brief The number of vertex weights in the graph file.
    */
    dim_t m_numVertexWeights;


    /**
    * @brief Whether or not the graph file stores edge weights.
    */
    bool m_hasEdgeWeights;


    /**
     * @brief Line buffer.
     */
    std::string m_line;


    /**
     * @brief The underlying text file.
     */
    TextFile m_file;


    /**
     * @brief Get the flags representing the weights associated with this
     * graph.
     *
     * @return 
     */
    int getWeightFlags();


    /**
     * @brief Determine if teh given line is a comment.
     *
     * @param line The line.
     *
     * @return True if the line is a comment.
     */
    virtual bool isComment(
        std::string const & line) const noexcept;


    /**
    * @brief Get the next non-comment line from the file.
    *
    * @param line The line to fill.
    *
    * @return True if the line was filled.
    */
    bool nextNoncommentLine(
        std::string & line);


    /**
     * @brief Read the header of this matrix file. Populates internal fields
     * with the header information.
     */
    virtual void readHeader();


    /**
     * @brief Write the header of this matrix file. The header consists of
     * internal fields set by "setInfo()".
     */
    virtual void writeHeader(); 


    /**
     * @brief Reset the current position in the graph to the first vertex.
     */
    virtual void firstVertex();


    /**
     * @brief Get the information of the next vertex.
     *
     * @param vertexWeights The vertex weight(s) (must be null or at least of
     * length equal to the number of constraints).
     * @param numEdges The number of edges incident to this vertex (output).
     * @param edgeDests The destination of each edge leaving this vertex (must
     * be of length equal to the number of edges of the vertex).
     * @param edgeWeights The weight of each edge leaving this vertex (must
     * null or be of length equal to the number of edges of the vertex).
     *
     * @return True if another vertex was found in the file.
     */
    virtual bool getNextVertex(
        val_t * vertexWeights,
        dim_t * numEdges,
        dim_t * edgeDests,
        val_t * edgeWeights);


    /**
     * @brief Set the adjacency list and vertex weight of the next vertex.
     *
     * @param vwgts The vertex weights for this vertex.
     * @param list The adjacecny list.
     */
    virtual void setNextVertex(
        std::vector<val_t> const & vwgts,
        std::vector<matrix_entry_struct> const & list);




};




}




#endif
