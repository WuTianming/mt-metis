/**
* @file SNAPFile.hpp
* @brief The SNAPFile class.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2018
* @version 1
* @date 2018-04-07
*/



#ifndef WILDRIVER_SNAPFILE_HPP
#define WILDRIVER_SNAPFILE_HPP

#include "IGraphReader.hpp"
#include "IGraphWriter.hpp"
#include "TextFile.hpp"


namespace WildRiver
{

/**
* @brief A class for reading and writing files in the SNAP graph text format.
* SNAP is an unweighted IJ graph format.
*
* ```
* # Directed graph (each unordered pair of nodes is saved once): example.txt 
* # Description of the graph.
* # Nodes: 6 Edges: 5
* # FromNodeId	ToNodeId
* 0	1
* 0	2
* 0	3
* 0	4
* 0	5
* ```
* The first line must either start with `# Undirected graph`, or
* `# Directed graph`. Undirected graphs get their edges counted in both
* directions, whereas directed graphs get the counted in only one.
* 
* @misc{snapnets,
*   author       = {Jure Leskovec and Andrej Krevl},
*   title        = {{SNAP Datasets}: {Stanford} Large Network Dataset Collection},
*   howpublished = {\url{http://snap.stanford.edu/data}},
*   month        = jun,
*   year         = 2014
* }
*/
class SNAPFile :
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
     * @brief Create a new MetisFile for reading and writing.
     *
     * @param fname The filename/path.
     */
    SNAPFile(
        std::string const & fname);


    /**
     * @brief Close file and free any memory.
     */
    virtual ~SNAPFile();


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
    * @brief Whether or not the graph file stores edge weights.
    */
    bool m_hasEdgeWeights;


    /**
    * @brief Whether or not the graph is directed.
    */
    bool m_directed;


    /**
     * @brief Line buffer.
     */
    std::string m_line;


    /**
     * @brief The underlying text file.
     */
    TextFile m_file;


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




};



}



#endif
