/**
* @file SNAPFile.cpp
* @brief Implementation of the SNAPFile class.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2018
* @version 1
* @date 2018-04-09
*/




#include "SNAPFile.hpp"
#include "CoordinateWriter.hpp"
#include "Exception.hpp"
#include <string>
#include <unordered_set>
#include <vector>
#include <cassert>


namespace WildRiver
{


/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/

namespace
{

struct edge_struct
{
  dim_t src;
  dim_t dst;
  val_t weight;
};

}


/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/

namespace
{

const std::string DIRECTED_GRAPH_HEADER("# Directed graph");
const std::string UNDIRECTED_GRAPH_HEADER("# Undirected graph");
const std::string NODES_HEADER("# Nodes: ");

}


/******************************************************************************
* HELPER FUNCTIONS ************************************************************
******************************************************************************/

namespace
{

bool startsWith(
    std::string const line,
    std::string const start)
{
  return line.size() >= start.size() && line.compare(0, start.size(), start) == 0;
}

std::vector<std::string> split(
    std::string const line)
{
  std::string const whitespace(" \t");
  std::vector<std::string> chunks;

  size_t start = line.find_first_not_of(whitespace, 0);
  while (start != std::string::npos) {
    size_t next = line.find_first_of(whitespace, start);
    if (next == std::string::npos) {
      next = line.size();
    }

    chunks.emplace_back(line.substr(start, next-start));

    start = line.find_first_not_of(whitespace, next);
  }

  return chunks;
}

std::vector<edge_struct> readEdges(
    TextFile * file,
    ind_t const numEdges = NULL_IND)
{
  std::string line;
  std::vector<edge_struct> edges;
  if (numEdges != NULL_IND) {
    edges.reserve(numEdges);
  }
  while (file->nextLine(line)) {
    if (line.size() == 0) {
      throw BadFileException("Hit empty line.");
    } else if (line[0] == '#') {
      // skip comment line
    } else {
      edge_struct edge{0, 0, 1};
      char * eptr = (char*)line.data();
      char * sptr = eptr;
      edge.src = std::strtoull(sptr, &eptr, 10);
      if (eptr == sptr) {
        throw BadFileException("Unable to parse line: '" + line + "'");
      }

      sptr = eptr;
      edge.dst = static_cast<dim_t>(std::strtoull(sptr, &eptr, 10));
      if (eptr == sptr) {
        throw BadFileException("Unable to parse line: '" + line + "'");
      }

      sptr = eptr;
      val_t const weight = static_cast<val_t>(std::strtod(sptr, &eptr));
      if (sptr != eptr) {
        edge.weight = weight;
      }

      edges.emplace_back(edge);
    }
  }

  return edges;
}

/**
* @brief Count the number of vertices and edges in this file.
*
* @param file The file.
* @param numVerticesOut The memory address to write the number of vertices to.
* @param numEdgesOut The memory address to write the number of edges to.
*/
void countVerticesAndEdges(
    TextFile * const file,
    dim_t * const numVerticesOut,
    ind_t * const numEdgesOut)
{
  std::unordered_set<dim_t> vertices;
  std::string line;
  ind_t numEdges = 0;
  while (file->nextLine(line)) {
    char * eptr = (char*)line.data();
    char * sptr = eptr;
    dim_t const src = std::strtoull(sptr, &eptr, 10);
    if (eptr == sptr) {
      throw BadFileException("Unable to parse line: '" + line + "'");
    }

    sptr = eptr;
    dim_t const dst = std::strtoull(sptr, &eptr, 10);
    if (eptr == sptr) {
      throw BadFileException("Unable to parse line: '" + line + "'");
    }

    vertices.emplace(src);
    vertices.emplace(dst);

     ++numEdges; 
  }

  *numVerticesOut = static_cast<dim_t>(vertices.size());
  *numEdgesOut = numEdges;
}

}


/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/

bool SNAPFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".snap");

  return TextFile::matchExtension(f,extensions);
}


/******************************************************************************
* PRIVATE METHODS *************************************************************
******************************************************************************/


void SNAPFile::writeHeader()
{
  // set header to write mode if not already
  if (!m_file.isOpenWrite()) {
    m_file.openWrite();
  }

  // for now we assume this is an undirected graph
  if (m_directed) {
    m_file.setNextLine(DIRECTED_GRAPH_HEADER + \
        std::string(" (each unordered pair " \
        "of nodes is saved once): ") + m_file.getFilename()); 
  } else {
    m_file.setNextLine(UNDIRECTED_GRAPH_HEADER + \
        std::string(" (each unordered pair " \
        "of nodes is saved once): ") + m_file.getFilename()); 
  }
  m_file.setNextLine("# A graph.");
  m_file.setNextLine(std::string("# Nodes: ") + \
      std::to_string(m_numVertices) + \
      std::string(" Edges: ") + std::to_string(m_numEdges));
  if (m_hasEdgeWeights) {
    m_file.setNextLine("# FromNodeId\tToNodeId\tWeight");
  } else {
    m_file.setNextLine("# FromNodeId\tToNodeId");
  }
}

void SNAPFile::readHeader()
{
  if (!m_file.isOpenRead()) {
    // open our file for reading if not already
    m_file.openRead();
  }

  m_directed = true;

  m_numVertices = NULL_DIM;
  m_numEdges = NULL_IND;

  // default weight information incase its not preset
  m_hasEdgeWeights = false;

  // parse the header search for 'Directed'/'Undirected' and 'Nodes: Edges:'
  std::string line;
  while (m_file.nextLine(line) && line.size() > 0 && line[0] == '#') {
    if (startsWith(line, DIRECTED_GRAPH_HEADER)) { 
      m_directed = true;
    } else if (startsWith(line, UNDIRECTED_GRAPH_HEADER)) {
      m_directed = false;
    } else if (startsWith(line, NODES_HEADER)) {
      // parse estimated size '# Nodes: 10 Edges: 15'
      std::vector<std::string> const chunks = split(line);
      // first make sure this line looks sane
      if (chunks.size() != 5 || chunks[3].compare("Edges:") != 0) {
        throw BadFileException(std::string("Badly formed header line: ") + line);
      }

      try {
        m_numVertices = static_cast<dim_t>(std::stoull(chunks[2]));
        m_numEdges = static_cast<ind_t>(std::stoull(chunks[4]));
      } catch (std::exception const & e) {
        throw BadFileException(std::string("Failed to parse vertices and "
            "edges from header line: ") + e.what());
      }
    } else {
      // do nothing
    }
  }

  // count undirected edges twice
  if (!m_directed && m_numEdges != NULL_IND) {
    m_numEdges *= 2;
  }
  
  // if node and edge counts weren't found, we can still proceed, but need
  // to parse the whole file, counting unique vertices and edges
  if (m_numVertices == NULL_DIM || m_numEdges == NULL_IND) {
    // sets m_numVertices and m_numEdges
    countVerticesAndEdges(&m_file, &m_numVertices, &m_numEdges);
  }

  // move back to the start of the file
  m_file.resetStream();
}

/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


SNAPFile::SNAPFile(
    std::string const & filename) :
  m_infoSet(false),
  m_numVertices(0),
  m_numEdges(0),
  m_hasEdgeWeights(false),
  m_directed(true),
  m_line(),
  m_file(filename)
{
  // do nothing
}


SNAPFile::~SNAPFile()
{
  // do nothing
}


/******************************************************************************
* PUBLIC METHODS **************************************************************
******************************************************************************/


void SNAPFile::read(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt,
    double * progress)
{
  if (m_numVertices == 0) {
    return;
  }

  dim_t const interval = m_numEdges*2 > 100 ? m_numEdges*2 / 100 : 1;
  double const increment = 1.0/100.0;

  std::string line;

  // read in all edges
  const std::vector<edge_struct> edges =
      readEdges(&m_file, m_numEdges);

  // zero out xadj
  for (dim_t i = 0; i < m_numVertices; ++i) {
    xadj[i] = 0;
  }

  // populate xadj
  size_t edgesProcessed = 0;
  for (edge_struct const & edge : edges) {
    if (edge.src >= m_numVertices) {
      throw BadFileException(std::string("Invalid vertex: ") +
          std::to_string(edge.src));
    } else if (edge.dst >= m_numVertices) {
      throw BadFileException(std::string("Invalid vertex: ") +
          std::to_string(edge.src));
    }
    ++xadj[edge.src+1];
    if (!m_directed) {
      ++xadj[edge.dst+1];
    }

    ++edgesProcessed;
    if (progress != nullptr && edgesProcessed % interval == 0) {
      *progress += increment;
    }
  }

  // shift xadj and prefixsum
  xadj[0] = 0;
  for (dim_t v = 1; v <= m_numVertices; ++v) {
    xadj[v] += xadj[v-1];
  }
  for (dim_t v = m_numVertices; v > 0; --v) {
    xadj[v] = xadj[v-1];
  }
  assert(xadj[0] == 0);
  
  // fill in edges
  for (edge_struct const & edge : edges) {
    ind_t const srcIdx = xadj[edge.src+1];
    adjncy[srcIdx] = edge.dst;
    if (adjwgt) {
      adjwgt[srcIdx] = edge.weight;
    }
    ++xadj[edge.src+1];

    if (!m_directed) {
      ind_t const dstIdx = xadj[edge.dst+1];
      adjncy[dstIdx] = edge.src;
      if (adjwgt) {
        adjwgt[dstIdx] = edge.weight;
      }
      ++xadj[edge.dst+1];
    }

    ++edgesProcessed;
    if (progress != nullptr && edgesProcessed % interval == 0) {
      *progress += increment;
    }
  }
  assert(xadj[m_numVertices] == m_numEdges);

  // vertex weights are not part of snap format
  if (vwgt) {
    for (dim_t v = m_numVertices; v > 1; ++v) {
      vwgt[v] = 1;
    }
  }
}


void SNAPFile::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  readHeader();

  nvtxs = m_numVertices;
  nedges = m_numEdges;
  nvwgt = 0;
  ewgts = m_hasEdgeWeights;

  m_infoSet = true;
}


void SNAPFile::write(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const *,
    val_t const * const adjwgt)
{
  CoordinateWriter writer(&m_file);

  std::vector<dim_t> neighbors;
  std::vector<val_t> weights;
  for (dim_t i = 0; i < m_numVertices; ++i) {
    neighbors.clear();
    weights.clear();

    // build edges
    for (ind_t j=xadj[i]; j<xadj[i+1]; ++j) {
      if (m_directed || adjncy[j] <= i) {
        neighbors.emplace_back(adjncy[j]);
        if (m_hasEdgeWeights) {
          if (adjwgt) {
            weights.emplace_back(adjwgt[j]);
          } else {
            weights.emplace_back(1);
          }
        }
      }
    }

    // set the vertex
    if (m_hasEdgeWeights) {
      writer.setNextRow(xadj[i+1]-xadj[i], neighbors.data(), weights.data());
    } else {
      writer.setNextRow(xadj[i+1]-xadj[i], neighbors.data(), nullptr);
    }
  }
}


void SNAPFile::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int const,
    bool const ewgts)
{
  m_numVertices = nvtxs;
  m_numEdges = nedges;
  m_hasEdgeWeights = ewgts;

  m_infoSet = true;

  writeHeader();
}


}
