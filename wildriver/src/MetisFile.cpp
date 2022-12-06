/**
 * @file MetisFile.cpp
 * @brief Class functions for reading and writing metis graph files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include <sstream>
#include "MetisFile.hpp"




namespace WildRiver
{


/******************************************************************************
* INTERNAL TYPES **************************************************************
******************************************************************************/


enum {
  HAS_NOWEIGHTS = 0,
  HAS_EDGEWEIGHTS = 1,
  HAS_VERTEXWEIGHTS = 10
};




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


namespace
{


static size_t const BUFFER_SIZE = 4096;


}




/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/


bool MetisFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".graph");
  extensions.push_back(".metis");
  extensions.push_back(".chaco");

  return TextFile::matchExtension(f,extensions);
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


int MetisFile::getWeightFlags()
{
  int flags = HAS_NOWEIGHTS;

  if (m_hasEdgeWeights) {
    flags |= HAS_EDGEWEIGHTS;
  }

  if (m_numVertexWeights > 0) {
    flags |= HAS_VERTEXWEIGHTS;
  }

  return flags;
}


bool MetisFile::isComment(
    std::string const & line) const noexcept
{
  // awful solution since I can't declare this statically in c++ -- at
  // somepoint generate all 256 entries using template programming
  bool comment_chars[256] = {false};
  comment_chars[static_cast<size_t>('#')] = true;
  comment_chars[static_cast<size_t>('%')] = true;
  comment_chars[static_cast<size_t>('"')] = true;
  comment_chars[static_cast<size_t>('/')] = true;

  return line.size() > 0 && comment_chars[static_cast<uint8_t>(line[0])];
}


void MetisFile::readHeader()
{
  if (!m_file.isOpenRead()) {
    // open our file for reading if not already
    m_file.openRead();
  }

  // get the first line
  std::string line;
  nextNoncommentLine(line);

  // parse out my header
  size_t offset;
  m_numVertices = std::stoull(line,&offset,10);
  line = line.substr(offset);

  m_numEdges = std::stoull(line,&offset,10)*2;
  line = line.substr(offset);

  // default weight information uncase its not preset
  m_hasEdgeWeights = false;
  m_numVertexWeights = 0;

  // handle weights
  if (line.find_first_not_of(" \t") != std::string::npos) {
    int flags = std::stoul(line, &offset, 10);
    line = line.substr(offset);

    if (flags & HAS_EDGEWEIGHTS) {
      m_hasEdgeWeights = true;
    }

    if (flags & HAS_VERTEXWEIGHTS) {
      m_numVertexWeights = std::stoul(line, &offset, 10);
    }
  }
}


void MetisFile::writeHeader()
{
  // make sure we have a valid number of edges/nnz
  if (m_numEdges % 2 != 0) {
    throw BadParameterException("Metis files are required to be symmetric: " \
        "odd number of non-zeroes specified.");
  }

  // set header to write mode if not already
  if (!m_file.isOpenWrite()) {
    m_file.openWrite();
  }

  std::stringstream stream;

  // write the header -- edges in metis files are undirected (symmetric).
  stream << m_numVertices << " " << (m_numEdges/2);

  int weightflag = getWeightFlags();

  if (weightflag != HAS_NOWEIGHTS) {
    // write weight flags
    stream << " " << weightflag;
    if (weightflag & HAS_VERTEXWEIGHTS) {
      // write number of vertex weights
      stream << " " << m_numVertexWeights;
    }
  }

  m_file.setNextLine(stream.str());
}


void MetisFile::firstVertex()
{
  m_file.resetStream();

  // go past all comments -- read and discard header
  nextNoncommentLine(m_line);
}


bool MetisFile::getNextVertex(
        val_t * const vertexWeights,
        dim_t * const numEdges,
        dim_t * const edgeDests,
        val_t * const edgeWeights)
{
  dim_t const ncon = m_numVertexWeights;

  // get my m_line
  if (!nextNoncommentLine(m_line)) {
    return false;
  }

  char * eptr = (char*)m_line.data();

  // read in vertex weights
  for (dim_t k=0; k<ncon; ++k) {
    char * const sptr = eptr;
    val_t const val = std::strtod(sptr,&eptr);
    if (sptr == eptr) {
      throw BadFileException(std::string("Failed to read vertex weight on " \
            "line ") + std::to_string(m_file.getCurrentLine()));
    }

    if (vertexWeights != nullptr) {
      vertexWeights[k] = val;
    }
  }

  dim_t degree = 0;
  // read in edges
  while (true) {
    char * sptr = eptr;
    dim_t const dst = static_cast<dim_t>(std::strtoull(sptr,&eptr,10))-1;
    if (sptr == eptr) {
      break;
    }

    // make sure this is a valid edge
    if (dst >= m_numVertices) {
      throw BadFileException(std::string("Edge with destination of ") + \
          std::to_string(dst) + std::string("/") + \
          std::to_string(m_numVertices));
    }

    if (edgeDests != nullptr) {
      edgeDests[degree] = dst;
    }

    val_t wgt = static_cast<val_t>(1);
    if (m_hasEdgeWeights) {
      sptr = eptr;
      wgt = static_cast<val_t>(std::strtod(sptr,&eptr));
      if (sptr == eptr) {
        throw BadFileException(std::string("Could not read edge weight at "
              "line ") + std::to_string(m_file.getCurrentLine()));
      }
    } 
    
    if (edgeWeights != nullptr) {
      edgeWeights[degree] = wgt;
    }

    ++degree;
  }

  *numEdges = degree;

  // indicate that we successfully found a vertex
  return true;
}


void MetisFile::setNextVertex(
    std::vector<val_t> const & vwgts,
    std::vector<matrix_entry_struct> const & list)
{
  dim_t const ncon = m_numVertexWeights;
  dim_t const nadj = list.size();

  std::stringstream stream;

  // set vertex weights
  for (dim_t k=0; k<ncon; ++k) {
    stream << vwgts[k];
    if (k < ncon-1 || nadj > 0) {
      stream << " ";
    }
  }

  // set adjacency list and edge weights
  bool const ewgts = m_hasEdgeWeights;

  for (dim_t j = 0; j < nadj; ++j) {
    matrix_entry_struct e = list[j];
    stream << (e.ind+1);
    if (ewgts) {
      stream << " " << e.val;
    }
    if (j < nadj-1) {
      // add space
      stream << " ";
    }
  }

  m_file.setNextLine(stream.str());

  ++m_currentVertex;
}


bool MetisFile::nextNoncommentLine(
    std::string & line)
{
  do {
    if (!m_file.nextLine(line)) {
      return false;
    }
  } while (isComment(line));

  return true;
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MetisFile::MetisFile(
    std::string const & fname) :
  m_infoSet(false),
  m_numVertices(NULL_DIM),
  m_numEdges(NULL_IND),
  m_currentVertex(0),
  m_numVertexWeights(0),
  m_hasEdgeWeights(false),
  m_line(BUFFER_SIZE,'\0'),
  m_file(fname)
{
  // do nothing
}


MetisFile::~MetisFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MetisFile::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  if (!m_infoSet) {
    readHeader();
    m_infoSet = true;
  }

  nvtxs = m_numVertices;
  nedges = m_numEdges;
  nvwgt = m_numVertexWeights;
  ewgts = m_hasEdgeWeights;
}


void MetisFile::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int nvwgt,
    bool ewgts)
{
  m_numVertices = nvtxs;
  m_numEdges = nedges;
  m_numVertexWeights = nvwgt;
  m_hasEdgeWeights = ewgts;

  m_infoSet = true;

  writeHeader();
}


void MetisFile::read(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt,
    double * progress)
{
  dim_t const interval = m_numVertices > 100 ? m_numVertices / 100 : 1;
  double const increment = 1/100.0;

  xadj[0] = 0;
  for (dim_t i = 0; i < m_numVertices; ++i) {
    dim_t degree;

    val_t * const vwgtStart = (m_numVertexWeights > 0 && vwgt) ? \
        vwgt+(i*m_numVertexWeights) : nullptr;
    dim_t * const adjncyStart = adjncy+xadj[i];
    val_t * const adjwgtStart = adjwgt ? adjwgt+xadj[i] : nullptr;

    // retrieve row
    if (!getNextVertex(vwgtStart,&degree,adjncyStart,adjwgtStart)) {
      throw BadFileException(std::string("Premature end of file: ") + \
          std::to_string(i) + std::string("/") + \
          std::to_string(m_numVertices) + std::string(" vertices found."));
    }

    // handle vertex weights
    if (m_numVertexWeights == 0 && vwgt) {
      // set unit vertex weights
      vwgt[i] = 1;
    }

    xadj[i+1] = xadj[i]+degree;

    if (progress != nullptr && i % interval == 0) {
      *progress += increment;
    }
  }
}


void MetisFile::write(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  for (dim_t i = 0; i < m_numVertices; ++i) {
    std::vector<val_t> vwgts;
    std::vector<matrix_entry_struct> list;

    // build vertex weight vector
    for (dim_t k = 0; k < m_numVertexWeights; ++k) {
      vwgts.push_back(vwgt[(i*m_numVertexWeights)+k]);
    }

    // build edges
    for (ind_t j=xadj[i]; j<xadj[i+1]; ++j) {
      matrix_entry_struct e;
      e.ind = adjncy[j];
      if (m_hasEdgeWeights) {
        if (adjwgt) {
          e.val = adjwgt[j];
        } else {
          e.val = 1;
        }
      }
      list.push_back(e);
    }

    // set the vertex
    setNextVertex(vwgts,list);
  }
}





}
