/**
* @file MatrixMarketFile.cpp
* @brief Implementation of the MatrixMarketFile class.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2017
* @version 1
* @date 2017-07-28
*/




#include "MatrixMarketFile.hpp"

#include "Exception.hpp"
#include "Util.hpp"

#include <map>
#include <cassert>





namespace WildRiver
{


/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


namespace
{


size_t const BUFFER_SIZE = 1024;

std::vector<std::string> const EXTENSIONS{".mm", ".mtx"};


enum matrix_market_attributes {
  MATRIX_MARKET_NULL,
  MATRIX_MARKET_MATRIX,
  MATRIX_MARKET_VECTOR,
  MATRIX_MARKET_GENERAL,
  MATRIX_MARKET_SYMMETRIC,
  MATRIX_MARKET_PATTERN,
  MATRIX_MARKET_REAL,
  MATRIX_MARKET_INTEGER,
  MATRIX_MARKET_COMPLEX,
  MATRIX_MARKET_ARRAY,
  MATRIX_MARKET_COORDINATE
};


enum matrix_market_symmetric_orientation {
  ORIENTATION_UNKNOWN,
  ORIENTATION_LOWER,
  ORIENTATION_UPPER
};


std::string const BASE_HEADER("%%MatrixMarket");


std::map<std::string, int> const ENTITY_MAPPING{
  {"matrix", MATRIX_MARKET_MATRIX},
  {"vector", MATRIX_MARKET_VECTOR}
};


std::map<std::string, int> const STORAGE_MAPPING{
  {"general", MATRIX_MARKET_GENERAL},
  {"symmetric", MATRIX_MARKET_SYMMETRIC}
};


std::map<std::string, int> const FORMAT_MAPPING{
  {"array", MATRIX_MARKET_ARRAY},
  {"coordinate", MATRIX_MARKET_COORDINATE}
};


std::map<std::string, int> const TYPE_MAPPING{
  {"real", MATRIX_MARKET_REAL},
  {"complex", MATRIX_MARKET_REAL},
  {"integer", MATRIX_MARKET_INTEGER},
  {"pattern", MATRIX_MARKET_PATTERN}
};




}




/******************************************************************************
* HELPER FUNCTIONS ************************************************************
******************************************************************************/


namespace
{


bool isComment(
    std::string const & line) noexcept
{
  return line.size() > 0 && line[0] == '%';
}


/**
 * @brief Parse a line containing a triplet (row, column, value).
 *
 * @tparam D The type of the 'row' and 'column'.
 * @tparam V The type of the 'value'.
 * @param line The line.
 * @param row The row (output).
 * @param col The column (output).
 * @param val The value (output).
 */
template<typename D, typename V>
void parseTriplet(
    std::string const * const line,
    D * const row,
    D * const col,
    V * const val)
{
  char * sptr;
  std::vector<char> stringBuffer(line->c_str(), line->c_str() + \
      line->length()+1);
  char * eptr = stringBuffer.data();
  sptr = eptr;
  *row = static_cast<D>(std::strtoull(sptr, &eptr, 10));
  if (eptr == sptr) {
    throw BadFileException(std::string("Unable to parse triplet row: ") + \
        *line);
  }

  sptr = eptr;
  *col = static_cast<D>(std::strtoull(sptr, &eptr, 10));
  if (eptr == sptr) {
    throw BadFileException(std::string("Unable to parse triplet " \
        "column: ") + *line);
  }

  V tmpVal;
  sptr = eptr;
  if (std::is_floating_point<V>::value) {
    tmpVal = static_cast<V>(std::strtod(sptr, &eptr));
  } else {
    tmpVal = static_cast<V>(std::strtoll(sptr, &eptr, 10));
  }

  if (val != nullptr) {
    if (eptr == sptr) {
      throw BadFileException(std::string("Unable to parse triplet " \
          "value: ") + *line);
    }

    *val = tmpVal;
  }
}




}




/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/


bool MatrixMarketFile::hasExtension(
    std::string const & f)
{
  return TextFile::matchExtension(f, EXTENSIONS);
}




/******************************************************************************
* PRIVATE SECTIONS ************************************************************
******************************************************************************/


bool MatrixMarketFile::nextNoncommentLine(
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
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


MatrixMarketFile::MatrixMarketFile(
    std::string const & filename) :
  m_infoSet(false),
  m_nrows(NULL_DIM),
  m_ncols(NULL_DIM),
  m_nnz(NULL_IND),
  m_line(BUFFER_SIZE,'\0'),
  m_file(filename),
  m_entity(MATRIX_MARKET_NULL),
  m_format(MATRIX_MARKET_NULL),
  m_type(MATRIX_MARKET_NULL),
  m_symmetric(false)
{
  // do nothing
}



MatrixMarketFile::~MatrixMarketFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixMarketFile::getInfo(
    dim_t & nrows,
    dim_t & ncols,
    ind_t & nnz)
{
  if (!m_infoSet) {
    readHeader();
  }

  nrows = m_nrows;
  ncols = m_ncols;
  nnz = m_nnz;
}


void MatrixMarketFile::read(
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval,
    double * const progress)
{
  if (!m_infoSet) {
    // read header
    throw UnsetInfoException("Cannot call read() before calling getInfo()");
  }

  if (m_format == MATRIX_MARKET_COORDINATE) {
    // sparse
    if (m_symmetric) {
      readSymmetricCoordinates(rowptr, rowind, rowval, progress);
    } else {
      readCoordinates(rowptr, rowind, rowval, progress);
    }
  } else if (m_format == MATRIX_MARKET_ARRAY) {
    // dense
    readArray();
  }
}


void MatrixMarketFile::setInfo(
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz)
{
  if (m_infoSet) {
    throw BadFileStateException("Info already set.");
  }

  m_nrows = nrows;
  m_ncols = ncols;
  m_nnz = nnz;

  writeHeader();

  m_infoSet = true;
}


void MatrixMarketFile::write(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  writeCoordinates(rowptr, rowind, rowval);
}


void MatrixMarketFile::readHeader()
{
  m_file.openRead();

  // read the special header
  if (!m_file.nextLine(m_line)) {
    throw BadFileException(std::string("Not a valid MatrixMarket file '") + \
        m_file.getFilename() + std::string("': is empty."));
  }
  if (m_line.find(BASE_HEADER + std::string(" ")) != 0) {
    throw BadFileException(std::string("Not a valid MatrixMarket file '") + \
        m_file.getFilename() + std::string("': first line is: '") + m_line + \
      std::string("'."));
  }

  // discard base header
  m_line = m_line.substr(BASE_HEADER.length()+1);

  // break rest of line up into whitespace delimited tokens
  const std::vector<std::string> chunks = Util::split(m_line);

  if (chunks.size() == 0) {
    throw BadFileException(std::string("Not a valid MatrixMarket file '") + \
        m_file.getFilename() + std::string("': first line is: '") + m_line + \
      std::string("'."));
  }

  // figure out the entity
  m_entity = ENTITY_MAPPING.at(chunks[0]);

  if (m_entity == MATRIX_MARKET_MATRIX) {
    // parse out matrix data
    if (chunks.size() < 2) {
      throw BadFileException(std::string("Missing 'coordinate' or 'array' " \
            "specifier for matrix '") + m_file.getFilename() + \
            std::string("'."));
    }
    m_format = FORMAT_MAPPING.at(chunks[1]);

    if (m_format == MATRIX_MARKET_COORDINATE) {
      // read past all comments until we get to the size line
      do {
        if (!m_file.nextLine(m_line)) {
          throw BadFileException(std::string("Failed to find header line " \
              "in '") + m_file.getFilename() + std::string("'."));
        }
      } while (isComment(m_line));

      parseTriplet(&m_line, &m_nrows, &m_ncols, &m_nnz);
    } else {
      throw BadFileException("Array matrices are not yet supported.");
    }

    if (chunks.size() < 3) {
      throw BadFileException(std::string("Missing 'real', " \
            "'integer', or 'pattern' specifier for matrix '") + \
          m_file.getFilename() + std::string("'."));
    }
    m_type = TYPE_MAPPING.at(chunks[2]);
    if (m_type == MATRIX_MARKET_COMPLEX) {
      throw BadFileException("Complex numbers are not supported");
    }

    if (chunks.size() < 4) {
      throw BadFileException(std::string("Missing 'general' or " \
          "'symmetric' specifier for matrix '") + m_file.getFilename() + \
          std::string("'."));
    }
    m_symmetric = STORAGE_MAPPING.at(chunks[3]) == MATRIX_MARKET_SYMMETRIC;
    if (m_symmetric) {
      // we need to pesimistically predict m_nnz -- there could be up to
      // twice as many as reported if there are no diagonal entries.
      m_nnz *= 2;
    }

  } else if (m_entity == MATRIX_MARKET_VECTOR) {
    // parse out vector data
    throw BadFileException("Vectors are currently unsupported.");
  } else {
    throw BadFileException(std::string("Unknown matrix type code ") + \
        std::to_string(m_type));
  }

  m_infoSet = true;
}


void MatrixMarketFile::writeHeader()
{
  m_file.openWrite();

  // TODO: this locks in writing CSRs
  m_file.setNextLine(BASE_HEADER + \
      std::string(" matrix coordinate real general"));
  m_file.setNextLine("%====================================================");
  m_file.setNextLine("%= Generated by wildriver. =");
  m_file.setNextLine("%====================================================");
  m_file.setNextLine(std::to_string(m_nrows) + std::string(" ") + \
    std::to_string(m_ncols) + std::string(" ") + std::to_string(m_nnz));
}


void MatrixMarketFile::readCoordinates(
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval,
    double * const progress)
{
  // make these large enough to hold whatever value is in the file
  int64_t row, col;
  val_t value;

  // TODO: I'd like to take advantage of cases where the triplets are in order:
  // if the triplets are in row major order, we can do it in a single pass
  // and be very happy. If they are out of order, we need to passes 1) to count
  // the size of each row and 2) to fill each row.

  // for now we'll just assume they are out of order.

  // the strategy is to allocate a new 'row' array, read in our coordinate
  // data, and sort it
  std::vector<dim_t> rows(m_nnz);

  // zero out rowptr
  for (size_t i = 0; i < m_nrows+1; ++i) {
    rowptr[i] = 0;
  }

  dim_t const interval = m_nnz > 100 ? m_nnz / 100 : 1;
  double const increment = 1.0/100.0;

  // count non-zeros per row
  for (ind_t nnz = 0; nnz < m_nnz; ++nnz) {
    if (!nextNoncommentLine(m_line)) {
      throw BadFileException(std::string("Only found ") + \
          std::to_string(nnz) + std::string("/") + std::to_string(m_nnz) + \
          std::string(" non-zeros."));
    }

    if (m_type == MATRIX_MARKET_PATTERN) {
      parseTriplet(&m_line, &row, &col, static_cast<val_t*>(nullptr));
      value = 1;
    } else if (m_type == MATRIX_MARKET_REAL || \
        m_type == MATRIX_MARKET_INTEGER) {
      parseTriplet(&m_line, &row, &col, &value);
    } else {
      throw BadFileException("Complex types are not supported.");
    }

    // handle 1-based rows
    if (row <= 0) {
      throw BadFileException(std::string("Invalid row ") + \
          std::to_string(row) + std::string(" must be 1-based indexing."));
    } else if (row > m_nrows) {
      throw BadFileException(std::string("Invalid row ") + \
          std::to_string(row) + std::string(" exceeds total rows ") + \
          std::to_string(m_nrows) + std::string("."));
    } 
    --row;

    // handle 1-based columns
    if (col <= 0) {
      throw BadFileException(std::string("Invalid column ") + \
          std::to_string(col) + std::string(" must be 1-based indexing."));
    } else if (row > m_nrows) {
      throw BadFileException(std::string("Invalid column ") + \
          std::to_string(col) + std::string(" exceeds total columns ") + \
          std::to_string(m_ncols) + std::string("."));
    }
    --col;

    rows[nnz] = static_cast<dim_t>(row);
    rowind[nnz] = static_cast<dim_t>(col);
    rowval[nnz] = value;

    ++rowptr[row+1];

    if (progress != nullptr && nnz % interval == 0) {
      *progress += increment;
    }
  }

  // prefix sum rows in the second row
  for (ind_t i = 1; i < m_nrows+1; ++i) {
    rowptr[i] += rowptr[i-1]; 
  }
  assert(rowptr[0] == 0);
  assert(rowptr[m_nrows] == m_nnz);

  // at this point we have a rowptr ready for insertion, and all triplets
  // stored in memory --  

  // determine the source of each nz in the 'rows' variable
  std::vector<ind_t> source(m_nnz);
  for (ind_t nnz = 0; nnz < m_nnz; ++nnz) {
    ind_t const dest = rowptr[rows[nnz]]++;
    source[dest] = nnz;
  }

  // we now need to sort the rowind and rowval arrays

  // re-use rows to swap out rowind -- copy rowind over and then use source to
  // copy it back
  rows.assign(rowind, rowind+m_nnz);
  for (ind_t nnz = 0; nnz < m_nnz; ++nnz) {
    ind_t const src = source[nnz];
    rowind[nnz] = rows[src];
  }
  // throw away rows
  std::vector<dim_t>().swap(rows);

  // copy over rowval
  std::vector<val_t> vals(rowval, rowval+m_nnz);
  for (ind_t nnz = 0; nnz < m_nnz; ++nnz) {
    ind_t const src = source[nnz];
    rowval[nnz] = vals[src];
  }

  // throw away vals
  std::vector<val_t>().swap(vals);

  // shift
  for (ind_t i = m_nrows; i > 0; --i) {
    rowptr[i] = rowptr[i-1]; 
  }
  rowptr[0] = 0;

  assert(rowptr[m_nrows] == m_nnz);
}


void MatrixMarketFile::readSymmetricCoordinates(
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval,
    double * const progress)
{
  // make these large enough to hold whatever value is in the file
  int64_t row, col;
  val_t value;
  int orientation = ORIENTATION_UNKNOWN;

  // TODO: Avoiding the excess memory is more tricky here than for the
  // 'general' case, as even if the nnz's are in order, the correspodning
  // upper diagonal entries won't be. Still would be nice to do it in place.

  // the strategy is to allocate a new 'row' array, read in our coordinate
  // data, and sort it
  std::vector<dim_t> rows(m_nnz);

  // zero out rowptr
  for (size_t i = 0; i < m_nrows+1; ++i) {
    rowptr[i] = 0;
  }

  // count non-zeros per row
  ind_t nnz = 0;
  // we doubled the nnz when reading the header, so halving it hear for the
  // number of lines is safe
  ind_t nlines = m_nnz / 2;

  dim_t const interval = nlines > 100 ? nlines / 100 : 1;
  double const increment = 1.0/100.0;

  for (ind_t line = 0; line < nlines; ++line) {
    if (!nextNoncommentLine(m_line)) {
      throw BadFileException(std::string("Only found ") + \
          std::to_string(nnz) + std::string("/") + std::to_string(m_nnz) + \
          std::string(" non-zeros."));
    }

    if (m_type == MATRIX_MARKET_PATTERN) {
      parseTriplet(&m_line, &row, &col, static_cast<val_t*>(nullptr));
      value = 1;
    } else if (m_type == MATRIX_MARKET_REAL || \
        m_type == MATRIX_MARKET_INTEGER) {
      parseTriplet(&m_line, &row, &col, &value);
    } else {
      throw BadFileException("Complex types are not supported.");
    }

    // handle 1-based rows
    if (row <= 0) {
      throw BadFileException(std::string("Invalid row ") + \
          std::to_string(row) + std::string(" must be 1-based indexing."));
    } else if (row > m_nrows) {
      throw BadFileException(std::string("Invalid row ") + \
          std::to_string(row) + std::string(" exceeds total rows ") + \
          std::to_string(m_nrows) + std::string("."));
    } 
    --row;

    // handle 1-based columns
    if (col <= 0) {
      throw BadFileException(std::string("Invalid column ") + \
          std::to_string(col) + std::string(" must be 1-based indexing."));
    } else if (col > m_ncols) {
      throw BadFileException(std::string("Invalid column ") + \
          std::to_string(col) + std::string(" exceeds total columns ") + \
          std::to_string(m_ncols) + std::string("."));
    }
    --col;

    // handle symmetry
    if (orientation == ORIENTATION_UNKNOWN) {
      // if we're off diagonal, set orienation
      if (col > row) {
        orientation = ORIENTATION_UPPER;
      } else if (col < row) {
        orientation = ORIENTATION_LOWER;
      }
    } else if (orientation == ORIENTATION_LOWER && col > row) {
      // invalid for symmetric
      throw BadFileException(std::string("Non-zero in upper triangle: (") + \
          std::to_string(row) + std::string(", ") + std::to_string(col) + \
          std::string(") when lower triangle non-zeros have been found."));
    } else if (orientation == ORIENTATION_UPPER && col < row) {
      // invalid for symmetric
      throw BadFileException(std::string("Non-zero in lower triangle: (") + \
          std::to_string(row) + std::string(", ") + std::to_string(col) + \
          std::string(") when upper triangle non-zeros have been found."));
    }

    rows[nnz] = static_cast<dim_t>(row);
    rowind[nnz] = static_cast<dim_t>(col);
    rowval[nnz] = value;

    ++rowptr[row+1];
    ++nnz;

    // add corresponding entry (if not diagonal)
    if (row != col) {
      rows[nnz] = static_cast<dim_t>(col);
      rowind[nnz] = static_cast<dim_t>(row);
      rowval[nnz] = value;

      ++rowptr[col+1];
      ++nnz;
    }

    if (progress != nullptr && line % interval == 0) {
      *progress += increment;
    }
  }
  // set proper nnz count
  m_nnz = nnz;

  // prefix sum rows in the second row
  for (ind_t i = 1; i < m_nrows+1; ++i) {
    rowptr[i] += rowptr[i-1]; 
  }
  assert(rowptr[0] == 0);
  assert(rowptr[m_nrows] == m_nnz);

  // at this point we have a rowptr ready for insertion, and all triplets
  // stored in memory --  

  // determine the source of each nz in the 'rows' variable
  std::vector<ind_t> source(m_nnz);
  for (ind_t nnz = 0; nnz < m_nnz; ++nnz) {
    ind_t const dest = rowptr[rows[nnz]]++;
    source[dest] = nnz;
  }

  // we now need to sort the rowind and rowval arrays

  // re-use rows to swap out rowind -- copy rowind over and then use source to
  // copy it back
  rows.assign(rowind, rowind+m_nnz);
  for (nnz = 0; nnz < m_nnz; ++nnz) {
    ind_t const src = source[nnz];
    rowind[nnz] = rows[src];
  }
  // throw away rows
  std::vector<dim_t>().swap(rows);

  // copy over rowval
  std::vector<val_t> vals(rowval, rowval+m_nnz);
  for (nnz = 0; nnz < m_nnz; ++nnz) {
    ind_t const src = source[nnz];
    rowval[nnz] = vals[src];
  }

  // throw away vals
  std::vector<val_t>().swap(vals);

  // shift
  for (ind_t i = m_nrows; i > 0; --i) {
    rowptr[i] = rowptr[i-1]; 
  }
  rowptr[0] = 0;

  assert(rowptr[m_nrows] == m_nnz);
}


void MatrixMarketFile::readArray()
{
  throw BadFileException("Reading arrays unimplemented.");
}


void MatrixMarketFile::writeCoordinates(
    ind_t const * rowptr,
    dim_t const * rowind,
    val_t const * rowval)
{
  for (dim_t row = 0; row < m_nrows; ++row) {
    for (ind_t nnz = rowptr[row]; nnz < rowptr[row+1]; ++nnz) {
      m_file.setNextLine(std::to_string(row+1) + std::string(" ") + \
          std::to_string(rowind[nnz]+1) + std::string(" ") + \
          std::to_string(rowval[nnz]));
    }
  }
}




}
