/**
 * @file CSRFile_test.cpp
 * @brief Test for reading and writing CSR formatted csrs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2017
 * @version 1
 *
 */




#include <iostream>
#include <fstream>
#include <memory>

#include "CSRFile.hpp"
#include "DomTest.hpp"




using namespace WildRiver;




namespace DomTest
{


static void writeTest(
    std::string const & testFile)
{
  CSRFile csr(testFile);

  csr.setInfo(6,6,14);

  wildriver_ind_t rowptr[] = {0,2,4,7,10,12,14};
  wildriver_dim_t rowind[] = {1,2,0,2,0,1,3,2,4,5,3,5,3,4};
  wildriver_val_t rowval[] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5};

  csr.write(rowptr,rowind,rowval);
}


static void readTestBaseZero(
    std::string const & testFile)
{
  std::ofstream fout(testFile, std::ofstream::trunc);
  fout << "1 1.0 2 2.0" << std::endl;
  fout << "0 3.0 2 4.0" << std::endl;
  fout << "0 5.0 1 6.0 3 7.0" << std::endl;
  fout << "2 8.0 4 9.0 5 1.0" << std::endl;
  fout << "3 2.0 5 3.0" << std::endl;
  fout << "3 4.0 4 5.0" << std::endl;

  fout.close();

  CSRFile csr(testFile);

  wildriver_dim_t nrows, ncols;
  wildriver_ind_t nnz;

  csr.getInfo(nrows,ncols,nnz);

  testEquals(nrows,6);
  testEquals(ncols,6);
  testEquals(nnz,14);

  std::unique_ptr<wildriver_ind_t[]> rowptr(new wildriver_ind_t[nrows+1]);
  std::unique_ptr<wildriver_dim_t[]> rowind(new wildriver_dim_t[nnz]);
  std::unique_ptr<wildriver_val_t[]> rowval(new wildriver_val_t[nnz]);

  csr.read(rowptr.get(),rowind.get(),rowval.get(),nullptr);

  // test rowptr
  testEquals(rowptr[0],0);
  testEquals(rowptr[1],2);
  testEquals(rowptr[2],4);
  testEquals(rowptr[3],7);
  testEquals(rowptr[4],10);
  testEquals(rowptr[5],12);
  testEquals(rowptr[6],14);

  // test rowind
  testEquals(rowind[0],1);
  testEquals(rowind[1],2);
  testEquals(rowval[0],1);
  testEquals(rowval[1],2);

  testEquals(rowind[2],0);
  testEquals(rowind[3],2);
  testEquals(rowval[2],3);
  testEquals(rowval[3],4);

  testEquals(rowind[4],0);
  testEquals(rowind[5],1);
  testEquals(rowind[6],3);
  testEquals(rowval[4],5);
  testEquals(rowval[5],6);
  testEquals(rowval[6],7);

  testEquals(rowind[7],2);
  testEquals(rowind[8],4);
  testEquals(rowind[9],5);
  testEquals(rowval[7],8);
  testEquals(rowval[8],9);
  testEquals(rowval[9],1);

  testEquals(rowind[10],3);
  testEquals(rowind[11],5);
  testEquals(rowval[10],2);
  testEquals(rowval[11],3);

  testEquals(rowind[12],3);
  testEquals(rowind[13],4);
  testEquals(rowval[12],4);
  testEquals(rowval[13],5);
}


static void readTestBaseOne(
    std::string const & testFile)
{
  std::ofstream fout(testFile, std::ofstream::trunc);
  fout << "2 1.0 3 2.0" << std::endl;
  fout << "1 3.0 3 4.0" << std::endl;
  fout << "1 5.0 2 6.0 4 7.0" << std::endl;
  fout << "3 8.0 5 9.0 6 1.0" << std::endl;
  fout << "4 2.0 6 3.0" << std::endl;
  fout << "4 4.0 5 5.0" << std::endl;

  fout.close();

  CSRFile csr(testFile);

  wildriver_dim_t nrows, ncols;
  wildriver_ind_t nnz;

  csr.getInfo(nrows,ncols,nnz);

  testEquals(nrows,6);
  testEquals(ncols,6);
  testEquals(nnz,14);

  std::unique_ptr<wildriver_ind_t[]> rowptr(new wildriver_ind_t[nrows+1]);
  std::unique_ptr<wildriver_dim_t[]> rowind(new wildriver_dim_t[nnz]);
  std::unique_ptr<wildriver_val_t[]> rowval(new wildriver_val_t[nnz]);

  csr.read(rowptr.get(),rowind.get(),rowval.get(),nullptr);

  // test rowptr
  testEquals(rowptr[0],0);
  testEquals(rowptr[1],2);
  testEquals(rowptr[2],4);
  testEquals(rowptr[3],7);
  testEquals(rowptr[4],10);
  testEquals(rowptr[5],12);
  testEquals(rowptr[6],14);

  // test rowind
  testEquals(rowind[0],1);
  testEquals(rowind[1],2);
  testEquals(rowval[0],1);
  testEquals(rowval[1],2);

  testEquals(rowind[2],0);
  testEquals(rowind[3],2);
  testEquals(rowval[2],3);
  testEquals(rowval[3],4);

  testEquals(rowind[4],0);
  testEquals(rowind[5],1);
  testEquals(rowind[6],3);
  testEquals(rowval[4],5);
  testEquals(rowval[5],6);
  testEquals(rowval[6],7);

  testEquals(rowind[7],2);
  testEquals(rowind[8],4);
  testEquals(rowind[9],5);
  testEquals(rowval[7],8);
  testEquals(rowval[8],9);
  testEquals(rowval[9],1);

  testEquals(rowind[10],3);
  testEquals(rowind[11],5);
  testEquals(rowval[10],2);
  testEquals(rowval[11],3);

  testEquals(rowind[12],3);
  testEquals(rowind[13],4);
  testEquals(rowval[12],4);
  testEquals(rowval[13],5);
}



void Test::run()
{
  std::string testFile("./CSRFile_test.csr");

  writeTest(testFile);

  remove(testFile.c_str());

  readTestBaseZero(testFile);

  remove(testFile.c_str());

  readTestBaseOne(testFile);

  remove(testFile.c_str());

}




}
