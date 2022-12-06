/**
 * @file GraphInHandle_test.cpp
 * @brief Test for reading matrices
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015
 * @version 1
 *
 */




#include <iostream>
#include <fstream>
#include <memory>

#include "GraphInHandle.hpp"
#include "DomTest.hpp"




using namespace WildRiver;




namespace DomTest
{


static void writeMetis(
    std::string const & testFile)
{
  std::fstream stream(testFile,std::fstream::out | std::fstream::trunc);

  stream << "6 7 1" << std::endl;

  stream << "2 1 3 2" << std::endl;
  stream << "1 3 3 4" << std::endl;
  stream << "1 5 2 6 4 7" << std::endl;
  stream << "3 8 5 9 6 1" << std::endl;
  stream << "4 2 6 3" << std::endl;
  stream << "4 4 5 5" << std::endl;
}


static void writeSparse(
    std::string const & testFile)
{
  std::fstream stream(testFile,std::fstream::out | std::fstream::trunc);

  stream << "1 1 2 2" << std::endl;
  stream << "0 3 2 4" << std::endl;
  stream << "0 5 1 6 3 7" << std::endl;
  stream << "2 8 4 9 5 1" << std::endl;
  stream << "3 2 5 3" << std::endl;
  stream << "3 4 4 5" << std::endl;
}


static void readSparse(
    std::string const & testFile)
{
  GraphInHandle handle(testFile);

  wildriver_dim_t nvtxs;
  wildriver_ind_t nedges;
  int nvwgt;
  bool ewgt;
  handle.getInfo(nvtxs, nedges, nvwgt, ewgt);

  testEquals(nvtxs,6);
  testEquals(nedges,14);
  testEquals(nvwgt,0);
  testEquals(ewgt,true);

  std::unique_ptr<wildriver_ind_t[]> xadj(new wildriver_ind_t[nvtxs+1]);
  std::unique_ptr<wildriver_dim_t[]> adjncy(new wildriver_dim_t[nedges]);
  std::unique_ptr<wildriver_val_t[]> vwgt(new wildriver_val_t[nvtxs]);
  std::unique_ptr<wildriver_val_t[]> adjwgt(new wildriver_val_t[nedges]);

  handle.readGraph(xadj.get(), adjncy.get(), vwgt.get(), adjwgt.get(), nullptr);

  // test xadj
  testEquals(xadj[0],0);
  testEquals(xadj[1],2);
  testEquals(xadj[2],4);
  testEquals(xadj[3],7);
  testEquals(xadj[4],10);
  testEquals(xadj[5],12);
  testEquals(xadj[6],14);

  // test vwgt
  testEquals(vwgt[0],1);
  testEquals(vwgt[1],1);
  testEquals(vwgt[2],1);
  testEquals(vwgt[3],1);
  testEquals(vwgt[4],1);
  testEquals(vwgt[5],1);

  // test adjncy
  testEquals(adjncy[0],1);
  testEquals(adjncy[1],2);
  testEquals(adjwgt[0],1);
  testEquals(adjwgt[1],2);

  testEquals(adjncy[2],0);
  testEquals(adjncy[3],2);
  testEquals(adjwgt[2],3);
  testEquals(adjwgt[3],4);

  testEquals(adjncy[4],0);
  testEquals(adjncy[5],1);
  testEquals(adjncy[6],3);
  testEquals(adjwgt[4],5);
  testEquals(adjwgt[5],6);
  testEquals(adjwgt[6],7);

  testEquals(adjncy[7],2);
  testEquals(adjncy[8],4);
  testEquals(adjncy[9],5);
  testEquals(adjwgt[7],8);
  testEquals(adjwgt[8],9);
  testEquals(adjwgt[9],1);

  testEquals(adjncy[10],3);
  testEquals(adjncy[11],5);
  testEquals(adjwgt[10],2);
  testEquals(adjwgt[11],3);

  testEquals(adjncy[12],3);
  testEquals(adjncy[13],4);
  testEquals(adjwgt[12],4);
  testEquals(adjwgt[13],5);
}


void Test::run()
{
  // generate test metis file
  std::string metisFile("./GraphInHandle_test.graph");
  writeMetis(metisFile);
  readSparse(metisFile);

  std::string csrFile("./GraphInHandle_test.csr");
  writeSparse(csrFile);
  readSparse(csrFile);
}




}
