/**
 * @file PlainVectorFile_test.cpp
 * @brief Unit tests for the TextFile class. 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2017
 * @version 1
 */




#include "PlainVectorFile.hpp"
#include "DomTest.hpp"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>




using namespace WildRiver;




namespace DomTest
{


namespace
{

const std::vector<std::string> LINES{
  "my first line",
  "my second line",
  "my third line"
};

}


static void writeTest(
    std::string const & testFile)
{
  TextFile f(testFile);

  testTrue(!f.isOpenWrite());
  testTrue(!f.isOpenRead());

  f.openWrite();

  testTrue(f.isOpenWrite());
  testTrue(!f.isOpenRead());

  for (std::string const & line : LINES) {
    f.setNextLine(line);
  }
}


static void readTest(
    std::string const & testFile)
{
  TextFile f(testFile);

  testTrue(!f.isOpenWrite());
  testTrue(!f.isOpenRead());

  f.openRead();

  testTrue(!f.isOpenWrite());
  testTrue(f.isOpenRead());

  size_t i = 0;
  std::string line;
  while (f.nextLine(line)) {
    testEquals(line, LINES[i]);
    ++i;
  }
  testEquals(i, 3);
}


void Test::run()
{
  std::string testFile("./test.txt");

  writeTest(testFile);
  readTest(testFile);

  Test::removeFile(testFile);
}




}
