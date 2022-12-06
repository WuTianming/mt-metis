/**
 * @file Util_test.cpp
 * @brief Test for Util class.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2017
 * @version 1
 *
 */




#include <iostream>
#include <fstream>
#include <memory>

#include "Util.hpp"
#include "DomTest.hpp"




using namespace WildRiver;




namespace DomTest
{


static void splitTest()
{
  std::string tripletLine = "9 8 1.0";
  std::vector<std::string> chunks = Util::split(tripletLine);
  testEquals(chunks.size(), 3);
  testEquals(chunks[0].compare("9"), 0);
  testEquals(chunks[1].compare("8"), 0);
  testEquals(chunks[2].compare("1.0"), 0);

  std::string spaceAndTabs = "hello\t world";
  chunks = Util::split(spaceAndTabs);
  testEquals(chunks.size(), 2);
  testEquals(chunks[0].compare("hello"), 0);
  testEquals(chunks[1].compare("world"), 0);

  std::string trailingAndLeading = "  hello\t world   ";
  chunks = Util::split(trailingAndLeading);
  testEquals(chunks.size(), 2);
  testEquals(chunks[0].compare("hello"), 0);
  testEquals(chunks[1].compare("world"), 0);
}



void Test::run()
{
  splitTest();
}




}
