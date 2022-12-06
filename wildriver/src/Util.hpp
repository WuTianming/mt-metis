/**
 * @file Util.hpp
 * @brief Utility functions.
 * @author Dominique LaSalle <dominique@solidlake.com>
 * Copyright 2017
 * @version 1
 * @date 2017-08-03
 */




#ifndef WILDIRVER_UTIL_HPP
#define WILDIRVER_UTIL_HPP




#include <cstring>
#include <string>
#include <vector>




namespace WildRiver
{



class Util
{
  public:
    /**
     * @brief Split a string into chunks separated by 'delims'.
     *
     * @param input The string to split.
     * @param delims The characters to split the string on.
     *
     * @return A vector non-empty strings.
     */
    static std::vector<std::string> split(
        std::string const & input,
        char const * const delims = " \t") noexcept
    {
      // create a character map of all 0's
      bool map[256] = {false};

      // mark my delimiters
      const size_t n = strlen(delims);
      for (size_t i = 0; i < n; ++i) {
        map[(int)delims[i]] = true;
      }

      // fill in chunks
      std::vector<std::string> chunks;
      size_t last = 0;
      for (size_t i = 0; i < input.length(); ++i) {
        if (map[static_cast<uint8_t>(input[i])]) {
          if (i > last) {
            chunks.emplace_back(input.substr(last, i-last));
          }
          last = i+1;
        }
      }
      // check for last chunk
      if (last < input.length()) {
        chunks.emplace_back(input.substr(last, input.length()-last));
      }

      return chunks;
    }




};




}




#endif
