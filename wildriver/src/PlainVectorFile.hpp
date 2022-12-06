/**
 * @file PlainVectorFile.hpp
 * @brief Concrete class text vector files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-02-08
 */





#ifndef WILDIRVER_PLAINVECTORFILE_HPP
#define WILDIRVER_PLAINVECTORFILE_HPP




#include <string>


#include "IVectorReader.hpp"
#include "IVectorWriter.hpp"
#include "TextFile.hpp"




namespace WildRiver
{


class PlainVectorFile :
    public IVectorReader,
    public IVectorWriter
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
     * @brief Open a plain vector file for reading or writing.
     *
     * @param name The name of the file.
     */
    PlainVectorFile(
        std::string const & name);


    /**
     * @brief Close file and associated memory.
     */
    ~PlainVectorFile();


    /**
    * @brief Set the size of the vector.
    *
    * @param size The new size of the vector.
    */
    virtual void setSize(
        const ind_t size) override;


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    virtual ind_t getSize() override;


    /**
     * @brief Read the values of the vector. 
     *
     * @param vals The values in the vector (output).
     * @param progress The variable to update as teh vector is loaded (can be
     * null).
     */
    virtual void read(
        val_t * vals,
        double * progress) override;



    /**
     * @brief Save the values of the vector.
     *
     * @param vals The values in the vector.
     * @param progress The variable to update as the vector is saved (can be
     * null).
     */
    virtual void write(
        val_t const * vals,
        double * progress) override;



  private:
    /**
    * @brief Underlying text file.
    */
    TextFile m_file;


    /**
    * @brief Size of the vector.
    */
    ind_t m_size;


    /**
     * @brief The string for buffer read lines into.
     */
    std::string m_buffer;


    /**
    * @brief Get the next non-comment line from the file.
    *
    * @param line The line buffer.
    *
    * @return True if a line filled the buffer.
    */
    bool nextNoncommentLine(
        std::string & line);




};




}




#endif
