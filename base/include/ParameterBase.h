/** \file ParameterBase.h
 *  \brief Implimentations for the Base Class of Parameters. Several classes in the
 *  code are built on top of this parameters class. That gives those classes a set of
 *  parameters that can be accessed "easily". A few typedefs are added here also, but
 *  many of them might be obsolete (TODO)
 **/

#ifndef SS_PARAMETERSBASE_H
#define SS_PARAMETERSBASE_H

#include "typedefs.h"

/**
 *  The class _SS_Parameters is the base class for
 *  the several of the classes in the solver selecter. It is basically
 *  a wrapper for the std::map< std::string, std::pair\< string, string \>.
 **/
class
    _SS_Parameters
{
public:

    /** Constructor: Build the _SS_Parameters */
    _SS_Parameters();

    /**
     * Get the keys that define the parameters set in this class
     **/
    _SS_ErrorFlag GetKeys( std::set<std::string> &keys /**< output, key for each parameter */) const ;

    /**
     * Add a new parameter to the set. OVerrites if already exists
     **/
    _SS_ErrorFlag Add ( const std::string &key, /**< input, parameter key */
                        const std::string &value, /**< input, value for parameter */
                        const std::string &description /**< input, description of the parameter*/
                      ) ;

    /**
     * Set a parameter in the parameter set. If the value does not exist in the set, nothing
     * is done.
     **/
    _SS_ErrorFlag Set ( const std::string &key,       /**< input, parameter key */
                        const std::string &value   /**< input, parameter value */
                      );

    /**
     * Get the values of the a parameter in the parameter set.
     **/
    _SS_ErrorFlag Get ( const std::string &key, /***< input, parameter key */
                        std::string &value    /**< output, parameter value1 */
                      ) const;
    /**
     * Clear the parameters
     **/
    virtual _SS_ErrorFlag Clear();

    /**
     * Print the parameters
     **/
    virtual _SS_ErrorFlag Print() const ;

private:
    std::map< std::string, std::pair<std::string, std::string > > parameters; /**< The main parameter set */
};

#endif

