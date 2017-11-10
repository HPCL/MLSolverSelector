/** \file ParameterBase.C
 *  \brief Implimentations for the Base Class of Parameters. Several classes in the
 *  code are built on top of this parameters class.
 **/

#include "ParameterBase.h"

_SS_Parameters::_SS_Parameters() { } ;

_SS_ErrorFlag _SS_Parameters::GetKeys( std::set<std::string> &keys /**< output, key for each parameter */) const
{
    for (auto it : parameters)
        keys.insert(it.first);
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Parameters::Add ( const std::string &key, /**< input, parameter key */
                                    const std::string &value, /**< input, value for parameter */
                                    const std::string &description /**< input, description of the parameter*/
                                  )
{
    parameters[key] = std::make_pair(value,description);
    return _SS_error_flag;
}


_SS_ErrorFlag _SS_Parameters::Set ( const std::string &key,       /**< input, parameter key */
                                    const std::string &value   /**< input, parameter value */
                                  )
{

    auto search = parameters.find(key);
    if ( search != parameters.end() )
    {
        search->second.first = value;
    }
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Parameters::Get ( const std::string &key, /***< input, parameter key */
                                    std::string &value    /**< output, parameter value1 */
                                  ) const
{
    auto search = parameters.find(key);
    if ( search != parameters.end() )
    {
        value = (search->second).first;
    }
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Parameters::Clear()
{
    parameters.clear();
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Parameters::Print() const
{
    for ( auto it : parameters )
        std::cout << it.first << " : " << it.second.first << " : " << it.second.second << std::endl;
    return _SS_error_flag;
}

