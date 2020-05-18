/*!
 * \file CRTSNError.h
 * \brief Class to handle Error Messages
 * \author S. Schotthoefer
 */

#ifndef CRTSNERROR_H
#define CRTSNERROR_H

#include <iostream>

#include "spdlog/spdlog.h"

class ErrorMessages
{
public:
  ErrorMessages();

  inline static void Error(std::string ErrorMsg, std::string FunctionName){
      auto log = spdlog::get( "event" );
      // if (Rank == 0){ //For MPI implementation later
      log->error( "\n" );
      log->error( "Error in \"{0}\": ", FunctionName );
      log->error( "-------------------------------------------------------------------------" );
      log->error( ErrorMsg );
      log->error( "------------------------------ Error Exit -------------------------------" );
      log->error( "\n" );
    //}
    exit(EXIT_FAILURE);
  }

  inline static void OptionNotSetError(std::string OptionName, std::string FunctionName){
    //if (Rank == 0){
      std::cout << std::endl << std::endl;
      std::cout << "Error in \"" << FunctionName << "\": " << std::endl;
      std::cout <<  "-------------------------------------------------------------------------" << std::endl;
      std::cout << "Option "<< OptionName << " not set. Please check your config file." << std::endl;
      std::cout <<  "------------------------------ Error Exit -------------------------------" << std::endl;
      std::cout << std::endl << std::endl;
    //}
    exit(EXIT_FAILURE);
  }
};

#endif // CRTSNERROR_H


/* Depending on the compiler, define the correct macro to get the current function name */

#if defined(__GNUC__) || (defined(__ICC) && (__ICC >= 600))
# define CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
# define CURRENT_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600))
# define CURRENT_FUNCTION __FUNCTION__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
# define CURRENT_FUNCTION __func__
#elif defined(__cplusplus) && (__cplusplus >= 201103)
# define CURRENT_FUNCTION __func__
#else
# define CURRENT_FUNCTION "(unknown)"
#endif
