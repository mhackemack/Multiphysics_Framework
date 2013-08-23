//===========================================================================//
/*!
 * \class Error
 * \brief
 *
 * Routine to print error statements and end program if necessary.
 */
//===========================================================================//

#ifndef ERROR_H
#define ERROR_H

#include "../Global/GlobalHeaders.h"

class Error {

  public:
  // Constructor
  Error() {};
  
  // Destructor
  ~Error() {};

  // Routines
  void warning(char message_);
  void fatal_error(char message_);

  private:

};


void Error::warning(char message) {
  
  std::cout << "===> Warning:  " << message;

  return;
};

void Error::fatal_error(char message) {
  
  std::cout << "===> Fatal Error:  " << message << "\n";
  std::cout << "\n";

  return;
};


#endif
