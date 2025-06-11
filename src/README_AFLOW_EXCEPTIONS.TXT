****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE: materials.duke.edu/AFLOW/README_AFLOW_EXCEPTIONS.TXT
****************************************************************************************************
Written by Marco Esters, Duke University (2018).
****************************************************************************************************

This document describes exceptions in AFLOW. It is primarily intended for developers, but users may
 find useful information in here as well.

Exceptions are thrown after an error in AFLOW occurs. They return a message with the function in
which the error occurred, the kind of error that occurred, and a detailed description of the error.
Exceptions in AFLOW are designed to terminate an AFLOW run and are preferable to hard exits. If an
error should not cause AFLOW to exit, please use warnings instead.

AFLOW exceptions are part of AUROSTD and can be accessed by including the aurostd or the
aurostd_exception header file.

----------------------------------------------------------------------------------------------------

Table of Contents:

1.     Using Exceptions in AFLOW
1.1.     Message Structure
1.2.     How to Throw an Exception
1.3.     Error Types and Error Codes

2.     Descriptions of the Exceptions
2.1.     List of Exceptions with Error Codes
2.2.     Detailed Descriptions of the Exceptions

3.     Changing and Adding New Exceptions

----------------------------------------------------------------------------------------------------

1. Using Exceptions in AFLOW
1.1. Message Structure
An exception contains:

  - The name of the function that throws the exception (f_name).
  - An error code (error_code).
  - An error type (error_type) with a more detailed error.
  - A message accompanying the exception.

The output is structured as follows:

"f_name ERROR error_code: error_type(error) - message"

After printing the message, AFLOW will terminate, using the error code as the exit code.

1.2. How to Throw an Exception
To use AFLOW exceptions, the aurostd_exception.h file must be included in the cpp file. The
exception can be thrown using:

    throw aurostd::xerror(function, message, error_code);

The parameters are:

  - function (std::string or std::stringstream): the name of the function from which the exception
        is thrown.
  - message (std::string or std::stringstream): a detailed message about the error that caused the
        exception.
  - error_code (int): the error code (defaults to 1).

For more information on error codes, see 1.3, and for more information on which exception to throw,
see 2.2.

1.3 Error Types and Error Codes
AFLOW exceptions are assigned to error types, and each type contains more specific errors (see
section 2. for details). The first digit of the error code corresponds to the error type, the second
digit to the error associated with the error type. For example, all errors associated with input in
aflow.in will have error codes between 10 and 19.

The following error codes are reserved:

  - Error code 0.
  - Error codes without error type.
  - Error codes for generic errors.

Error code 0:
  Since the error codes are also used as exit codes, 0 cannot be used for error codes as it implies
  that AFLOW terminated normally.

Error codes without error types:
  Error codes without error types refer to general AFLOW errors that cannot be classified with an
  error type. The error codes reserved for these types are 1 - 9.

Generic errors:
  Generic errors are errors that can be described using an error type, but cannot be specified with
  a more detailed category and do not warrant the creation of a new error classification. Their
  error codes end with 0. For example, a method requiring a uniform supercell expansion, but
  receiving a structure with a non-uniform expansion would be such a generic runtime error.

When using error codes, it is advisable to use the names of the error code constants shown in 2.1.
instead of the numbers themselves. It will make the code more legible and changes much easier to
implement.

2. Description of the Exceptions
2.1. List of Exceptions with Error Codes

This list can be output by running aflow -e|--errors|--exceptions.

--------------------------------------------------------------------------------
Code     Error Type          Error                     Name of Constant
--------------------------------------------------------------------------------
   1     N/A                 Generic error             _GENERIC_ERROR_
   2                         Illegal error code        _ILLEGAL_CODE_
  10     Input Error         generic                   _INPUT_ERROR_
  11                         unknown flag              _INPUT_UNKNOWN_
  12                         missing flag              _INPUT_MISSING_
  13                         input ambiguous           _INPUT_AMBIGUOUS_
  14                         illegal parameter         _INPUT_ILLEGAL_
  15                         number of parameters      _INPUT_NUMBER_
  20     File Error          generic                   _FILE_ERROR_
  21                         file not found            _FILE_NOT_FOUND_
  22                         wrong format              _FILE_WRONG_FORMAT_
  23                         file corrupt              _FILE_CORRUPT_
  30     Value Error         generic                   _VALUE_ERROR_
  31                         illegal value             _VALUE_ILLEGAL_
  32                         out of range              _VALUE_RANGE_
  40     Index Error         generic                   _INDEX_ERROR_
  41                         illegal value             _INDEX_ILLEGAL_
  42                         out of bounds             _INDEX_BOUNDS_
  43                         mismatch                  _INDEX_MISMATCH_
  50     Runtime Error       generic                   _RUNTIME_ERROR_
  51                         not initialized           _RUNTIME_INIT_
  52                         SQL error                 _RUNTIME_SQL_
  53                         busy                      _RUNTIME_BUSY_
  54                         external command missing  _RUNTIME_EXTERNAL_MISS_
  55                         external command failed   _RUNTIME_EXTERNAL_FAIL_
  56                         HTTP error                _RUNTIME_HTTP_
  60     Allocation Error    generic                   _ALLOC_ERROR_
  61                         could not allocate        _ALLOC_ALLOCATE_
  62                         insufficient memory       _ALLOC_INSUFFICIENT_
-------------------------------------------------------------------------------

2.2. Detailed Descriptions of the Exceptions

[1-9] No error type: Errors that cannot be classified under an error type.

   1. Generic error (default): Generic AFLOW error. Throw when nothing else fits.

   2. Illegal error code: Thrown when the exception is supplied with an invalid error code.

[10 - 19] Input Error: Errors within the user input, e.g. in aflow.in or in command line arguments.
  This error type should only be used for illegal input and not for errors that result from bad
  parameters.

  For example, tolerance = 1E-18 is legal input, but some programs may require this value to be
  higher. In that case, a Value Error (out of range) should be raised instead. On the other hand,
  tolerance = dog is illegal input and should thus result in an Input Error (illegal parameter).

  11. Unknown flag: Used when aflow.in contains a flag that is not a valid AFLOW input flag.

  12. Missing flag: Used when aflow.in is missing a required flag.

  13. Input ambiguous: Used when aflow.in contains two contradictory or mutually exclusive entries.

  14. Illegal parameter: Used when a flag in aflow.in is supplied with an invalid parameter. For
                         example, a flag that requires a positive number, but is supplied with a
                         negative number or with non-numeric characters, would require this
                         exception to be thrown.

  15. Number of parameters: Used when the number of parameters for input flags are inconsistent.
                            For example, when the number of entries in MAGMOM is not the same as
                            the number of atoms in the structure, this exception should be thrown.

[20 - 29] File Error: Errors associated with files.

  21. File not found: Used when the requested file is not found.

  22. Wrong format: Used when the requested file is not in the correct format, e.g. when AFLOW
                    tries to read a text file, but the file is binary.

  23. File corrupt: Used when the requested file is incomplete or corrupt, e.g. an incomplete APL
                    hibernation (APL.xml) file.

[30 - 39] Value Error:

  31. Illegal value: Used when a supplied value is illegal, e.g. a negative number for a parameter
                     that requires positive numbers. This exception should not be thrown when the
                     illegal value is detected while user input is processed. Use Input
                      Error (illegal parameter) instead.

  32. Out of range: Used when a supplied value is too high or too low.

[40 - 49] Index Error:

  41. Illegal value: Used when the value for an index or for the size of an array, vector, etc. is
                     illegal (e.g. a negative number).

  42. Out of bounds: Used when an index is larger than the size of an array, vector, etc.

  43: Mismatch: Used when two objects have different shapes even though the same shape is required,
                e.g. when adding two xmatrix objects.

[50 - 59] Runtime Error:

  51. Not initialized: Used when an uninitialized object is used. The increment operator in xcombos
                       may throw such an error.

  52. SQL Error: Used when SQLite returns an error.

  53. Busy: Used when a required process is busy or when a file is already accessed by another
            process. Can also be used to prevent concurrent execution of a function or process.
  
  54. external command not found: Used when an external subsystem call cannot be found.

  55. external command failed: Used when an external subsystem call fails.

  56. HTTP error: Used when an http request fails.

[61 - 69] Allocation Error: Errors associated with memory allocation.

  61. Could not allocate: Used when requested memory could not be allocated.

  62. Insufficient memory: Used when more memory is requested than available.


3. Changing and Adding New Exceptions
Two std::string arrays found in aurostd_exception.cpp contain the exception descriptors:

  - "error_types" contains the error types.
  - "errors" contains the errors associated with each error type.

To change the name of an existing exception, simply replace the names in these
arrays.

****************************************************************************************************
*                                             WARNING!                                             *
*                                                                                                  *
*               Avoid renaming existing constants! Doing so requires that they need                *
*               to be changed in EVERY file that uses them.                                        *
*                                                                                                  *
****************************************************************************************************

The more exceptions exist, the more difficult it is to keep track of them. This reduces the
legibility of the code and makes implementing meaningful error messages, and thus debugging, more
difficult. Therefore, before adding a new exception, please double-check that the error cannot be
described with an existing exception. If a new exception needs to be added, it should be as generic
as possible and should be applicable to many different scenarios.

To add an error to an existing error type, the following steps are necessary:

  - Replace one of the empty strings in the errors array. Make sure that the subarray maps to the
    correct error type.
  - Add a shorthand name for the error code into the definitions at the top of aurostd_exception.h.
  - Update the table in section 2.1. of this README.
  - Update the banner in aflow.cpp.

To add a new error type, do the following:

  - Add the new error type in the error_types array.
  - Update the _AFLOW_NUM_ERR_TYPES_ at the top of aurostd_exception.cpp.
  - Add a new line with nine (9) empty strings to the errors array. 
  - Update the table in section 2.1. of this README.
  - Update the banner in aflow.cpp.

****************************************************************************************************
*                                             WARNING!                                             *
*                                                                                                  *
*     It is inadvisable to change the order of the error types and errors. If the order is         *
*     changed, all exceptions that use the affected error codes need to be updated. So, only       *
*     append new error types and errors to the existing arrays instead of inserting them.          *
*     However, if the order absolutely has to be changed, make sure that the error codes in        *
*     ALL cpp files are updated as well, especially the definitions in aurostd_exception.h.        *
*     While the named constants are designed to make these changes easier, there is no             *
*     guarantee that everyone uses them!                                                           *
*                                                                                                  *
****************************************************************************************************
****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
