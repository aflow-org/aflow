// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
/// @file
/// @brief This file contains the preprocessor macros to ensure a proper instantiation of all aurostd functions.
/// @authors
/// @mod{HE,20230622,auto template system}
/// @note the macro iterator logic is based on a blog post by Jonathan Heathcote
/// @see
/// @xlink{Blog about macro iterator,http://jhnet.co.uk/articles/cpp_magic#turning-recursion-into-an-iterator}
/// source code examples next to the macros #AST_GEN_1 #AST_GEN_2 #AST_GEN_3

#ifndef _AUROSTD_AUTOMATIC_TEMPLATE_H_
#define _AUROSTD_AUTOMATIC_TEMPLATE_H_

// utype lists
#define AST_UTYPE_BASE_NUM int, uint, float, double
#define AST_UTYPE_NUM int, unsigned int, long int, long unsigned int, long long int, long long unsigned int, float, double, long double
#define AST_UTYPE_FLOAT float, double, long double
#define AST_UTYPE_INT int, long int, long long int
#define AST_UTYPE_AINT int, unsigned int, long int, unsigned long int, long long int, unsigned long long int
#define AST_UTYPE_UINT unsigned int, unsigned long int, unsigned long long int

// Nested evaluations
#define AST_FIRST(a, ...) a
#define AST_SECOND(a, b, ...) b

#define AST_EMPTY()

#define AST_EVAL(...) AST_EVAL1024(__VA_ARGS__)
#define AST_EVAL1024(...) AST_EVAL512(AST_EVAL512(__VA_ARGS__))
#define AST_EVAL512(...) AST_EVAL256(AST_EVAL256(__VA_ARGS__))
#define AST_EVAL256(...) AST_EVAL128(AST_EVAL128(__VA_ARGS__))
#define AST_EVAL128(...) AST_EVAL64(AST_EVAL64(__VA_ARGS__))
#define AST_EVAL64(...) AST_EVAL32(AST_EVAL32(__VA_ARGS__))
#define AST_EVAL32(...) AST_EVAL16(AST_EVAL16(__VA_ARGS__))
#define AST_EVAL16(...) AST_EVAL8(AST_EVAL8(__VA_ARGS__))
#define AST_EVAL8(...) AST_EVAL4(AST_EVAL4(__VA_ARGS__))
#define AST_EVAL4(...) AST_EVAL2(AST_EVAL2(__VA_ARGS__))
#define AST_EVAL2(...) AST_EVAL1(AST_EVAL1(__VA_ARGS__))
#define AST_EVAL1(...) __VA_ARGS__

#define AST_DEFER1(m) m AST_EMPTY()
#define AST_DEFER2(m) m AST_EMPTY AST_EMPTY()()
#define AST_DEFER3(m) m AST_EMPTY AST_EMPTY AST_EMPTY()()()
#define AST_DEFER4(m) m AST_EMPTY AST_EMPTY AST_EMPTY AST_EMPTY()()()()

#define AST_IS_AST_PROBE(...) AST_SECOND(__VA_ARGS__, 0)
#define AST_PROBE() ~, 1

#define AST_CAT(a,b) a ## b

#define AST_NOT(x) AST_IS_AST_PROBE(AST_CAT(_AST_NOT_, x))
#define _AST_NOT_0 AST_PROBE()

#define AST_BOOL(x) AST_NOT(AST_NOT(x))

#define AST_IF_ELSE(condition) _AST_IF_ELSE(AST_BOOL(condition))
#define _AST_IF_ELSE(condition) AST_CAT(_IF_, condition)

#define _IF_1(...) __VA_ARGS__ _IF_1_ELSE
#define _IF_0(...)             _IF_0_ELSE

#define _IF_1_ELSE(...)
#define _IF_0_ELSE(...) __VA_ARGS__

#define AST_HAS_ARGS(...) AST_BOOL(AST_FIRST(_AST_END_OF_ARGUMENTS_ __VA_ARGS__)())
#define AST_HAS_ARGS_2(skip1, ...) AST_BOOL(AST_FIRST(_AST_END_OF_ARGUMENTS_ __VA_ARGS__)())
#define AST_HAS_ARGS_3(skip1, skip2, ...) AST_BOOL(AST_FIRST(_AST_END_OF_ARGUMENTS_ __VA_ARGS__)())
#define _AST_END_OF_ARGUMENTS_() 0
#define AST_INVOKE(macro, ...) macro(__VA_ARGS__)

// logic for iterations up to 3 types
// single utype
#define AST_MAP(m, first, ...)             \
  m(first)                                 \
  AST_IF_ELSE(AST_HAS_ARGS(__VA_ARGS__))(  \
    AST_DEFER2(_AST_MAP)()(m, __VA_ARGS__) \
  )(                                       \
    /* Do nothing, just terminate */       \
  )
#define _AST_MAP() AST_MAP

// dual utype
#define AST_MAP2_1(m, first, ...)             \
  m(first, first)                             \
  AST_IF_ELSE(AST_HAS_ARGS(__VA_ARGS__))(     \
    AST_DEFER2(_AST_MAP2_1)()(m, __VA_ARGS__) \
  )(                                          \
    /* Do nothing, just terminate */          \
  )
#define _AST_MAP2_1() AST_MAP2_1

#define AST_MAP2_2(m, first,...)                    \
  AST_DEFER3(AST_MAP2_2_SUB(m, first, __VA_ARGS__)) \
  AST_IF_ELSE(AST_HAS_ARGS_2(__VA_ARGS__))(         \
    AST_DEFER2(_AST_MAP2_2)()(m, __VA_ARGS__)       \
  )(                                                \
    /* Do nothing, just terminate */                \
  )
#define _AST_MAP2_2() AST_MAP2_2

#define AST_MAP2_2_SUB(m, first, second, ...)           \
  m(first, second)                                      \
  m(second, first)                                      \
  AST_IF_ELSE(AST_HAS_ARGS(__VA_ARGS__))(               \
    AST_DEFER2(_AST_MAP2_2_SUB)()(m,first, __VA_ARGS__) \
  )(                                                    \
    /* Do nothing, just terminate */                    \
  )
#define _AST_MAP2_2_SUB() AST_MAP2_2_SUB

// Three types
#define AST_MAP3_1(m, first, ...)             \
  m(first, first, first)                      \
  AST_IF_ELSE(AST_HAS_ARGS(__VA_ARGS__))(     \
    AST_DEFER2(_AST_MAP3_1)()(m, __VA_ARGS__) \
  )(                                          \
    /* Do nothing, just terminate */          \
  )
#define _AST_MAP3_1() AST_MAP3_1

#define AST_MAP3_2(m, first, ...)                   \
  AST_DEFER4(AST_MAP3_2_SUB(m, first, __VA_ARGS__)) \
  AST_IF_ELSE(AST_HAS_ARGS_3(__VA_ARGS__))(         \
    AST_DEFER3(_AST_MAP3_2)()(m, __VA_ARGS__)       \
  )(                                                \
    /* Do nothing, just terminate */                \
  )
#define _AST_MAP3_2() AST_MAP3_2


#define AST_MAP3_2_SUB(m, first, second, ...)                \
  m(first, second, second)                                   \
  m(second, first, second)                                   \
  m(first, second, second)                                   \
  m(second, first, first)                                    \
  m(first, second, first)                                    \
   AST_IF_ELSE(AST_HAS_ARGS(__VA_ARGS__))(                   \
    AST_DEFER3(AST_MAP3_3_SUB(m, first, second,__VA_ARGS__)) \
    AST_DEFER2(_AST_MAP3_2_SUB)()(m,first, __VA_ARGS__)      \
  )(                                                         \
    /* Do nothing, just terminate */                         \
  )
#define _AST_MAP3_2_SUB() AST_MAP3_2_SUB


#define AST_MAP3_3_SUB(m, first, second, third, ...)             \
    m(first, second, third)                                      \
    m(first, third, second)                                      \
    m(second, first, third)                                      \
    m(second, third, first)                                      \
    m(third, first, second)                                      \
    m(third, second, first)                                      \
    AST_IF_ELSE(AST_HAS_ARGS(__VA_ARGS__))(                      \
    AST_DEFER2(_AST_MAP3_3_SUB)()(m, first, second, __VA_ARGS__) \
  )(                                                             \
    /* Do nothing, just terminate */                             \
  )
#define _AST_MAP3_3_SUB() AST_MAP3_3_SUB


/// @brief autogenerate 1D code based on #AST_TEMPLATE
/// @param type_selection use one of the `AST_UTYPE_*`
/// @note always define the #AST_TEMPLATE before using this macro function
/// @code{cpp}
///  #define AST_TEMPLATE(utype) template class xvector<utype>;
///    AST_GEN_1(AST_UTYPE_NUM)
///  #undef AST_TEMPLATE
/// @endcode
/// the preprocessor will generate (without newlines):
/// @code{cpp}
/// template class xvector<int>;
/// template class xvector<unsigned int>;
/// template class xvector<long int>;
/// template class xvector<long unsigned int>;
/// template class xvector<long long int>;
/// template class xvector<long long unsigned int>;
/// template class xvector<float>;
/// template class xvector<double>;
/// template class xvector<long double>;
/// @endcode
#define AST_GEN_1(type_selection) \
      AST_EVAL(AST_INVOKE(AST_MAP, AST_TEMPLATE, type_selection))

/// @brief autogenerate 2D code based on #AST_TEMPLATE
/// @param type_selection use one of the `AST_UTYPE_*`
/// @note always define the #AST_TEMPLATE before using this macro function
/// @code{cpp}
///  #define AST_TEMPLATE(atype, btype) template xvector<atype> xvector2utype(const xvector<btype>& a);
///    AST_GEN_2(AST_UTYPE_FLOAT)
///  #undef AST_TEMPLATE
/// @endcode
/// the preprocessor will generate (without newlines):
/// @code{cpp}
/// template xvector<float> xvector2utype(const xvector<float>& a); 
/// template xvector<double> xvector2utype(const xvector<double>& a);
/// template xvector<long double> xvector2utype(const xvector<long double>& a);
/// template xvector<float> xvector2utype(const xvector<double>& a);
/// template xvector<double> xvector2utype(const xvector<float>& a);
/// template xvector<float> xvector2utype(const xvector<long double>& a);
/// template xvector<long double> xvector2utype(const xvector<float>& a);
/// template xvector<double> xvector2utype(const xvector<long double>& a);
/// template xvector<long double> xvector2utype(const xvector<double>& a);
/// @endcode
#define AST_GEN_2(type_selection) \
      AST_EVAL(AST_INVOKE(AST_MAP2_1, AST_TEMPLATE, type_selection)) \
      AST_EVAL(AST_INVOKE(AST_MAP2_2, AST_TEMPLATE, type_selection))

/// @brief autogenerate 3D code based on #AST_TEMPLATE
/// @param type_selection use one of the `AST_UTYPE_*`
/// @note The maximal amount of generated lines should not be bigger than 1024
/// (10 types for AST_GEN_3), while larger system work they are very slow
///
/// @note always define the #AST_TEMPLATE before using this macro function
/// @code{cpp}
///  #define AST_TEMPLATE(atype, btype, ctype) template void quicksort3(unsigned long n, xvector<atype>&, xvector<btype>&, xvector<ctype>&);
///    AST_GEN_3(AST_UTYPE_FLOAT)
///  #undef AST_TEMPLATE
/// @endcode
/// the preprocessor will generate (without newlines):
/// @code
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<double>&, xvector<double>&);
/// template void quicksort3(unsigned long n, xvector<double>&, xvector<float>&, xvector<double>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<double>&, xvector<double>&);
/// template void quicksort3(unsigned long n, xvector<double>&, xvector<float>&, xvector<float>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<double>&, xvector<float>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<double>&, xvector<long double>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<long double>&, xvector<double>&);
/// template void quicksort3(unsigned long n, xvector<double>&, xvector<float>&, xvector<long double>&);
/// template void quicksort3(unsigned long n, xvector<double>&, xvector<long double>&, xvector<float>&);
/// template void quicksort3(unsigned long n, xvector<long double>&, xvector<float>&, xvector<double>&);
/// template void quicksort3(unsigned long n, xvector<long double>&, xvector<double>&, xvector<float>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<long double>&, xvector<long double>&);
/// template void quicksort3(unsigned long n, xvector<long double>&, xvector<float>&, xvector<long double>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<long double>&, xvector<long double>&);
/// template void quicksort3(unsigned long n, xvector<long double>&, xvector<float>&, xvector<float>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<long double>&, xvector<float>&);
/// template void quicksort3(unsigned long n, xvector<float>&, xvector<float>&, xvector<float>&);
/// template void quicksort3(unsigned long n, xvector<double>&, xvector<double>&, xvector<double>&);
/// template void quicksort3(unsigned long n, xvector<long double>&, xvector<long double>&, xvector<long double>&);
/// @endcode
#define AST_GEN_3(type_selection) \
      AST_EVAL(AST_INVOKE(AST_MAP3_2, AST_TEMPLATE, type_selection)) \
      AST_EVAL(AST_INVOKE(AST_MAP3_1, AST_TEMPLATE, type_selection)) \


#endif //_AUROSTD_AUTOMATIC_TEMPLATE_H_
