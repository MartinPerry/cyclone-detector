#ifndef MY_MACROS_H
#define MY_MACROS_H

#ifdef _MSC_VER
#	ifndef my_fopen 
#		define my_fopen(a, b, c) fopen_s(a, b, c)	
#	endif	
#	ifndef my_fseek 
#		define my_fseek(a, b, c) _fseeki64(a, b, c)	
#	endif
#	ifndef my_ftell
#		define my_ftell(a) _ftelli64(a)	
#	endif
#else
#	ifndef my_fopen 
#		define my_fopen(a, b, c) (*a = fopen(b, c))
#	endif	
#	ifndef my_fseek 
#		define my_fseek(a, b, c) fseeko64(a, b, c)	
#	endif
#	ifndef my_ftell
#		define my_ftell(a) ftello64(a)	
#	endif
#endif


#ifndef SAFE_DELETE
#   define SAFE_DELETE(a) {if (a != nullptr) { delete   a; a = nullptr; }};
#endif

#ifndef SAFE_DELETE_ARRAY
#   define SAFE_DELETE_ARRAY(a) {if (a != nullptr) { delete[] a; a = nullptr; }};
#endif


#endif
