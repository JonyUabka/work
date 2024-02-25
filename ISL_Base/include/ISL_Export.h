#ifndef ISLExport_H
#define ISLExport_H

/**
* @brief	Linux 与 windows 都能够通用的库设置
* @param	No
* @return	No
*/
//#ifdef WIN32
#if defined(WIN32) || defined(WIN64) || defined(_WINDOWS)
    #ifdef MAKE_ISL_LIB
        #define ISL_EXPORT __declspec(dllexport)
    #else
        #define ISL_EXPORT __declspec(dllimport)
    #endif
#else
    #define ISL_EXPORT
#endif


#endif