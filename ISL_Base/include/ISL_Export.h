#ifndef ISLExport_H
#define ISLExport_H

/**
* @brief	Linux �� windows ���ܹ�ͨ�õĿ�����
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