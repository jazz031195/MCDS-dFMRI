//!  Constants ==============================================================================/
/*!
*   \details   Defines constants, macros, and values to be used through the implementation.
*   \author    Jonathan Rafael
*   \date      November 2016
*
============================================================================================*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdio.h>
#include <string>
#include <limits>
#define PRINT_EXPETED_TIME 1


#define SH_BG_RED               "\033[1;41m"
#define SH_BG_LIGHT_YELLOW      "\033[1;43m"
#define SH_FG_LIGHT_RED         "\033[1;31m"
#define SH_FG_LIGHT_YELLOW      "\033[1;33m"
#define SH_FG_PURPLE            "\033[1;35m"
#define SH_FG_GRAY              "\033[0;47m"
#define SH_FG_GREEN             "\033[0;32m"
#define SH_DEFAULT              "\033[0m"

#define VERSION_ID              "1.50.000"

constexpr double m_to_mm = 1e3;                     /*!< meters to milimeters constant                              */
constexpr double s_to_ms = 1e3;                     /*!< seconds to milisecodns constant                            */
constexpr double EPS_VAL = 1e-12;                   /*!< numerical epsilon value                                    */
constexpr double m2_to_mm2 = 1e6;                   /*!< squared meters to squared milimeters                       */
constexpr double giro = 267.51525e3;                /*!< Gyromagnetic radio given in rad/(ms*T)                     */
constexpr double DIFF_CONST = 2.02e-7;              /*!< Default diffusion coeficient                               */
constexpr double barrier_tickness = 1e-6;           /*!< Defines the defaul tickness of a obstacle barrier          */
constexpr double max_number_bouncings = 10000.0;    /*!< Defines the maximum number of bouncing a particle can make.*/
constexpr double triangle_eps = 1e-10;              /*!< Extra area for the PLY triangles. Help to numerical erros  */
constexpr unsigned max_rejections =25;              /*!< Max number of tries to unstuck a particle in a single step */
constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();  /*!< numerical infinity value           */

#ifdef _WIN64
typedef unsigned int ulong;
typedef unsigned int uint;
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#elif __APPLE__
typedef unsigned long ulong;
typedef unsigned int uint;
#endif

/*!< Colision optimization parameters                                                                           */
#define PRECISE_T_MIN_D 0

/*!< outher collision sphere relative size      DEPRECATED                                                      */
//const double outher_col_dist_factor   = 5.0;

/*!< inner collision sphere relative size                                                                       */
constexpr double inner_col_dist_factor      = 0.25;


#endif // CONSTANTS_H
