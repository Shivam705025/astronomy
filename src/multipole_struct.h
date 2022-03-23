/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_MULTIPOLE_STRUCT_H
#define SWIFT_MULTIPOLE_STRUCT_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "timeline.h"

/**
 * @brief Alignment of the structure for allocation
 */
#define multipole_align 128

/**
 * @brief Field tensor components at the location of the multipole.
 */
struct grav_tensor {

  /* 0th order terms */
  MyFloat F_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* 1st order terms */
  MyFloat F_100, F_010, F_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 2nd order terms */
  MyFloat F_200, F_020, F_002;
  MyFloat F_110, F_101, F_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 3rd order terms */
  MyFloat F_300, F_030, F_003;
  MyFloat F_210, F_201;
  MyFloat F_120, F_021;
  MyFloat F_102, F_012;
  MyFloat F_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 4th order terms */
  MyFloat F_400, F_040, F_004;
  MyFloat F_310, F_301;
  MyFloat F_130, F_031;
  MyFloat F_103, F_013;
  MyFloat F_220, F_202, F_022;
  MyFloat F_211, F_121, F_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* 5th order terms */
  MyFloat F_005, F_014, F_023;
  MyFloat F_032, F_041, F_050;
  MyFloat F_104, F_113, F_122;
  MyFloat F_131, F_140, F_203;
  MyFloat F_212, F_221, F_230;
  MyFloat F_302, F_311, F_320;
  MyFloat F_401, F_410, F_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Number of gparts interacted through the tree. */
  long long num_interacted_tree;

  /* Number of gparts interacted through the FFT mesh */
  long long num_interacted_pm;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* Total number of gpart this field tensor interacted with */
  long long num_interacted;

  /* Last time this tensor was zeroed */
  integertime_t ti_init;

#endif

  /* Has this tensor received any contribution? */
  char interacted;
};

/**
 * @brief The properties of the gravity multipole.
 */
struct multipole {

  /*! Bulk velocity */
  MyFloat vel[3];

  /*! Maximal velocity along each axis of all #gpart */
  float max_delta_vel[3];

  /*! Minimal velocity along each axis of all #gpart */
  float min_delta_vel[3];

  /*! Maximal co-moving softening of all the #gpart in the mulipole */
  float max_softening;

  /*! Minimal acceleration norm of all the #gpart in the mulipole */
  MyFloat min_old_a_grav_norm;

  /*! Mulipole power for the different orders */
  MyFloat power[SELF_GRAVITY_MULTIPOLE_ORDER + 1];

  /* 0th order term */
  MyFloat M_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* 1st order terms (all 0 since we expand around CoM) */
  // MyFloat M_100, M_010, M_001;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 2nd order terms */
  MyFloat M_200, M_020, M_002;
  MyFloat M_110, M_101, M_011;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 3rd order terms */
  MyFloat M_300, M_030, M_003;
  MyFloat M_210, M_201;
  MyFloat M_120, M_021;
  MyFloat M_102, M_012;
  MyFloat M_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 4th order terms */
  MyFloat M_400, M_040, M_004;
  MyFloat M_310, M_301;
  MyFloat M_130, M_031;
  MyFloat M_103, M_013;
  MyFloat M_220, M_202, M_022;
  MyFloat M_211, M_121, M_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* 5th order terms */
  MyFloat M_005, M_014, M_023;
  MyFloat M_032, M_041, M_050;
  MyFloat M_104, M_113, M_122;
  MyFloat M_131, M_140, M_203;
  MyFloat M_212, M_221, M_230;
  MyFloat M_302, M_311, M_320;
  MyFloat M_401, M_410, M_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

  /* Total number of gpart in this multipole */
  long long num_gpart;

#endif
};

/**
 * @brief The multipole expansion of a mass distribution.
 */
struct gravity_tensors {

  union {

    /*! Linking pointer for "memory management". */
    struct gravity_tensors *next;

    /*! The actual content */
    struct {

      /*! Field tensor for the potential */
      struct grav_tensor pot;

      /*! Multipole mass */
      struct multipole m_pole;

      /*! Centre of mass of the matter dsitribution */
      double CoM[3];

      /*! Centre of mass of the matter dsitribution at the last rebuild */
      double CoM_rebuild[3];

      /*! Upper limit of the CoM<->gpart distance */
      double r_max;

      /*! Upper limit of the CoM<->gpart distance at the last rebuild */
      double r_max_rebuild;
    };
  };
} SWIFT_STRUCT_ALIGN;

/**
 * @brief Values returned by the M2P kernel.
 */
struct reduced_grav_tensor {

  /* 0th order terms */
  MyFloat F_000;

  /* 1st order terms */
  MyFloat F_100;
  MyFloat F_010;
  MyFloat F_001;
};

#ifdef WITH_MPI
/* MPI datatypes for transfers */
extern MPI_Datatype multipole_mpi_type;
extern MPI_Op multipole_mpi_reduce_op;

void multipole_create_mpi_types(void);
void multipole_free_mpi_types(void);
#endif

#endif /* SWIFT_MULTIPOLE_STRUCT_H */
