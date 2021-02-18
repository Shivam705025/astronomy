/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_SPHENIX_LOGGER_HYDRO_H
#define SWIFT_SPHENIX_LOGGER_HYDRO_H

#include "../config.h"

/* local includes */
#include "hydro.h"
#include "hydro_logger.h"
#include "logger_interpolation.h"
#include "logger_loader_io.h"
#include "logger_python_tools.h"

/* Index of the mask in the header mask array */
extern int hydro_logger_local_to_global[hydro_logger_field_count];

/* Size for each mask */
extern const int hydro_logger_field_size[hydro_logger_field_count];

/**
 * @brief Create the link between the fields and their derivatives.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void
hydro_logger_reader_link_derivatives(struct header *head) {

  /* Set the first and second derivatives */
  const int pos_id =
      hydro_logger_local_to_global[hydro_logger_field_coordinates];
  const int vel_id =
      hydro_logger_local_to_global[hydro_logger_field_velocities];
  const int acc_id =
      hydro_logger_local_to_global[hydro_logger_field_accelerations];

  /* Coordinates */
  header_set_first_derivative(head, pos_id, vel_id);
  header_set_second_derivative(head, pos_id, acc_id);

  /* Velocities */
  header_set_first_derivative(head, vel_id, acc_id);
}

/**
 * @brief Interpolate a field of the particle at the given time.
 * Here we use a linear interpolation for most of the fields.
 * For the position (velocity), we use a quintic (cubic) hermite interpolation
 * based on the positions, velocities and accelerations at the time of the two
 * particles.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 * @param field The field to reconstruct (follows the order of
 * #hydro_logger_fields).
 */
__attribute__((always_inline)) INLINE static void
hydro_logger_interpolate_field(const double t_before,
                               const struct logger_field *restrict before,
                               const double t_after,
                               const struct logger_field *restrict after,
                               void *restrict output, const double t,
                               const int field) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the times */
  if (t_before > t || t_after < t) {
    error(
        "The times for the interpolation are not "
        "monotonically increasing %g < %g < %g.",
        t_before, t, t_after);
  }
#endif

  switch (field) {
    /* Do the position */
    case hydro_logger_field_coordinates:
      interpolate_quintic_double_float_ND(t_before, before, t_after, after,
                                          output, t, /* dimension= */ 3);
      break;
      /* Do the velocity */
    case hydro_logger_field_velocities:
      interpolate_cubic_float_ND(t_before, before, t_after, after, output, t,
                                 /* dimension= */ 3);
      break;
    case hydro_logger_field_accelerations:
      interpolate_linear_float_ND(t_before, before, t_after, after, output, t,
                                  /* dimension= */ 3);
      break;
      /* Do the linear interpolation of float. */
    case hydro_logger_field_masses:
    case hydro_logger_field_smoothing_lengths:
    case hydro_logger_field_internal_energies:
    case hydro_logger_field_densities:
      interpolate_linear_float(t_before, before, t_after, after, output, t);
      break;
      /* Check the ids */
    case hydro_logger_field_particle_ids:
      interpolate_ids(t_before, before, t_after, after, output, t);
      break;

    case hydro_logger_field_secondary_fields: {
      /* Get some variables */
      const float wa = (t - t_before) / (t_after - t_before);
      const float wb = 1. - wa;
      float *x = (float *)output;
      const float *bef = (float *)before->field;
      const float *aft = (float *)after->field;

      /* Entropy + pressure + Viscosity + Diffusion + Laplacian*/
      const int n_linear = 5;
      for (int i = 0; i < n_linear; i++) {
        x[i] = wa * aft[i] + wb * bef[i];
      }

      /* Div v */
      x[n_linear] = interpolate_cubic_hermite_spline(
          t_before, bef[n_linear], bef[n_linear + 1], t_after, aft[n_linear],
          aft[n_linear + 1], t);

      /* d Div v / dt */
      x[n_linear + 1] = wa * aft[n_linear + 1] + wb * bef[n_linear + 1];
      break;
    }

    default:
      error("Not implemented");
  }
}

#ifdef HAVE_PYTHON
__attribute__((always_inline)) INLINE static void hydro_logger_generate_python(
    struct logger_python_field *fields) {

  fields[hydro_logger_field_coordinates] =
      logger_loader_python_field(/* Dimension */ 3, NPY_DOUBLE);
  fields[hydro_logger_field_velocities] =
      logger_loader_python_field(/* Dimension */ 3, NPY_FLOAT32);
  fields[hydro_logger_field_accelerations] =
      logger_loader_python_field(/* Dimension */ 3, NPY_FLOAT32);
  fields[hydro_logger_field_masses] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);
  fields[hydro_logger_field_smoothing_lengths] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);
  fields[hydro_logger_field_internal_energies] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);
  fields[hydro_logger_field_particle_ids] =
      logger_loader_python_field(/* Dimension */ 1, NPY_LONGLONG);
  fields[hydro_logger_field_densities] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);

  /* Now separate the secondary fields */
  fields[hydro_logger_field_secondary_fields] =
      logger_loader_python_field(/* Dimension */ 7, CUSTOM_NPY_TYPE);

  logger_loader_python_field_add_subfield(
      &fields[hydro_logger_field_secondary_fields], "Entropies", /* offset */ 0,
      /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[hydro_logger_field_secondary_fields], "Pressures",
      /* offset */ sizeof(float), /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[hydro_logger_field_secondary_fields], "ViscosityParameters",
      /* offset */ 2 * sizeof(float), /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[hydro_logger_field_secondary_fields], "DiffusionParameters",
      /* offset */ 3 * sizeof(float), /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[hydro_logger_field_secondary_fields], "LaplacianInternalEnergies",
      /* offset */ 4 * sizeof(float), /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[hydro_logger_field_secondary_fields], "VelocityDivergences",
      /* offset */ 5 * sizeof(float), /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[hydro_logger_field_secondary_fields],
      "VelocityDivergenceTimeDifferentials", /* offset */ 6 * sizeof(float),
      /* Dimension */ 1, NPY_FLOAT32);
}

#endif  // HAVE_PYTHON
#endif  /* SWIFT_SPHENIX_LOGGER_HYDRO_H */
