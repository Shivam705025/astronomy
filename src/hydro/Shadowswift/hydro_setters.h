/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_SHADOWSWIFT_HYDRO_SETTERS_H
#define SWIFT_SHADOWSWIFT_HYDRO_SETTERS_H

/**
 * @brief Set the primitive variables for the given particle to the given
 * values.
 *
 * @param p Particle.
 * @param W Primitive variables.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_set_primitive_variables(struct part* restrict p, const float* W) {

  p->rho = W[0];
  p->v[0] = W[1];
  p->v[1] = W[2];
  p->v[2] = W[3];
  p->P = W[4];
  p->A = W[5];
}

/**
 * @brief Set the conserved variables for the given particle to the given
 * values.
 *
 * @param p Particle.
 * @param Q Conserved variables.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_set_conserved_variables(struct part* restrict p, const float* Q) {

  p->conserved.mass = Q[0];
  p->conserved.momentum[0] = Q[1];
  p->conserved.momentum[1] = Q[2];
  p->conserved.momentum[2] = Q[3];
  p->conserved.energy = Q[4];
  p->conserved.entropy = Q[5];
}

/**
 * @brief Sets the mass of a particle
 *
 * @param p The particle of interest
 * @param m The mass to set.
 */
__attribute__((always_inline)) INLINE static void hydro_set_mass(
    struct part* restrict p, float m) {

  p->conserved.mass = m;
}

/**
 * @brief Sets the time derivative of the co-moving internal energy of a
 * particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the comoving internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy_dt(struct part* restrict p,
                                      const float du_dt) {
  error("Needs implementing");
}

/**
 * @brief Sets the time derivative of the physical internal energy of a particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param du_dt The time derivative of the physical internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_dt(struct part* restrict p,
                                      const struct cosmology* restrict cosmo,
                                      const float du_dt) {
  error("Needs implementing");
}
/**
 * @brief Sets the physical entropy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param entropy The physical entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_physical_entropy(
    struct part* p, struct xpart* xp, const struct cosmology* cosmo,
    const float entropy) {

  error("Needs implementing");
}

/**
 * @brief Sets the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy(struct part* p, struct xpart* xp,
                                   const struct cosmology* cosmo,
                                   const float u) {
  error("Need implementing");
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(struct part* p,
                                           const struct cosmology* cosmo,
                                           const float u) {
  error("Need implementing");
}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    struct part* restrict p, float alpha) {
  /* Purposefully left empty */
}

/**
 * @brief Update the value of the viscosity alpha to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_diffusive_feedback_reset(struct part* restrict p) {
  /* Purposefully left empty */
}

/**
 * @brief Modifies the thermal state of a particle to the imposed internal
 * energy
 *
 * This overrides the current state of the particle but does *not* change its
 * time-derivatives
 * NOTE: This function may violate energy conservation.
 *
 * @param p The particle
 * @param u The new internal energy
 */
__attribute__((always_inline)) INLINE static void hydro_set_internal_energy(
    struct part* restrict p, float u) {

  /* conserved.energy is NOT the specific energy (u), but the total thermal
     energy (u*m) */
  p->thermal_energy = u * p->conserved.mass;

  /* Update the total energy */
  p->conserved.energy =
      p->thermal_energy + 0.5f * p->conserved.mass *
                              (p->conserved.momentum[0] * p->v[0] +
                               p->conserved.momentum[1] * p->v[1] +
                               p->conserved.momentum[2] * p->v[2]);

  p->P = gas_pressure_from_internal_energy(p->rho, u);
  p->A = gas_entropy_from_internal_energy(p->rho, u);
}

/**
 * @brief Modifies the thermal state of a particle to the imposed entropy
 *
 * This overrides the current state of the particle but does *not* change its
 * time-derivatives
 *
 * @param p The particle
 * @param A The new entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_entropy(
    struct part* restrict p, float A) {

  p->thermal_energy = A * pow_gamma_minus_one(p->rho) *
                      hydro_one_over_gamma_minus_one * p->conserved.mass;
  /* add the kinetic energy */
  p->conserved.energy =
      p->thermal_energy + 0.5f * p->conserved.mass *
                              (p->conserved.momentum[0] * p->v[0] +
                               p->conserved.momentum[1] * p->v[1] +
                               p->conserved.momentum[2] * p->v[2]);

  p->P = gas_pressure_from_entropy(p->rho, A);
  p->A = A;
}

/**
 * @brief Overwrite the initial internal energy of a particle.
 *
 * Note that in the cases where the thermodynamic variable is not
 * internal energy but gets converted later, we must overwrite that
 * field. The conversion to the actual variable happens later after
 * the initial fake time-step.
 *
 * @param p The #part to write to.
 * @param u_init The new initial internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_init_internal_energy(struct part* p, float u_init) {

  /* We store the initial energy per unit mass in the energy
   * variable as the conversion to energy will be done later,
   * in hydro_first_init_part(). */
  p->conserved.energy = u_init;
}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_SETTERS_H */