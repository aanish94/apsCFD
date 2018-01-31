#!/usr/bin/env python

""" Defines fluids and provides an interface to access properties.

"""

import CoolProp.CoolProp as CP

# Fluids (http://www.coolprop.org/dev/fluid_properties/PurePseudoPure.html#list-of-fluids)
AIR = "Air"
HELIUM = "Helium"
HYDROGEN = "Hydrogen"
METHANE = "Methane"
NITROGEN = "Nitrogen"
OXYGEN = "Oxygen"
XENON = 'Xenon'
VALID_FLUIDS = [AIR, HELIUM, HYDROGEN, METHANE, NITROGEN, OXYGEN, XENON]


class Fluid(object):
    """Define the properties of a generic fluid mixture.

    Fluid properties (http://www.coolprop.org/dev/coolprop/HighLevelAPI.html#id18)
    """

    def __init__(self, names, mole_fractions=None):
        """
        :param <list> name: Fluid Names
        :param <list> mole_fractions: Mole fractions
        """

        # Check fluid names are valid
        for name in names:
            if name not in VALID_FLUIDS:
                raise ValueError("Unsupported fluid : {}".format(name))

        # Single-fluid case
        if not mole_fractions:
            if len(names) == 1:
                mole_fractions = [1]
            else:
                raise ValueError("No mole fractions specified")

        if len(names) != len(mole_fractions):
            raise ValueError("Mis-matching size of input arrays!")

        # Check mole fractions are valid - below 1
        mole_fraction_sum = 0
        for mole_fraction in mole_fractions:
            if mole_fraction > 1:
                raise ValueError("Mole fraction must be less than!")

            mole_fraction_sum += mole_fraction

        # Check mole fractions sum to 1 - by definition
        if round(mole_fraction_sum, 1) != 1.0:
            raise ValueError(
                "Mole fractions must sum to 1! Instead : {}".format(mole_fraction_sum))

        self.names = names
        self.mole_fractions = mole_fractions
        self.name = self.construct_fluid_name(names, mole_fractions)

    def construct_fluid_name(self, names, mole_fractions):
        """Encode mixture components & mole fractions into CoolProp format.

        :param <list> name: Fluid Names
        :param <list> mole_fractions: Mole fractions

        :return <str> mixture_name: CoolProp mixture name format
        """

        # Mixture name format must be like: HEOS::R32[0.697615]&R125[0.302385]
        mixture_name = ""

        # Single-fluid case
        if len(names) == 1:
            return names[0]

        # Iterate through each fluid
        for idx, name in enumerate(names):
            # Construct string
            mixture_name += name
            mole_fraction = str(mole_fractions[idx])
            mixture_name += '[{}]'.format(mole_fraction)

            # Add & everytime except for the last element
            if idx < len(names) - 1:
                mixture_name += '&'

        return mixture_name

    def density(self, T, P):
        """Calculate mass density.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <float> D: Mass density [kg/m^3]
        """

        return CP.PropsSI('D', 'T', T, 'P', P, self.name)

    def enthalpy(self, T, P):
        """Calculate mass specific enthalpy.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <float> H: Enthalpy [J/kg]
        """

        return CP.PropsSI('H', 'T', T, 'P', P, self.name)

    def entropy(self, T, P):
        """Calculate mass specific entropy.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <int> S: Mass specific entropy [J/kg]
        """

        return CP.PropsSI('S', 'T', T, 'P', P, self.name)

    def internal_energy(self, T, P):
        """Calculate internal energy.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <int> U: Mass specific internal energy [J/kg]
        """

        return CP.PropsSI('U', 'T', T, 'P', P, self.name)

    def molar_mass(self, T, P):
        """Calculate molar mass.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <int> M: Molar Mass [kg/mol]
        """

        return CP.PropsSI('M', 'T', T, 'P', P, self.name)

    def pressure(self, D, T):
        """Calculate pressure.

        :param <int> D: Density [kg/m^3]
        :param <int> T: Temperature [K]

        :return <int> P: Pressure [Pa]
        """

        return CP.PropsSI('P', 'D', D, 'T', T, self.name)

    def specific_heat_cp(self, T, P):
        """Calculate mass specific constant pressure specific heat.
        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <int> CP: Specific heat constant pressure [J/kg/K]
        """

        return CP.PropsSI('C', 'T', T, 'P', P, self.name)

    def specific_heat_cv(self, T, P):
        """Calculate mass specific constant volume specific heat.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <int> CV: Specific heat constant volume [J/kg/K]
        """

        return CP.PropsSI('CVMASS', 'T', T, 'P', P, self.name)

    def speed_of_sound(self, T, P):
        """Calculate the speed of sound.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <float> A: Speed of Sound [m/s]
        """

        return CP.PropsSI('A', 'T', T, 'P', P, self.name)

    def temperature(self, D, P):
        """Calculate temperature.

        :param <float> D: Density [K]
        :param <float> P: Pressure [Pa]

        :return <float> T: Temperature [K]
        """

        return CP.PropsSI('T', 'D', D, 'P', P, self.name)

    def thermal_conductivity(self, T, P):
        """Calculate thermal conductivity.

        :param <float> T: Temperature [K]
        :param <float> P: Pressure [Pa]

        :return <float> K: Thermal conductivity [W/m/K]
        """

        return CP.PropsSI('CONDUCTIVITY', 'T', T, 'P', P, self.name)

    def viscosity(self, T, P):
        """Calculate viscosity.

        :param <int> T: Temperature [K]
        :param <int> P: Pressure [Pa]

        :return <int> V: Viscosity [Pa s]
        """

        return CP.PropsSI('V', 'T', T, 'P', P, self.name)
