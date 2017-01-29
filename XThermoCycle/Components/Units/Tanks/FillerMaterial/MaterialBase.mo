within Exergy.XThermoCycle.Components.Units.Tanks.FillerMaterial;
record MaterialBase
  "Solid media properties of sensible heat storage materials from {Atear OE. Storage of Thermal Energy. Encyclopedia of Life support systems (EOLSS); 2008}"
parameter Modelica.SIunits.Density rho_sol = 2707
    "Density of the filler material";
parameter Modelica.SIunits.SpecificHeatCapacity cp_sol = 896
    "Specific heat capacity of the filler material";
parameter Modelica.SIunits.ThermalConductivity k_sol = 204
    "Thermal conductivity of the filler material";
end MaterialBase;
