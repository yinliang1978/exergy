within Exergy.Utilities.Types;
type EnergyTypes = enumeration(
    Enthalpy,
    InternalEnergy,
    Power,
    HeatTransfer)
  "Enumeration defining Reference State Of Water when evaluating exergy"
  annotation (Evaluate=true);
