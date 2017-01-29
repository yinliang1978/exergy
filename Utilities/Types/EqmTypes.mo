within Exergy.Utilities.Types;
type EqmTypes = enumeration(
    TM " thermo–mechanical equilibrium",
    TMC "thermo-mechanical-chemical equilibrium")
  "Enumeration defining equilibrium when evaluating exergy"
annotation (Evaluate=true);
