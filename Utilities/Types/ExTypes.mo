within Exergy.Utilities.Types;
type ExTypes = enumeration(
    Flow "flow EnthalpyExergy",
    Vol "volume InternalEnergyExergy")
  "Enumeration defining what exergy is referring to when evaluating exergy"
annotation (Evaluate=true);
