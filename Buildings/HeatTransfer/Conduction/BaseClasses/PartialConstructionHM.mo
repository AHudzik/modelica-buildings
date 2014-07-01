within Buildings.HeatTransfer.Conduction.BaseClasses;
model PartialConstructionHM "Partial model for multi-layer constructions"
  extends Buildings.BaseClasses.BaseIcon;
  parameter Modelica.SIunits.Area A "Heat transfer area";

  replaceable parameter
    Buildings.HeatTransfer.Data.OpaqueConstructions.GenericHM
    layers "Construction definition from Data.OpaqueConstructions"
    annotation (Evaluate=true, choicesAllMatching=true, Placement(transformation(extent={{60,60},
            {80,80}})));

  final parameter Integer nLay(min=1, fixed=true) = layers.nLay
    "Number of layers";
  final parameter Integer nSta[nLay](min=1)={layers.material[i].nSta for i in 1:nLay}
    "Number of states"  annotation(Evaluate=true);

  parameter Modelica.SIunits.Temperature T_a_start=293.15
    "Initial temperature at port_a, used if steadyStateInitial = false"
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  parameter Modelica.SIunits.Temperature T_b_start=293.15
    "Initial temperature at port_b, used if steadyStateInitial = false"
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics));

end PartialConstructionHM;
