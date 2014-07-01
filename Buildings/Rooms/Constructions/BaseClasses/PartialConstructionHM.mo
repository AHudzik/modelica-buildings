within Buildings.Rooms.Constructions.BaseClasses;
partial model PartialConstructionHM
  "Partial model for exterior construction that has no window"

 parameter Modelica.SIunits.Area A "Heat and mass transfer area";
 parameter Modelica.SIunits.Area AOpa
    "Heat and mass transfer area of opaque construction"
   annotation (Dialog(group="Opaque construction"));

  replaceable parameter
    Buildings.HeatTransfer.Data.OpaqueConstructions.GenericHM
    layers "Material properties of opaque construction"
    annotation(Dialog(group="Opaque construction"),
               Evaluate=true, choicesAllMatching=true, Placement(transformation(extent={{146,258},
            {166,278}})));
  parameter Modelica.SIunits.Angle til "Surface tilt";

  final parameter Boolean isFloor=til > 2.74889125 and til < 3.53428875
    "Flag, true if construction is a floor" annotation (Evaluate=true);
  final parameter Boolean isCeiling=til > -0.392699 and til < 0.392699
    "Flag, true if construction is a floor" annotation (Evaluate=true);

  HeatTransfer.Interfaces.HeatMassPort_a opa_a
    annotation (Placement(transformation(extent={{-108,-6},{-94,8}})));
  HeatTransfer.Interfaces.HeatMassPort_b opa_b
    annotation (Placement(transformation(extent={{92,-6},{106,8}})));

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

  HeatTransfer.Conduction.MultiLayerHM opa(
    final A=AOpa,
    final layers=layers,
    final T_a_start=T_a_start,
    final T_b_start=T_b_start)
    "Model for heat transfer through opaque construction"
    annotation (Placement(transformation(extent={{-30,-26},{24,28}})));

equation
  connect(opa_a, opa.heatMassPort_a) annotation (Line(
      points={{-101,1},{-30,1}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(opa.heatMassPort_b, opa_b) annotation (Line(
      points={{24,1},{99,1}},
      color={127,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-28,40},{-20,-38}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{20,40},{28,-40}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-100,2},{96,0}},
          lineColor={0,0,0},
          fillColor={127,0,127},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{12,14},{14,14},{16,10},{18,4},{18,0},{16,-6},{10,-10},{4,-12},
              {-2,-16},{-6,-12},{-12,-6},{-14,4},{-12,10},{-10,14},{-8,16},{-6,
              18},{-2,20},{2,20},{8,18},{12,14}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-20,30},{20,30}},
          color={0,0,255},
          pattern=LinePattern.Dash,
          smooth=Smooth.None),
        Line(
          points={{-20,-20},{20,-20}},
          color={0,0,255},
          pattern=LinePattern.Dash,
          smooth=Smooth.None)}));
end PartialConstructionHM;
