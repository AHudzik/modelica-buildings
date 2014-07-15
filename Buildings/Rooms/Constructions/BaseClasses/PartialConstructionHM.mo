within Buildings.Rooms.Constructions.BaseClasses;
partial model PartialConstructionHM
  "Partial model for exterior construction that has no window"

 parameter Modelica.SIunits.Area A "Heat and mass transfer area";
 parameter Modelica.SIunits.Area AOpa
    "Heat and mass transfer area of opaque construction"
   annotation (Dialog(group="Opaque construction"));

  parameter Integer switch_w "switch for the water content"
    annotation (Dialog(tab="HM"));
  parameter Integer switch_lamb "switch for the heat conductivity"
    annotation (Dialog(tab="HM"));
  parameter Integer switch_dw
    "switch for the liquid transport coefficient for suction "
    annotation (Dialog(tab="HM"));

 parameter Boolean activatesuction;

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
    annotation (Placement(transformation(extent={{-308,194},{-294,208}}),
        iconTransformation(extent={{-322,184},{-284,222}})));
  HeatTransfer.Interfaces.HeatMassPort_b opa_b
    annotation (Placement(transformation(extent={{292,194},{306,208}}),
        iconTransformation(extent={{282,182},{322,220}})));

   final parameter Integer nLay(min=1, fixed=true) = layers.nLay
    "Number of layers";
  final parameter Integer nSta[nLay](min=1)={layers.material[i].nSta for i in 1:nLay}
    "Number of states"  annotation(Evaluate=true);

  HeatTransfer.Conduction.MultiLayerHM opa(
    final A=AOpa,
    final layers=layers,
    final switch_w=switch_w,
    final switch_lamb=switch_lamb,
    final switch_dw=switch_dw,
    final activatesuction=activatesuction)
    "Model for heat transfer through opaque construction"
    annotation (Placement(transformation(extent={{-64,142},{54,260}})));

equation
  connect(opa_a, opa.heatMassPort_a) annotation (Line(
      points={{-301,201},{-64,201}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(opa.heatMassPort_b, opa_b) annotation (Line(
      points={{54,201},{299,201}},
      color={127,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-300,
            -300},{300,300}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-300,-300},{300,300}}), graphics={
        Rectangle(
          extent={{-64,258},{-48,144}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{46,258},{62,144}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-294,200},{296,204}},
          lineColor={0,0,0},
          fillColor={127,0,127},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{12,208},{14,208},{16,204},{18,198},{18,194},{16,188},{10,184},
              {4,182},{-2,178},{-6,182},{-12,188},{-14,198},{-12,204},{-10,208},
              {-8,210},{-6,212},{-2,214},{2,214},{8,212},{12,208}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-50,228},{46,228}},
          color={0,0,255},
          pattern=LinePattern.Dash,
          smooth=Smooth.None),
        Line(
          points={{-50,172},{48,172}},
          color={0,0,255},
          pattern=LinePattern.Dash,
          smooth=Smooth.None)}));
end PartialConstructionHM;
