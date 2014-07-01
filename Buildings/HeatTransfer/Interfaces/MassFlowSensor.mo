within Buildings.HeatTransfer.Interfaces;
model MassFlowSensor "Mass flow rate sensor"
  extends Modelica.Icons.RotationalSensor;

Modelica.Blocks.Interfaces.RealOutput m_flow(unit="kg/s")
    "Mass flow from port_a to port_b as output signal"
                                                     annotation (Placement(
        transformation(
        origin={0,-100},
        extent={{-10,-10},{10,10}},
        rotation=270)));

  HeatMassPort_a heatMassPort_a annotation (Placement(transformation(extent={{-102,
            -10},{-82,10}}), iconTransformation(extent={{-102,-10},{-82,10}})));
  HeatMassPort_b heatMassPort_b annotation (Placement(transformation(extent={{80,
            -10},{100,10}}), iconTransformation(extent={{80,-10},{100,10}})));
equation
    heatMassPort_a.heatPort.T = heatMassPort_b.heatPort.T;
    heatMassPort_a.massPort.Xi_outflow =heatMassPort_b.massPort.Xi_outflow;
  heatMassPort_a.massPort.m_flow + heatMassPort_b.massPort.m_flow = 0;
  m_flow = heatMassPort_a.massPort.m_flow;
  heatMassPort_a.heatPort.Q_flow = heatMassPort_b.heatPort.Q_flow;
 annotation (
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Line(points={{-70,0},{-95,0}}, color={170,170,255}),
        Line(points={{0,-70},{0,-90}}, color={0,0,127}),
        Line(points={{94,0},{69,0}}, color={170,170,255})}),
    Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                        graphics={
        Text(
          extent={{5,-86},{116,-110}},
          lineColor={0,0,0},
          textString="m_flow"),
        Line(points={{-70,0},{-90,0}}, color={170,170,255}),
        Line(points={{69,0},{90,0}}, color={170,170,255}),
        Line(points={{0,-70},{0,-90}}, color={0,0,127}),
        Text(
          extent={{-150,125},{150,85}},
          textString="%name",
          lineColor={0,0,255})}),
    Documentation(info="<HTML>
<p>
This model is capable of monitoring the mass flow rate flowing through
this component. The sensed value of mass flow rate is the amount that
passes through this sensor while keeping the humidity drop across the
sensor zero.  This is an ideal model.
The output signal is positive, if the mass flows from port_a to port_b.
</p>
</html>"));
end MassFlowSensor;
