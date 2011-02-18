within Buildings.Fluid.Sensors.BaseClasses;
partial model PartialFlowSensor
  "Partial component to model sensors that measure flow properties"
  extends Modelica.Fluid.Interfaces.PartialTwoPort;

  parameter Medium.MassFlowRate m_flow_nominal(min=0)
    "Nominal mass flow rate, used for regularization near zero flow"
    annotation(Dialog(group = "Nominal condition"));
  parameter Medium.MassFlowRate m_flow_small(min=0) = 1E-4*m_flow_nominal
    "For bi-directional flow, temperature is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
    annotation(Dialog(group="Advanced"));

equation
  // mass balance
  0 = port_a.m_flow + port_b.m_flow;

  // momentum equation (no pressure loss)
  port_a.p = port_b.p;

  // isenthalpic state transformation (no storage and no loss of energy)
  port_a.h_outflow = inStream(port_b.h_outflow);
  port_b.h_outflow = inStream(port_a.h_outflow);

  port_a.Xi_outflow = inStream(port_b.Xi_outflow);
  port_b.Xi_outflow = inStream(port_a.Xi_outflow);

  port_a.C_outflow = inStream(port_b.C_outflow);
  port_b.C_outflow = inStream(port_a.C_outflow);
  annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.<br>
The model includes zero-volume balance equations. Sensor models inheriting from
this partial class should add a medium instance to calculate the measured property.
</p>
</html>"),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics),
    Icon(graphics));
end PartialFlowSensor;