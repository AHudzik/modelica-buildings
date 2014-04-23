within Buildings.HeatTransfer.Interfaces;
connector MassPort
    Modelica.SIunits.MassFraction Xi_outflow
    "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";

  flow Modelica.SIunits.MassFlowRate m_flow
    "Mass flow rate from the connection point into the component";

end MassPort;
